"""End-to-end design: feed -> kinetics -> CSTR -> mixing -> heating."""
from __future__ import annotations

from dataclasses import dataclass, field, asdict
import numpy as np

from .adm1 import ADM1Params, IX, LIQ
from .cstr import CSTRSpec, steady_state, history_series, cod_balance
from .data import load_feedstocks, load_kinetics, clone
from .feed import FeedInputs, blend_properties, slurry_cp
from .geometry import tank_geometry
from .mixing import size_agitator
from .thermal import heat_load, size_coil


@dataclass
class DesignInputs:
    # --- the five user inputs ---
    Q_tpd: float = 100.0
    T_feed_C: float = 25.0
    VS_pct_TS: float = 85.0
    TS_pct: float = 8.0
    blend: dict = field(default_factory=lambda: {'Napier grass': 40, 'Cow dung': 40, 'Chicken litter': 20})
    # --- design basis (defaults agreed) ---
    HRT_d: float = 20.0
    T_op_C: float = 37.0
    T_amb_C: float = 15.0
    insul_t_m: float = 0.10
    n_tanks: int = 1
    H_over_D: float = 1.0
    alk_cat_kmol_m3: float = 0.045
    S_IC_kmol_m3: float = 0.02
    kinetics_set: str = 'ADM1_BSM2_mesophilic'
    theta: float | None = None
    eps_W_m3: float = 6.0
    mu_override_Pa_s: float | None = None
    T_w_in: float = 70.0
    T_w_out: float = 55.0
    U_coil_override: float | None = None
    k_overrides: dict = field(default_factory=dict)   # {'uptake.ac.k_m': 9.0, ...}


def apply_overrides(kset: dict, overrides: dict) -> dict:
    ks = clone(kset)
    for path, val in overrides.items():
        cur = ks
        parts = path.split('.')
        for q in parts[:-1]:
            cur = cur[q]
        cur[parts[-1]] = val
    return ks


def run_design(di: DesignInputs, kinetics: dict | None = None, feedstocks: dict | None = None,
               full_history: bool = True) -> dict:
    feedstocks = feedstocks or load_feedstocks()
    kinetics = kinetics or load_kinetics()
    kset = apply_overrides(kinetics[di.kinetics_set], di.k_overrides)

    fi = FeedInputs(di.Q_tpd, di.T_feed_C, di.TS_pct, di.VS_pct_TS, di.blend,
                    di.alk_cat_kmol_m3, di.S_IC_kmol_m3)
    fp = blend_properties(fi, feedstocks)
    V_liq = fp['Q_m3d'] * di.HRT_d
    geo = tank_geometry(V_liq, di.H_over_D, n_tanks=di.n_tanks)

    p = ADM1Params(kset, di.T_op_C, hyd=fp['hyd'], theta=di.theta)
    spec = CSTRSpec(V_liq=V_liq, V_gas=geo['V_gas_each'] * di.n_tanks, Q=fp['Q_m3d'], T_C=di.T_op_C)
    res = steady_state(spec, p, fp['feed'])
    cb = cod_balance(res, spec, fp['feed'])

    y = res['y']
    VFA = y[IX['S_ac']] + y[IX['S_pro']] + y[IX['S_bu']] + y[IX['S_va']]
    gas = res['gas']
    perf = dict(
        pH=res['diag']['pH'], VFA_kgCOD_m3=VFA, S_ac=y[IX['S_ac']], S_IN_mgN_L=y[IX['S_IN']] * 14e3,
        NH3_mgN_L=res['diag']['nh3'] * 14e3, I_nh3=res['diag']['I_nh3'], I_ac=res['diag']['I_ac'],
        biogas_Nm3_d=gas['total'], CH4_Nm3_d=gas['ch4'], CH4_pct=100 * gas['ch4_frac'],
        CO2_Nm3_d=gas['co2'], SMY_Nm3_kgVS=gas['ch4'] / fp['VS_load_kg_d'] if fp['VS_load_kg_d'] else 0,
        BMP_blend=fp['BMP'], conversion_pct=100 * gas['ch4'] / fp['theoretical_CH4_Nm3d'] if fp['theoretical_CH4_Nm3d'] else 0,
        OLR_kgVS_m3_d=fp['VS_load_kg_d'] / V_liq, VS_load_kg_d=fp['VS_load_kg_d'],
        COD_removal_pct=100 * (1 - (cb['cod_effluent']) / cb['cod_in']) if cb['cod_in'] else 0,
        CH4_energy_kW=gas['ch4'] * 35.8 / 3.6 / 24,       # LHV 35.8 MJ/Nm3 -> kW
        converged=res['converged'], cod_closure=cb['closure'],
        X_ac=y[IX['X_ac']], X_h2=y[IX['X_h2']],
    )
    flags = []
    if not res['converged']:
        flags.append('Steady state not reached within the integration horizon; treat results with caution.')
    if perf['pH'] < 6.6:
        flags.append(f"pH {perf['pH']:.2f} is below the methanogenic range: acidification / failure risk.")
    if VFA > 2.0:
        flags.append(f'VFA {VFA:.1f} kg COD/m3 is high: process unstable.')
    if perf['NH3_mgN_L'] > 80:
        flags.append(f"Free ammonia {perf['NH3_mgN_L']:.0f} mgN/L: acetoclastic inhibition risk (consider less litter).")
    if perf['OLR_kgVS_m3_d'] > 4.5:
        flags.append(f"OLR {perf['OLR_kgVS_m3_d']:.1f} kg VS/m3/d is high for a CSTR (typical 1-4).")
    if y[IX['X_ac']] < 0.05:
        flags.append('Acetoclastic biomass has washed out at this HRT.')
    if perf['S_IN_mgN_L'] < 30:
        flags.append(f"Inorganic N is only {perf['S_IN_mgN_L']:.0f} mgN/L: nitrogen-limited. Add an N-rich co-substrate (litter, dung) or urea.")
    if not 30 <= di.T_op_C <= 42:
        flags.append('Operating temperature is outside the mesophilic range the ADM1 constants were calibrated for; '
                     'temperature correction is an extrapolation. Calibrate before relying on these numbers.')
    if di.HRT_d < 15:
        flags.append('HRT < 15 d is short for lignocellulosic blends.')

    rho, cp = fp['rho'], slurry_cp(di.TS_pct)
    mix = size_agitator(geo['V_liq_each'], geo['D'], geo['H_liq'], di.TS_pct, rho,
                        mu=di.mu_override_Pa_s, eps_W_m3=di.eps_W_m3)
    heat = heat_load(fp['Q_m3d'] / di.n_tanks, rho, cp, di.T_feed_C, di.T_op_C, geo['D'], geo['H_liq'],
                     geo['H_free'], di.T_amb_C, t_ins=di.insul_t_m)
    T_w_in, T_w_out = di.T_w_in, di.T_w_out
    if T_w_out < di.T_op_C + 8.0:      # keep a workable approach temperature
        shift = di.T_op_C + 8.0 - T_w_out
        T_w_in, T_w_out = T_w_in + shift, T_w_out + shift
        flags.append(f'Hot-water temperatures raised to {T_w_in:.0f}/{T_w_out:.0f} C to keep an '
                     f'8 K approach to the digester temperature.')
    coil = size_coil(heat['Q_design_kW'], di.T_op_C, di.TS_pct, geo['D'], geo['H_liq'],
                     T_w_in, T_w_out, U=di.U_coil_override)

    energy = dict(
        biogas_kW=perf['CH4_energy_kW'], heat_kW=heat['Q_design_kW'] * di.n_tanks,
        mixing_kW=mix['P_shaft_kW'] * mix['n_agitators'] * di.n_tanks / 0.95)
    energy['net_kW'] = energy['biogas_kW'] - energy['heat_kW'] - energy['mixing_kW']
    energy['parasitic_pct'] = 100 * (energy['heat_kW'] + energy['mixing_kW']) / max(energy['biogas_kW'], 1e-9)

    out = dict(inputs=asdict(di), feed=fp, geometry=geo, performance=perf, state=dict(zip(LIQ, y[:len(LIQ)])),
               mixing=mix, heat=heat, coil=coil, energy=energy, flags=flags, kinetics=kset,
               cod_balance=cb, spec=spec, _params=p, _res=res)
    if full_history:
        out['history'] = history_series(res, p, spec, stride=max(1, len(res['t']) // 400))
    return out
