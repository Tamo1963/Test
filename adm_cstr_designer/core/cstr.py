"""Anaerobic CSTR: mass balances, headspace, transient + steady-state solve."""
from __future__ import annotations

from dataclasses import dataclass, field
import numpy as np
from scipy.integrate import solve_ivp

from .adm1 import (ADM1Params, IX, LIQ, NL, NG, STATES, rates, NM3_PER_KMOL, P_ATM, R_GAS)


@dataclass
class CSTRSpec:
    V_liq: float                 # m3
    V_gas: float                 # m3
    Q: float                     # m3/d
    T_C: float                   # degC
    kp: float | None = None      # gas pipe resistance, m3/d/bar (auto-scaled if None)


def default_initial_state(feed: np.ndarray, seed_biomass: float = 1.0, p: ADM1Params | None = None) -> np.ndarray:
    """Seeded digester: biomass present, substrates low, headspace mostly CH4/CO2."""
    y = np.zeros(len(STATES))
    y[:NL] = feed
    for n in ('S_su', 'S_aa', 'S_fa', 'S_va', 'S_bu', 'S_pro', 'S_h2', 'X_c', 'X_ch', 'X_pr', 'X_li'):
        y[IX[n]] = 0.0
    y[IX['S_ac']] = 0.05
    for n in ('X_su', 'X_aa', 'X_fa', 'X_c4', 'X_pro', 'X_ac', 'X_h2'):
        y[IX[n]] = seed_biomass * (1.0 if n != 'X_ac' else 1.5)
    y[IX['S_IC']] = max(feed[IX['S_IC']], 0.08)
    y[IX['S_IN']] = max(feed[IX['S_IN']], 0.02)
    y[IX['G_ch4']] = 0.06
    y[IX['G_co2']] = 0.015
    y[IX['G_h2']] = 1e-5
    if p is not None:      # dissolved gases in equilibrium with the headspace (no start-up transfer artefact)
        y[IX['S_ch4']] = 64 * p.KH_ch4 * (y[IX['G_ch4']] * R_GAS * p.T / 64)
        y[IX['S_h2']] = 16 * p.KH_h2 * (y[IX['G_h2']] * R_GAS * p.T / 16)
    return y


def make_rhs(spec: CSTRSpec, p: ADM1Params, feed: np.ndarray):
    kp = spec.kp if spec.kp is not None else 5.0e4 * spec.V_liq / 3400.0
    M = p.M
    D = spec.Q / spec.V_liq
    vr = spec.V_liq / spec.V_gas

    def rhs(t, y):
        rho, d = rates(y, p)
        dy = np.zeros_like(y)
        dy[:NL] = D * (feed - y[:NL]) + M.T @ rho
        q_gas = max(kp * (d['p_gas'] - P_ATM), 0.0)
        dy[NL:] = -y[NL:] * q_gas / spec.V_gas + rho[19:22] * vr
        return dy

    return rhs, kp


def gas_flows(rho: np.ndarray, V_liq: float) -> dict:
    """Biogas flows from interphase transfer rates, Nm3/d."""
    h2 = rho[19] * V_liq / 16 * NM3_PER_KMOL
    ch4 = rho[20] * V_liq / 64 * NM3_PER_KMOL
    co2 = rho[21] * V_liq * NM3_PER_KMOL
    tot = h2 + ch4 + co2
    return dict(h2=h2, ch4=ch4, co2=co2, total=tot, ch4_frac=ch4 / tot if tot > 0 else 0.0)


def simulate(spec: CSTRSpec, p: ADM1Params, feed: np.ndarray, y0: np.ndarray | None = None,
             t_end: float | None = None, n_out: int = 300, rtol=1e-6, atol=1e-9):
    y0 = default_initial_state(feed, p=p) if y0 is None else y0
    rhs, kp = make_rhs(spec, p, feed)
    t_end = t_end or max(400.0, 12 * spec.V_liq / spec.Q)
    sol = solve_ivp(rhs, (0, t_end), y0, method='BDF', rtol=rtol, atol=atol,
                    t_eval=np.linspace(0, t_end, n_out))
    return sol, rhs


def steady_state(spec: CSTRSpec, p: ADM1Params, feed: np.ndarray, y0=None, max_rounds=4):
    """Integrate until d(y)/dt is negligible; return state, diagnostics, history."""
    y0 = default_initial_state(feed, p=p) if y0 is None else y0
    t_end = max(400.0, 12 * spec.V_liq / spec.Q)
    ts, ys = [], []
    toff = 0.0
    y = y0
    for _ in range(max_rounds):
        sol, rhs = simulate(spec, p, feed, y, t_end=t_end)
        ts.append(sol.t + toff)
        ys.append(sol.y.T)
        y = sol.y[:, -1]
        toff += t_end
        dy = rhs(0, y)
        scale = np.maximum(np.abs(y), 1e-3)
        if np.max(np.abs(dy) / scale) < 1e-4 / 1.0:
            converged = True
            break
    else:
        converged = False
    t = np.concatenate(ts)
    Y = np.vstack(ys)
    rho, diag = rates(y, p)
    return dict(y=y, rho=rho, diag=diag, converged=converged, t=t, Y=Y,
                gas=gas_flows(rho, spec.V_liq), rhs=rhs)


def history_series(res: dict, p: ADM1Params, spec: CSTRSpec, stride: int = 1) -> dict:
    """pH, CH4 flow, VFA time series from the stored trajectory."""
    t, Y = res['t'][::stride], res['Y'][::stride]
    pH, ch4, vfa = [], [], []
    for row in Y:
        rho, d = rates(row, p)
        pH.append(d['pH'])
        ch4.append(rho[20] * spec.V_liq / 64 * NM3_PER_KMOL)
        vfa.append(row[IX['S_ac']] + row[IX['S_pro']] + row[IX['S_bu']] + row[IX['S_va']])
    return dict(t=t, pH=np.array(pH), ch4=np.array(ch4), vfa=np.array(vfa))


def cod_balance(res: dict, spec: CSTRSpec, feed: np.ndarray) -> dict:
    """COD in vs COD out (effluent organics + CH4 + H2 to gas), kg COD/d.

    Effluent includes dissolved CH4/H2; gas CH4/H2 leave via transfer.
    """
    cod_ix = [IX[n] for n in LIQ if n.startswith(('S_', 'X_')) and n not in ('S_IC', 'S_IN')]
    cod_in = spec.Q * feed[cod_ix].sum()
    cod_eff = spec.Q * res['y'][cod_ix].sum()
    g = res['gas']
    cod_gas = (g['ch4'] / NM3_PER_KMOL * 64) + (g['h2'] / NM3_PER_KMOL * 16)
    return dict(cod_in=cod_in, cod_effluent=cod_eff, cod_gas=cod_gas,
                closure=(cod_eff + cod_gas) / cod_in if cod_in else float('nan'))
