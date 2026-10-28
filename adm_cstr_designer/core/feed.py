"""Substrate blend -> ADM1 influent vector."""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from .adm1 import IX, LIQ

COD_PER_NM3_CH4 = 64.0 / 22.414  # kg COD per Nm3 CH4 = 2.855


@dataclass
class FeedInputs:
    Q_tpd: float                 # wet substrate mass flow, t/d (incl. any dilution water)
    T_feed_C: float
    TS_pct: float                # % of wet mass, blend as fed to the digester
    VS_pct_TS: float             # % of TS
    blend: dict                  # {feedstock name: % of the blend (by VS)}
    alk_cat_kmol_m3: float = 0.045   # net cations (S_cat - S_an) in feed water/buffer
    S_IC_kmol_m3: float = 0.02       # inorganic carbon in feed


def normalise_blend(blend: dict) -> dict:
    tot = sum(max(v, 0.0) for v in blend.values())
    if tot <= 0:
        raise ValueError('Blend must contain at least one feedstock > 0 %')
    return {k: max(v, 0.0) / tot for k, v in blend.items() if v > 0}


def slurry_density(TS_pct: float) -> float:
    return 1000.0 + 6.0 * TS_pct         # kg/m3, indicative


def slurry_cp(TS_pct: float) -> float:
    return 4.18 * (1 - TS_pct / 100) + 1.6 * TS_pct / 100   # kJ/kg/K


def blend_properties(fi: FeedInputs, lib: dict) -> dict:
    """Blend composition, COD split, hydrolysis constants and ADM1 feed vector."""
    w = normalise_blend(fi.blend)                # VS-based fractions
    codvs = sum(w[k] * lib[k]['COD_per_VS'] for k in w)
    bmp = sum(w[k] * lib[k]['BMP_Nm3CH4_per_kgVS'] for k in w)      # per kg VS
    # fraction of each feedstock's COD that is biodegradable (BMP-consistent, 5 % to biomass)
    fbd = {k: min(0.95, lib[k]['BMP_Nm3CH4_per_kgVS'] * COD_PER_NM3_CH4 / (lib[k]['COD_per_VS'] * 0.95))
           for k in w}
    cod_w = {k: w[k] * lib[k]['COD_per_VS'] for k in w}               # COD share basis
    cod_tot = sum(cod_w.values())
    cod_share = {k: cod_w[k] / cod_tot for k in w}
    # COD split of the blend
    f_ch = sum(cod_share[k] * fbd[k] * lib[k]['frac_ch_pr_li'][0] for k in w)
    f_pr = sum(cod_share[k] * fbd[k] * lib[k]['frac_ch_pr_li'][1] for k in w)
    f_li = sum(cod_share[k] * fbd[k] * lib[k]['frac_ch_pr_li'][2] for k in w)
    vfa = sum(cod_share[k] * lib[k]['vfa_frac_COD'] for k in w)
    f_inert = 1.0 - (f_ch + f_pr + f_li)
    # blend hydrolysis constants weighted by the COD of the substrate hydrolysed
    def kw(idx, key):
        num = sum(cod_share[k] * fbd[k] * lib[k]['frac_ch_pr_li'][idx] * lib[k][key] for k in w)
        den = sum(cod_share[k] * fbd[k] * lib[k]['frac_ch_pr_li'][idx] for k in w)
        return num / den
    hyd = dict(k_hyd_ch=kw(0, 'k_hyd_ch'), k_hyd_pr=kw(1, 'k_hyd_pr'), k_hyd_li=kw(2, 'k_hyd_li'))
    tan = sum(w[k] * lib[k]['TAN_pct_TS'] for k in w) / 100.0          # kg N / kg TS (approx)
    ts_lib = sum(w[k] * lib[k]['TS_pct_wet'] for k in w)
    vs_lib = sum(w[k] * lib[k]['VS_pct_TS'] for k in w)

    rho = slurry_density(fi.TS_pct)
    Q_m3d = fi.Q_tpd * 1000.0 / rho
    VS_kg_m3 = rho * fi.TS_pct / 100 * fi.VS_pct_TS / 100
    COD_kg_m3 = VS_kg_m3 * codvs
    TAN_kg_N_m3 = rho * fi.TS_pct / 100 * tan

    feed = np.zeros(len(LIQ))
    s_ac = vfa * COD_kg_m3
    part = COD_kg_m3 - s_ac
    feed[IX['S_ac']] = s_ac
    scale = part / max(f_ch + f_pr + f_li + f_inert, 1e-12)
    feed[IX['X_ch']] = f_ch * scale
    feed[IX['X_pr']] = f_pr * scale
    feed[IX['X_li']] = f_li * scale
    feed[IX['X_I']] = f_inert * scale
    feed[IX['S_IN']] = TAN_kg_N_m3 / 14.0
    feed[IX['S_IC']] = fi.S_IC_kmol_m3
    feed[IX['S_cat']] = fi.alk_cat_kmol_m3
    feed[IX['S_an']] = 0.0
    return dict(feed=feed, hyd=hyd, Q_m3d=Q_m3d, rho=rho, VS_kg_m3=VS_kg_m3,
                COD_kg_m3=COD_kg_m3, COD_per_VS=codvs, BMP=bmp, weights=w,
                fractions=dict(ch=f_ch, pr=f_pr, li=f_li, inert=f_inert, vfa=vfa),
                TS_lib=ts_lib, VS_lib=vs_lib, TAN_kg_m3=TAN_kg_N_m3,
                VS_load_kg_d=Q_m3d * VS_kg_m3,
                theoretical_CH4_Nm3d=Q_m3d * VS_kg_m3 * bmp)
