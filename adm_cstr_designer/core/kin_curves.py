"""Data for kinetic plots/tables shared by the app and the report."""
from __future__ import annotations

import numpy as np
import pandas as pd


def kinetics_table(kset: dict, T_C: float, theta: float | None = None) -> pd.DataFrame:
    th = kset.get('theta', 1.05) if theta is None else theta
    fT = th ** (T_C - kset['T_ref_C'])
    rows = []
    for key, v in kset['uptake'].items():
        rows.append({'Group': v['label'], 'k_m @ref (1/d)': v['k_m'], f'k_m @ {T_C:g} C (1/d)': v['k_m'] * fT,
                     'K_S (kg COD/m3)': v['K_S'], 'Y (kg/kg)': v['Y'], 'b (1/d)': v['b'] * fT,
                     'mu_max = Y*k_m (1/d)': v['Y'] * v['k_m'] * fT,
                     'Extra inhibition': ', '.join(f"{k}={v[k]:g}" for k in ('KI_h2', 'KI_nh3') if k in v)})
    return pd.DataFrame(rows)


def monod_curves(kset: dict, T_C: float, theta: float | None = None, n: int = 120) -> dict:
    th = kset.get('theta', 1.05) if theta is None else theta
    fT = th ** (T_C - kset['T_ref_C'])
    out = {}
    for key, v in kset['uptake'].items():
        S = np.linspace(0, max(5 * v['K_S'], 1e-6), n)
        out[v['label']] = (S, v['k_m'] * fT * S / (v['K_S'] + S))
    return out


def inhibition_curves(kset: dict) -> dict:
    pH = np.linspace(4, 9, 200)
    res = {}
    for name, (ll, ul) in kset['pH_limits'].items():
        res[f'pH inhibition ({name})'] = (pH, np.where(pH >= ul, 1.0, np.exp(-3 * ((pH - ul) / (ul - ll)) ** 2)))
    u = kset['uptake']
    h2 = np.logspace(-8, -3, 200)
    for k in ('fa', 'c4', 'pro'):
        res[f'H2 inhibition ({k})'] = (h2, 1 / (1 + h2 / u[k]['KI_h2']))
    nh3 = np.linspace(0, 0.01, 200)
    res['NH3 inhibition (acetoclastic)'] = (nh3, 1 / (1 + nh3 / u['ac']['KI_nh3']))
    return res
