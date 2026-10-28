"""Scenario sweeps: how the system responds to modification."""
from __future__ import annotations

import dataclasses
import pandas as pd

from .design import DesignInputs, run_design


def _row(label, d):
    p = d['performance']
    return dict(Scenario=label, V_liq_m3=d['spec'].V_liq, HRT_d=d['inputs']['HRT_d'],
                T_op_C=d['inputs']['T_op_C'], pH=round(p['pH'], 2), VFA=round(p['VFA_kgCOD_m3'], 2),
                NH3_mgN_L=round(p['NH3_mgN_L']), OLR=round(p['OLR_kgVS_m3_d'], 2),
                CH4_Nm3_d=round(p['CH4_Nm3_d']), CH4_pct=round(p['CH4_pct'], 1),
                SMY=round(p['SMY_Nm3_kgVS'], 3), Net_kW=round(d['energy']['net_kW'], 1),
                Parasitic_pct=round(d['energy']['parasitic_pct'], 1), No_flags=not d['flags'])


def sweep(di: DesignInputs, kind: str, values, **kw) -> pd.DataFrame:
    rows = []
    for v in values:
        if kind == 'HRT':
            nd = dataclasses.replace(di, HRT_d=v)
            label = f'HRT {v:g} d'
        elif kind == 'Temperature':
            nd = dataclasses.replace(di, T_op_C=v)
            label = f'T {v:g} C'
        elif kind == 'TS':
            nd = dataclasses.replace(di, TS_pct=v)
            label = f'TS {v:g} %'
        elif kind == 'Blend':
            nd = dataclasses.replace(di, blend=v[1])
            label = v[0]
        else:
            raise ValueError(kind)
        try:
            d = run_design(nd, full_history=False, **kw)
            rows.append(_row(label, d))
        except Exception as e:     # keep sweep alive, show failure
            rows.append(dict(Scenario=label, Error=str(e)))
    return pd.DataFrame(rows)
