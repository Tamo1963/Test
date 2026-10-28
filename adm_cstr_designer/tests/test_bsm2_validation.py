"""Validate the ADM1 core against the BSM2 benchmark steady state.

Reference (MATLAB BSM2, as quoted in QSDsan tests/test_junctions.py):
Q = 178.4674 m3/d, V_liq = 3400, V_gas = 300, T = 35 C.
"""
import json, pathlib
import numpy as np
import pytest
from core.adm1 import ADM1Params, IX, LIQ, check_balances
from core.cstr import CSTRSpec, steady_state

DATA = pathlib.Path(__file__).resolve().parents[1] / 'data'
KSET = json.loads((DATA / 'kinetics_db.json').read_text())['sets']['ADM1_BSM2_mesophilic']

PRE = dict(S_aa=0.04388, S_IC=0.0079326, S_IN=0.0019721, S_I=0.028067, X_ch=3.7236,
           X_pr=15.9235, X_li=8.047, X_I=17.0106, S_an=0.0052101)
POST = dict(S_su=0.012394, S_aa=0.0055432, S_fa=0.10741, S_va=0.012333, S_bu=0.014003,
            S_pro=0.017584, S_ac=0.089315, S_h2=2.5055e-07, S_ch4=0.05549,
            S_IC=0.095149, S_IN=1.3226 / 14, S_I=0.13087, X_c=0.10792, X_ch=0.020517,
            X_pr=0.08422, X_li=0.043629, X_su=0.31222, X_aa=0.93167, X_fa=0.33839,
            X_c4=0.33577, X_pro=0.10112, X_ac=0.67724, X_h2=0.28484, X_I=17.2162)


def test_stoichiometry_closes_C_and_N():
    p = ADM1Params(KSET, 35.0)
    bal = check_balances(p.M)
    assert bal['C'] < 1e-9 and bal['N'] < 1e-9


@pytest.fixture(scope='module')
def bsm2():
    p = ADM1Params(KSET, 35.0)
    feed = np.zeros(len(LIQ))
    for k, v in PRE.items():
        feed[IX[k]] = v
    spec = CSTRSpec(V_liq=3400, V_gas=300, Q=178.4674, T_C=35.0)
    res = steady_state(spec, p, feed)
    return res, p


def test_converged(bsm2):
    assert bsm2[0]['converged']


@pytest.mark.parametrize('name,rtol', [
    ('S_ac', 0.10), ('S_ch4', 0.10), ('S_IC', 0.05), ('S_IN', 0.05), ('X_ac', 0.10),
    ('X_h2', 0.15), ('X_I', 0.02), ('S_fa', 0.15), ('S_pro', 0.20)])
def test_state_vs_matlab(bsm2, name, rtol):
    res, _ = bsm2
    assert res['y'][IX[name]] == pytest.approx(POST[name], rel=rtol)


def test_pH_and_gas(bsm2):
    res, _ = bsm2
    assert 7.1 < res['diag']['pH'] < 7.5
    g = res['gas']
    assert 0.55 < g['ch4_frac'] < 0.70
    assert g['total'] > 0


def test_cod_balance_closes(bsm2):
    from core.cstr import cod_balance
    res, _ = bsm2
    spec = CSTRSpec(V_liq=3400, V_gas=300, Q=178.4674, T_C=35.0)
    feed = np.zeros(len(LIQ))
    for k, v in PRE.items():
        feed[IX[k]] = v
    assert cod_balance(res, spec, feed)['closure'] == pytest.approx(1.0, abs=0.01)
