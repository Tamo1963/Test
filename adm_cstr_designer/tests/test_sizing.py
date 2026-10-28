import math
import pytest
from core.data import load_feedstocks
from core.feed import FeedInputs, blend_properties, normalise_blend
from core.mixing import size_agitator
from core.thermal import heat_load, size_coil, lmtd
from core.design import DesignInputs, run_design

LIB = load_feedstocks()


def test_blend_normalises_and_rejects_empty():
    assert sum(normalise_blend({'a': 30, 'b': 10}).values()) == pytest.approx(1.0)
    with pytest.raises(ValueError):
        normalise_blend({'a': 0})


def test_feed_cod_and_vs_mass_balance():
    fp = blend_properties(FeedInputs(100, 25, 8, 85, {'Napier grass': 50, 'Cow dung': 50}), LIB)
    from core.adm1 import IX
    f = fp['feed']
    cod = sum(f[IX[n]] for n in ('X_ch', 'X_pr', 'X_li', 'X_I', 'S_ac'))
    assert cod == pytest.approx(fp['COD_kg_m3'], rel=1e-9)
    assert fp['VS_load_kg_d'] == pytest.approx(100 * 1000 * 0.08 * 0.85, rel=0.01)


def test_agitator_power_identity():
    m = size_agitator(1000, 12.0, 12.0, 8.0, 1050.0)
    N = m['N_rpm'] / 60
    P = m['Np'] * 1050.0 * N ** 3 * m['D_impeller_m'] ** 5
    assert P / 1000 == pytest.approx(m['P_hydraulic_kW'], rel=1e-9)
    assert m['motor_kW'] >= m['P_shaft_kW'] * 1.15


def test_coil_area_matches_duty():
    c = size_coil(100.0, 37.0, 8.0, 13.0, 13.0)
    assert c['area_m2'] * c['U_W_m2K'] * c['LMTD_K'] == pytest.approx(100e3, rel=1e-9)
    assert lmtd(30, 18) == pytest.approx(12 / math.log(30 / 18))


def test_heat_load_scales_with_delta_T():
    a = heat_load(100, 1050, 4.0, 25, 37, 13, 13, 0.5, 15)
    b = heat_load(100, 1050, 4.0, 15, 37, 13, 13, 0.5, 15)
    assert b['Q_feed_kW'] > a['Q_feed_kW']


def test_design_end_to_end_default():
    d = run_design(DesignInputs(), full_history=False)
    p = d['performance']
    assert p['converged'] and 0.99 < p['cod_closure'] < 1.01
    assert p['SMY_Nm3_kgVS'] <= p['BMP_blend'] + 1e-9      # cannot exceed the blend BMP
    assert 6.5 < p['pH'] < 8.0
