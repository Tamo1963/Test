"""ADM1 (BSM2 formulation) with algebraic pH, liquid + headspace states.

Units follow BSM2: soluble/particulate organics kg COD/m3, S_IC and S_IN
kmol/m3, S_cat/S_an kmol/m3, gas states G_h2/G_ch4 kg COD/m3 and G_co2 kmol/m3.
"""
from __future__ import annotations

import copy
import numpy as np
from scipy.optimize import brentq

LIQ = ['S_su', 'S_aa', 'S_fa', 'S_va', 'S_bu', 'S_pro', 'S_ac', 'S_h2', 'S_ch4',
       'S_IC', 'S_IN', 'S_I', 'X_c', 'X_ch', 'X_pr', 'X_li', 'X_su', 'X_aa',
       'X_fa', 'X_c4', 'X_pro', 'X_ac', 'X_h2', 'X_I', 'S_cat', 'S_an']
GAS = ['G_h2', 'G_ch4', 'G_co2']
STATES = LIQ + GAS
NL, NG = len(LIQ), len(GAS)
IX = {n: i for i, n in enumerate(STATES)}

R_GAS = 0.083145      # bar m3 / kmol / K
T_BASE = 298.15
P_ATM = 1.013         # bar
NM3_PER_KMOL = 22.414

# carbon (kmol C / kg COD) and nitrogen (kmol N / kg COD) content
N_XC, N_I, N_AA, N_BAC = 0.0376 / 14, 0.06 / 14, 0.007, 0.08 / 14
C = dict(S_su=.0313, S_aa=.03, S_fa=.0217, S_va=.024, S_bu=.025, S_pro=.0268,
         S_ac=.0313, S_h2=0., S_ch4=.0156, S_I=.03, X_c=.02786, X_ch=.0313,
         X_pr=.03, X_li=.022, X_I=.03)
for b in ('X_su', 'X_aa', 'X_fa', 'X_c4', 'X_pro', 'X_ac', 'X_h2'):
    C[b] = .0313
N = dict(S_aa=N_AA, S_I=N_I, X_c=N_XC, X_pr=N_AA, X_I=N_I)
for b in ('X_su', 'X_aa', 'X_fa', 'X_c4', 'X_pro', 'X_ac', 'X_h2'):
    N[b] = N_BAC

PROCS = ['disint', 'hyd_ch', 'hyd_pr', 'hyd_li', 'up_su', 'up_aa', 'up_fa',
         'up_va', 'up_bu', 'up_pro', 'up_ac', 'up_h2', 'dec_su', 'dec_aa',
         'dec_fa', 'dec_c4', 'dec_pro', 'dec_ac', 'dec_h2',
         'T_h2', 'T_ch4', 'T_co2']
BIOMASS = ['X_su', 'X_aa', 'X_fa', 'X_c4', 'X_pro', 'X_ac', 'X_h2']


def petersen(st: dict) -> np.ndarray:
    """Stoichiometric matrix (processes x liquid states). C and N conserved."""
    Y = st['Y']
    M = np.zeros((len(PROCS), NL))
    li = {n: i for i, n in enumerate(LIQ)}

    def row(p, d):
        for k, v in d.items():
            M[PROCS.index(p), li[k]] += v

    row('disint', {'X_c': -1, 'S_I': st['f_sI_xc'], 'X_ch': st['f_ch_xc'],
                   'X_pr': st['f_pr_xc'], 'X_li': st['f_li_xc'], 'X_I': st['f_xI_xc']})
    row('hyd_ch', {'X_ch': -1, 'S_su': 1})
    row('hyd_pr', {'X_pr': -1, 'S_aa': 1})
    row('hyd_li', {'X_li': -1, 'S_su': 1 - st['f_fa_li'], 'S_fa': st['f_fa_li']})
    f = st
    h2_su = 1 - f['f_bu_su'] - f['f_pro_su'] - f['f_ac_su']
    row('up_su', {'S_su': -1, 'S_bu': (1 - Y['su']) * f['f_bu_su'],
                  'S_pro': (1 - Y['su']) * f['f_pro_su'], 'S_ac': (1 - Y['su']) * f['f_ac_su'],
                  'S_h2': (1 - Y['su']) * h2_su, 'X_su': Y['su']})
    h2_aa = 1 - f['f_va_aa'] - f['f_bu_aa'] - f['f_pro_aa'] - f['f_ac_aa']
    row('up_aa', {'S_aa': -1, 'S_va': (1 - Y['aa']) * f['f_va_aa'],
                  'S_bu': (1 - Y['aa']) * f['f_bu_aa'], 'S_pro': (1 - Y['aa']) * f['f_pro_aa'],
                  'S_ac': (1 - Y['aa']) * f['f_ac_aa'], 'S_h2': (1 - Y['aa']) * h2_aa,
                  'X_aa': Y['aa']})
    row('up_fa', {'S_fa': -1, 'S_ac': (1 - Y['fa']) * f['f_ac_fa'],
                  'S_h2': (1 - Y['fa']) * (1 - f['f_ac_fa']), 'X_fa': Y['fa']})
    row('up_va', {'S_va': -1, 'S_pro': (1 - Y['c4']) * f['f_pro_va'],
                  'S_ac': (1 - Y['c4']) * f['f_ac_va'],
                  'S_h2': (1 - Y['c4']) * (1 - f['f_pro_va'] - f['f_ac_va']), 'X_c4': Y['c4']})
    row('up_bu', {'S_bu': -1, 'S_ac': (1 - Y['c4']) * f['f_ac_bu'],
                  'S_h2': (1 - Y['c4']) * (1 - f['f_ac_bu']), 'X_c4': Y['c4']})
    row('up_pro', {'S_pro': -1, 'S_ac': (1 - Y['pro']) * f['f_ac_pro'],
                   'S_h2': (1 - Y['pro']) * (1 - f['f_ac_pro']), 'X_pro': Y['pro']})
    row('up_ac', {'S_ac': -1, 'S_ch4': 1 - Y['ac'], 'X_ac': Y['ac']})
    row('up_h2', {'S_h2': -1, 'S_ch4': 1 - Y['h2'], 'X_h2': Y['h2']})
    for k, x in zip(('su', 'aa', 'fa', 'c4', 'pro', 'ac', 'h2'), BIOMASS):
        row('dec_' + k, {x: -1, 'X_c': 1})

    # carbon / nitrogen balance closure via S_IC and S_IN (biochemical rows only)
    iC, iN = li['S_IC'], li['S_IN']
    for p in range(19):
        sc = sn = 0.0
        for j, name in enumerate(LIQ):
            if name in C and M[p, j]:
                sc += M[p, j] * C[name]
            if name in N and M[p, j]:
                sn += M[p, j] * N[name]
        M[p, iC] = -sc
        M[p, iN] = -sn
    # gas transfer: liquid loses what goes to the headspace
    M[PROCS.index('T_h2'), li['S_h2']] = -1
    M[PROCS.index('T_ch4'), li['S_ch4']] = -1
    M[PROCS.index('T_co2'), iC] = -1
    return M


def check_balances(M: np.ndarray, tol=1e-9) -> dict:
    """C and N residual per biochemical process (should be ~0)."""
    cvec = np.array([C.get(n, 0.) if n != 'S_IC' else 1. for n in LIQ])
    nvec = np.array([N.get(n, 0.) if n != 'S_IN' else 1. for n in LIQ])
    return {'C': float(np.abs(M[:19] @ cvec).max()), 'N': float(np.abs(M[:19] @ nvec).max())}


class ADM1Params:
    """Holds a kinetic set (dict from kinetics_db.json) plus feed-specific hydrolysis."""

    def __init__(self, kset: dict, T_C: float, hyd: dict | None = None, theta: float | None = None):
        self.k = copy.deepcopy(kset)
        self.T = T_C + 273.15
        theta = self.k.get('theta', 1.05) if theta is None else theta
        self.fT = theta ** (T_C - self.k['T_ref_C'])
        self.uptake = self.k['uptake']
        self.hyd = dict(self.k['hydrolysis'])
        if hyd:
            self.hyd.update(hyd)
        self.Y = {s: v['Y'] for s, v in self.uptake.items()}
        self.st = dict(self.k['stoich'], Y=self.Y)
        self.M = petersen(self.st)
        T = self.T
        a = 1.0 / T_BASE - 1.0 / T
        self.Kw = 1e-14 * np.exp(55900 / (100 * R_GAS) * a)
        self.Ka_co2 = 10 ** -6.35 * np.exp(7646 / (100 * R_GAS) * a)
        self.Ka_IN = 10 ** -9.25 * np.exp(51965 / (100 * R_GAS) * a)
        self.Ka_ac, self.Ka_pro = 10 ** -4.76, 10 ** -4.88
        self.Ka_bu, self.Ka_va = 10 ** -4.82, 10 ** -4.86
        self.KH_co2 = 0.035 * np.exp(-19410 / (100 * R_GAS) * a)
        self.KH_ch4 = 0.0014 * np.exp(-14240 / (100 * R_GAS) * a)
        self.KH_h2 = 7.8e-4 * np.exp(-4180 / (100 * R_GAS) * a)
        self.p_h2o = 0.0313 * np.exp(5290 * (1 / 298.15 - 1 / T))
        self.kLa = 200.0


def solve_pH(x: np.ndarray, p: ADM1Params) -> dict:
    """Charge-balance pH. Returns h, pH and speciation."""
    S_ac, S_pro, S_bu, S_va = x[IX['S_ac']], x[IX['S_pro']], x[IX['S_bu']], x[IX['S_va']]
    S_IC, S_IN, S_cat, S_an = x[IX['S_IC']], x[IX['S_IN']], x[IX['S_cat']], x[IX['S_an']]

    def charge(pH):
        h = 10.0 ** (-pH)
        nh3 = S_IN * p.Ka_IN / (p.Ka_IN + h)
        hco3 = S_IC * p.Ka_co2 / (p.Ka_co2 + h)
        vfa = (S_ac / 64 * p.Ka_ac / (p.Ka_ac + h) + S_pro / 112 * p.Ka_pro / (p.Ka_pro + h)
               + S_bu / 160 * p.Ka_bu / (p.Ka_bu + h) + S_va / 208 * p.Ka_va / (p.Ka_va + h))
        return S_cat + (S_IN - nh3) + h - hco3 - vfa - p.Kw / h - S_an

    try:
        pH = brentq(charge, 2.0, 12.0, xtol=1e-9)
    except ValueError:
        pH = 2.0 if charge(2.0) < 0 else 12.0
    h = 10.0 ** (-pH)
    nh3 = S_IN * p.Ka_IN / (p.Ka_IN + h)
    hco3 = S_IC * p.Ka_co2 / (p.Ka_co2 + h)
    return dict(pH=pH, h=h, nh3=nh3, hco3=hco3, co2=S_IC - hco3)


def rates(x: np.ndarray, p: ADM1Params, pH_override: float | None = None) -> tuple[np.ndarray, dict]:
    """Process rates rho (22) and diagnostics. x = full state vector."""
    ph = solve_pH(x, p)
    if pH_override is not None:
        ph = dict(ph, pH=pH_override, h=10 ** -pH_override)
    h = ph['h']
    g = lambda n: max(x[IX[n]], 0.0)
    u, fT = p.uptake, p.fT

    def pHinh(ul_ll):
        ll, ul = ul_ll
        pH = ph['pH']
        return 1.0 if pH >= ul else np.exp(-3.0 * ((pH - ul) / (ul - ll)) ** 2)

    pl = p.k['pH_limits']
    I_aa, I_ac, I_h2p = pHinh(pl['aa']), pHinh(pl['ac']), pHinh(pl['h2'])
    S_IN = g('S_IN')
    I_IN = S_IN / (S_IN + p.k['K_S_IN'])
    I_nh3 = 1.0 / (1.0 + ph['nh3'] / u['ac']['KI_nh3'])
    inhH = lambda s: 1.0 / (1.0 + g('S_h2') / u[s]['KI_h2'])

    def mon(s, S):
        return u[s]['k_m'] * fT * S / (u[s]['K_S'] + S)

    rho = np.zeros(len(PROCS))
    rho[0] = p.hyd['k_dis'] * fT * g('X_c')
    rho[1] = p.hyd['k_hyd_ch'] * fT * g('X_ch')
    rho[2] = p.hyd['k_hyd_pr'] * fT * g('X_pr')
    rho[3] = p.hyd['k_hyd_li'] * fT * g('X_li')
    rho[4] = mon('su', g('S_su')) * g('X_su') * I_aa * I_IN
    rho[5] = mon('aa', g('S_aa')) * g('X_aa') * I_aa * I_IN
    rho[6] = mon('fa', g('S_fa')) * g('X_fa') * I_aa * I_IN * inhH('fa')
    sva, sbu = g('S_va'), g('S_bu')
    rho[7] = mon('c4', sva) * g('X_c4') * (sva / (sva + sbu + 1e-12)) * I_aa * I_IN * inhH('c4')
    rho[8] = mon('c4', sbu) * g('X_c4') * (sbu / (sva + sbu + 1e-12)) * I_aa * I_IN * inhH('c4')
    rho[9] = mon('pro', g('S_pro')) * g('X_pro') * I_aa * I_IN * inhH('pro')
    rho[10] = mon('ac', g('S_ac')) * g('X_ac') * I_ac * I_IN * I_nh3
    rho[11] = mon('h2', g('S_h2')) * g('X_h2') * I_h2p * I_IN
    for j, (s, xn) in enumerate(zip(('su', 'aa', 'fa', 'c4', 'pro', 'ac', 'h2'), BIOMASS)):
        rho[12 + j] = u[s]['b'] * fT * g(xn)
    # gas transfer
    pg_h2 = x[IX['G_h2']] * R_GAS * p.T / 16
    pg_ch4 = x[IX['G_ch4']] * R_GAS * p.T / 64
    pg_co2 = x[IX['G_co2']] * R_GAS * p.T
    rho[19] = p.kLa * (x[IX['S_h2']] - 16 * p.KH_h2 * pg_h2)
    rho[20] = p.kLa * (x[IX['S_ch4']] - 64 * p.KH_ch4 * pg_ch4)
    rho[21] = p.kLa * (ph['co2'] - p.KH_co2 * pg_co2)
    diag = dict(ph, I_aa=I_aa, I_ac=I_ac, I_h2=I_h2p, I_IN=I_IN, I_nh3=I_nh3,
                p_gas=pg_h2 + pg_ch4 + pg_co2 + p.p_h2o)
    return rho, diag
