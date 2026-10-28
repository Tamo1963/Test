"""Agitator sizing for a top-entry axial-flow (hydrofoil) mixer in a digester."""
from __future__ import annotations

import math

IEC_KW = [0.37, 0.55, 0.75, 1.1, 1.5, 2.2, 3.0, 4.0, 5.5, 7.5, 11, 15, 18.5, 22, 30,
          37, 45, 55, 75, 90, 110, 132, 160, 200]


def apparent_viscosity(TS_pct: float) -> float:
    """Indicative apparent viscosity of digestate, Pa.s (user may override)."""
    return 0.001 * math.exp(0.45 * TS_pct)


def size_agitator(V_liq, T_diam, H_liq, TS_pct, rho, mu=None, eps_W_m3=6.0, DT=0.35,
                  Np=0.30, Nq=0.55, n_agit=None, eta_mech=0.90, service=1.15,
                  tau_allow_MPa=40.0, peak_torque_factor=2.5):
    """Size by specific power input; check Re, tip speed, G and turnover."""
    mu = apparent_viscosity(TS_pct) if mu is None else mu
    n = n_agit or max(1, math.ceil(T_diam / 12.0))
    V_each = V_liq / n
    P_hyd = eps_W_m3 * V_each                                  # W hydraulic per agitator
    # one agitator serves a tank sector: take impeller D from the tank diameter
    D = DT * (T_diam if n == 1 else T_diam / math.sqrt(n))
    N = (P_hyd / (Np * rho * D ** 5)) ** (1.0 / 3.0)             # rev/s
    Re = rho * N * D ** 2 / mu
    tip = math.pi * N * D
    Qp = Nq * N * D ** 3                                         # m3/s pumping
    turnover_min = V_each / Qp / 60.0
    G = math.sqrt(P_hyd / (mu * V_each))                         # 1/s
    P_shaft = P_hyd / eta_mech
    P_motor_req = P_shaft * service / 1000.0                     # kW
    motor = next((m for m in IEC_KW if m >= P_motor_req), math.ceil(P_motor_req))
    torque = P_shaft / (2 * math.pi * N)
    d_shaft = (16 * torque * peak_torque_factor / (math.pi * tau_allow_MPa * 1e6)) ** (1.0 / 3.0)
    flags = []
    if Re < 1e4:
        flags.append('Re < 10,000: transitional flow; hydrofoil Np rises, check with vendor CFD.')
    if tip > 4.0:
        flags.append('Tip speed > 4 m/s: reduce speed / use more agitators.')
    if tip < 0.3:
        flags.append('Tip speed < 0.3 m/s: low; confirm solids suspension.')
    if turnover_min > 60:
        flags.append('Turnover time > 60 min: mixing may be poor; raise power or add agitators.')
    if not 5 <= G <= 100:
        flags.append(f'G = {G:.0f} 1/s outside the 5-100 1/s digester guide range.')
    return dict(n_agitators=n, mu_Pa_s=mu, eps_W_m3=eps_W_m3, D_impeller_m=D,
                D_over_T=D / T_diam, N_rpm=N * 60, Np=Np, Nq=Nq, Re=Re,
                tip_speed_m_s=tip, pumping_m3_h=Qp * 3600, turnover_min=turnover_min,
                G_1_s=G, P_hydraulic_kW=P_hyd / 1000, P_shaft_kW=P_shaft / 1000,
                motor_kW=motor, torque_Nm=torque, shaft_diameter_mm=d_shaft * 1000,
                specific_power_total_kW=n * P_hyd / 1000, flags=flags)
