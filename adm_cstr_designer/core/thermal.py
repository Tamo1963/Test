"""Heat load and heating-coil sizing."""
from __future__ import annotations

import math


def heat_load(Q_m3d, rho, cp, T_feed, T_op, D, H_liq, V_gas_dome_h, T_amb,
              t_conc=0.30, k_conc=1.7, t_ins=0.10, k_ins=0.04, U_roof=0.60,
              U_floor=0.50, T_ground=10.0, margin=0.15):
    m_dot = Q_m3d * rho / 86400.0                                # kg/s
    Q_feed = m_dot * cp * (T_op - T_feed)                        # kW (cp in kJ/kg/K)
    R_wall = t_conc / k_conc + t_ins / k_ins + 0.13 + 0.04
    U_wall = 1.0 / R_wall
    A_wall = math.pi * D * (H_liq + V_gas_dome_h)
    A_roof = math.pi * D ** 2 / 4 * 1.15                         # dome surface ~ 1.15 x plan
    A_floor = math.pi * D ** 2 / 4
    Q_wall = U_wall * A_wall * (T_op - T_amb) / 1000
    Q_roof = U_roof * A_roof * (T_op - T_amb) / 1000
    Q_floor = U_floor * A_floor * (T_op - T_ground) / 1000
    loss = Q_wall + Q_roof + Q_floor
    design = (Q_feed + loss) * (1 + margin)
    return dict(m_dot_kg_s=m_dot, Q_feed_kW=Q_feed, Q_wall_kW=Q_wall, Q_roof_kW=Q_roof,
                Q_floor_kW=Q_floor, Q_loss_kW=loss, margin=margin, Q_design_kW=design,
                U_wall=U_wall, A_wall=A_wall, A_roof=A_roof, A_floor=A_floor)


def lmtd(dT1, dT2):
    if abs(dT1 - dT2) < 1e-6:
        return dT1
    return (dT1 - dT2) / math.log(dT1 / dT2)


def size_coil(Q_kW, T_op, TS_pct, D_tank, H_liq, T_w_in=70.0, T_w_out=55.0,
              OD_m=0.0603, ID_m=0.0525, U=None, v_max=1.5, L_circuit_max=150.0,
              wall_clearance=0.5):
    """Wall-mounted helical coil (internal), hot water inside, digestate outside."""
    if T_w_out <= T_op:
        raise ValueError('Hot-water return temperature must exceed the digester temperature.')
    U = U if U else 200.0 / (1 + 0.08 * TS_pct)                    # W/m2K indicative
    dT = lmtd(T_w_in - T_op, T_w_out - T_op)
    A = Q_kW * 1000 / (U * dT)
    L = A / (math.pi * OD_m)
    cp_w = 4.18
    m_w = Q_kW / (cp_w * (T_w_in - T_w_out))                      # kg/s
    q_w = m_w / 980.0                                            # m3/s
    a_pipe = math.pi * ID_m ** 2 / 4
    n_par = max(math.ceil(q_w / (v_max * a_pipe)), math.ceil(L / L_circuit_max), 1)
    v = q_w / n_par / a_pipe
    L_c = L / n_par
    D_coil = D_tank - 2 * wall_clearance
    turns_total = L / (math.pi * D_coil)
    # one circuit per vertical bank; stack total turns over the liquid height
    H_avail = 0.75 * H_liq
    pitch = H_avail / max(turns_total, 1e-9)
    flags = []
    if pitch < 2.5 * OD_m:
        flags.append('Coil pitch < 2.5 x OD; add a second concentric coil, increase water temperature or tube size.')
    if v < 0.4:
        flags.append('Water velocity < 0.4 m/s; fouling / laminar risk.')
    return dict(U_W_m2K=U, LMTD_K=dT, area_m2=A, length_m=L, n_circuits=n_par,
                length_per_circuit_m=L_c, OD_mm=OD_m * 1000, ID_mm=ID_m * 1000,
                water_flow_kg_s=m_w, water_flow_m3_h=q_w * 3600, water_velocity_m_s=v,
                coil_diameter_m=D_coil, turns_total=turns_total, pitch_m=pitch,
                H_coil_band_m=H_avail, T_w_in=T_w_in, T_w_out=T_w_out, flags=flags)
