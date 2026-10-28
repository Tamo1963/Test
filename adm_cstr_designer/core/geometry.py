"""Tank geometry from working volume."""
from __future__ import annotations

import math


def tank_geometry(V_liq, H_over_D=1.0, gas_frac=0.10, n_tanks=1, wall_t=0.30, floor_t=0.40):
    V_each = V_liq / n_tanks
    D = (4 * V_each / (math.pi * H_over_D)) ** (1 / 3)
    H = H_over_D * D
    V_gas = gas_frac * V_each
    # spherical-cap dome of base diameter D: V = (pi h/6)(3 (D/2)^2 + h^2)
    h = 0.0
    for _ in range(60):
        f = math.pi * h / 6 * (3 * (D / 2) ** 2 + h ** 2) - V_gas
        df = math.pi / 6 * (3 * (D / 2) ** 2 + 3 * h ** 2)
        h -= f / df
    cyl_free = 0.5
    return dict(n_tanks=n_tanks, V_liq_each=V_each, D=D, H_liq=H, H_free=cyl_free, h_dome=h,
                V_gas_each=V_gas + math.pi * D ** 2 / 4 * cyl_free,
                H_total=H + cyl_free + h, wall_t=wall_t, floor_t=floor_t)
