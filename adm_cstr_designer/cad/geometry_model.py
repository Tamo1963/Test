"""Pure-geometry description of the digester, used by both the trimesh preview
and the FreeCAD script (single source of dimensions)."""
from __future__ import annotations

import math


def cad_params(design: dict) -> dict:
    g, m, c = design['geometry'], design['mixing'], design['coil']
    D = g['D']
    n = m['n_agitators']
    return dict(
        D=D, H_liq=g['H_liq'], H_free=g['H_free'], h_dome=g['h_dome'],
        wall_t=g['wall_t'], floor_t=g['floor_t'],
        n_agit=n, D_imp=m['D_impeller_m'], shaft_d=m['shaft_diameter_mm'] / 1000.0,
        imp_z=0.30 * g['H_liq'], n_blades=3,
        agit_xy=[(0.0, 0.0)] if n == 1 else [
            (0.5 * D / 2 * math.cos(2 * math.pi * i / n), 0.5 * D / 2 * math.sin(2 * math.pi * i / n))
            for i in range(n)],
        coil_D=c['coil_diameter_m'], coil_OD=c['OD_mm'] / 1000.0, coil_turns=max(c['turns_total'], 1.0),
        coil_z0=0.15 * g['H_liq'], coil_H=c['H_coil_band_m'] / 1.0, n_circuits=c['n_circuits'],
        coil_pitch=c['pitch_m'],
    )
