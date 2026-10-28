# -*- coding: utf-8 -*-
# FreeCAD build script for the anaerobic CSTR designed by adm_cstr_designer.
# Run with:  freecadcmd this_file.py      (or paste into FreeCAD's Python console / an MCP execute-code tool)
# PARAMS and OUT_DIR are defined at the top of the generated file.
import math
import os
import FreeCAD as App
import Part

V = App.Vector
P = PARAMS  # noqa: F821  (injected by the generator)
R = P['D'] / 2.0
zt = P['H_liq'] + P['H_free']          # top of cylindrical shell
doc = App.newDocument('AnaerobicCSTR')
objs = []


def add(name, shape, color=None):
    o = doc.addObject('Part::Feature', name)
    o.Shape = shape
    if color is not None and hasattr(o, 'ViewObject') and o.ViewObject is not None:
        o.ViewObject.ShapeColor = color
    objs.append(o)
    return o


# --- tank: wall, floor slab, domed roof ---------------------------------------
outer = Part.makeCylinder(R, zt, V(0, 0, 0), V(0, 0, 1))
inner = Part.makeCylinder(R - P['wall_t'], zt + 0.01, V(0, 0, 0), V(0, 0, 1))
add('TankWall', outer.cut(inner), (0.62, 0.65, 0.69))
add('FloorSlab', Part.makeCylinder(R, P['floor_t'], V(0, 0, -P['floor_t']), V(0, 0, 1)), (0.42, 0.46, 0.52))

hd = P['h_dome']
Rs = (R ** 2 + hd ** 2) / (2.0 * hd)             # sphere radius of the cap
zc = zt + hd - Rs                                 # sphere centre height
cap = Part.makeSphere(Rs, V(0, 0, zc)).common(Part.makeCylinder(R, hd + 0.01, V(0, 0, zt), V(0, 0, 1)))
roof_t = min(0.20, P['wall_t'])
cap = cap.cut(Part.makeSphere(Rs - roof_t, V(0, 0, zc)))
add('RoofDome', cap, (0.80, 0.82, 0.85))

# --- nozzles: feed (low), digestate overflow (at liquid level), gas outlet (roof) -------------
nz_r = 0.15
add('FeedNozzle', Part.makeCylinder(nz_r, 1.2, V(R - P['wall_t'] - 0.1, 0, 1.0), V(1, 0, 0)))
add('OutletNozzle', Part.makeCylinder(nz_r, 1.2, V(-(R + 1.0), 0, P['H_liq'] - 0.3), V(1, 0, 0)))
add('GasOutlet', Part.makeCylinder(0.20, 1.0, V(0, 0, zt + hd - 0.05), V(0, 0, 1)))

# --- heating coil: helical pipe(s) --------------------------------------------
nc = max(int(P['n_circuits']), 1)
turns_c = P['coil_turns'] / nc
rc = P['coil_D'] / 2.0
pitch = P['coil_H'] / max(turns_c, 1e-6)
for c in range(nc):
    helix = Part.makeHelix(pitch, P['coil_H'], rc)
    # profile circle at helix start (rc,0,0), normal along the helix tangent
    tang = V(0, 2 * math.pi * rc, pitch)
    tang.normalize()
    prof = Part.Wire(Part.makeCircle(P['coil_OD'] / 2.0, V(rc, 0, 0), tang))
    pipe = Part.Wire(helix).makePipeShell([prof], True, True)
    pipe.rotate(V(0, 0, 0), V(0, 0, 1), 360.0 * c / nc)
    pipe.translate(V(0, 0, P['coil_z0']))
    add('HeatingCoil_%d' % (c + 1), pipe, (0.85, 0.33, 0.31))

# --- agitators: shaft, hub, pitched blades ------------------------------------
Lb, Wb, Tb = P['D_imp'] / 2.0, 0.16 * P['D_imp'], 0.03 * P['D_imp']
for i, (x, y) in enumerate(P['agit_xy']):
    shaft = Part.makeCylinder(P['shaft_d'] / 2.0, zt + 0.8 - P['imp_z'], V(x, y, P['imp_z']), V(0, 0, 1))
    add('AgitatorShaft_%d' % (i + 1), shaft, (0.17, 0.18, 0.21))
    add('Hub_%d' % (i + 1), Part.makeCylinder(0.12 * P['D_imp'] / 2.0 + 0.05, 0.25, V(x, y, P['imp_z'] - 0.12), V(0, 0, 1)))
    for k in range(int(P['n_blades'])):
        blade = Part.makeBox(Lb, Wb, Tb, V(0, -Wb / 2.0, -Tb / 2.0), V(1, 0, 0))
        blade.rotate(V(0, 0, 0), V(1, 0, 0), 20.0)                       # blade pitch
        blade.rotate(V(0, 0, 0), V(0, 0, 1), 360.0 * k / P['n_blades'])
        blade.translate(V(x, y, P['imp_z']))
        add('Blade_%d_%d' % (i + 1, k + 1), blade, (0.12, 0.47, 0.71))

doc.recompute()

# --- export -----------------------------------------------------------------
os.makedirs(OUT_DIR, exist_ok=True)  # noqa: F821
fc = os.path.join(OUT_DIR, 'anaerobic_cstr.FCStd')  # noqa: F821
doc.saveAs(fc)
try:
    import Import
    Import.export(objs, os.path.join(OUT_DIR, 'anaerobic_cstr.step'))  # noqa: F821
except Exception as e:  # STEP export needs the Import module
    print('STEP export skipped:', e)
try:
    import Mesh
    Mesh.export(objs, os.path.join(OUT_DIR, 'anaerobic_cstr.stl'))  # noqa: F821
except Exception as e:
    print('STL export skipped:', e)
print('FreeCAD model written to', OUT_DIR)  # noqa: F821
