"""Lightweight 3D model (numpy meshes) for the in-app preview and STL export.

This is NOT the FreeCAD model: it is an independent visualisation built from the
same dimensions (cad.geometry_model.cad_params) so the app works without FreeCAD.
"""
from __future__ import annotations

import io
import math
import numpy as np

PART_COLORS = {'Tank wall': '#9aa5b1', 'Floor': '#6b7785', 'Roof dome': '#c9d1da',
               'Heating coil': '#d9534f', 'Agitator shaft': '#2b2f36', 'Impeller': '#1f78b4',
               'Liquid level': '#4fb3d9'}


def _tube(path: np.ndarray, r: float, n: int = 10):
    """Swept circle along a polyline -> (verts, faces)."""
    path = np.asarray(path, float)
    T = np.gradient(path, axis=0)
    T /= np.linalg.norm(T, axis=1, keepdims=True) + 1e-12
    ref = np.array([0, 0, 1.0])
    verts, faces = [], []
    th = np.linspace(0, 2 * np.pi, n, endpoint=False)
    for i, (p, t) in enumerate(zip(path, T)):
        a = np.cross(t, ref)
        if np.linalg.norm(a) < 1e-6:
            a = np.cross(t, [1, 0, 0])
        a /= np.linalg.norm(a)
        b = np.cross(t, a)
        verts.extend(p + r * (np.cos(k) * a + np.sin(k) * b) for k in th)
    for i in range(len(path) - 1):
        for k in range(n):
            a0, a1 = i * n + k, i * n + (k + 1) % n
            b0, b1 = (i + 1) * n + k, (i + 1) * n + (k + 1) % n
            faces += [(a0, b0, a1), (a1, b0, b1)]
    return np.array(verts), np.array(faces)


def _revolve(profile, n=48):
    """Revolve (r,z) profile about z."""
    profile = np.asarray(profile, float)
    th = np.linspace(0, 2 * np.pi, n, endpoint=False)
    verts = np.array([[r * math.cos(t), r * math.sin(t), z] for r, z in profile for t in th])
    faces = []
    m = len(profile)
    for i in range(m - 1):
        for k in range(n):
            a0, a1 = i * n + k, i * n + (k + 1) % n
            b0, b1 = (i + 1) * n + k, (i + 1) * n + (k + 1) % n
            faces += [(a0, b0, a1), (a1, b0, b1)]
    return verts, np.array(faces)


def _box(c, size, rot=None):
    sx, sy, sz = [s / 2 for s in size]
    v = np.array([[x, y, z] for x in (-sx, sx) for y in (-sy, sy) for z in (-sz, sz)], float)
    if rot is not None:
        v = v @ rot.T
    f = np.array([[0, 1, 3], [0, 3, 2], [4, 6, 7], [4, 7, 5], [0, 4, 5], [0, 5, 1],
                  [2, 3, 7], [2, 7, 6], [0, 2, 6], [0, 6, 4], [1, 5, 7], [1, 7, 3]])
    return v + np.asarray(c), f


def build_parts(cp: dict) -> dict:
    D, R = cp['D'], cp['D'] / 2
    H, Hf, hd = cp['H_liq'], cp['H_free'], cp['h_dome']
    t, tf = cp['wall_t'], cp['floor_t']
    zt = H + Hf
    parts = {}
    # wall: outer & inner cylinder profile revolve (closed ring section)
    prof = [(R, 0), (R, zt), (R - t, zt), (R - t, 0), (R, 0)]
    parts['Tank wall'] = _revolve(prof)
    parts['Floor'] = _revolve([(0, -tf), (R, -tf), (R, 0), (0, 0)])
    # spherical cap roof
    Rs = (R ** 2 + hd ** 2) / (2 * hd)
    zc = zt + hd - Rs
    ang = np.linspace(0, math.asin(R / Rs), 14)
    parts['Roof dome'] = _revolve([(Rs * math.sin(a), zc + Rs * math.cos(a)) for a in ang])
    # coil
    nc = max(int(cp['n_circuits']), 1)
    turns = cp['coil_turns'] / nc
    cv, cf, off = [], [], 0
    for c in range(nc):
        s = np.linspace(0, 1, max(int(turns * 40), 30))
        ph = 2 * math.pi * c / nc
        a = 2 * math.pi * turns * s + ph
        path = np.c_[cp['coil_D'] / 2 * np.cos(a), cp['coil_D'] / 2 * np.sin(a),
                     cp['coil_z0'] + cp['coil_H'] * s]
        v, f = _tube(path, cp['coil_OD'] / 2, 8)
        cv.append(v); cf.append(f + off); off += len(v)
    parts['Heating coil'] = (np.vstack(cv), np.vstack(cf))
    # agitators
    sv, sf, iv, if_, so, io_ = [], [], [], [], 0, 0
    for (x, y) in cp['agit_xy']:
        v, f = _revolve([(0, zt + 0.8), (cp['shaft_d'] / 2, zt + 0.8), (cp['shaft_d'] / 2, cp['imp_z']),
                         (0, cp['imp_z'])], 16)
        sv.append(v + [x, y, 0]); sf.append(f + so); so += len(v)
        for k in range(cp['n_blades']):
            th = 2 * math.pi * k / cp['n_blades']
            rot = np.array([[math.cos(th), -math.sin(th), 0], [math.sin(th), math.cos(th), 0], [0, 0, 1]])
            pitch = np.radians(20)
            tilt = np.array([[1, 0, 0], [0, math.cos(pitch), -math.sin(pitch)], [0, math.sin(pitch), math.cos(pitch)]])
            L, W, Th = cp['D_imp'] / 2, 0.16 * cp['D_imp'], 0.03 * cp['D_imp']
            v, f = _box([L / 2, 0, 0], [L, W, Th], tilt)
            v = v @ rot.T + [x, y, cp['imp_z']]
            iv.append(v); if_.append(f + io_); io_ += len(v)
    parts['Agitator shaft'] = (np.vstack(sv), np.vstack(sf))
    parts['Impeller'] = (np.vstack(iv), np.vstack(if_))
    return parts


def to_stl_bytes(parts: dict) -> bytes:
    import struct
    tris = []
    for v, f in parts.values():
        tris.extend(v[f])
    tris = np.asarray(tris, dtype=np.float32)
    n = np.cross(tris[:, 1] - tris[:, 0], tris[:, 2] - tris[:, 0])
    n /= (np.linalg.norm(n, axis=1, keepdims=True) + 1e-12)
    buf = io.BytesIO()
    buf.write(b'adm_cstr_designer preview'.ljust(80, b' '))
    buf.write(struct.pack('<I', len(tris)))
    for nn, tr in zip(n.astype(np.float32), tris):
        buf.write(nn.tobytes() + tr.tobytes() + b'\x00\x00')
    return buf.getvalue()


def plotly_figure(parts: dict, H_liq: float | None = None, D: float | None = None, show_liquid=False):
    import plotly.graph_objects as go
    fig = go.Figure()
    for name, (v, f) in parts.items():
        fig.add_trace(go.Mesh3d(x=v[:, 0], y=v[:, 1], z=v[:, 2], i=f[:, 0], j=f[:, 1], k=f[:, 2],
                                name=name, color=PART_COLORS.get(name, '#888'),
                                opacity=0.35 if name in ('Tank wall', 'Roof dome', 'Floor') else 1.0,
                                flatshading=True, showlegend=True))
    fig.update_layout(scene=dict(aspectmode='data', xaxis_title='x (m)', yaxis_title='y (m)', zaxis_title='z (m)'),
                      margin=dict(l=0, r=0, t=0, b=0), legend=dict(orientation='h'))
    return fig


def render_png(parts: dict, path: str, elev=22, azim=-55):
    """Matplotlib snapshot for the Word report."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    fig = plt.figure(figsize=(6.5, 6.5))
    ax = fig.add_subplot(111, projection='3d')
    allv = np.vstack([v for v, _ in parts.values()])
    for name, (v, f) in parts.items():
        pc = Poly3DCollection(v[f], facecolor=PART_COLORS.get(name, '#888'),
                              alpha=0.18 if name in ('Tank wall', 'Roof dome', 'Floor') else 1.0,
                              edgecolor='none')
        ax.add_collection3d(pc)
    mid = allv.mean(0)
    span = (allv.max(0) - allv.min(0)).max() / 2
    ax.set_xlim(mid[0] - span, mid[0] + span); ax.set_ylim(mid[1] - span, mid[1] + span)
    ax.set_zlim(allv[:, 2].min(), allv[:, 2].min() + 2 * span)
    ax.view_init(elev, azim)
    ax.set_axis_off()
    fig.savefig(path, dpi=160, bbox_inches='tight')
    plt.close(fig)
