#!/usr/bin/env python3
"""Extract one cube face and one advection step of a C24 trajectory (devkit traj/<case>) into raw
r8 files for the single-face TL/AD driver (model/adjoint/driver).

Everything is written in FV3 index order (Fortran order, little-endian r8, no record markers):

  dxa, dya, area      (isd:ied, jsd:jed)            30 x 30       grid metrics with halo
  dx                  (isd:ied, jsd:jed+1)          30 x 31
  dy                  (isd:ied+1, jsd:jed)          31 x 30
  sin_sg              (isd:ied, jsd:jed, 9)
  cx, mfx             (is:ie+1, jsd:jed, npz)       25 x 30 x 72  x-face values incl. y halo
  cy, mfy             (isd:ied, js:je+1, npz)       30 x 25 x 72  y-face values incl. x halo
  ple0, ple1, pleadv  (is:ie, js:je, npz+1)         dry edges, TOA first (pleadv: GCHP check value)
  q_in                (isd:ied, jsd:jed, npz)       CO2 before the step, TOA first; halo from the
                                                    neighbouring faces, cube-corner ghosts = 0
  q_out_gchp          (is:ie, js:je, npz)           CO2 after the step from GCHP (check value)
  ak, bk              (npz+1)                       fitted from the merra2 dry edges (GEOS 72L)

The halo values (q, CX/CY, MFX/MFY, metrics) are taken from the neighbouring faces by matching
grid points, so the rotation across face edges is handled geometrically: a flux through an edge is
expressed in this face's +i / +j direction. Grid metrics follow the FV3 formulas (fv_grid_tools,
fv_grid_utils) on the halo-extended grid. They are close to FV3's but not bit-identical in the
halo (FV3 fills cube-corner ghosts with fill_corners and patches sin_sg at the edges); that does
not matter for the transpose tests, which only need a well-defined operator.

usage: make_face_inputs.py --traj <devkit>/traj/8ms --face 4 --step 1 --out <dir>
       (face 1..6 as in GCHP, step 1..48; step 1 starts from restart_IC.nc4)
"""
import argparse, os
import numpy as np
import netCDF4 as nc

NG, N, NPZ = 3, 24, 72
NPX = N + 1
RADIUS = 6371.0e3          # MAPL_RADIUS


def xyz(lon, lat):
    return np.stack([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)], -1)


def key(p):
    return tuple(np.round(np.asarray(p) * 1e5).astype(np.int64))


def gc_dist(p, q):
    return RADIUS * np.arctan2(np.linalg.norm(np.cross(p, q), axis=-1), np.sum(p * q, axis=-1))


def mid(p, q):
    m = p + q
    return m / np.linalg.norm(m, axis=-1, keepdims=True)


def cos_angle(p1, p2, p3):
    """FV3 cos_angle: cos of the angle at p1 between the great circles p1->p2 and p1->p3."""
    pv = np.cross(p1, p2)
    qv = np.cross(p1, p3)
    ddd = np.sqrt(np.sum(pv * pv, -1) * np.sum(qv * qv, -1))
    return np.where(ddd > 0, np.sum(pv * qv, -1) / np.where(ddd > 0, ddd, 1.0), 1.0)


def sph_area(p1, p2, p3, p4):
    """Area of the spherical quadrilateral p1-p2-p3-p4 (two triangles, l'Huilier-free formula)."""
    def tri(a, b, c):
        num = np.abs(np.sum(a * np.cross(b, c), -1))
        den = 1 + np.sum(a * b, -1) + np.sum(b * c, -1) + np.sum(c * a, -1)
        return 2 * np.arctan2(num, den)
    return RADIUS ** 2 * (tri(p1, p2, p3) + tri(p1, p3, p4))


class Cube:
    """Corner points of the 6 faces, plus lookup tables for points, cells and stored C-grid edges."""

    def __init__(self, ds):
        lon = np.radians(np.asarray(ds['corner_lons'][:], 'f8'))
        lat = np.radians(np.asarray(ds['corner_lats'][:], 'f8'))
        self.C = xyz(lon, lat)                    # (face, J, I, 3), J/I = 0..24
        self.pt = {}                              # point key -> list of (g, I, J)
        self.cell = {}                            # frozenset of 4 point keys -> (g, I, J)
        self.edge = {}                            # frozenset of 2 point keys -> ('W'|'S', g, I, J)
        C = self.C
        for g in range(6):
            for J in range(NPX):
                for I in range(NPX):
                    self.pt.setdefault(key(C[g, J, I]), []).append((g, I, J))
            for J in range(N):
                for I in range(N):
                    k = [key(C[g, J, I]), key(C[g, J, I + 1]), key(C[g, J + 1, I]), key(C[g, J + 1, I + 1])]
                    self.cell[frozenset(k)] = (g, I, J)
                    self.edge[frozenset([k[0], k[2]])] = ('W', g, I, J)
                    self.edge[frozenset([k[0], k[1]])] = ('S', g, I, J)


def extended_grid(cube, f):
    """Corner points G[i, j] of face f, FV3 indices i, j = 1-NG .. NPX+NG (offset o = NG).
    Cube-corner ghost points are NaN."""
    o = NG
    G = np.full((NPX + 2 * NG, NPX + 2 * NG, 3), np.nan)
    G[o:o + NPX, o:o + NPX] = np.transpose(cube.C[f], (1, 0, 2))     # G[i-1+o, j-1+o] = C[f, j-1, i-1]
    sides = {   # side: (edge points as FV3 (i, j) list, outward unit step in (i, j))
        'W': ([(1, j) for j in range(1, NPX + 1)], (-1, 0)),
        'E': ([(NPX, j) for j in range(1, NPX + 1)], (1, 0)),
        'S': ([(i, 1) for i in range(1, NPX + 1)], (0, -1)),
        'N': ([(i, NPX) for i in range(1, NPX + 1)], (0, 1)),
    }
    for side, (pts, (di, dj)) in sides.items():
        cand = None
        for (i, j) in pts:
            gs = {g for (g, I, J) in cube.pt[key(G[i - 1 + o, j - 1 + o])] if g != f}
            cand = gs if cand is None else cand & gs
        assert len(cand) == 1, (side, cand)
        g = cand.pop()
        idx = [next((I, J) for (gg, I, J) in cube.pt[key(G[i - 1 + o, j - 1 + o])] if gg == g) for (i, j) in pts]
        Is = {I for I, J in idx}
        if len(Is) == 1:             # edge is a constant-I line of face g
            I0 = Is.pop(); inward = (1 if I0 == 0 else -1, 0)
        else:
            J0 = {J for I, J in idx}.pop(); inward = (0, 1 if J0 == 0 else -1)
        for (i, j), (I, J) in zip(pts, idx):
            for d in range(1, NG + 1):
                G[i + d * di - 1 + o, j + d * dj - 1 + o] = cube.C[g, J + d * inward[1], I + d * inward[0]]
    return G


def face_fields(cube, f, G, glob):
    """glob: dict name -> array (6, N, N[, ...]) of cell / west-edge / south-edge values (face, J, I, ...).
    Returns cell values with halo (q), x-edge values (cx, mfx) and y-edge values (cy, mfy)."""
    o = NG
    nd = NPX - 1 + 2 * NG          # 30 cells with halo

    def gp(i, j):
        return G[i - 1 + o, j - 1 + o]

    def cellkey(i, j):
        ps = [gp(i, j), gp(i + 1, j), gp(i, j + 1), gp(i + 1, j + 1)]
        if any(np.isnan(p).any() for p in ps):
            return None
        return frozenset(key(p) for p in ps)

    def cell_of(i, j):
        k = cellkey(i, j)
        return None if k is None else cube.cell[k]

    out = {}
    extra = glob['q'].shape[3:]
    q = np.zeros((nd, nd) + extra)
    for j in range(1 - NG, N + NG + 1):
        for i in range(1 - NG, N + NG + 1):
            c = cell_of(i, j)
            if c is not None:
                g, I, J = c
                q[i - 1 + o, j - 1 + o] = glob['q'][g, J, I]
    out['q'] = q

    def edge_value(p, r, plus, minus, names):
        """Value of the flux through the edge p-r, positive towards the cell 'plus' of face f."""
        if np.isnan(p).any() or np.isnan(r).any():
            return None
        typ, g, I, J = cube.edge[frozenset([key(p), key(r)])]
        stored_plus = (g, I, J)
        if plus is not None:
            sgn = 1.0 if plus == stored_plus else -1.0
        else:
            sgn = -1.0 if minus == stored_plus else 1.0
        src = names[0] if typ == 'W' else names[1]
        return sgn * glob[src][g, J, I]

    lev = extra
    for pair in (('cx', 'cy'), ('mfx', 'mfy')):
        ex = np.zeros((N + 1, nd) + lev)            # (is:ie+1, jsd:jed)
        for j in range(1 - NG, N + NG + 1):
            for i in range(1, N + 2):
                v = edge_value(gp(i, j), gp(i, j + 1), cell_of(i, j), cell_of(i - 1, j), pair)
                if v is not None:
                    ex[i - 1, j - 1 + o] = v
        ey = np.zeros((nd, N + 1) + lev)            # (isd:ied, js:je+1)
        for j in range(1, N + 2):
            for i in range(1 - NG, N + NG + 1):
                v = edge_value(gp(i, j), gp(i + 1, j), cell_of(i, j), cell_of(i, j - 1), pair)
                if v is not None:
                    ey[i - 1 + o, j - 1] = v
        out[pair[0]], out[pair[1]] = ex, ey
    return out


def metrics(G):
    """FV3 grid metrics on the extended grid. Cube-corner ghost cells get the values of the cell
    that copy_corners(dir=1) maps them to (dxa and dya swapped), so they are finite and positive."""
    o, nd = NG, N + 2 * NG
    P = lambda i, j: G[i - 1 + o, j - 1 + o]
    m = {k: np.full((nd, nd), np.nan) for k in ('dxa', 'dya', 'area')}
    m['dx'] = np.full((nd, nd + 1), np.nan)
    m['dy'] = np.full((nd + 1, nd), np.nan)
    sg = np.full((nd, nd, 9), np.nan)
    for j in range(1 - NG, N + NG + 2):
        for i in range(1 - NG, N + NG + 2):
            if i <= N + NG and j <= N + NG + 1 and not (np.isnan(P(i, j)).any() or np.isnan(P(i + 1, j)).any()):
                m['dx'][i - 1 + o, j - 1 + o] = gc_dist(P(i, j), P(i + 1, j))
            if j <= N + NG and i <= N + NG + 1 and not (np.isnan(P(i, j)).any() or np.isnan(P(i, j + 1)).any()):
                m['dy'][i - 1 + o, j - 1 + o] = gc_dist(P(i, j), P(i, j + 1))
    for j in range(1 - NG, N + NG + 1):
        for i in range(1 - NG, N + NG + 1):
            p1, p2, p3, p4 = P(i, j), P(i + 1, j), P(i + 1, j + 1), P(i, j + 1)
            if any(np.isnan(p).any() for p in (p1, p2, p3, p4)):
                continue
            c = p1 + p2 + p3 + p4; c /= np.linalg.norm(c)
            a = (i - 1 + o, j - 1 + o)
            m['dxa'][a] = gc_dist(mid(p1, p4), mid(p2, p3))
            m['dya'][a] = gc_dist(mid(p1, p2), mid(p4, p3))
            m['area'][a] = sph_area(p1, p2, p3, p4)
            cs = np.empty(9)
            cs[5] = cos_angle(p1, p2, p4)
            cs[6] = -cos_angle(p2, p1, p3)
            cs[7] = cos_angle(p3, p2, p4)
            cs[8] = -cos_angle(p4, p1, p3)
            cs[0] = cos_angle(mid(p1, p4), c, p4)
            cs[1] = cos_angle(mid(p1, p2), p2, c)
            cs[2] = cos_angle(mid(p2, p3), c, p2)
            cs[3] = cos_angle(mid(p4, p3), p4, c)
            e1 = mid(p2, p3) - mid(p1, p4); e1 -= c * (e1 @ c); e1 /= np.linalg.norm(e1)
            e2 = mid(p4, p3) - mid(p1, p2); e2 -= c * (e2 @ c); e2 /= np.linalg.norm(e2)
            cs[4] = e1 @ e2
            sg[a] = np.minimum(1.0, np.sqrt(np.maximum(0.0, 1.0 - cs ** 2)))
    # cube-corner ghost cells: copy_corners(dir=1) mapping
    for j in range(1 - NG, N + NG + 1):
        for i in range(1 - NG, N + NG + 1):
            a = (i - 1 + o, j - 1 + o)
            if not np.isnan(m['area'][a]):
                continue
            if i < 1 and j < 1:
                si, sj = j, 1 - i
            elif i > N and j < 1:
                si, sj = NPX - j, i - NPX + 1
            elif i > N and j > N:
                si, sj = j, 2 * NPX - 1 - i
            else:
                si, sj = NPX - j, i - 1 + NPX
            b = (si - 1 + o, sj - 1 + o)
            m['dxa'][a], m['dya'][a], m['area'][a] = m['dya'][b], m['dxa'][b], m['area'][b]
            sg[a] = sg[b]
    for k in ('dx', 'dy'):
        v = m[k]; v[np.isnan(v)] = np.nanmean(v)
    m['sin_sg'] = sg
    assert all(np.isfinite(v).all() for v in m.values())
    return m


def write(outdir, name, arr):
    np.asarray(arr, '<f8').ravel(order='F').tofile(os.path.join(outdir, name + '.r8'))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--traj', required=True)
    ap.add_argument('--face', type=int, default=4)
    ap.add_argument('--step', type=int, default=1)
    ap.add_argument('--out', required=True)
    ap.add_argument('--merra2', default=None, help='traj/merra2 dir for the ak/bk fit (default: ../merra2)')
    a = ap.parse_args()
    f, n = a.face - 1, a.step - 1
    os.makedirs(a.out, exist_ok=True)
    dC = nc.Dataset(os.path.join(a.traj, 'AdvTrajC.nc4'))
    dE = nc.Dataset(os.path.join(a.traj, 'AdvTrajE.nc4'))
    cube = Cube(dC)
    G = extended_grid(cube, f)

    rd = lambda ds, v, t: np.moveaxis(np.asarray(ds[v][t], 'f8'), 0, -1)      # (face, Y, X, lev)
    glob = {v.lower(): rd(dC, v + '_R4', n) for v in ('CX', 'CY', 'MFX', 'MFY')}
    if n == 0:
        r = nc.Dataset(os.path.join(a.traj, 'restart_IC.nc4'))
        q0 = np.asarray(r['SPC_CO2'][:], 'f8').reshape(NPZ, 6, N, N)          # TOA first, faces along y
        glob['q'] = np.moveaxis(q0, 0, -1)
    else:
        glob['q'] = rd(dC, 'SpeciesConcVV_CO2', n - 1)[..., ::-1]            # surface first -> TOA first
    fl = face_fields(cube, f, G, glob)
    m = metrics(G)

    for k in ('dxa', 'dya', 'area', 'dx', 'dy', 'sin_sg'):
        write(a.out, k, m[k])
    for k in ('cx', 'cy', 'mfx', 'mfy'):
        write(a.out, k, fl[k])
    write(a.out, 'q_in', fl['q'])
    T = lambda x: np.transpose(x, (1, 0, 2))                                   # (Y, X, lev) -> (i, j, lev)
    for v, nm in (('DryPLE0_R4', 'ple0'), ('DryPLE1_R4', 'ple1'), ('PLEadv_R4', 'pleadv')):
        write(a.out, nm, T(rd(dE, v, n)[f]))
    write(a.out, 'q_out_gchp', T(rd(dC, 'SpeciesConcVV_CO2', n)[f][..., ::-1]))

    mdir = a.merra2 or os.path.join(os.path.dirname(os.path.abspath(a.traj)), 'merra2')
    p = np.asarray(nc.Dataset(os.path.join(mdir, 'AdvTrajE.nc4'))['DryPLE1_R4'][0], 'f8').reshape(NPZ + 1, -1)
    A = np.vstack([np.ones(p.shape[1]), p[-1]]).T
    ak, bk = np.linalg.lstsq(A, p.T, rcond=None)[0]
    bk[0] = 0.0
    write(a.out, 'ak', ak); write(a.out, 'bk', bk)

    with open(os.path.join(a.out, 'meta.txt'), 'w') as fh:
        fh.write(f'face {a.face}\nstep {a.step}\nnpx {NPX}\nnpz {NPZ}\nng {NG}\nptop {ak[0]:.6f}\n')
        fh.write(f'traj {os.path.abspath(a.traj)}\n')
    # sanity: interior of q_in and the halo coverage
    qi = fl['q'][NG:NG + N, NG:NG + N]
    print(f'face {a.face} step {a.step}: q_in interior {qi.min():.6e}..{qi.max():.6e}, '
          f'halo cells set {(fl["q"][..., 0] != 0).sum() - N * N} of {(N + 2 * NG) ** 2 - N * N - 4 * NG * NG}')
    print(f'  max|cx| {np.abs(fl["cx"]).max():.3e}  max|cy| {np.abs(fl["cy"]).max():.3e}  '
          f'area {m["area"][NG:NG + N, NG:NG + N].sum() / (4 * np.pi * RADIUS ** 2) * 6:.6f} of a face')
    print(f'  ptop {ak[0]}  ak/bk fit residual {np.abs(A @ np.vstack([ak, bk]) - p.T).max():.3e} Pa')


if __name__ == '__main__':
    main()
