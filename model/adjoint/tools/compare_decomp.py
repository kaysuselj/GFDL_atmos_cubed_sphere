#!/usr/bin/env python3
"""Compare fields written per rank by test_adv (dump_field: <tag>_n<ranks>_t<tile>_i<is>_j<js>.r8, compute
domain, npz levels) between two decompositions: assemble each face from its pieces and report the
largest difference.  usage: compare_decomp.py <data dir> <tag> [<tag> ...]"""
import glob, os, re, sys
import numpy as np

N, NPZ = 24, 72


def faces(d, tag):
    """{layout piece size: {tile: (N, N, NPZ) array}} from the files of one tag"""
    out = {}
    for f in glob.glob(os.path.join(d, f'{tag}_n*_t*_i*_j*.r8')):
        t, i, j = map(int, re.search(r'_n\d+_t(\d+)_i(\d+)_j(\d+)\.r8$', f).groups())
        a = np.fromfile(f, '<f8')
        n = int(round((a.size / NPZ) ** 0.5))
        out.setdefault(n, {}).setdefault(t, np.full((N, N, NPZ), np.nan))[i-1:i-1+n, j-1:j-1+n] = \
            a.reshape((n, n, NPZ), order='F')
    return out


d = sys.argv[1]
for tag in sys.argv[2:]:
    f = faces(d, tag)
    sizes = sorted(f, reverse=True)
    if len(sizes) < 2:
        print(f'{tag}: only one decomposition found ({sizes})'); continue
    a, b = f[sizes[0]], f[sizes[1]]
    A = np.stack([a[t] for t in range(1, 7)]); B = np.stack([b[t] for t in range(1, 7)])
    assert not np.isnan(A).any() and not np.isnan(B).any()
    print(f'{tag}: {6*(N//sizes[0])**2} vs {6*(N//sizes[1])**2} ranks: bit-identical {np.array_equal(A, B)}, '
          f'max |diff| / max |field| = {np.abs(A-B).max()/np.abs(A).max():.2e}')
