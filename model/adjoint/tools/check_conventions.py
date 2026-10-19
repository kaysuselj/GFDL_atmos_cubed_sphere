#!/usr/bin/env python3
"""Find how the advection inputs of a GCHP adjoint run relate to those of the forward run, step by
step: for each record of the adjoint run's AdvTrajC/E files (MFX_R4, MFY_R4, CX_R4, CY_R4,
DryPLE0_R4, DryPLE1_R4, written by the same History collections as traj/run_traj.pbs), find the
forward record and the sign of each flux, and which forward edges DryPLE0/1 of the adjoint run are.
The adapter (gchp/adv_adj_discrete.F90) needs the forward-step inputs: the AdvCore patch assumes
mfx_fwd = -MFX_adj, cx_fwd = -CX_adj, ple0_fwd = DryPLE0_adj import, ple1_fwd = DryPLE1_adj import
(ADV_ADJ_MASSFLUX: forward). This script checks that on real output.

usage: check_conventions.py <forward dir> <adjoint dir>     (each with AdvTrajC.nc4, AdvTrajE.nc4)
       check_conventions.py <dir> <same dir>                  (self-test: identity, sign +1)
"""
import sys
import numpy as np
import netCDF4 as nc


def load(d):
    c = nc.Dataset(f'{d}/AdvTrajC.nc4'); e = nc.Dataset(f'{d}/AdvTrajE.nc4')
    f = {v: np.asarray(c[v + '_R4'][:], 'f8') for v in ('MFX', 'MFY', 'CX', 'CY')}
    f.update({v: np.asarray(e[v + '_R4'][:], 'f8') for v in ('DryPLE0', 'DryPLE1')})
    return f


def rel(a, b):
    return np.abs(a - b).max() / max(np.abs(b).max(), 1e-300)


fw, ad = load(sys.argv[1]), load(sys.argv[2])
nf, na = fw['MFX'].shape[0], ad['MFX'].shape[0]
print(f'forward records: {nf}, adjoint records: {na}')
print('adj rec -> fwd rec  sign(MFX) rel.diff   sign(CX) rel.diff   DryPLE0_adj ~ fwd   DryPLE1_adj ~ fwd')
for r in range(na):
    best = min(((rel(s*ad['MFX'][r], fw['MFX'][n]), n, s) for n in range(nf) for s in (1., -1.)))
    _, n, s = best
    sc = min(((rel(t*ad['CX'][r], fw['CX'][n]), t) for t in (1., -1.)))
    pe = []
    for v in ('DryPLE0', 'DryPLE1'):
        pe.append(min(((rel(ad[v][r], fw[w][n]), w) for w in ('DryPLE0', 'DryPLE1'))))
    print(f'{r+1:7d} -> {n+1:7d}  {s:+.0f} {best[0]:9.1e}   {sc[1]:+.0f} {sc[0]:9.1e}   '
          f'{pe[0][1]} {pe[0][0]:8.1e}   {pe[1][1]} {pe[1][0]:8.1e}')
