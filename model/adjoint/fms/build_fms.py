#!/usr/bin/env python3
"""Build the part of FMS (geoschem/FMS 259759d1) that fv_mp_mod.F90 needs, standalone, with the
CPP definitions GCHP uses for its fms_r4 library (FMS/CMakeLists.txt: OVERLOAD_R4 OVERLOAD_R8 MAPL_MODE SPMD
use_libMPI use_netCDF PORTABLE_KINDS HAVE_GETTID INTERNAL_FILE_NML, -fcray-pointer). FMS's own
CMakeLists.txt needs the GEOS ESMA build macros, so it is not used here.

usage: build_fms.py <FMS source dir> <build dir> [root modules ...]
Writes <build dir>/libfms.a and the .mod files. Default roots: what fv_mp_mod and the adjoint
driver use (fms_mod, fms_io_mod, mpp_mod, mpp_domains_mod, mpp_parameter_mod, ensemble_manager_mod).
"""
import os, re, subprocess, sys

FC = os.environ.get('MPIFC', 'mpif90')
DEFS = ['-DOVERLOAD_R4', '-DOVERLOAD_R8', '-DMAPL_MODE', '-DSPMD', '-Duse_libMPI', '-Duse_netCDF', '-DPORTABLE_KINDS',
        '-DHAVE_GETTID', '-DINTERNAL_FILE_NML', '-D__APPLE__']
FFLAGS = ['-O2', '-g', '-fcray-pointer', '-ffree-line-length-none', '-fallow-argument-mismatch',
          '-fallow-invalid-boz']


def nf_flags():
    inc = subprocess.run(['nf-config', '--includedir'], capture_output=True, text=True).stdout.strip()
    return ['-I' + inc]


def main():
    src, bld = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2])
    roots = sys.argv[3:] or ['fms_mod', 'fms_io_mod', 'mpp_mod', 'mpp_domains_mod', 'mpp_parameter_mod',
                             'ensemble_manager_mod']
    os.makedirs(bld, exist_ok=True)
    incs = ['-I' + os.path.join(src, d) for d in ('include', 'mpp/include', 'fms', 'fft', 'mpp')] + nf_flags()
    files = []
    for d, _, fs in os.walk(src):
        if '/test_' in d or '/.git' in d:
            continue
        files += [os.path.join(d, f) for f in fs if f.endswith(('.F90', '.f90'))]
    provides, needs = {}, {}
    for f in files:
        pp = subprocess.run([FC, '-E', '-cpp', '-P'] + DEFS + incs + [f], capture_output=True, text=True).stdout
        mods = set(m.lower() for m in re.findall(r'(?im)^\s*module\s+(?!procedure\b)(\w+)\s*$', pp))
        uses = set(u.lower() for u in re.findall(r'(?im)^\s*use\s*(?:,\s*\w+\s*::)?\s*(\w+)', pp))
        for m in mods:
            provides[m] = f
        needs[f] = uses - mods
    order, seen = [], set()

    def visit(m):
        f = provides.get(m)
        if f is None or f in seen:
            return
        seen.add(f)
        for u in sorted(needs[f]):
            visit(u)
        order.append(f)
    for r in roots:
        if r not in provides:
            sys.exit(f'module {r} not found')
        visit(r)
    objs = []
    for f in order:
        o = os.path.join(bld, os.path.basename(f).rsplit('.', 1)[0] + '.o')
        cmd = [FC, '-c', '-cpp'] + DEFS + FFLAGS + incs + ['-J' + bld, f, '-o', o]
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            print(r.stdout[-3000:], r.stderr[-3000:])
            sys.exit(f'failed: {f}')
        objs.append(o)
        print('compiled', os.path.relpath(f, src))
    # C helpers used by mpp and memutils (FMS/CMakeLists.txt)
    for c in ('memutils/memuse.c', 'mpp/threadloc.c', 'mpp/nsclock.c', 'affinity/affinity.c'):
        o = os.path.join(bld, os.path.basename(c)[:-2] + '_c.o')
        r = subprocess.run([os.environ.get('MPICC', 'mpicc'), '-c', '-O2'] + DEFS + incs + [os.path.join(src, c), '-o', o],
                           capture_output=True, text=True)
        if r.returncode != 0:
            print(r.stderr[-3000:])
            sys.exit(f'failed: {c}')
        objs.append(o)
        print('compiled', c)
    lib = os.path.join(bld, 'libfms.a')
    if os.path.exists(lib):
        os.remove(lib)
    subprocess.run(['ar', 'rcs', lib] + objs, check=True)
    print(f'{lib}: {len(objs)} objects')


if __name__ == '__main__':
    main()
