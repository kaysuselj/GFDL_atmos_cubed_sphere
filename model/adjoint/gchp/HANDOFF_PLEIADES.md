# Handoff to Pleiades: the discrete advection adjoint in GCHP (task 4)

Written 2026-10-19 on the laptop. The option-B code of `FULL_ADVECTION_ADJOINT_SCOPE.md` is on the
branch `adv_adjoint_B` of `kaysuselj/GFDL_atmos_cubed_sphere` (the fvdycore fork), under
`model/adjoint/`. The forward FV3 files are unchanged. Laptop status and all decisions:
`model/adjoint/STATUS_laptop.md` (a copy of the devkit's STATUS.md of 2026-10-19); code layout:
`model/adjoint/README.md`.

## 1. What is done and validated on the laptop

- TAF (6.8.11, FastOpt server) TL/AD of one offline advection step, generated from the unmodified
  `tp_core.F90`, `fv_tracer2d.F90`, `fv_mapz.F90`, `fv_fill.F90`. Covered: the q halo exchange
  (GCHP's group update; adjoint = FMS `mpp_update_domains_ad`), `tracer_2d` with sub-cycling,
  `fv_tp_2d` (`copy_corners`, `xppm`/`yppm`), `map1_q2`/`scalar_profile`, and the tropospheric
  mass fixer (frozen factor or linearized s(q)).
- Exact transpose of the linear scheme (`hord_tr 2`, `kord_tr 17`) on six faces, 6 and 24 ranks:
  dot tests ≤ 2e-14; 48-step receptor tests (8ms / P2a receptors, tiltedrot / P2c, merra2 / T2):
  adj/FD − 1 between 3e-14 and 1e-11 (2e-6 for one receptor with J ≈ 1e-8, FD round-off).
- The adapter `gchp/adv_adj_discrete.F90` (this directory), tested in a GCHP-like build with the
  real `fv_arrays.F90` and `fv_mp_mod.F90` (`make gchp`, `build_gchp/test_gchp_adapter`): it gives
  bit-identical results to the direct calls (forward, frozen-s0 adjoint, linearized adjoint).

## 2. Files for the GCHP build

| file | what |
|---|---|
| `model/adjoint/taf/gen_gchp/*_ad_gchp.F90` | the generated TL/AD as GCHP compiles it: all modules renamed (`stb_*` for the stand-ins of FMS/FV3 modules, `adj_tp_core_mod`, `adj_fv_tracer2d_mod`, `adj_fv_mapz_mod`, `adj_fv_fill_mod`), lines ≤ 132 columns. Do not edit; regenerate with `taf/run_taf.sh` (needs the TAF server, laptop). |
| `model/adjoint/comm/adj_comm_fms.F90` | communication back-end: GCHP's group halo update, `mpp_update_domains_ad`, `mpp_sum`, `mp_reduce_max`; `backend_attach` uses GCHP's domain |
| `model/adjoint/gchp/adv_adj_discrete.F90` | the adapter: `adv_adj_discrete_init` (once), `adv_adj_discrete_adj` (one adjoint step of one tracer), `adv_adj_discrete_fwd` (one forward step; checks), `adv_adj_discrete_conv_check` |
| `model/adjoint/gchp/sources.cmake` | adds these 9 files to the fvdycore library when `ADJOINT` is on |
| `model/adjoint/gchp/patches/fvdycore_CMakeLists.patch` | fvdycore `CMakeLists.txt`: `include` of `sources.cmake` (5 lines) |
| `model/adjoint/gchp/patches/FVdycoreCubed_AdvCore_discrete_linear.patch` | `AdvCore_GridCompMod.F90` of `kaysuselj/FVdycoreCubed_GridComp` (branch `adjoint_port_14.8.0`): `ADV_ADJ_METHOD: discrete_linear` |
| `model/adjoint/tools/check_conventions.py` | compares the advection inputs of an adjoint run with those of the forward run (see 4.2) |

## 3. Apply and build (Pleiades)

```bash
cd /home5/ksuselj1/nobackup/gchp_14.8.0_adj/code/src/GCHP_GridComp/FVdycoreCubed_GridComp/fvdycore
git fetch origin && git checkout adv_adjoint_B          # the laptop branch, pushed 2026-10-19
git apply model/adjoint/gchp/patches/fvdycore_CMakeLists.patch
cd ..                                                    # FVdycoreCubed_GridComp
git checkout -b adv_adjoint_B
git apply fvdycore/model/adjoint/gchp/patches/FVdycoreCubed_AdvCore_discrete_linear.patch
# then the usual adjoint build (ADJOINT=ON) into a new build directory, e.g. build_disc
```

Things that may need a fix in the first build (not testable on the laptop: no MAPL, gfortran only):
- `MAPL_CommsAllReduceMax` in the check block of the AdvCore patch: if MAPL has no such interface,
  use `MPI_ALLREDUCE(..., MPI_MAX, ...)`.
- Intel: the generated code is plain Fortran 2008 and compiled with gfortran `-ffree-line-length-132`;
  `stubs_ad_gchp.F90` uses `error stop`.
- The adapter assumes `real(FVPRC)` is the default real (R8 FV build, `-r8`), as the rest of fvdycore.
- FMS: the laptop used the fms_r4 library with OVERLOAD_R8 (r8 interfaces); GCHP's R8 FV links
  `fms_r8`; `mpp_update_domains_ad` exists in both.

## 4. Validation plan

### 4.1 Regression
With `ADV_ADJ_METHOD: reverse_wind_pure` everything must be bit-identical to the current binary
(the patch only adds a branch). Run P1a and one P2a adjoint, compare with `tests_final`.

### 4.2 Input convention (do this before trusting any adjoint number)
The adapter needs the inputs of the forward step that the adjoint step reverses. The AdvCore patch
assumes (with `ADV_ADJ_MASSFLUX: forward`):
`mfx_fwd = -MFX`, `mfy_fwd = -MFY`, `cx_fwd = -CX`, `cy_fwd = -CY` (as AdvCore gets them in the
adjoint run), `ple0_fwd = DryPLE1`, `ple1_fwd = DryPLE0` after AdvCore's own swap (i.e. the
imports DryPLE0 / DryPLE1 unchanged). Check it on real output:
1. Forward run of 8ms and merra2 with the AdvTrajC/E History collections (as `traj/run_traj.pbs`).
2. Adjoint run of the same window with the same collections (MFX_R4, MFY_R4, CX_R4, CY_R4,
   DryPLE0_R4, DryPLE1_R4 are GCHPctmEnv exports, i.e. what AdvCore imports).
3. `python3 model/adjoint/tools/check_conventions.py <forward dir> <adjoint dir>`: for every adjoint
   record it prints the matching forward step, the sign of MFX and CX, and which forward edges
   DryPLE0/1 of the adjoint run are. Expected: sign −1 for MFX and CX, DryPLE0 ↔ forward DryPLE0,
   DryPLE1 ↔ forward DryPLE1, and each adjoint record matching the forward step it reverses (this
   also confirms the order of the backward steps). If not, change the six assignment lines in the
   patch (`adjMfxF = -MFX`, ..., `adjPle0F = DryPLE1; adjPle1F = DryPLE0`).
   The adapter's `adv_adj_discrete_conv_check` gives a second check if the forward PLEadv of the
   step is at hand: the inputs must reproduce it to ~0.01 Pa (laptop: 0.009–0.011 Pa; with the
   winds flipped 1.6e3–4.0e3 Pa).

### 4.3 Consistency inside GCHP
`ADV_ADJ_DISCRETE_CHECK: 1`: every adjoint step also runs a forward step of a test field with the
real `offline_tracer_advection` and with the discrete code on the same inputs. Expected
`max |offline_tracer_advection - discrete forward|` ~1e-16 relative (round-off; the laptop gave
bit-identical results). This checks the grid copy, the halos and the build in GCHP.

### 4.4 The target: the linear-scheme dot test
`baselines/run_P2lin.pbs` (tests_final/P2_advection) with `METHOD=discrete_linear` (set inline in the
script, next to `DAYS`; it writes `ADV_ADJ_METHOD`; the `ADV_ADJ_WEIGHT` key it also writes is read
but not used by this method), `HORD=2 KORD=17` (the defaults), stages `base`, `fd`, `analyze`. Today (reverse_wind_pure): mean |adj/FD − 1| 48–87 % at the edges, 195–1000 % at
the corners. Expected with discrete_linear: ≲1e-6 (FD with eps 2e-5; the remaining mass-fixer
difference: the forward uses s(q), the adjoint `ADV_ADJ_DISCRETE_FIXER: none`, |s − 1| ~ 5e-12 per
step in the idealized cases).
Then P2a/P2c/T2/P2e with the linear scheme in both forward and adjoint (exact), and with the
operational forward (hord 12 / kord 8) against the linear adjoint (not exact: limiter; for
comparison with reverse_wind_pure).

### 4.5 Still open
- `ADV_ADJ_DISCRETE_FIXER: linearized` needs the q before each forward advection (the open task-0
  item: which checkpointed SPC field). Then call `adv_adj_discrete_adj(..., 'linearized', q_traj=q)`.
- Option A (operational scheme): the AD of the kord 8 remap is not exact yet (dot test 7e-6,
  STATUS_laptop.md); hord 12 horizontal is consistent.
- Timing (task 7): not measured yet. The generated AD keeps no tape; it recomputes forward values
  inside the adjoint step (taf/gen/taf.log lists the recomputations).

## 5. Laptop-only pieces (not needed on Pleiades)
`stubs/`, `comm/adj_comm_single.F90`, `comm/fv_arrays_fms.F90`, `driver/` (tests), `prep/`,
`fms/build_fms.py`, `taf/gen/` (the un-renamed TAF output), `results/`.
