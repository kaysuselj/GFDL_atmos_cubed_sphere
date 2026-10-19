# Status: advection adjoint, option B (laptop)

Last updated 2026-10-19. Plan: `docs/FULL_ADVECTION_ADJOINT_SCOPE.md`. Kit contents: `README.md`, `HANDOFF.md`.

## 1. Decided by you

| date | decision |
|---|---|
| 2026-10-16 | Option B first: the exact transpose (TL + AD) of a *linear* version of the GCHP 14.8.0 FV3 tracer advection. Option A (trajectory-based adjoint of hord 12 / kord 8) comes later and reuses the harness. |
| 2026-10-16 | TL/AD are generated with TAF from our own 14.8.0 source. `fv3jedi_lm/` is only a reference: we do not ship it, which also avoids its GPL v2 licence. |
| 2026-10-16 | Start from the linear scheme in `baselines/run_P2lin.pbs` (`hord_tr = 2`, `kord_tr = 17`). |
| 2026-10-16 | Precision r8 (`-fdefault-real-8 -fdefault-double-8`), CPU only (no CUDA). |
| 2026-10-16 | Work on git branch `adv_adjoint_B` (from `adjoint_port_14.8.0` of the `kaysuselj/GFDL_atmos_cubed_sphere` fork). New code goes in a new directory (e.g. `model/adjoint/`). **The forward FV3 files are not modified.** |
| 2026-10-16 | Commits and pushes go to the kaysuselj forks only, and only after your approval. |
| 2026-10-16 | Local commits on `adv_adjoint_B` (in `advection_adjoint_14.8.0/fvdycore`) are OK without asking; **no push** without approval. |
| 2026-10-16 | Mass fixer in M_lin: s frozen at its trajectory value s0 (default, exactly linear), plus a switch for the full rank-1 TL/AD of s(q) for FD tests against the GCHP forward. |
| 2026-10-16 | TAF is the server version `/Users/ksuselj/TAF/staf` (FastOpt server over ssh). Uploading our FV3 source to the server is approved. |
| 2026-10-16 | Budget: **$800** of Claude usage for the laptop tasks (0–3). Claude estimates the spend from the session and stops at the limit, with the state written here. Rough estimate used by the end of task 3 (2026-10-18): ~$250–400 (check with `/cost`). |
| 2026-10-18 | Continue option B (task 2 onwards). Installing `netcdf-fortran` (Homebrew) and building FMS on the laptop is approved. |
| 2026-10-19 | Push `adv_adjoint_B` to the kaysuselj fork; write the GCHP adapter (task 4) on the laptop plus a handoff for Pleiades; keep the test data in `work/`. |
| 2026-10-16 | Single-face tests in two variants: (a) halo q as an active input (tests the edge/corner stencils and `copy_corners`), (b) halo q = 0 (a closed one-face operator). |

## 1b. Decided by Claude while you were away (please review; all reversible)

| decision | why |
|---|---|
| The unmodified forward files are compiled against single-face **stubs** of FMS and the FV3-side modules (`model/adjoint/stubs/stubs.F90`), and TAF gets these exact files. The GMAO `tp-core-driver` was used only as the model; its stubs were too thin for `tracer_2d`/`fv_mapz`. | TAF must differentiate GCHP's own code. |
| TAF differentiates **`top_step`** (new, `model/adjoint/driver/adj_tops.F90`), not `offline_tracer_advection`. `top_step` repeats the ~60 lines of glue code (pe1/pe2, kStart, mass fixer) and calls the unmodified `tracer_2d` and `map1_q2`. It is checked to be bit-identical to `offline_tracer_advection` (at -O0; 1e-16 at -O2 for hord 12). | TAF's reverse mode cannot handle the xL/yL buffer reuse and the GOTOs of the unused `mapn_tracer` inside `offline_tracer_advection`. The frozen-s0 switch is explicit there. |
| In the **TAF input copy** of `fv_mapz.F90` only: `map1_q2` is replaced by a GOTO-free version (same arithmetic), and `ppm_limiters` declares its assumed-size dummies with their length. The repository file is unchanged, and the TAF build's forward is bit-identical to the original. | The GOTOs make TAF lose the q dependence (irreducible CFG); TAF cannot store assumed-size arrays. |
| The generated TL/AD (`taf/gen/`) is committed. | It can then be built on Pleiades without TAF. |
| The grid metrics (dxa, dya, dx, dy, area, sin_sg) are recomputed in Python from the file's cell corners, with the halo taken from the neighbouring faces. They are not bit-identical to FV3's in the halo (fill_corners and sin_sg edge patches are approximated). | Not needed for the transpose tests. The forward still matches GCHP to 4e-8 – 2e-6. In GCHP (task 4) the real gridstruct is used anyway. |
| Task 2: the communication goes through three external **black-box routines** that TAF gets as flow directives (`stubs.F90`): `adj_halo_q` (q halo update; TL = itself, AD = `adj_halo_q_ad`), `adj_global_sum` (mass-fixer sums; self-adjoint), `adj_reduce_max` (Courant max; passive). The six-face back-end calls the real `fv_mp_mod` (`domain_decomp`, GCHP's group halo update) and uses FMS `mpp_update_domains_ad` as the adjoint. | TAF cannot see FMS. The same generated TL/AD then serves one face (no MPI) and six faces (6/24 ranks), and in GCHP (task 4) the black boxes map directly onto fv_mp_mod / FMS. |
| In the MPI build the stub modules are renamed by the preprocessor (`-Dmpp_mod=stb_mpp_mod`, ...), so the TAF code and the real FMS / fv_mp_mod link together. FMS 259759d1 is built standalone by `fms/build_fms.py` with GCHP's fms_r4 definitions (its CMakeLists needs the GEOS ESMA macros). | Smallest change that keeps the forward files and the generated code unchanged. |
| Two small post-processing steps in `taf/run_taf.sh`, both checked there: q_ad is made `intent(inout)` in the generated `start_halo_4d_ad` (TAF declares it `intent(in)` although the black-box adjoint overwrites it), and the generated black-box calls are checked against the back-ends' argument lists (implicit interfaces, so a mismatch would not be caught by the compiler; this caught one bug). | Generated code must match the hand-written back-ends. |
| The frozen mass-fixer factor `adj_s0` is only read in the differentiated code, s(q) goes to a separate `adj_s_last`. | Otherwise TAF makes `adj_s0` active and its TL copy carries over between calls. |
| GCHP's fvdycore CPP definitions (`-DMAPL_MODE -DSPMD -DTIMING`, hydrostatic) are used for all builds and the TAF input. | Fidelity. They do not affect the advection path (only `Lagrangian_to_Eulerian`). |

## 2. Done

### Task 0 audit (part 1): fv3jedi_lm vs 14.8.0 (`src/fvdycore/model`)

**tp_core** (`tp_core.F90` vs `fv3jedi_lm/tp_core_nlm.F90`)
- `hord = 2` (the linear scheme): `xppm`/`yppm` are **identical**. This includes the interior 4th-order edge values, the face-edge formulas `al(0), al(1), al(2)` / `al(npx-1..npx+1)` and the flux formula. In q (winds fixed), hord 2 is exactly linear.
- `hord >= 8` branch (8/10/11/12/13, with the `pmp/lac` Huynh constraint used by hord 12): identical.
- `fv_tp_2d`, `copy_corners`, `pert_ppm`, `deln_flux`: identical except the `lim_fac` argument.
- Differences:
  1. 14.8.0 adds `lim_fac`. It is used only in `mord == 1`, so hord 2 and hord 12 do not depend on it.
  2. 14.8.0 uses `mord = abs(iord)`. Negative `iord` clips `al >= 0`, and a new `mord == 1` scheme was added. JEDI has `iord == 7` (positivity fix of `al`) instead.
  3. ord 3/4: the `hi5/hi6` switch logic differs (14.8.0: `smt5(i-1).and.smt5(i)`, `smt6(i-1).or.smt6(i)`; JEDI: an upwind-biased combination). This does not matter for us.
  4. JEDI has `hord = 333`, a linear 3rd-order upwind scheme (Holdaway & Kent 2015) with no special face-edge treatment. 14.8.0 does not have it.
- JEDI **TLM** `XPPM_TLM`/`YPPM_TLM` contain only the linear branches (hord 1, 2, 333). The JEDI ADM contains all branches.

**fv_tracer2d** (`tracer_2d`, `offline_tracer_advection`)
- `tracer_2d`: the same algorithm (xfx/yfx, `cmax`, `nsplt`, per-sub-step `dp2`, `fv_tp_2d`, flux-form update), with two differences:
  1. 14.8.0 updates `dp1` while `it /= ksplt(k)` and returns `dpA = dp2` of the last sub-step. JEDI uses `it /= nsplt` and sets `dpA = dp1`. With `nsplt = 1`, JEDI's `dpA` is therefore the *pre*-advection dp, which is wrong for the remap.
  2. `lim_fac` is passed through in 14.8.0.
- `offline_tracer_advection`:
  - 14.8.0 zeroes `q3` and `dpL` including the halo (needed because the far corner ghosts are read by `al(0)`).
  - 14.8.0 has the optional `pleadv`, the `mapn_tracer` path for nq > 5, `rescale_trop_only` (`calcScalingFactorTrop`, kStart = first level with bk > 0) and no r4 rounding of the scaling factor.
  - JEDI rescales the full column with `calcScalingFactor`.
  - **JEDI has no TL/AD of `offline_tracer_advection` or of the mass fixer at all.**

**fv_mapz** (`map1_q2`, `scalar_profile`)
- `kord = 17` (`abs(kord) > 16`): the "perfectly linear" branch. It does the tridiagonal edge solve and returns before any limiter, so it is linear in q for given pressures. `qs` is not used for `iv = 0`. JEDI has the same branch.
- `map1_q2` is the same in both. It uses GOTO-based search loops; this needs a check that TAF accepts them, otherwise we write a structured copy.

### Consequence for the linear operator (scheme proposal, pending your OK)

With winds and pressures fixed, every part of one step is linear in q except the global mass fixer:

| part | `hord_tr = 2`, `kord_tr = 17` |
|---|---|
| halo update of q, `copy_corners` | linear |
| `fv_tp_2d` hord 2 (including the face-edge `al` formulas) | linear |
| `tracer_2d` dp-weighted update, sub-cycling (nsplt fixed from winds) | linear |
| `map1_q2` kord 17 | linear |
| mass fixer `q * s(q)`, where s is a ratio of two global sums of q | **nonlinear** (degree-1 homogeneous) |

Proposal: **M_lin = hord 2 + kord 17 + the mass fixer with s frozen at its trajectory value s0** (a scalar per step, so the transpose is trivially `s0 * λ`). Optionally, a switch for the full rank-1 TL/AD of s(q) for the FD tests against the GCHP forward. In the idealized cases s0 ≈ 1.

### Task 1 (single face): TAF TL/AD + tests — **done, all pass** (commits `7da74be`, `d0ef0dc` on `adv_adjoint_B`, local only)

Code: `advection_adjoint_14.8.0/fvdycore/model/adjoint/` (see its `README.md`). Data: `advection_adjoint_14.8.0/work/<case>_f<face>_s<step>/`.

- **Forward check** (`build/test_fwd`, hord 12 / kord 8, q halo from the neighbouring faces) against GCHP's trajectory: q to 3.7e-8 relative on 8ms face 4 step 1 (r4 file precision; the step changes q by 2.6e-6). PLEadv to 0.008 Pa. Other cases are in the table below. The surface layer is excluded, because GCHP's output there includes the surface fluxes. The mass-fixer factor is divided out, because one face's sums are not GCHP's global sums.
- **TL/AD** (`build_ad/test_adv`, hord 2 / kord 17), from TAF on the unmodified forward code. Results:

| operator | dot test, random x, y | dot test, point seeds y = e_k (3 cube corners, 2 face edges, 2nd ring, interior) | TL vs centred FD |
|---|---|---|---|
| `fv_tp_2d` (level 40) | 1.6e-15 (halo active), 4.9e-16 (halo 0) | – | 1.1e-14 |
| `tracer_2d` (72 levels) | 1.3e-14 / 1.8e-14 | ≤ 3.8e-16 | 3.2e-14 |
| full step, frozen s0 | 1.3e-14 | ≤ 8.4e-16 | 1.8e-13 |
| full step, s(q) linearised | 1.3e-14 | ≤ 9.4e-15 | 3.8e-13 |

  `run_tests.sh` on 6 cases (8ms f4 s1, f1 s24, f6 s48; tiltedrot f3 s10; merra2 f2 s30, f5 s12): all pass, worst dot-test difference 6.2e-14. Forward vs GCHP: 3.7e-8, 1.7e-7, 1.0e-7, 7.9e-8, 3.2e-7, 2.2e-6.

### Task 2 (six faces, FMS halo exchange): **done, all pass** (commit `03fb801`, local only)

Build: `make mpi` in `model/adjoint` (FMS from `fms/build_fms.py`, MPICH, netcdf-fortran). Run in a directory with an `input.nml` (fms_nml), e.g. `work/run_mpi`.

| test | result |
|---|---|
| q halo update alone (`build_mpi/test_halo`): ⟨H x, y⟩ vs ⟨x, Hᵀ y⟩, H = GCHP's group update, Hᵀ = FMS `mpp_update_domains_ad` | 5.0e-15 (6 ranks), 3.0e-15 (24 ranks) |
| FMS halos of q, CX/CY, MFX/MFY vs the geometric halos of the prep script (rotation across face edges) | identical (difference 0) |
| `build_mpi/test_adv`, 8ms step 1, all six faces, hord 2 / kord 17: fv_tp_2d, tracer_2d, full step (with the halo exchange and the global mass fixer, frozen s0 and s(q)); point seeds at cube corners on faces 1, 4, 6, edges, interior | all pass on 6 and 24 ranks, worst 1.8e-14; TL vs FD 4e-13 |
| 6 vs 24 ranks: forward, TL and AD fields | equal to 1.4e-16 / 2.0e-16 / 8.3e-16 (not bitwise: the global sums of the mass fixer depend on the decomposition, as in GCHP) |
| sub-cycling: merra2 step 12, winds × 3 (max \|C\| = 2.1, nsplt = 3), 24 ranks | all pass, worst 2.6e-14 (unscaled merra2: 7.1e-15) |

### Found on the way (for option A, not B)
- hord 12 / kord 8: the TL/AD of `fv_tp_2d` and `tracer_2d` are consistent (dot test 2e-15 / 3e-14). The full step **fails at 6.8e-6**, so the AD of the kord 8 remap (`scalar_profile` limiters) is not exact. TAF reports recomputation warnings (`taf/gen/taf.log`). Remedy for option A: TAF store directives in the TAF input copy, or restructuring.
- hord 12 TL vs FD differs by 1–30 % at every ε, because the limiter branches at the nearly uniform q0 differ from those taken with a random perturbation. This is the expected non-differentiability.

## 3. To do

### Task 0 (remaining)
- [x] Scheme agreed: hord_tr 2 + kord_tr 17 + mass fixer with frozen s0 (+ optional rank-1 TL).
- [x] Cloned `kaysuselj/GFDL_atmos_cubed_sphere` `adjoint_port_14.8.0` (ca57320) into `advection_adjoint_14.8.0/fvdycore`, branch `adv_adjoint_B` created. Model/tools files are byte-identical to the kit's `src/`.
- [ ] Confirm which checkpointed SPC field equals q before advection (needed in task 4).
- [x] TAF (server, `staf`, TAF 6.8.11) works. Derived types with allocatable components and pointer aliases are fine with `-f08`. Top-levels must be external subroutines.

### Task 1 (single face)
- [x] Driver, stubs, prep script, TAF TL/AD, dot and FD tests (see section 2).
- [ ] Look at the merra2 face 5 step 12 forward difference (2.2e-6, against ≤ 3e-7 elsewhere): r4 input precision, or the approximate halo metrics?
- [x] Sub-cycling (nsplt > 1): tested with merra2 winds × 3 (task 2).

### Task 2 — done (see section 2)
- [x] FMS built (MPICH 5.0.1, netcdf-fortran 4.6.4; Homebrew also updated gcc to 16.2).
- [x] Six faces on 6 and 24 ranks, halo exchange adjoint, global dot tests, cube corners, sub-cycling.

### Task 3 (standalone multi-step tests) — **done, all pass** (commits `b1f1fe5`, `5a44f8d`)
`build_mpi/test_chain`, 1 day (48 steps) of each trajectory, 24 ranks, J = q at receptor cells after the day, x = random field on the receptor level; FD vs ⟨λ, x⟩ (λ = adjoint run backwards through the 48 steps) vs ⟨TL x, e_r⟩; mass fixer frozen (s0 per step) and s(q). Outputs in `model/adjoint/results/`.

| case | receptors | max \|adj/FD − 1\| | max \|adj − TL\| / (‖x‖‖λ‖) |
|---|---|---|---|
| 8ms | P2a: westofseam, onseam, eastofseam, corner_west, corner_on, corner_east | 1.4e-12 | 4.5e-17 |
| tiltedrot | P2c: inside_f1, cross_f2_f1 (crosses the face 1/2 edge) | 1.1e-11 | 2.2e-17 |
| merra2 | T2 points tropical_10N, jet_30N; face-4 seam and two cube corners | 2.2e-6 (jet_30N only: J ≈ 1e-8 because x is on another level, FD round-off); others ≤ 8.8e-11 | 3.9e-17 |

For comparison, the current reversed-wind adjoint in GCHP: P2a_lin 48–1000 % at these receptors (5 days, inside GCHP; the standalone test is 1 day, with approximate halo metrics).

### Task 4 preparation on the laptop — **done** (2026-10-19; Pleiades steps: `fvdycore/model/adjoint/gchp/HANDOFF_PLEIADES.md`)
- Adapter `model/adjoint/gchp/adv_adj_discrete.F90`: init from GCHP's real `gridstruct`/`flagstruct`/`bd`/`domain` (copies the metrics into the stub-typed data of the generated code; the FMS back-end keeps a copy of GCHP's domain), one adjoint step per tracer, one forward step, and a convention check against the forward PLEadv.
- GCHP-like laptop build (`make gchp`: real `fv_arrays.F90` and `fv_mp_mod.F90`, the generated code as GCHP compiles it): the adapter is **bit-identical** to the direct calls (forward, frozen-s0 adjoint, linearized adjoint; 6 and 24 ranks, merra2). Convention check: 0.009–0.011 Pa against GCHP's forward PLEadv with the right inputs, 1.6e3–4.0e3 Pa with the winds flipped.
- `taf/gen_gchp/`: the generated code with all modules renamed (`stb_*`, `adj_*`; GCHP has the real `tp_core_mod` etc.) and lines wrapped at 132 columns (Intel), so GCHP needs no preprocessor renaming. Compiled on the laptop with `-ffree-line-length-132`.
- `gchp/sources.cmake` + `gchp/patches/fvdycore_CMakeLists.patch` (add the 9 sources to fvdycore when ADJOINT) and `gchp/patches/FVdycoreCubed_AdvCore_discrete_linear.patch` (`ADV_ADJ_METHOD: discrete_linear`, `ADV_ADJ_DISCRETE_FIXER: none`, `ADV_ADJ_DISCRETE_CHECK`). Both patches apply cleanly (`git apply --check`); **not compiled** (needs MAPL/GCHP).
- `tools/check_conventions.py`: finds, from History output of a forward and an adjoint run, the record mapping, the sign of MFX/CX and which forward edges DryPLE0/1 are. The patch assumes forward-step inputs = (−MFX, −MFY, −CX, −CY) and (DryPLE1, DryPLE0) after AdvCore's swap; to be confirmed on Pleiades first.
- Laptop-only quirk: the MPI test programs end with a barrier and `exit` (MPI_Finalize hangs after repeated loads, MPICH 5.0.1); once in ~10 runs mpirun then printed a termination message after all results were written. Not relevant for GCHP.
- Open: which checkpointed SPC field is q before advection (needed for `ADV_ADJ_DISCRETE_FIXER: linearized`; task 0 item, needs the GCHP checkpoint code on Pleiades).

### Later (scope doc §5)
- [ ] Task 2: halo / corner adjoints with FMS (6 and 24 ranks).
- [x] Task 3: full-step tests on traj/8ms, tiltedrot, merra2 (see above).
- [ ] Task 4: integration in AdvCore (`ADV_ADJ_METHOD: discrete_linear`), on Pleiades.
- [ ] Tasks 5–7: validation against `baselines/` (target: beat P2a_lin, 48–1000 % → ~1e-12), C48, timing.
