!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the FV3 dynamical core.
!*
!* The FV3 dynamical core is free software: you can redistribute it
!* and/or modify it under the terms of the
!* GNU Lesser General Public License as published by the
!* Free Software Foundation, either version 3 of the License, or
!* (at your option) any later version.
!*
!* The FV3 dynamical core is distributed in the hope that it will be
!* useful, but WITHOUT ANYWARRANTY; without even the implied warranty
!* of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!* See the GNU General Public License for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with the FV3 dynamical core.
!* If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************

!>@brief The module 'adv_adj_discrete_mod' is the interface between GCHP (AdvCore) and the
!! TAF-generated discrete TL/AD of one offline tracer-advection step (model/adjoint).
!>@details The generated code works on the data of adj_face_mod (stub types, see stubs/stubs.F90).
!! This module fills adj_face_mod from GCHP's real FV3 objects and calls it:
!!  - adv_adj_discrete_init: once; copies the grid metrics of gridstruct, the bounds, the tracer
!!    flags, ak/bk, and attaches the FMS back-end to GCHP's domain;
!!  - adv_adj_discrete_fwd:  one forward step of one tracer, q(is:ie,js:je,npz) -> q (the same
!!    operator as offline_tracer_advection; for checks and for the frozen mass-fixer factor);
!!  - adv_adj_discrete_adj:  one adjoint step of one tracer, lambda after -> lambda before.
!! All inputs are those of the FORWARD step (as offline_tracer_advection gets them in a forward
!! run): mfx, mfy, cx, cy on (is:ie, js:je, npz) with the forward sign, ple0 the dry pressure edges
!! at the start and ple1 at the end of the forward step, (is:ie, js:je, npz+1), TOA first.
!! The mass fixer: 'none' (factor 1), 'frozen' (the given factor s0, e.g. from
!! adv_adj_discrete_fwd at the forward trajectory) or 'linearized' (s(q) linearized around the
!! trajectory q_traj, the q before the forward step).
!!  - adv_adj_discrete_conv_check: checks the convention of the inputs (wind sign, which pressure
!!    edges are ple0 / ple1) against the forward run's PLEadv of the same step.
!! Only q is active; for the linear scheme (hord_tr 2, kord_tr 17) the adjoint is the exact
!! transpose of the forward step (standalone tests: STATUS.md of the devkit).
module adv_adj_discrete_mod

  use fv_arrays_mod,   only: fv_grid_type, fv_flags_type, fv_grid_bounds_type
  use mpp_domains_mod, only: domain2d
  use adj_face_mod,    only: af_gs => gridstruct, af_fl => flagstruct, af_bd => bd, &
                             af_is => is, af_ie => ie, af_js => js, af_je => je, &
                             af_isd => isd, af_ied => ied, af_jsd => jsd, af_jed => jed, &
                             af_tile => tile, af_npx => npx, af_npy => npy, af_npz => npz, af_dt => dt, &
                             af_ak => ak, af_bk => bk, af_ptop => ptop, &
                             cx_e, cy_e, mfx_e, mfy_e, af_ple0 => ple0, af_ple1 => ple1, &
                             pleadv_out, dpa_out, adj_freeze_s, adj_s0, adj_s_last
  use adj_backend_mod, only: backend_attach, backend_cgrid_halo

  implicit none
  private
  public :: adv_adj_discrete_init, adv_adj_discrete_fwd, adv_adj_discrete_adj, adv_adj_discrete_conv_check

  logical, save :: initialized = .false.

contains

  subroutine adv_adj_discrete_init(gridstruct, flagstruct, bd, domain, tile, npx, npy, npz, ak, bk, ptop, dt)
    type(fv_grid_type), intent(in) :: gridstruct
    type(fv_flags_type), intent(in) :: flagstruct
    type(fv_grid_bounds_type), intent(in) :: bd
    type(domain2d), intent(in) :: domain
    integer, intent(in) :: tile, npx, npy, npz
    real, intent(in) :: ak(npz+1), bk(npz+1), ptop, dt
    integer :: isd, ied, jsd, jed

    if (flagstruct%nord_tr /= 0 .or. flagstruct%trdm2 > 1.e-4) &
       call fail('tracer damping (nord_tr, trdm2) is not covered by the discrete adjoint')
    if (flagstruct%z_tracer) call fail('z_tracer = .true. (tracer_2d_1L) is not covered')
    if (gridstruct%nested) call fail('nested grids are not covered')

    call free()
    af_npx = npx; af_npy = npy; af_npz = npz; af_dt = dt; af_tile = tile
    af_is = bd%is; af_ie = bd%ie; af_js = bd%js; af_je = bd%je
    af_isd = bd%isd; af_ied = bd%ied; af_jsd = bd%jsd; af_jed = bd%jed
    af_bd%is = bd%is; af_bd%ie = bd%ie; af_bd%js = bd%js; af_bd%je = bd%je
    af_bd%isd = bd%isd; af_bd%ied = bd%ied; af_bd%jsd = bd%jsd; af_bd%jed = bd%jed
    isd = bd%isd; ied = bd%ied; jsd = bd%jsd; jed = bd%jed

    ! grid metrics used by tracer_2d / fv_tp_2d (nord_tr = 0: no deln_flux)
    allocate(af_gs%area(isd:ied,jsd:jed), af_gs%rarea(isd:ied,jsd:jed))
    allocate(af_gs%dxa(isd:ied,jsd:jed), af_gs%dya(isd:ied,jsd:jed))
    allocate(af_gs%dx(isd:ied,jsd:jed+1), af_gs%dy(isd:ied+1,jsd:jed))
    allocate(af_gs%sin_sg(isd:ied,jsd:jed,9), af_gs%area_64(isd:ied,jsd:jed))
    allocate(af_gs%rdxc(isd:ied+1,jsd:jed), af_gs%rdyc(isd:ied,jsd:jed+1))
    allocate(af_gs%del6_u(isd:ied,jsd:jed+1), af_gs%del6_v(isd:ied+1,jsd:jed))
    allocate(af_gs%cosa_s(isd:ied,jsd:jed), af_gs%rsin2(isd:ied,jsd:jed))
    af_gs%area    = gridstruct%area(isd:ied, jsd:jed)
    af_gs%rarea   = gridstruct%rarea(isd:ied, jsd:jed)
    af_gs%dxa     = gridstruct%dxa(isd:ied, jsd:jed)
    af_gs%dya     = gridstruct%dya(isd:ied, jsd:jed)
    af_gs%dx      = gridstruct%dx(isd:ied, jsd:jed+1)
    af_gs%dy      = gridstruct%dy(isd:ied+1, jsd:jed)
    af_gs%sin_sg  = gridstruct%sin_sg(isd:ied, jsd:jed, 1:9)
    af_gs%area_64 = gridstruct%area_64(isd:ied, jsd:jed)
    af_gs%rdxc = 0.; af_gs%rdyc = 0.; af_gs%del6_u = 0.; af_gs%del6_v = 0.
    af_gs%cosa_s = 0.; af_gs%rsin2 = 1.
    af_gs%da_min = gridstruct%da_min
    af_gs%nested = .false.
    af_gs%grid_type = gridstruct%grid_type
    af_gs%sw_corner = gridstruct%sw_corner; af_gs%se_corner = gridstruct%se_corner
    af_gs%ne_corner = gridstruct%ne_corner; af_gs%nw_corner = gridstruct%nw_corner

    af_fl%hord_tr = flagstruct%hord_tr; af_fl%kord_tr = flagstruct%kord_tr
    af_fl%nord_tr = flagstruct%nord_tr; af_fl%q_split = flagstruct%q_split
    af_fl%trdm2 = flagstruct%trdm2;     af_fl%lim_fac = flagstruct%lim_fac
    af_fl%fill = .false.;               af_fl%z_tracer = .false.

    allocate(af_ak(npz+1), af_bk(npz+1))
    af_ak = ak; af_bk = bk; af_ptop = ptop

    allocate(cx_e(af_is:af_ie+1, jsd:jed, npz), mfx_e(af_is:af_ie+1, jsd:jed, npz))
    allocate(cy_e(isd:ied, af_js:af_je+1, npz), mfy_e(isd:ied, af_js:af_je+1, npz))
    allocate(af_ple0(af_is:af_ie, af_js:af_je, npz+1), af_ple1(af_is:af_ie, af_js:af_je, npz+1))
    allocate(pleadv_out(af_is:af_ie, af_js:af_je, npz+1), dpa_out(af_is:af_ie, af_js:af_je, npz))

    call backend_attach(domain, bd%is, bd%ie, bd%js, bd%je, tile)
    initialized = .true.
  end subroutine adv_adj_discrete_init

  !> winds and pressure edges of one forward step (interior arrays, forward sign); the C-grid
  !! halos come from the FMS calls of offline_tracer_advection
  subroutine set_step(mfx, mfy, cx, cy, ple0, ple1)
    real, intent(in) :: mfx(af_is:af_ie, af_js:af_je, af_npz), mfy(af_is:af_ie, af_js:af_je, af_npz)
    real, intent(in) :: cx(af_is:af_ie, af_js:af_je, af_npz), cy(af_is:af_ie, af_js:af_je, af_npz)
    real, intent(in) :: ple0(af_is:af_ie, af_js:af_je, af_npz+1), ple1(af_is:af_ie, af_js:af_je, af_npz+1)
    if (.not. initialized) call fail('adv_adj_discrete_init has not been called')
    call backend_cgrid_halo(cx, cy, cx_e, cy_e)
    call backend_cgrid_halo(mfx, mfy, mfx_e, mfy_e)
    af_ple0 = ple0
    af_ple1 = ple1
  end subroutine set_step

  subroutine set_fixer(fixer, s0)
    character(len=*), intent(in) :: fixer
    real, intent(in), optional :: s0
    select case (trim(fixer))
    case ('none')
       adj_freeze_s = .true.; adj_s0 = 1.
    case ('frozen')
       if (.not. present(s0)) call fail('fixer = frozen needs s0')
       adj_freeze_s = .true.; adj_s0 = s0
    case ('linearized', 'full')
       adj_freeze_s = .false.
    case default
       call fail('fixer must be none, frozen or linearized (full in adv_adj_discrete_fwd)')
    end select
  end subroutine set_fixer

  !> one forward step of one tracer: q (before) -> q (after); s_out = the mass-fixer factor used
  !! (fixer = 'full': s(q) as in offline_tracer_advection; 'none' / 'frozen' as in set_fixer)
  subroutine adv_adj_discrete_fwd(q, mfx, mfy, cx, cy, ple0, ple1, fixer, s0, s_out)
    real, intent(inout) :: q(af_is:af_ie, af_js:af_je, af_npz)
    real, intent(in) :: mfx(af_is:af_ie, af_js:af_je, af_npz), mfy(af_is:af_ie, af_js:af_je, af_npz)
    real, intent(in) :: cx(af_is:af_ie, af_js:af_je, af_npz), cy(af_is:af_ie, af_js:af_je, af_npz)
    real, intent(in) :: ple0(af_is:af_ie, af_js:af_je, af_npz+1), ple1(af_is:af_ie, af_js:af_je, af_npz+1)
    character(len=*), intent(in) :: fixer
    real, intent(in), optional :: s0
    real, intent(out), optional :: s_out
    real :: qh(af_isd:af_ied, af_jsd:af_jed, af_npz, 1), fx(af_is:af_ie+1, af_js:af_je), fy(af_is:af_ie, af_js:af_je+1)
    call set_step(mfx, mfy, cx, cy, ple0, ple1)
    call set_fixer(fixer, s0)
    qh = 0.; fx = 0.; fy = 0.
    call top_adv(3, qh, q, fx, fy)
    if (present(s_out)) then
       if (adj_freeze_s) then
          s_out = adj_s0
       else
          s_out = adj_s_last
       endif
    endif
  end subroutine adv_adj_discrete_fwd

  !> one adjoint step of one tracer: lam (adjoint of q after the forward step) -> lam (adjoint of
  !! q before it). q_traj = q before the forward step; needed for fixer = 'linearized' only.
  subroutine adv_adj_discrete_adj(lam, mfx, mfy, cx, cy, ple0, ple1, fixer, s0, q_traj)
    real, intent(inout) :: lam(af_is:af_ie, af_js:af_je, af_npz)
    real, intent(in) :: mfx(af_is:af_ie, af_js:af_je, af_npz), mfy(af_is:af_ie, af_js:af_je, af_npz)
    real, intent(in) :: cx(af_is:af_ie, af_js:af_je, af_npz), cy(af_is:af_ie, af_js:af_je, af_npz)
    real, intent(in) :: ple0(af_is:af_ie, af_js:af_je, af_npz+1), ple1(af_is:af_ie, af_js:af_je, af_npz+1)
    character(len=*), intent(in) :: fixer
    real, intent(in), optional :: s0
    real, intent(in), optional :: q_traj(af_is:af_ie, af_js:af_je, af_npz)
    real :: q(af_is:af_ie, af_js:af_je, af_npz), lam4(af_is:af_ie, af_js:af_je, af_npz, 1)
    real :: qh(af_isd:af_ied, af_jsd:af_jed, af_npz, 1), qh_ad(af_isd:af_ied, af_jsd:af_jed, af_npz, 1)
    real :: fx(af_is:af_ie+1, af_js:af_je), fy(af_is:af_ie, af_js:af_je+1)
    real :: fx_ad(af_is:af_ie+1, af_js:af_je), fy_ad(af_is:af_ie, af_js:af_je+1)
    call set_step(mfx, mfy, cx, cy, ple0, ple1)
    call set_fixer(fixer, s0)
    if (.not. adj_freeze_s .and. .not. present(q_traj)) call fail('fixer = linearized needs q_traj')
    q = 0.
    if (present(q_traj)) q = q_traj
    qh = 0.; qh_ad = 0.; fx = 0.; fy = 0.; fx_ad = 0.; fy_ad = 0.
    lam4(:, :, :, 1) = lam
    call top_adv_ad(3, qh, qh_ad, q, lam4, fx, fx_ad, fy, fy_ad)
    lam = lam4(:, :, :, 1)
  end subroutine adv_adj_discrete_adj

  !> release what adj_face_mod holds (from an earlier init, or from the standalone loader)
  !> Convention check of the inputs (wind sign; which edges are ple0 / ple1): the horizontal step
  !! turns ple0 into the edges after advection (pleadv), which must equal the forward run's PLEadv
  !! export of the same step (pleadv_ref; 0.008 Pa in the standalone test, r4 output). Returns the
  !! largest |pleadv - pleadv_ref| (Pa, all levels and ranks) for the inputs as given (err) and with
  !! the winds' sign flipped (err_flip; the C24 tests give 1e2 - 1e3 Pa). Note: pleadv is not close
  !! to ple1 (GCHP's mass fluxes do not close the pressure budget; the remap and the mass fixer
  !! absorb the difference), so ple1 cannot serve as the reference.
  subroutine adv_adj_discrete_conv_check(mfx, mfy, cx, cy, ple0, ple1, pleadv_ref, err, err_flip)
    use adj_backend_mod, only: backend_allmax
    real, intent(in) :: mfx(af_is:af_ie, af_js:af_je, af_npz), mfy(af_is:af_ie, af_js:af_je, af_npz)
    real, intent(in) :: cx(af_is:af_ie, af_js:af_je, af_npz), cy(af_is:af_ie, af_js:af_je, af_npz)
    real, intent(in) :: ple0(af_is:af_ie, af_js:af_je, af_npz+1), ple1(af_is:af_ie, af_js:af_je, af_npz+1)
    real, intent(in) :: pleadv_ref(af_is:af_ie, af_js:af_je, af_npz+1)
    real, intent(out) :: err, err_flip
    real :: q(af_is:af_ie, af_js:af_je, af_npz)
    q = 0.
    call adv_adj_discrete_fwd(q, mfx, mfy, cx, cy, ple0, ple1, 'none')
    err = backend_allmax(maxval(abs(pleadv_out - pleadv_ref)))
    q = 0.
    call adv_adj_discrete_fwd(q, -mfx, -mfy, -cx, -cy, ple0, ple1, 'none')
    err_flip = backend_allmax(maxval(abs(pleadv_out - pleadv_ref)))
  end subroutine adv_adj_discrete_conv_check

  subroutine free()
    if (allocated(af_gs%area))    deallocate(af_gs%area)
    if (allocated(af_gs%rarea))   deallocate(af_gs%rarea)
    if (allocated(af_gs%dxa))     deallocate(af_gs%dxa)
    if (allocated(af_gs%dya))     deallocate(af_gs%dya)
    if (allocated(af_gs%dx))      deallocate(af_gs%dx)
    if (allocated(af_gs%dy))      deallocate(af_gs%dy)
    if (allocated(af_gs%sin_sg))  deallocate(af_gs%sin_sg)
    if (allocated(af_gs%area_64)) deallocate(af_gs%area_64)
    if (allocated(af_gs%rdxc))    deallocate(af_gs%rdxc)
    if (allocated(af_gs%rdyc))    deallocate(af_gs%rdyc)
    if (allocated(af_gs%del6_u))  deallocate(af_gs%del6_u)
    if (allocated(af_gs%del6_v))  deallocate(af_gs%del6_v)
    if (allocated(af_gs%cosa_s))  deallocate(af_gs%cosa_s)
    if (allocated(af_gs%rsin2))   deallocate(af_gs%rsin2)
    if (allocated(af_ak))         deallocate(af_ak)
    if (allocated(af_bk))         deallocate(af_bk)
    if (allocated(cx_e))          deallocate(cx_e)
    if (allocated(cy_e))          deallocate(cy_e)
    if (allocated(mfx_e))         deallocate(mfx_e)
    if (allocated(mfy_e))         deallocate(mfy_e)
    if (allocated(af_ple0))       deallocate(af_ple0)
    if (allocated(af_ple1))       deallocate(af_ple1)
    if (allocated(pleadv_out))    deallocate(pleadv_out)
    if (allocated(dpa_out))       deallocate(dpa_out)
  end subroutine free

  subroutine fail(msg)
    character(len=*), intent(in) :: msg
    write(*,*) 'adv_adj_discrete_mod: ', trim(msg)
    error stop 'adv_adj_discrete_mod'
  end subroutine fail

end module adv_adj_discrete_mod
