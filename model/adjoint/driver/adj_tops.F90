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

!>@brief Top-level routines that TAF differentiates (external subroutines, as TAF requires).
!! The tracer q is the only active variable; winds, pressures and grid come from adj_face_mod and are copied to locals on every call, because
!! the forward routines modify some of them in place (tracer_2d rescales cx/mfx when nsplt > 1 and
!! updates dp1).
!>@details
!! - top_tp2d:      fv_tp_2d on level k_test, q(isd:ied,jsd:jed) -> fluxes fx, fy
!! - top_tracer2d:  tracer_2d, q(isd:ied,jsd:jed,npz,nq) with halo -> q (interior is the result)
!! - top_offline:   offline_tracer_advection (horizontal step, remap, mass fixer),
!!                  q(is:ie,js:je,npz,nq) -> q; the q halo is zero inside, as in GCHP
!!                  (forward reference only; TAF differentiates top_step instead)
!! - top_step:      the same step, bit-identical, in a form TAF can differentiate (see there)
subroutine top_tp2d(q, fx, fy)
  use adj_face_mod
  use tp_core_mod, only: fv_tp_2d
  implicit none
  real, intent(inout) :: q(isd:ied, jsd:jed)
  real, intent(out)   :: fx(is:ie+1, js:je), fy(is:ie, js:je+1)
  real :: crx(is:ie+1, jsd:jed), xfx(is:ie+1, jsd:jed)
  real :: cry(isd:ied, js:je+1), yfx(isd:ied, js:je+1)
  real :: ra_x(is:ie, jsd:jed), ra_y(isd:ied, js:je)
  real :: mfx(is:ie+1, js:je), mfy(is:ie, js:je+1)
  integer :: i, j, k

  k = k_test
  ! as in tracer_2d (nsplt = 1)
  do j = jsd, jed
     do i = is, ie+1
        crx(i,j) = cx_e(i,j,k)
        if (crx(i,j) > 0.) then
           xfx(i,j) = crx(i,j)*gridstruct%dxa(i-1,j)*gridstruct%dy(i,j)*gridstruct%sin_sg(i-1,j,3)
        else
           xfx(i,j) = crx(i,j)*gridstruct%dxa(i,j)*gridstruct%dy(i,j)*gridstruct%sin_sg(i,j,1)
        endif
     enddo
  enddo
  do j = js, je+1
     do i = isd, ied
        cry(i,j) = cy_e(i,j,k)
        if (cry(i,j) > 0.) then
           yfx(i,j) = cry(i,j)*gridstruct%dya(i,j-1)*gridstruct%dx(i,j)*gridstruct%sin_sg(i,j-1,4)
        else
           yfx(i,j) = cry(i,j)*gridstruct%dya(i,j)*gridstruct%dx(i,j)*gridstruct%sin_sg(i,j,2)
        endif
     enddo
  enddo
  do j = jsd, jed
     do i = is, ie
        ra_x(i,j) = gridstruct%area(i,j) + (xfx(i,j) - xfx(i+1,j))
     enddo
  enddo
  do j = js, je
     do i = isd, ied
        ra_y(i,j) = gridstruct%area(i,j) + (yfx(i,j) - yfx(i,j+1))
     enddo
  enddo
  mfx = mfx_e(is:ie+1, js:je, k)
  mfy = mfy_e(is:ie, js:je+1, k)
  call fv_tp_2d(q, crx, cry, npx, npy, flagstruct%hord_tr, fx, fy, xfx, yfx, &
                gridstruct, bd, ra_x, ra_y, flagstruct%lim_fac, mfx=mfx, mfy=mfy)
end subroutine top_tp2d

subroutine top_tracer2d(q)
  use adj_face_mod
  use fv_tracer2d_mod, only: tracer_2d
  implicit none
  real, intent(inout) :: q(isd:ied, jsd:jed, npz, nq)
  real :: dp1(isd:ied, jsd:jed, npz)
  real :: mfxl(is:ie+1, js:je, npz), mfyl(is:ie, js:je+1, npz)
  real :: cxl(is:ie+1, jsd:jed, npz), cyl(isd:ied, js:je+1, npz)

  dp1 = 0.
  dp1(is:ie, js:je, :) = ple0(:, :, 2:npz+1) - ple0(:, :, 1:npz)
  mfxl = mfx_e(is:ie+1, js:je, :)
  mfyl = mfy_e(is:ie, js:je+1, :)
  cxl = cx_e
  cyl = cy_e
  call tracer_2d(q, dp1, mfxl, mfyl, cxl, cyl, gridstruct, bd, domain, npx, npy, npz, nq, &
                 flagstruct%hord_tr, flagstruct%q_split, dt, 0, q_pack, &
                 flagstruct%nord_tr, flagstruct%trdm2, flagstruct%lim_fac, dpA=dpa_out)
end subroutine top_tracer2d

subroutine top_offline(q)
  use adj_face_mod
  use fv_tracer2d_mod, only: offline_tracer_advection
  implicit none
  real, intent(inout) :: q(is:ie, js:je, npz, nq)
  real :: cx(is:ie, js:je, npz), cy(is:ie, js:je, npz)
  real :: mfx(is:ie, js:je, npz), mfy(is:ie, js:je, npz)
  real :: ple1l(is:ie, js:je, npz+1)

  cx  = cx_e(is:ie, js:je, :)
  cy  = cy_e(is:ie, js:je, :)
  mfx = mfx_e(is:ie, js:je, :)
  mfy = mfy_e(is:ie, js:je, :)
  ple1l = ple1
  call offline_tracer_advection(q, ple0, ple1l, mfx, mfy, cx, cy, gridstruct, flagstruct, bd, domain, &
                                ak, bk, ptop, npx, npy, npz, nq, dt, pleadv=pleadv_out)
end subroutine top_offline

!> One advection step as offline_tracer_advection (fv_tracer2d.F90) does it for nq = 1 and
!! z_tracer = .false.: tracer_2d, map1_q2 per j row, and the tropospheric mass fixer
!! (calcScalingFactorTrop, rescale_trop_only = .true.). It calls the unmodified tracer_2d and
!! map1_q2, and repeats the glue code of offline_tracer_advection statement by statement, so the
!! result is the same (test_fwd: bit-identical at -O0; at -O2 within 1e-16 for hord 12 / kord 8,
!! from code generation, and bit-identical for hord 2 / kord 17). Differences in form only:
!!  - the C-grid wind halos are taken from adj_face_mod directly (the values that the
!!    mpp_get_boundary / mpp_update_domains calls deliver); they are passive;
!!  - the mass-fixer factor is s(q) (adj_freeze_s = .false.; written to adj_s_last) or the given s0
!!    (adj_freeze_s = .true., s0 = adj_s0, e.g. adj_s_last of a forward call at the trajectory).
!! TAF cannot differentiate offline_tracer_advection itself in reverse mode (reuse of the xL/yL
!! buffers for CX and MFX, and the GOTOs in the unused mapn_tracer path).
subroutine top_step(q)
  use adj_face_mod
  use fv_mp_mod,       only: start_group_halo_update
  use fv_tracer2d_mod, only: tracer_2d
  use fv_mapz_mod,     only: map1_q2
  use mpp_mod,         only: mpp_sum
  implicit none
  real, intent(inout) :: q(is:ie, js:je, npz, nq)
  real :: cxl(is:ie+1, jsd:jed, npz), cyl(isd:ied, js:je+1, npz)
  real :: mfxl(is:ie+1, js:je, npz), mfyl(is:ie, js:je+1, npz)
  real :: dpl(isd:ied, jsd:jed, npz), dpa(is:ie, js:je, npz)
  real :: q1(is:ie, js:je, npz), q3(isd:ied, jsd:jed, npz, nq)
  real :: pe1(is:ie, npz+1), pe2(is:ie, npz+1), dp2(is:ie, js:je, npz)
  real :: partialsums(3, npz), globalsums(3), tempsum, scalingfactor
  real, parameter :: tiny_denominator = tiny(1.0)
  integer :: i, j, k, iq, kstart
  logical :: found

  cxl  = cx_e
  cyl  = cy_e
  mfxl = mfx_e(is:ie+1, js:je, :)
  mfyl = mfy_e(is:ie, js:je+1, :)
  dpl = 0.0
  dpl(is:ie, js:je, :) = ple0(:, :, 2:npz+1) - ple0(:, :, 1:npz)
  q3 = 0.0
  q3(is:ie, js:je, :, :) = q(is:ie, js:je, :, :)
  call start_group_halo_update(q_pack, q3, domain)
  call tracer_2d(q3, dpl, mfxl, mfyl, cxl, cyl, gridstruct, bd, domain, npx, npy, npz, nq, &
                 flagstruct%hord_tr, flagstruct%q_split, dt, 0, q_pack, &
                 flagstruct%nord_tr, flagstruct%trdm2, flagstruct%lim_fac, dpA=dpa)

  pleadv_out(:, :, 1) = ptop
  do k = 2, npz+1
     pleadv_out(:, :, k) = pleadv_out(:, :, k-1) + dpa(:, :, k-1)
  end do

  kstart = 0
  found = .false.
  do k = 1, npz
     if (.not. found .and. bk(k+1) > 0.d0) then
        kstart = k
        found = .true.
     end if
  end do

  do iq = 1, nq
     do j = js, je
        pe1(:, 1) = ptop
        do k = 2, npz+1
           pe1(:, k) = pe1(:, k-1) + dpa(:, j, k-1)
        enddo
        pe2(:, 1) = ptop
        pe2(:, npz+1) = pe1(:, npz+1)
        do k = 2, npz
           pe2(:, k) = ak(k) + bk(k)*pe1(:, npz+1)
        enddo
        do k = 1, npz
           dp2(:, j, k) = pe2(:, k+1) - pe2(:, k)
        enddo
        call map1_q2(npz, pe1, q3(isd, jsd, 1, iq), npz, pe2, q1(:, j, :), dp2(:, j, :), &
                     is, ie, 0, flagstruct%kord_tr, j, isd, ied, jsd, jed, 0.)
        q3(is:ie, j, 1:npz, iq) = q1(:, j, :)
     enddo
  enddo

  do iq = 1, nq
     if (adj_freeze_s) then
        scalingfactor = adj_s0
     else
        ! calcScalingFactorTrop
        partialsums = 0.0d0
        do k = 1, npz
           partialsums(1, k) = sum(q3(is:ie, js:je, k, iq)*dp2(:, :, k)*gridstruct%area(is:ie, js:je))
           tempsum = sum(q3(is:ie, js:je, k, iq)*(ple1(:, :, k+1) - ple1(:, :, k))*gridstruct%area(is:ie, js:je))
           partialsums(2, k) = tempsum
           if (k .lt. kstart) then
              partialsums(3, k) = tempsum
           end if
        end do
        globalsums(1) = sum(partialsums(1, :))
        globalsums(2) = sum(partialsums(2, :))
        globalsums(3) = sum(partialsums(3, :))
        call mpp_sum(globalsums, 3)
        globalsums(1) = globalsums(1) - globalsums(3)
        globalsums(2) = globalsums(2) - globalsums(3)
        if (globalsums(2) > tiny_denominator) then
           scalingfactor = globalsums(1) / globalsums(2)
        else
           scalingfactor = 1.d0
        end if
        adj_s_last = scalingfactor
     endif
     q(is:ie, js:je, 1:(kstart-1), iq) = q3(is:ie, js:je, 1:(kstart-1), iq)
     q(is:ie, js:je, kstart:npz, iq) = q3(is:ie, js:je, kstart:npz, iq) * scalingfactor
  enddo
end subroutine top_step

!> Single entry point for TAF (one generated TL/AD code for all test operators):
!!   mode 1: top_tp2d     (qh(:,:,k_test) -> fx, fy)
!!   mode 2: top_tracer2d (qh with halo -> qh)
!!   mode 3: top_step     (qc -> qc, full step as offline_tracer_advection)
subroutine top_adv(mode, qh, qc, fx, fy)
  use adj_face_mod
  implicit none
  integer, intent(in) :: mode
  real, intent(inout) :: qh(isd:ied, jsd:jed, npz, nq)
  real, intent(inout) :: qc(is:ie, js:je, npz, nq)
  real, intent(inout) :: fx(is:ie+1, js:je), fy(is:ie, js:je+1)
  if (mode == 1) then
     call top_tp2d(qh(isd, jsd, k_test, 1), fx, fy)
  elseif (mode == 2) then
     call top_tracer2d(qh)
  elseif (mode == 3) then
     call top_step(qc)
  endif
end subroutine top_adv
