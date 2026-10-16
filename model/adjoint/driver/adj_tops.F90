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

!>@brief Top-level routines that TAF differentiates. The tracer q is the only active variable;
!! winds, pressures and grid come from adj_face_mod and are copied to locals on every call, because
!! the forward routines modify some of them in place (tracer_2d rescales cx/mfx when nsplt > 1 and
!! updates dp1).
!>@details
!! - top_tp2d:      fv_tp_2d on level k_test, q(isd:ied,jsd:jed) -> fluxes fx, fy
!! - top_tracer2d:  tracer_2d, q(isd:ied,jsd:jed,npz,nq) with halo -> q (interior is the result)
!! - top_offline:   offline_tracer_advection (horizontal step, remap, mass fixer),
!!                  q(is:ie,js:je,npz,nq) -> q; the q halo is zero inside, as in GCHP
module adj_tops_mod

  use adj_face_mod
  use tp_core_mod,      only: fv_tp_2d
  use fv_tracer2d_mod,  only: tracer_2d, offline_tracer_advection

  implicit none
  private
  public :: top_tp2d, top_tracer2d, top_offline, pleadv_out, dpa_out

  real, save :: pleadv_out(is:ie, js:je, npz+1)
  real, save :: dpa_out(is:ie, js:je, npz)

contains

  subroutine top_tp2d(q, fx, fy)
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

end module adj_tops_mod
