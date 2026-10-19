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

!>@brief TL-vs-FD and dot-product tests of the TAF-generated TL/AD (taf/gen), on one face (single
!! back-end) or on the six faces (FMS back-end, 6 or 24 ranks).
!>@details usage: test_adv <data dir> [hord_tr kord_tr [layout [wind_scale]]]
!!   defaults 2 17 (the linear scheme), layout 1 (ranks per face = layout**2, FMS build only),
!!   wind_scale 1 (scales CX/CY/MFX/MFY; > 1 forces sub-cycling, nsplt = int(1 + max|C|)).
!! Operators (top_adv mode):
!!   1  fv_tp_2d on level k_test:  q (with halo) -> fx, fy
!!   2  tracer_2d:                 q (with halo) -> q (interior)
!!   3  full step (top_step):      q -> q, with the q halo exchange, mass fixer frozen (s0) or s(q)
!! Dot test: <TL x, y> against <x, AD y> (sums over all ranks), for x, y pseudo-random functions of
!! the global cell index (the same fields for any decomposition) and for y = e_k at cube corners,
!! face edges and interior cells. In mode 1 and 2, x is tested with the halo active (variant a) and
!! with the halo zero (variant b); on six faces mode 3 includes the halo exchange and its adjoint.
!! TL vs FD: centred differences (F(q0 + eps x) - F(q0 - eps x)) / (2 eps).
!! Mode 3 writes TL x and AD y per rank (tl_step_*, ad_step_*) for the decomposition comparison.
program test_adv
  use adj_face_mod
  use adj_load_mod,    only: load_step, dump_field
  use adj_backend_mod, only: multiface, backend_init, backend_finalize, backend_rank, backend_nranks, &
                             backend_allsum, backend_allmax
  implicit none
  character(len=512) :: dir, arg
  real, allocatable :: qh0(:,:,:,:), qc0(:,:,:,:)
  integer :: nfail, layout
  real :: wscale, cmax
  logical :: root

  layout = 1; wscale = 1.
  if (command_argument_count() >= 4) then
     call get_command_argument(4, arg); read(arg, *) layout
  endif
  if (command_argument_count() >= 5) then
     call get_command_argument(5, arg); read(arg, *) wscale
  endif
  call backend_init(layout)
  root = backend_rank() == 0
  call get_command_argument(1, dir)
  call load_step(trim(dir))
  flagstruct%hord_tr = 2; flagstruct%kord_tr = 17
  if (command_argument_count() >= 3) then
     call get_command_argument(2, arg); read(arg, *) flagstruct%hord_tr
     call get_command_argument(3, arg); read(arg, *) flagstruct%kord_tr
  endif
  if (wscale /= 1.) then
     cx_e = wscale*cx_e; cy_e = wscale*cy_e; mfx_e = wscale*mfx_e; mfy_e = wscale*mfy_e
  endif
  cmax = backend_allmax(max(maxval(abs(cx_e(is:ie, js:je, :))), maxval(abs(cy_e(is:ie, js:je, :)))))
  if (root) then
     write(*,'(a,a)') 'data: ', trim(dir)
     write(*,'(a,2i4,a,i3,a,i3,a,l2)') 'hord_tr, kord_tr =', flagstruct%hord_tr, flagstruct%kord_tr, &
          '   k_test =', k_test, '   ranks =', backend_nranks(), '   six faces:', multiface
     write(*,'(a,f6.2,a,f6.3,a,i2)') 'wind scale', wscale, ', max |C| =', cmax, ', nsplt =', int(1. + cmax)
  endif

  allocate(qh0(isd:ied, jsd:jed, npz, nq), qc0(is:ie, js:je, npz, nq))
  qh0(:, :, :, 1) = q_in
  qc0(:, :, :, 1) = q_in(is:ie, js:je, :)
  nfail = 0

  call test_mode(1, 'fv_tp_2d (level k_test)')
  call test_mode(2, 'tracer_2d')
  ! full step: first a forward with s(q) to store s0, then frozen s0 (linear) and full s(q)
  adj_freeze_s = .false.
  call fwd_dump()
  adj_s0 = adj_s_last
  if (root) write(*,'(a,es24.16)') 'mass-fixer factor s0 =', adj_s0
  adj_freeze_s = .true.
  call test_mode(3, 'full step, frozen s0')
  adj_freeze_s = .false.
  call test_mode(3, 'full step, s(q) linearized')

  if (root) then
     write(*,*)
     if (nfail == 0) then
        write(*,'(a)') 'ALL TESTS PASSED'
     else
        write(*,'(a,i3,a)') 'FAILED: ', nfail, ' test(s)'
     endif
  endif
  call backend_finalize()

contains

  !> forward copies (the top routines overwrite their arguments)
  subroutine fwd(mode, qh, qc, fx, fy)
    integer, intent(in) :: mode
    real, intent(inout) :: qh(isd:ied, jsd:jed, npz, nq), qc(is:ie, js:je, npz, nq)
    real, intent(out) :: fx(is:ie+1, js:je), fy(is:ie, js:je+1)
    fx = 0.; fy = 0.
    call top_adv(mode, qh, qc, fx, fy)
  end subroutine fwd

  !> one forward step from q_in (with s(q)), written for the comparison with the original forward
  subroutine fwd_dump()
    real :: qh(isd:ied, jsd:jed, npz, nq), qc(is:ie, js:je, npz, nq), fx(is:ie+1, js:je), fy(is:ie, js:je+1)
    character(len=32) :: tag
    qh = qh0; qc = qc0
    call fwd(3, qh, qc, fx, fy)
    write(tag, '(a,i0,a,i0)') 'step_taf_h', flagstruct%hord_tr, '_k', flagstruct%kord_tr
    call dump_field(trim(dir), trim(tag), qc, size(qc))
  end subroutine fwd_dump

  real function rel(a, b)
    real, intent(in) :: a, b
    rel = abs(a - b) / max(abs(a), abs(b), tiny(1.0))
  end function rel

  !> output-space inner product over all ranks (mode 1: fluxes; 2: q interior; 3: q)
  real function dot_out(mode, qh, qc, fx, fy, yh, yc, yfx, yfy)
    integer, intent(in) :: mode
    real, intent(in) :: qh(isd:ied, jsd:jed, npz, nq), qc(is:ie, js:je, npz, nq)
    real, intent(in) :: fx(is:ie+1, js:je), fy(is:ie, js:je+1)
    real, intent(in) :: yh(isd:ied, jsd:jed, npz, nq), yc(is:ie, js:je, npz, nq)
    real, intent(in) :: yfx(is:ie+1, js:je), yfy(is:ie, js:je+1)
    real :: s
    if (mode == 1) then
       s = sum(fx*yfx) + sum(fy*yfy)
    elseif (mode == 2) then
       s = sum(qh(is:ie, js:je, :, :)*yh(is:ie, js:je, :, :))
    else
       s = sum(qc*yc)
    endif
    dot_out = backend_allsum(s)
  end function dot_out

  !> input-space inner product over all ranks: with halo (modes 1, 2) or interior (mode 3)
  real function dot_in(mode, xh, xc, ah, ac)
    integer, intent(in) :: mode
    real, intent(in) :: xh(isd:ied, jsd:jed, npz, nq), xc(is:ie, js:je, npz, nq)
    real, intent(in) :: ah(isd:ied, jsd:jed, npz, nq), ac(is:ie, js:je, npz, nq)
    if (mode < 3) then
       dot_in = backend_allsum(sum(xh*ah))
    else
       dot_in = backend_allsum(sum(xc*ac))
    endif
  end function dot_in

  subroutine test_mode(mode, title)
    integer, intent(in) :: mode
    character(len=*), intent(in) :: title
    real :: xh(isd:ied, jsd:jed, npz, nq), xc(is:ie, js:je, npz, nq)
    real :: yh(isd:ied, jsd:jed, npz, nq), yc(is:ie, js:je, npz, nq)
    real :: yfx(is:ie+1, js:je), yfy(is:ie, js:je+1)
    real :: qh_tl(isd:ied, jsd:jed, npz, nq), qc_tl(is:ie, js:je, npz, nq)
    real :: qh_ad(isd:ied, jsd:jed, npz, nq), qc_ad(is:ie, js:je, npz, nq)
    real :: fx_tl(is:ie+1, js:je), fy_tl(is:ie, js:je+1)
    real :: qhp(isd:ied, jsd:jed, npz, nq), qcp(is:ie, js:je, npz, nq), fxp(is:ie+1, js:je), fyp(is:ie, js:je+1)
    real :: qhm(isd:ied, jsd:jed, npz, nq), qcm(is:ie, js:je, npz, nq), fxm(is:ie+1, js:je), fym(is:ie, js:je+1)
    real :: lhs, rhs, eps, fd, tl, scale
    integer :: ivar, ie_, ip, kk
    integer, parameter :: npt = 7
    integer :: pts(4, npt)
    character(len=12) :: ptname(npt)
    character(len=32) :: tag
    real, parameter :: tol_dot = 1.e-12
    integer, parameter :: nf = npx - 1

    if (root) write(*,*)
    if (root) write(*,'(a,i2,2a)') '=== mode ', mode, ': ', title
    kk = k_test
    ! point seeds (tile, i, j, k); on one face the tile is 1
    pts(:,1) = (/ 1, 1, 1, kk /);         ptname(1) = 'SW corner'
    pts(:,2) = (/ 1, nf, nf, kk /);       ptname(2) = 'NE corner'
    pts(:,3) = (/ 1, 1, nf/2, kk /);      ptname(3) = 'W edge'
    pts(:,4) = (/ 1, nf/2, nf, kk /);     ptname(4) = 'N edge'
    pts(:,5) = (/ 1, 2, 2, kk /);         ptname(5) = '2nd ring'
    pts(:,6) = (/ 1, nf/2, nf/2, kk /);   ptname(6) = 'interior'
    pts(:,7) = (/ 1, nf, 1, kk /);        ptname(7) = 'SE corner'
    if (multiface) then
       pts(1,2) = 4; pts(1,3) = 4; pts(1,7) = 6      ! other faces, other cube corners
       pts(:,6) = (/ 2, nf/2+1, nf/2+1, kk /)        ! on the ranks' boundary for layout 2
    endif

    ! ---- dot tests: x (variant a: halo active; b: halo zero), y
    do ivar = 1, merge(2, 1, mode < 3)
       call field_h(xh, 11 + ivar); call field_c(xc, 13)
       if (ivar == 2) call zero_halo(xh)
       if (mode == 1) call one_level(xh, kk)
       call field_h(yh, 21); call field_c(yc, 23); call field_fx(yfx, 25); call field_fy(yfy, 27)
       call run_tl(mode, xh, xc, qh_tl, qc_tl, fx_tl, fy_tl)
       lhs = dot_out(mode, qh_tl, qc_tl, fx_tl, fy_tl, yh, yc, yfx, yfy)
       call run_ad(mode, yh, yc, yfx, yfy, qh_ad, qc_ad)
       rhs = dot_in(mode, xh, xc, qh_ad, qc_ad)
       call report('dot, x'//merge(' (halo active)', ' (halo zero)  ', ivar == 1)//', y', lhs, rhs, tol_dot)
       if (mode == 3 .and. ivar == 1) then
          write(tag, '(a,l1)') 'tl_step_frozen', adj_freeze_s
          call dump_field(trim(dir), trim(tag), qc_tl, size(qc_tl))
          write(tag, '(a,l1)') 'ad_step_frozen', adj_freeze_s
          call dump_field(trim(dir), trim(tag), qc_ad, size(qc_ad))
       endif
    enddo

    ! ---- dot tests with point seeds y = e_k (q output modes), x with halo active
    if (mode >= 2) then
       call field_h(xh, 31); call field_c(xc, 33)
       call run_tl(mode, xh, xc, qh_tl, qc_tl, fx_tl, fy_tl)
       do ip = 1, npt
          yh = 0.; yc = 0.; yfx = 0.; yfy = 0.
          if (owns(pts(:,ip))) then
             yh(pts(2,ip), pts(3,ip), pts(4,ip), 1) = 1.
             yc(pts(2,ip), pts(3,ip), pts(4,ip), 1) = 1.
          endif
          lhs = dot_out(mode, qh_tl, qc_tl, fx_tl, fy_tl, yh, yc, yfx, yfy)
          call run_ad(mode, yh, yc, yfx, yfy, qh_ad, qc_ad)
          rhs = dot_in(mode, xh, xc, qh_ad, qc_ad)
          tag = ''
          if (multiface) write(tag, '(a,i0)') ' tile ', pts(1,ip)
          call report('dot, y = e_k at '//trim(ptname(ip))//trim(tag), lhs, rhs, tol_dot)
       enddo
    endif

    ! ---- TL vs centred FD, eps scan
    call field_h(xh, 41); call field_c(xc, 43)
    if (mode == 1) call one_level(xh, kk)
    call run_tl(mode, xh, xc, qh_tl, qc_tl, fx_tl, fy_tl)
    call field_h(yh, 51); call field_c(yc, 53); call field_fx(yfx, 55); call field_fy(yfy, 57)
    tl = dot_out(mode, qh_tl, qc_tl, fx_tl, fy_tl, yh, yc, yfx, yfy)
    scale = backend_allmax(maxval(abs(q_in(is:ie, js:je, :))))
    if (root) write(*,'(a)') '  TL vs FD  <(F(q0+eps x) - F(q0-eps x))/(2 eps), y>  vs  <TL x, y>'
    do ie_ = 1, 5
       eps = 10.**(-ie_)
       ! x is O(1), scaled by max q0, so eps is a relative perturbation
       qhp = qh0 + eps*scale*xh; qcp = qc0 + eps*scale*xc
       qhm = qh0 - eps*scale*xh; qcm = qc0 - eps*scale*xc
       call fwd(mode, qhp, qcp, fxp, fyp); call fwd(mode, qhm, qcm, fxm, fym)
       fd = (dot_out(mode, qhp, qcp, fxp, fyp, yh, yc, yfx, yfy) - &
             dot_out(mode, qhm, qcm, fxm, fym, yh, yc, yfx, yfy)) / (2.*eps*scale)
       if (root) write(*,'(a,es8.1,a,es24.16,a,es10.2)') '    eps ', eps, '  FD ', fd, '  rel. diff ', rel(fd, tl)
    enddo
    if (root) write(*,'(a,es24.16)') '    TL            ', tl
  end subroutine test_mode

  subroutine run_tl(mode, xh, xc, qh_tl, qc_tl, fx_tl, fy_tl)
    integer, intent(in) :: mode
    real, intent(in) :: xh(isd:ied, jsd:jed, npz, nq), xc(is:ie, js:je, npz, nq)
    real, intent(out) :: qh_tl(isd:ied, jsd:jed, npz, nq), qc_tl(is:ie, js:je, npz, nq)
    real, intent(out) :: fx_tl(is:ie+1, js:je), fy_tl(is:ie, js:je+1)
    real :: qh(isd:ied, jsd:jed, npz, nq), qc(is:ie, js:je, npz, nq)
    real :: fx(is:ie+1, js:je), fy(is:ie, js:je+1)
    qh = qh0; qc = qc0; fx = 0.; fy = 0.
    qh_tl = xh; qc_tl = xc; fx_tl = 0.; fy_tl = 0.
    call top_adv_tl(mode, qh, qh_tl, qc, qc_tl, fx, fx_tl, fy, fy_tl)
  end subroutine run_tl

  subroutine run_ad(mode, yh, yc, yfx, yfy, qh_ad, qc_ad)
    integer, intent(in) :: mode
    real, intent(in) :: yh(isd:ied, jsd:jed, npz, nq), yc(is:ie, js:je, npz, nq)
    real, intent(in) :: yfx(is:ie+1, js:je), yfy(is:ie, js:je+1)
    real, intent(out) :: qh_ad(isd:ied, jsd:jed, npz, nq), qc_ad(is:ie, js:je, npz, nq)
    real :: qh(isd:ied, jsd:jed, npz, nq), qc(is:ie, js:je, npz, nq)
    real :: fx(is:ie+1, js:je), fy(is:ie, js:je+1), fx_ad(is:ie+1, js:je), fy_ad(is:ie, js:je+1)
    qh = qh0; qc = qc0; fx = 0.; fy = 0.
    qh_ad = 0.; qc_ad = 0.; fx_ad = 0.; fy_ad = 0.
    if (mode == 1) then
       fx_ad = yfx; fy_ad = yfy
    elseif (mode == 2) then
       qh_ad(is:ie, js:je, :, :) = yh(is:ie, js:je, :, :)
    else
       qc_ad = yc
    endif
    call top_adv_ad(mode, qh, qh_ad, qc, qc_ad, fx, fx_ad, fy, fy_ad)
  end subroutine run_ad

  subroutine report(what, lhs, rhs, tol)
    character(len=*), intent(in) :: what
    real, intent(in) :: lhs, rhs, tol
    character(len=4) :: verdict
    verdict = 'ok'
    if (rel(lhs, rhs) > tol) then
       verdict = 'FAIL'
       nfail = nfail + 1
    endif
    if (root) write(*,'(2x,a,t52,2es24.16,es10.2,2x,a)') what, lhs, rhs, rel(lhs, rhs), verdict
  end subroutine report

  logical function owns(p)
    integer, intent(in) :: p(4)
    owns = p(1) == tile .and. p(2) >= is .and. p(2) <= ie .and. p(3) >= js .and. p(3) <= je
  end function owns

  !> pseudo-random value in (-0.5, 0.5) of a global index, the same for any decomposition
  real function prand(t, i, j, k, seed)
    integer, intent(in) :: t, i, j, k, seed
    real :: a
    a = sin(12.9898*i + 78.233*j + 37.719*k + 4.581*t + 0.917*seed)*43758.5453
    prand = a - floor(a) - 0.5
  end function prand

  subroutine field_h(a, seed)
    real, intent(out) :: a(isd:ied, jsd:jed, npz, nq)
    integer, intent(in) :: seed
    integer :: i, j, k
    do k = 1, npz
       do j = jsd, jed
          do i = isd, ied
             a(i,j,k,1) = prand(tile, i, j, k, seed)
          enddo
       enddo
    enddo
  end subroutine field_h

  subroutine field_c(a, seed)
    real, intent(out) :: a(is:ie, js:je, npz, nq)
    integer, intent(in) :: seed
    integer :: i, j, k
    do k = 1, npz
       do j = js, je
          do i = is, ie
             a(i,j,k,1) = prand(tile, i, j, k, seed)
          enddo
       enddo
    enddo
  end subroutine field_c

  subroutine field_fx(a, seed)
    real, intent(out) :: a(is:ie+1, js:je)
    integer, intent(in) :: seed
    integer :: i, j
    do j = js, je
       do i = is, ie+1
          a(i,j) = prand(tile, i, j, 0, seed)
       enddo
    enddo
  end subroutine field_fx

  subroutine field_fy(a, seed)
    real, intent(out) :: a(is:ie, js:je+1)
    integer, intent(in) :: seed
    integer :: i, j
    do j = js, je+1
       do i = is, ie
          a(i,j) = prand(tile, i, j, 0, seed)
       enddo
    enddo
  end subroutine field_fy

  subroutine zero_halo(a)
    real, intent(inout) :: a(isd:ied, jsd:jed, npz, nq)
    real :: t(is:ie, js:je, npz, nq)
    t = a(is:ie, js:je, :, :)
    a = 0.
    a(is:ie, js:je, :, :) = t
  end subroutine zero_halo

  subroutine one_level(a, k)
    real, intent(inout) :: a(isd:ied, jsd:jed, npz, nq)
    integer, intent(in) :: k
    a(:, :, 1:k-1, :) = 0.; a(:, :, k+1:, :) = 0.
  end subroutine one_level

end program test_adv
