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

!>@brief TL-vs-FD and dot-product tests of the TAF-generated TL/AD (taf/gen) on one face.
!>@details usage: test_adv <face data dir> [hord_tr kord_tr]   (default 2 17, the linear scheme)
!! Operators (top_adv mode):
!!   1  fv_tp_2d on level k_test:  q (with halo) -> fx, fy
!!   2  tracer_2d:                 q (with halo) -> q (interior)
!!   3  full step (top_step):      q -> q, mass fixer with frozen s0 or with s(q)
!! Dot test: <TL x, y> against <x, AD y>, for random x, y and for y = e_k at face corner, edge and
!! interior cells. In mode 1 and 2 x is tested with the halo active (variant a) and with the halo
!! zero (variant b). TL vs FD: centred differences (F(q0 + eps x) - F(q0 - eps x)) / (2 eps).
program test_adv
  use adj_face_mod
  implicit none
  character(len=512) :: dir, arg
  integer, parameter :: nh = (ied-isd+1)*(jed-jsd+1)*npz*nq, nc = n*n*npz*nq
  integer, parameter :: nfx = (n+1)*n, nfy = n*(n+1)
  real :: qh0(isd:ied, jsd:jed, npz, nq), qc0(is:ie, js:je, npz, nq)
  integer :: nfail

  call get_command_argument(1, dir)
  call load_face(trim(dir))
  flagstruct%hord_tr = 2; flagstruct%kord_tr = 17
  if (command_argument_count() >= 3) then
     call get_command_argument(2, arg); read(arg, *) flagstruct%hord_tr
     call get_command_argument(3, arg); read(arg, *) flagstruct%kord_tr
  endif
  write(*,'(a,a)') 'face data: ', trim(dir)
  write(*,'(a,2i4,a,i3)') 'hord_tr, kord_tr =', flagstruct%hord_tr, flagstruct%kord_tr, '   k_test =', k_test

  qh0(:, :, :, 1) = q_in
  qc0(:, :, :, 1) = q_in(is:ie, js:je, :)
  nfail = 0

  call test_mode(1, 'fv_tp_2d (level k_test)', .false.)
  call test_mode(2, 'tracer_2d', .false.)
  ! full step: first a forward with s(q) to store s0, then frozen s0 (linear) and full s(q)
  adj_freeze_s = .false.
  call fwd(3, qh0, qc0)
  call dump_step(trim(dir), 'taf', qc0)
  qc0(:, :, :, 1) = q_in(is:ie, js:je, :)
  write(*,'(a,es24.16)') 'mass-fixer factor s0 (this face) =', adj_s0
  adj_freeze_s = .true.
  call test_mode(3, 'full step, frozen s0', .false.)
  adj_freeze_s = .false.
  call test_mode(3, 'full step, s(q) linearized', .true.)

  write(*,*)
  if (nfail == 0) then
     write(*,'(a)') 'ALL TESTS PASSED'
  else
     write(*,'(a,i3,a)') 'FAILED: ', nfail, ' test(s)'
  endif

contains

  !> forward copies (the top routines overwrite their arguments)
  subroutine fwd(mode, qh, qc, fx, fy)
    integer, intent(in) :: mode
    real, intent(inout) :: qh(isd:ied, jsd:jed, npz, nq), qc(is:ie, js:je, npz, nq)
    real, intent(out), optional :: fx(is:ie+1, js:je), fy(is:ie, js:je+1)
    real :: fxl(is:ie+1, js:je), fyl(is:ie, js:je+1)
    fxl = 0.; fyl = 0.
    call top_adv(mode, qh, qc, fxl, fyl)
    if (present(fx)) fx = fxl
    if (present(fy)) fy = fyl
  end subroutine fwd

  real function rel(a, b)
    real, intent(in) :: a, b
    rel = abs(a - b) / max(abs(a), abs(b), tiny(1.0))
  end function rel

  !> output-space inner product (mode 1: fluxes; 2: q interior; 3: q)
  real function dot_out(mode, qh, qc, fx, fy, yh, yc, yfx, yfy)
    integer, intent(in) :: mode
    real, intent(in) :: qh(isd:ied, jsd:jed, npz, nq), qc(is:ie, js:je, npz, nq)
    real, intent(in) :: fx(is:ie+1, js:je), fy(is:ie, js:je+1)
    real, intent(in) :: yh(isd:ied, jsd:jed, npz, nq), yc(is:ie, js:je, npz, nq)
    real, intent(in) :: yfx(is:ie+1, js:je), yfy(is:ie, js:je+1)
    if (mode == 1) then
       dot_out = sum(fx*yfx) + sum(fy*yfy)
    elseif (mode == 2) then
       dot_out = sum(qh(is:ie, js:je, :, :)*yh(is:ie, js:je, :, :))
    else
       dot_out = sum(qc*yc)
    endif
  end function dot_out

  subroutine test_mode(mode, title, nonlinear)
    integer, intent(in) :: mode
    character(len=*), intent(in) :: title
    logical, intent(in) :: nonlinear
    real :: xh(isd:ied, jsd:jed, npz, nq), xc(is:ie, js:je, npz, nq)
    real :: yh(isd:ied, jsd:jed, npz, nq), yc(is:ie, js:je, npz, nq)
    real :: yfx(is:ie+1, js:je), yfy(is:ie, js:je+1)
    real :: qh(isd:ied, jsd:jed, npz, nq), qc(is:ie, js:je, npz, nq)
    real :: qh_tl(isd:ied, jsd:jed, npz, nq), qc_tl(is:ie, js:je, npz, nq)
    real :: qh_ad(isd:ied, jsd:jed, npz, nq), qc_ad(is:ie, js:je, npz, nq)
    real :: fx(is:ie+1, js:je), fy(is:ie, js:je+1), fx_tl(is:ie+1, js:je), fy_tl(is:ie, js:je+1)
    real :: fx_ad(is:ie+1, js:je), fy_ad(is:ie, js:je+1)
    real :: qhp(isd:ied, jsd:jed, npz, nq), qcp(is:ie, js:je, npz, nq), fxp(is:ie+1, js:je), fyp(is:ie, js:je+1)
    real :: qhm(isd:ied, jsd:jed, npz, nq), qcm(is:ie, js:je, npz, nq), fxm(is:ie+1, js:je), fym(is:ie, js:je+1)
    real :: lhs, rhs, eps, fd, tl, worst
    integer :: ivar, ie_, ip, kk, ks
    integer, parameter :: npt = 7
    integer :: pts(3, npt)
    character(len=12) :: ptname(npt)
    real, parameter :: tol_dot = 1.e-12

    write(*,*)
    write(*,'(a,i2,2a)') '=== mode ', mode, ': ', title
    kk = k_test
    pts(:,1) = (/ 1, 1, kk /);   ptname(1) = 'SW corner'
    pts(:,2) = (/ n, n, kk /);   ptname(2) = 'NE corner'
    pts(:,3) = (/ 1, n/2, kk /); ptname(3) = 'W edge'
    pts(:,4) = (/ n/2, n, kk /); ptname(4) = 'N edge'
    pts(:,5) = (/ 2, 2, kk /);   ptname(5) = '2nd ring'
    pts(:,6) = (/ n/2, n/2, kk /); ptname(6) = 'interior'
    pts(:,7) = (/ n, 1, kk /);   ptname(7) = 'SE corner'

    ! ---- dot tests: random x (variant a: halo active; b: halo zero), random y
    do ivar = 1, merge(2, 1, mode < 3)
       call random_fill(xh, 11 + ivar); call random_fill(xc, 13)
       if (ivar == 2) call zero_halo(xh)
       if (mode == 1) then
          xh(:, :, 1:kk-1, :) = 0.; xh(:, :, kk+1:, :) = 0.
       endif
       call random_fill(yh, 21); call random_fill(yc, 23); call random_fill(yfx, 25); call random_fill(yfy, 27)
       call run_tl(mode, xh, xc, qh_tl, qc_tl, fx_tl, fy_tl)
       lhs = dot_out(mode, qh_tl, qc_tl, fx_tl, fy_tl, yh, yc, yfx, yfy)
       call run_ad(mode, yh, yc, yfx, yfy, qh_ad, qc_ad)
       if (mode < 3) then
          rhs = sum(xh*qh_ad)
       else
          rhs = sum(xc*qc_ad)
       endif
       call report('dot, random x'//merge(' (halo active)', ' (halo zero)  ', ivar == 1)//', random y', lhs, rhs, tol_dot)
    enddo

    ! ---- dot tests with point seeds y = e_k (q output modes), random x with halo active
    if (mode >= 2) then
       call random_fill(xh, 31); call random_fill(xc, 33)
       call run_tl(mode, xh, xc, qh_tl, qc_tl, fx_tl, fy_tl)
       worst = 0.
       do ip = 1, npt
          yh = 0.; yc = 0.; yfx = 0.; yfy = 0.
          yh(pts(1,ip), pts(2,ip), pts(3,ip), 1) = 1.
          yc(pts(1,ip), pts(2,ip), pts(3,ip), 1) = 1.
          lhs = dot_out(mode, qh_tl, qc_tl, fx_tl, fy_tl, yh, yc, yfx, yfy)
          call run_ad(mode, yh, yc, yfx, yfy, qh_ad, qc_ad)
          if (mode == 2) then
             rhs = sum(xh*qh_ad)
          else
             rhs = sum(xc*qc_ad)
          endif
          call report('dot, y = e_k at '//ptname(ip), lhs, rhs, tol_dot)
       enddo
    endif

    ! ---- TL vs centred FD, eps scan, random x (halo active)
    call random_fill(xh, 41); call random_fill(xc, 43)
    if (mode == 1) then
       xh(:, :, 1:kk-1, :) = 0.; xh(:, :, kk+1:, :) = 0.
    endif
    call run_tl(mode, xh, xc, qh_tl, qc_tl, fx_tl, fy_tl)
    call random_fill(yh, 51); call random_fill(yc, 53); call random_fill(yfx, 55); call random_fill(yfy, 57)
    tl = dot_out(mode, qh_tl, qc_tl, fx_tl, fy_tl, yh, yc, yfx, yfy)
    write(*,'(a)') '  TL vs FD  <(F(q0+eps x) - F(q0-eps x))/(2 eps), y>  vs  <TL x, y>'
    do ie_ = 1, 5
       eps = 10.**(-ie_)
       ! x has the size of q0 (random factor), so eps is a relative perturbation
       qhp = qh0 + eps*xh*maxval(abs(q_in)); qcp = qc0 + eps*xc*maxval(abs(q_in))
       qhm = qh0 - eps*xh*maxval(abs(q_in)); qcm = qc0 - eps*xc*maxval(abs(q_in))
       call fwd(mode, qhp, qcp, fxp, fyp); call fwd(mode, qhm, qcm, fxm, fym)
       fd = (dot_out(mode, qhp, qcp, fxp, fyp, yh, yc, yfx, yfy) - &
             dot_out(mode, qhm, qcm, fxm, fym, yh, yc, yfx, yfy)) / (2.*eps*maxval(abs(q_in)))
       write(*,'(a,es8.1,a,es24.16,a,es10.2)') '    eps ', eps, '  FD ', fd, '  rel. diff ', rel(fd, tl)
    enddo
    write(*,'(a,es24.16)') '    TL            ', tl
  end subroutine test_mode

  subroutine run_tl(mode, xh, xc, qh_tl, qc_tl, fx_tl, fy_tl)
    integer, intent(in) :: mode
    real, intent(in) :: xh(isd:ied, jsd:jed, npz, nq), xc(is:ie, js:je, npz, nq)
    real, intent(out) :: qh_tl(isd:ied, jsd:jed, npz, nq), qc_tl(is:ie, js:je, npz, nq)
    real, intent(out) :: fx_tl(is:ie+1, js:je), fy_tl(is:ie, js:je+1)
    real :: qh(isd:ied, jsd:jed, npz, nq), qc(is:ie, js:je, npz, nq)
    real :: fx(is:ie+1, js:je), fy(is:ie, js:je+1)
    logical :: frz
    frz = adj_freeze_s
    qh = qh0; qc = qc0; fx = 0.; fy = 0.
    qh_tl = xh; qc_tl = xc; fx_tl = 0.; fy_tl = 0.
    call top_adv_tl(mode, qh, qh_tl, qc, qc_tl, fx, fx_tl, fy, fy_tl)
    adj_freeze_s = frz
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
    write(*,'(2x,a,t50,2es24.16,es10.2,2x,a)') what, lhs, rhs, rel(lhs, rhs), verdict
  end subroutine report

  subroutine random_fill(a, seed)
    real, intent(out) :: a(..)
    integer, intent(in) :: seed
    integer :: s(64), ns, ip_
    call random_seed(size=ns)
    s(1:ns) = seed + 37*(/ (ns-ip_, ip_ = 1, ns) /)
    call random_seed(put=s(1:ns))
    select rank (a)
    rank (2)
       call random_number(a); a = a - 0.5
    rank (4)
       call random_number(a); a = a - 0.5
    end select
  end subroutine random_fill

  subroutine zero_halo(a)
    real, intent(inout) :: a(isd:ied, jsd:jed, npz, nq)
    real :: t(is:ie, js:je, npz, nq)
    t = a(is:ie, js:je, :, :)
    a = 0.
    a(is:ie, js:je, :, :) = t
  end subroutine zero_halo

end program test_adv
