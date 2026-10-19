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

!>@brief Multi-step test of the full-step TL/AD (top_step) on the six faces: the P2a-type
!! receptor test of the linear scheme, standalone.
!>@details usage: test_chain <chain dir> <nsteps> <receptors.txt> [layout [hord kord]]
!!   <chain dir>/sNN/f1..f6 from prep/make_face_inputs.py (step NN of the trajectory).
!!   receptors.txt as in baselines/P2a: "name F,j,i L" (face, 0-based j and i, GEOS-Chem level L,
!!   surface first; FV3 level k = 73 - L).
!! J_r = q at receptor r after nsteps steps, starting from q_in of step 1. x = pseudo-random field
!! on the receptor level (all columns, all faces), as the rnd perturbations of run_P2lin.pbs.
!! For each receptor: FD = (J_r(q0 + eps x) - J_r(q0 - eps x)) / (2 eps), <lambda, x> with
!! lambda = AD e_r (the adjoint run backwards through all steps), and <TL x, e_r>.
!! Mass fixer: frozen s0 (per step, from the base trajectory; then the step operator is linear)
!! and s(q) linearized around the base trajectory.
program test_chain
  use adj_face_mod
  use adj_load_mod,    only: load_step
  use adj_backend_mod, only: backend_init, backend_finalize, backend_rank, backend_allsum, backend_allmax
  implicit none
  character(len=512) :: dir, arg, rfile
  integer :: nsteps, layout, n, ir, nrec, u, ios, f, jj, ii, lev, nfail
  integer, parameter :: maxrec = 32
  character(len=32) :: rname(maxrec)
  integer :: rcell(4, maxrec)                      !< tile, i, j, k (FV3 indices)
  real, allocatable :: qtraj(:,:,:,:,:), s0(:), x(:,:,:,:), tlx(:,:,:,:), lam(:,:,:,:), qp(:,:,:,:), qm(:,:,:,:)
  real :: scale, eps, fd, ad, tl, jp, jm, xnorm, lnorm, err
  logical :: root, frozen
  character(len=256) :: line
  integer :: mloop

  call get_command_argument(1, dir)
  call get_command_argument(2, arg); read(arg, *) nsteps
  call get_command_argument(3, rfile)
  layout = 1
  if (command_argument_count() >= 4) then
     call get_command_argument(4, arg); read(arg, *) layout
  endif
  call backend_init(layout)
  root = backend_rank() == 0
  call load_step(step_dir(1), quiet=.false.)
  flagstruct%hord_tr = 2; flagstruct%kord_tr = 17
  if (command_argument_count() >= 6) then
     call get_command_argument(5, arg); read(arg, *) flagstruct%hord_tr
     call get_command_argument(6, arg); read(arg, *) flagstruct%kord_tr
  endif

  ! receptors
  nrec = 0
  open(newunit=u, file=trim(rfile), status='old', action='read')
  do
     read(u, '(a)', iostat=ios) line
     if (ios /= 0) exit
     line = adjustl(line)
     if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
     nrec = nrec + 1
     ii = index(line, ' ')
     rname(nrec) = line(1:ii-1)
     line = adjustl(line(ii:))
     do jj = 1, len_trim(line)
        if (line(jj:jj) == ',') line(jj:jj) = ' '
     enddo
     read(line, *) f, jj, ii, lev
     rcell(:, nrec) = (/ f, ii + 1, jj + 1, 73 - lev /)
  enddo
  close(u)

  allocate(qtraj(is:ie, js:je, npz, nq, 0:nsteps), s0(nsteps))
  allocate(x(is:ie, js:je, npz, nq), tlx(is:ie, js:je, npz, nq), lam(is:ie, js:je, npz, nq))
  allocate(qp(is:ie, js:je, npz, nq), qm(is:ie, js:je, npz, nq))
  if (root) then
     write(*,'(a,a,a,i3,a,2i4)') 'chain ', trim(dir), ', steps', nsteps, ', hord_tr kord_tr', &
          flagstruct%hord_tr, flagstruct%kord_tr
     write(*,'(a,i3,a,i3)') 'receptors:', nrec, ', level k =', rcell(4, 1)
  endif

  ! base trajectory with s(q); s0(n) for the frozen operator
  qtraj(:, :, :, 1, 0) = q_in(is:ie, js:je, :)
  adj_freeze_s = .false.
  do n = 1, nsteps
     call load_step(step_dir(n), quiet=.true.)
     qtraj(:, :, :, :, n) = qtraj(:, :, :, :, n-1)
     call step_fwd(qtraj(:, :, :, :, n))
     s0(n) = adj_s_last
  enddo
  if (root) write(*,'(a,2es24.16)') 'mass-fixer factor s0: min, max over the steps', minval(s0), maxval(s0)

  call field_level(x, rcell(4, 1))
  scale = backend_allmax(maxval(abs(qtraj(:, :, :, :, 0))))
  eps = 1.e-2
  nfail = 0
  do mloop = 1, 2
     frozen = mloop == 1
     adj_freeze_s = frozen
     if (root) then
        write(*,*)
        if (frozen) then
           write(*,'(a)') '=== mass fixer frozen (s0 per step): linear operator'
        else
           write(*,'(a)') '=== mass fixer s(q), linearized around the base trajectory'
        endif
        write(*,'(a)') 'receptor       tile   i   j   k   FD (J+ - J-)/2eps        <lambda, x>              ' // &
                       '<TL x, e_r>          adj/FD - 1   adj vs TL   |adj-TL|/(|x||lambda|)'
     endif
     ! TL through all steps
     tlx = x
     do n = 1, nsteps
        call load_step(step_dir(n), quiet=.true.)
        if (frozen) adj_s0 = s0(n)
        call step_tl(qtraj(:, :, :, :, n-1), tlx)
     enddo
     ! FD: forward chains from q0 +- eps*scale*x
     qp = qtraj(:, :, :, :, 0) + eps*scale*x
     qm = qtraj(:, :, :, :, 0) - eps*scale*x
     do n = 1, nsteps
        call load_step(step_dir(n), quiet=.true.)
        if (frozen) adj_s0 = s0(n)
        call step_fwd(qp)
        if (frozen) adj_s0 = s0(n)
        call step_fwd(qm)
     enddo
     do ir = 1, nrec
        lam = 0.
        if (owns(rcell(:, ir))) lam(rcell(2, ir), rcell(3, ir), rcell(4, ir), 1) = 1.
        do n = nsteps, 1, -1
           call load_step(step_dir(n), quiet=.true.)
           if (frozen) adj_s0 = s0(n)
           call step_ad(qtraj(:, :, :, :, n-1), lam)
        enddo
        ad = backend_allsum(sum(lam*x))
        tl = backend_allsum(cell_value(tlx, rcell(:, ir)))
        jp = backend_allsum(cell_value(qp, rcell(:, ir)))
        jm = backend_allsum(cell_value(qm, rcell(:, ir)))
        fd = (jp - jm) / (2.*eps*scale)
        ! dot-test error on the scale of the inner product (Cauchy-Schwarz bound |<x, lambda>| <= |x||lambda|):
        ! relative to |J| it is not meaningful when J is a small difference of large terms
        xnorm = sqrt(backend_allsum(sum(x*x)))
        lnorm = sqrt(backend_allsum(sum(lam*lam)))
        err = abs(ad - tl) / max(xnorm*lnorm, tiny(1.))
        if (err > 1.e-13) nfail = nfail + 1
        if (root) write(*,'(a14,4i4,3es24.16,3es12.2)') rname(ir), rcell(:, ir), fd, ad, tl, &
             ad/fd - 1., abs(ad - tl)/max(abs(ad), abs(tl)), err
     enddo
  enddo
  if (root) then
     write(*,*)
     if (nfail == 0) then
        write(*,'(a)') 'ALL TESTS PASSED (|<TL x, e_r> - <x, AD e_r>| / (|x| |AD e_r|) < 1e-13 for every receptor)'
     else
        write(*,'(a,i3,a)') 'FAILED: ', nfail, ' receptor(s)'
     endif
  endif
  call backend_finalize()

contains

  function step_dir(n) result(d)
    integer, intent(in) :: n
    character(len=600) :: d
    write(d, '(2a,i2.2)') trim(dir), '/s', n
  end function step_dir

  subroutine step_fwd(q)
    real, intent(inout) :: q(is:ie, js:je, npz, nq)
    real :: qh(isd:ied, jsd:jed, npz, nq), fx(is:ie+1, js:je), fy(is:ie, js:je+1)
    qh = 0.; fx = 0.; fy = 0.
    call top_adv(3, qh, q, fx, fy)
  end subroutine step_fwd

  !> TL of one step at the trajectory state q0: q_tl <- M'(q0) q_tl
  subroutine step_tl(q0, q_tl)
    real, intent(in) :: q0(is:ie, js:je, npz, nq)
    real, intent(inout) :: q_tl(is:ie, js:je, npz, nq)
    real :: q(is:ie, js:je, npz, nq), qh(isd:ied, jsd:jed, npz, nq), qh_tl(isd:ied, jsd:jed, npz, nq)
    real :: fx(is:ie+1, js:je), fy(is:ie, js:je+1), fx_tl(is:ie+1, js:je), fy_tl(is:ie, js:je+1)
    q = q0; qh = 0.; qh_tl = 0.; fx = 0.; fy = 0.; fx_tl = 0.; fy_tl = 0.
    call top_adv_tl(3, qh, qh_tl, q, q_tl, fx, fx_tl, fy, fy_tl)
  end subroutine step_tl

  !> AD of one step at the trajectory state q0: lam <- M'(q0)^T lam
  subroutine step_ad(q0, lam)
    real, intent(in) :: q0(is:ie, js:je, npz, nq)
    real, intent(inout) :: lam(is:ie, js:je, npz, nq)
    real :: q(is:ie, js:je, npz, nq), qh(isd:ied, jsd:jed, npz, nq), qh_ad(isd:ied, jsd:jed, npz, nq)
    real :: fx(is:ie+1, js:je), fy(is:ie, js:je+1), fx_ad(is:ie+1, js:je), fy_ad(is:ie, js:je+1)
    q = q0; qh = 0.; qh_ad = 0.; fx = 0.; fy = 0.; fx_ad = 0.; fy_ad = 0.
    call top_adv_ad(3, qh, qh_ad, q, lam, fx, fx_ad, fy, fy_ad)
  end subroutine step_ad

  logical function owns(p)
    integer, intent(in) :: p(4)
    owns = p(1) == tile .and. p(2) >= is .and. p(2) <= ie .and. p(3) >= js .and. p(3) <= je
  end function owns

  real function cell_value(a, p)
    real, intent(in) :: a(is:ie, js:je, npz, nq)
    integer, intent(in) :: p(4)
    cell_value = 0.
    if (owns(p)) cell_value = a(p(2), p(3), p(4), 1)
  end function cell_value

  !> pseudo-random field in (-0.5, 0.5) on level k (a function of the global index), zero elsewhere
  subroutine field_level(a, k)
    real, intent(out) :: a(is:ie, js:je, npz, nq)
    integer, intent(in) :: k
    integer :: i, j
    real :: v
    a = 0.
    do j = js, je
       do i = is, ie
          v = sin(12.9898*i + 78.233*j + 37.719*k + 4.581*tile + 0.917*61)*43758.5453
          a(i, j, k, 1) = v - floor(v) - 0.5
       enddo
    enddo
  end subroutine field_level

end program test_chain
