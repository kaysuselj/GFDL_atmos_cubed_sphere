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

!>@brief Test of the GCHP adapter (gchp/adv_adj_discrete.F90) in a GCHP-like build: the real
!! fv_arrays.F90 and fv_mp_mod.F90, and the generated code with all its modules renamed (adj_*,
!! stb_*), as it is linked into GCHP.
!>@details usage: mpirun -n 6*layout**2 test_gchp_adapter <data dir (f1..f6)> [layout]
!! The step is loaded by the standalone loader; a real fv_grid_type / fv_flags_type is filled from
!! it as GCHP has them, and the adapter must give the same forward step, frozen-s0 adjoint and
!! linearized adjoint as the direct calls of the generated code.
program test_gchp_adapter
  use fv_arrays_mod,        only: fv_grid_type, fv_flags_type
  use adj_face_mod,         only: gs_stub => gridstruct, is, ie, js, je, isd, ied, jsd, jed, tile, &
                                  npx, npy, npz, nq, dt, ak, bk, ptop, q_in, cx_e, cy_e, mfx_e, mfy_e, &
                                  ple0, ple1, adj_freeze_s, adj_s0, adj_s_last, flagstruct_stub => flagstruct, &
                                  pleadv_gchp
  use adj_load_mod,         only: load_step
  use adj_backend_mod,      only: backend_init, backend_finalize, backend_rank, backend_allsum, backend_allmax, Atm
  use adv_adj_discrete_mod, only: adv_adj_discrete_init, adv_adj_discrete_fwd, adv_adj_discrete_adj, &
                                  adv_adj_discrete_conv_check
  implicit none
  character(len=512) :: dir, arg
  integer :: layout, i, j, k, nfail
  type(fv_grid_type), target :: gs
  type(fv_flags_type) :: fl
  integer, target :: grid_type_t
  logical, target :: nested_t
  real, allocatable :: cxi(:,:,:), cyi(:,:,:), mfxi(:,:,:), mfyi(:,:,:), p0(:,:,:), p1(:,:,:)
  real, allocatable :: q0(:,:,:), qref(:,:,:), qadp(:,:,:), y(:,:,:), lref(:,:,:), ladp(:,:,:)
  real, allocatable :: qh(:,:,:,:), qh_ad(:,:,:,:), lam4(:,:,:,:), q4(:,:,:,:)
  real, allocatable :: fx(:,:), fy(:,:), fx_ad(:,:), fy_ad(:,:), akl(:), bkl(:), pref(:,:,:)
  real :: sref, sadp, a, e1, e2

  call get_command_argument(1, dir)
  layout = 1
  if (command_argument_count() >= 2) then
     call get_command_argument(2, arg); read(arg, *) layout
  endif
  call backend_init(layout)
  call load_step(trim(dir), quiet=.true.)
  flagstruct_stub%hord_tr = 2; flagstruct_stub%kord_tr = 17
  nfail = 0

  ! forward-step inputs as GCHP gives them to offline_tracer_advection (interior arrays)
  allocate(cxi(is:ie,js:je,npz), cyi(is:ie,js:je,npz), mfxi(is:ie,js:je,npz), mfyi(is:ie,js:je,npz))
  allocate(p0(is:ie,js:je,npz+1), p1(is:ie,js:je,npz+1), akl(npz+1), bkl(npz+1))
  cxi = cx_e(is:ie, js:je, :); cyi = cy_e(is:ie, js:je, :)
  mfxi = mfx_e(is:ie, js:je, :); mfyi = mfy_e(is:ie, js:je, :)
  p0 = ple0; p1 = ple1; akl = ak; bkl = bk
  allocate(pref(is:ie,js:je,npz+1)); pref = pleadv_gchp      ! GCHP's forward PLEadv (r4)
  allocate(q0(is:ie,js:je,npz), qref(is:ie,js:je,npz), qadp(is:ie,js:je,npz), y(is:ie,js:je,npz))
  allocate(lref(is:ie,js:je,npz), ladp(is:ie,js:je,npz))
  allocate(qh(isd:ied,jsd:jed,npz,1), qh_ad(isd:ied,jsd:jed,npz,1), lam4(is:ie,js:je,npz,1), q4(is:ie,js:je,npz,1))
  allocate(fx(is:ie+1,js:je), fy(is:ie,js:je+1), fx_ad(is:ie+1,js:je), fy_ad(is:ie,js:je+1))
  q0 = q_in(is:ie, js:je, :)
  do k = 1, npz
     do j = js, je
        do i = is, ie
           a = sin(12.9898*i + 78.233*j + 37.719*k + 4.581*tile + 7.7)*43758.5453
           y(i,j,k) = a - floor(a) - 0.5
        enddo
     enddo
  enddo

  ! references: direct calls of the generated code on the loader's data
  q4(:,:,:,1) = q0; qh = 0.; fx = 0.; fy = 0.
  adj_freeze_s = .false.
  call top_adv(3, qh, q4, fx, fy)
  qref = q4(:,:,:,1); sref = adj_s_last
  call ad_ref(.true., sref, lref)

  ! a real gridstruct / flagstruct as GCHP has them
  allocate(gs%area(isd:ied,jsd:jed), gs%rarea(isd:ied,jsd:jed), gs%dxa(isd:ied,jsd:jed), gs%dya(isd:ied,jsd:jed))
  allocate(gs%dx(isd:ied,jsd:jed+1), gs%dy(isd:ied+1,jsd:jed), gs%sin_sg(isd:ied,jsd:jed,9), gs%area_64(isd:ied,jsd:jed))
  gs%area = gs_stub%area; gs%rarea = gs_stub%rarea; gs%dxa = gs_stub%dxa; gs%dya = gs_stub%dya
  gs%dx = gs_stub%dx; gs%dy = gs_stub%dy; gs%sin_sg = gs_stub%sin_sg; gs%area_64 = gs_stub%area_64
  gs%da_min = gs_stub%da_min
  gs%sw_corner = gs_stub%sw_corner; gs%se_corner = gs_stub%se_corner
  gs%ne_corner = gs_stub%ne_corner; gs%nw_corner = gs_stub%nw_corner
  grid_type_t = 0; nested_t = .false.
  gs%grid_type => grid_type_t; gs%nested => nested_t
  fl%hord_tr = 2; fl%kord_tr = 17; fl%nord_tr = 0; fl%q_split = 0; fl%trdm2 = 0.; fl%lim_fac = 1.
  fl%z_tracer = .false.; fl%fill = .false.

  call adv_adj_discrete_init(gs, fl, Atm%bd, Atm%domain, Atm%tile, npx, npy, npz, akl, bkl, ptop, dt)

  ! forward step through the adapter
  qadp = q0
  call adv_adj_discrete_fwd(qadp, mfxi, mfyi, cxi, cyi, p0, p1, 'full', s_out=sadp)
  call check('forward step (s(q))', qadp, qref)
  if (backend_rank() == 0) write(*,'(a,2es24.16)') '  mass-fixer factor: direct, adapter', sref, sadp
  ! adjoint, frozen s0
  ladp = y
  call adv_adj_discrete_adj(ladp, mfxi, mfyi, cxi, cyi, p0, p1, 'frozen', s0=sref)
  call check('adjoint, frozen s0', ladp, lref)
  ! adjoint, linearized s(q)
  call ad_ref(.false., sref, lref)
  ladp = y
  call adv_adj_discrete_adj(ladp, mfxi, mfyi, cxi, cyi, p0, p1, 'linearized', q_traj=q0)
  call check('adjoint, linearized s(q)', ladp, lref)

  ! convention check (as the AdvCore patch runs it with ADV_ADJ_DISCRETE_CHECK)
  call adv_adj_discrete_conv_check(mfxi, mfyi, cxi, cyi, p0, p1, pref, e1, e2)
  if (backend_rank() == 0) write(*,'(a,es10.2,a,es10.2,a)') '  convention check: max |pleadv - PLEadv_gchp| =', e1, &
       ' Pa as given,', e2, ' Pa with the winds flipped'
  if (e1 > 0.1 .or. e2 < 10.*e1) nfail = nfail + 1

  if (backend_rank() == 0) then
     if (nfail == 0) then
        write(*,'(a)') 'ALL TESTS PASSED'
     else
        write(*,'(a,i3,a)') 'FAILED: ', nfail, ' test(s)'
     endif
  endif
  call backend_finalize()

contains

  subroutine ad_ref(frozen, s, lam)
    logical, intent(in) :: frozen
    real, intent(in) :: s
    real, intent(out) :: lam(is:ie,js:je,npz)
    adj_freeze_s = frozen; adj_s0 = s
    q4(:,:,:,1) = q0; qh = 0.; qh_ad = 0.; fx = 0.; fy = 0.; fx_ad = 0.; fy_ad = 0.
    lam4(:,:,:,1) = y
    call top_adv_ad(3, qh, qh_ad, q4, lam4, fx, fx_ad, fy, fy_ad)
    lam = lam4(:,:,:,1)
  end subroutine ad_ref

  subroutine check(what, a, b)
    character(len=*), intent(in) :: what
    real, intent(in) :: a(:,:,:), b(:,:,:)
    real :: d, m
    d = backend_allmax(maxval(abs(a - b)))
    m = backend_allmax(maxval(abs(b)))
    if (d > 1.e-15*m) nfail = nfail + 1
    if (backend_rank() == 0) write(*,'(2x,a,t36,a,es10.2,a,l2)') what, 'adapter vs direct: max |diff| / max =', &
         d/max(m, tiny(1.)), '   bit-identical:', d == 0.
  end subroutine check

end program test_gchp_adapter
