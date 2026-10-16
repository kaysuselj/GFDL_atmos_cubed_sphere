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

!>@brief Single-face stand-ins for FMS and the FV3 modules that tp_core.F90, fv_tracer2d.F90,
!! fv_mapz.F90 and fv_fill.F90 use, so that these forward files compile unchanged for the
!! standalone TL/AD tests (model/adjoint). One process, one cube face, no MPI.
!>@details Halo exchanges of q are no-ops: the q halo is either an input of the test operator
!! (tracer_2d level) or zero (offline_tracer_advection zeroes it). The C-grid winds are not active
!! variables; their "exchange" copies the halo values that the driver loaded from the global field
!! (see mpp_domains_mod). mpp_sum is the identity on one face, or returns frozen sums (mass fixer
!! with a frozen scaling factor, see mpp_mod).

module platform_mod
  implicit none
  integer, parameter :: r8_kind = 8, r4_kind = 4, i4_kind = 4, i8_kind = 8
end module platform_mod

module constants_mod
  implicit none
  real, parameter :: radius = 6371.0e3
  real(kind=8), parameter :: pi_8 = 3.14159265358979323846d0
  real, parameter :: pi = pi_8
  real, parameter :: grav = 9.80665, rdgas = 287.04, rvgas = 461.6
  real, parameter :: cp_air = 1004.64, cp_vapor = 4.0*rvgas
  real, parameter :: hlv = 2.501e6, hlf = 3.50e5, hls = hlv + hlf
end module constants_mod

module field_manager_mod
  implicit none
  integer, parameter :: MODEL_ATMOS = 1
end module field_manager_mod

module tracer_manager_mod
  implicit none
contains
  integer function get_tracer_index(model, name)
    integer, intent(in) :: model
    character(len=*), intent(in) :: name
    get_tracer_index = -1
  end function get_tracer_index
end module tracer_manager_mod

module fv_timing_mod
  implicit none
contains
  subroutine timing_on(blk)
    character(len=*), intent(in) :: blk
  end subroutine timing_on
  subroutine timing_off(blk)
    character(len=*), intent(in) :: blk
  end subroutine timing_off
end module fv_timing_mod

module mpp_mod
  implicit none
  integer, parameter :: FATAL = 2, WARNING = 1, NOTE = 0
  !> Mass-fixer control: adj_sum_mode = 0: mpp_sum is the identity (one face);
  !! 1: also record the sums; 2: return the recorded sums (frozen scaling factor).
  integer :: adj_sum_mode = 0
  integer :: adj_sum_n = 0
  real :: adj_sum_saved(8) = 0.
contains
  subroutine mpp_error(level, msg)
    integer, intent(in) :: level
    character(len=*), intent(in) :: msg
    write(*,*) 'mpp_error: ', trim(msg)
    if (level == FATAL) error stop 'mpp_error FATAL'
  end subroutine mpp_error
  subroutine mpp_sum(a, n)
    integer, intent(in) :: n
    real, intent(inout) :: a(n)
    integer :: i
    if (adj_sum_mode == 1) then
       adj_sum_n = n
       do i = 1, n
          adj_sum_saved(i) = a(i)
       enddo
    elseif (adj_sum_mode == 2) then
       do i = 1, n
          a(i) = adj_sum_saved(i)
       enddo
    endif
  end subroutine mpp_sum
  subroutine mpp_max(a)
    real, intent(inout) :: a
  end subroutine mpp_max
  subroutine mpp_broadcast(a, n, pe)
    integer, intent(in) :: n, pe
    real, intent(inout) :: a(n)
  end subroutine mpp_broadcast
  subroutine mpp_send(a, n, pe)
    integer, intent(in) :: n, pe
    real, intent(in) :: a(n)
  end subroutine mpp_send
  subroutine mpp_recv(a, n, pe)
    integer, intent(in) :: n, pe
    real, intent(inout) :: a(n)
  end subroutine mpp_recv
  integer function get_unit()
    get_unit = 77
  end function get_unit
  integer function mpp_pe()
    mpp_pe = 0
  end function mpp_pe
  integer function mpp_root_pe()
    mpp_root_pe = 0
  end function mpp_root_pe
end module mpp_mod

module mpp_domains_mod
  implicit none
  private
  public :: domain2d, CGRID_NE, mpp_update_domains, mpp_get_boundary, adj_set_cgrid_halo

  integer, parameter :: CGRID_NE = 41

  type domain2d
     integer :: id = 0
  end type domain2d

  !> C-grid wind halos (passive), loaded by the driver:
  !! x-face fields on (is:ie+1, jsd:jed, npz), y-face fields on (isd:ied, js:je+1, npz).
  integer, save :: is, ie, js, je, isd, ied, jsd, jed, npz
  real, allocatable, save :: hcx(:,:,:), hcy(:,:,:), hmfx(:,:,:), hmfy(:,:,:)

  interface mpp_update_domains
     module procedure update_2d, update_3d, update_vec_3d
  end interface mpp_update_domains

contains

  subroutine adj_set_cgrid_halo(is_, ie_, js_, je_, ng, npz_, cx, cy, mfx, mfy)
    integer, intent(in) :: is_, ie_, js_, je_, ng, npz_
    real, intent(in) :: cx(is_:ie_+1, js_-ng:je_+ng, npz_), mfx(is_:ie_+1, js_-ng:je_+ng, npz_)
    real, intent(in) :: cy(is_-ng:ie_+ng, js_:je_+1, npz_), mfy(is_-ng:ie_+ng, js_:je_+1, npz_)
    is = is_; ie = ie_; js = js_; je = je_; npz = npz_
    isd = is - ng; ied = ie + ng; jsd = js - ng; jed = je + ng
    if (allocated(hcx)) deallocate(hcx, hcy, hmfx, hmfy)
    allocate(hcx(is:ie+1, jsd:jed, npz), hmfx(is:ie+1, jsd:jed, npz))
    allocate(hcy(isd:ied, js:je+1, npz), hmfy(isd:ied, js:je+1, npz))
    hcx = cx; hcy = cy; hmfx = mfx; hmfy = mfy
  end subroutine adj_set_cgrid_halo

  !> 1: the interior of x (is:ie, js:je) is the stored CX, 2: MFX, 0: neither
  integer function which_field(x)
    real, intent(in) :: x(isd:ied+1, jsd:jed, npz)
    which_field = 0
    if (all(x(is:ie, js:je, :) == hcx(is:ie, js:je, :))) then
       which_field = 1
    elseif (all(x(is:ie, js:je, :) == hmfx(is:ie, js:je, :))) then
       which_field = 2
    endif
  end function which_field

  subroutine mpp_get_boundary(x, y, domain, wbufferx, ebufferx, sbuffery, nbuffery, gridtype)
    real, intent(in) :: x(isd:ied+1, jsd:jed, npz), y(isd:ied, jsd:jed+1, npz)
    type(domain2d), intent(in) :: domain
    real, intent(inout) :: wbufferx(js:je, npz), ebufferx(js:je, npz)
    real, intent(inout) :: sbuffery(is:ie, npz), nbuffery(is:ie, npz)
    integer, intent(in) :: gridtype
    select case (which_field(x))
    case (1)
       wbufferx = hcx(is, js:je, :);    ebufferx = hcx(ie+1, js:je, :)
       sbuffery = hcy(is:ie, js, :);    nbuffery = hcy(is:ie, je+1, :)
    case (2)
       wbufferx = hmfx(is, js:je, :);   ebufferx = hmfx(ie+1, js:je, :)
       sbuffery = hmfy(is:ie, js, :);   nbuffery = hmfy(is:ie, je+1, :)
    case default
       error stop 'mpp_get_boundary stub: unknown field'
    end select
  end subroutine mpp_get_boundary

  subroutine update_vec_3d(x, y, domain, gridtype, complete)
    real, intent(inout) :: x(isd:ied+1, jsd:jed, npz), y(isd:ied, jsd:jed+1, npz)
    type(domain2d), intent(inout) :: domain
    integer, intent(in), optional :: gridtype
    logical, intent(in), optional :: complete
    integer :: i, j, k
    if (which_field(x) /= 1) error stop 'mpp_update_domains stub: only CX/CY are exchanged'
    do k = 1, npz
       do j = jsd, jed
          do i = is, ie+1
             if (j < js .or. j > je) x(i,j,k) = hcx(i,j,k)
          enddo
       enddo
       do j = js, je+1
          do i = isd, ied
             if (i < is .or. i > ie) y(i,j,k) = hcy(i,j,k)
          enddo
       enddo
    enddo
  end subroutine update_vec_3d

  subroutine update_3d(q, domain, whalo, ehalo, shalo, nhalo, complete)
    real, intent(inout) :: q(:,:,:)
    type(domain2d), intent(inout) :: domain
    integer, intent(in), optional :: whalo, ehalo, shalo, nhalo
    logical, intent(in), optional :: complete
  end subroutine update_3d

  subroutine update_2d(q, domain, whalo, ehalo, shalo, nhalo, complete)
    real, intent(inout) :: q(:,:)
    type(domain2d), intent(inout) :: domain
    integer, intent(in), optional :: whalo, ehalo, shalo, nhalo
    logical, intent(in), optional :: complete
  end subroutine update_2d

end module mpp_domains_mod

module fv_mp_mod
  use mpp_domains_mod, only: domain2d
  implicit none
  private
  public :: ng, is_master, mp_reduce_max, mp_gather, group_halo_update_type, &
            start_group_halo_update, complete_group_halo_update, adj_qhalo_fill, adj_qhalo

  integer, parameter :: ng = 3

  !> If adj_qhalo_fill, start_group_halo_update copies the q halo from adj_qhalo (passive data,
  !! e.g. the neighbouring faces' q for a forward check against GCHP); otherwise it is a no-op.
  !! adj_qhalo has the shape of the exchanged array (lower bounds 1, halo width ng).
  logical, save :: adj_qhalo_fill = .false.
  real, allocatable, save :: adj_qhalo(:,:,:,:)

  type group_halo_update_type
     integer :: id = 0
  end type group_halo_update_type

  interface start_group_halo_update
     module procedure start_halo_3d, start_halo_4d
  end interface start_group_halo_update

contains

  logical function is_master()
    is_master = .true.
  end function is_master

  subroutine mp_reduce_max(a, n)
    integer, intent(in) :: n
    real, intent(inout) :: a(n)
  end subroutine mp_reduce_max

  subroutine mp_gather(a, i1, i2, j1, j2, idim, jdim)
    integer, intent(in) :: i1, i2, j1, j2, idim, jdim
    real, intent(inout) :: a(idim, jdim)
  end subroutine mp_gather

  subroutine start_halo_3d(pack, q, domain)
    type(group_halo_update_type), intent(inout) :: pack
    real, intent(inout) :: q(:,:,:)
    type(domain2d), intent(inout) :: domain
  end subroutine start_halo_3d

  subroutine start_halo_4d(pack, q, domain)
    type(group_halo_update_type), intent(inout) :: pack
    real, intent(inout) :: q(:,:,:,:)
    type(domain2d), intent(inout) :: domain
    integer :: i, j, k, l, ni, nj
    if (.not. adj_qhalo_fill) return
    ni = size(q, 1); nj = size(q, 2)
    do l = 1, size(q, 4)
       do k = 1, size(q, 3)
          do j = 1, nj
             do i = 1, ni
                if (i <= ng .or. i > ni-ng .or. j <= ng .or. j > nj-ng) q(i,j,k,l) = adj_qhalo(i,j,k,l)
             enddo
          enddo
       enddo
    enddo
  end subroutine start_halo_4d

  subroutine complete_group_halo_update(pack, domain)
    type(group_halo_update_type), intent(inout) :: pack
    type(domain2d), intent(inout) :: domain
  end subroutine complete_group_halo_update

end module fv_mp_mod

module fv_arrays_mod
  implicit none
  private
  public :: R_GRID, fv_grid_bounds_type, fv_grid_type, fv_flags_type, fv_nest_type, &
            fv_atmos_type, fv_nest_BC_type_3D

  integer, parameter :: R_GRID = 8

  type fv_grid_bounds_type
     integer :: is, ie, js, je
     integer :: isd, ied, jsd, jed
  end type fv_grid_bounds_type

  type fv_grid_type
     real, allocatable, dimension(:,:) :: area, rarea, dxa, dya, dx, dy, rdxc, rdyc
     real, allocatable, dimension(:,:) :: del6_u, del6_v, cosa_s, rsin2
     real, allocatable, dimension(:,:,:) :: sin_sg
     real(kind=R_GRID), allocatable, dimension(:,:) :: area_64
     real(kind=R_GRID) :: da_min = 0.
     logical :: nested = .false.
     integer :: grid_type = 0
     logical :: sw_corner = .true., se_corner = .true., ne_corner = .true., nw_corner = .true.
  end type fv_grid_type

  type fv_flags_type
     integer :: hord_tr = 2, kord_tr = 17, nord_tr = 0, q_split = 0
     real :: trdm2 = 0., lim_fac = 1.
     logical :: fill = .false., z_tracer = .false., fv_debug = .false.
  end type fv_flags_type

  type fv_nest_BC_type_3D
     real, allocatable, dimension(:,:,:) :: west_t1, east_t1, south_t1, north_t1
  end type fv_nest_BC_type_3D

  type fv_atmos_type
     type(fv_flags_type) :: flagstruct
  end type fv_atmos_type

  type fv_nest_type
     integer :: nestbctype = 1, tracer_nest_timestep = 0
     type(fv_nest_BC_type_3D), allocatable :: q_BC(:)
  end type fv_nest_type

end module fv_arrays_mod

module boundary_mod
  use fv_arrays_mod, only: fv_grid_bounds_type, fv_nest_BC_type_3D
  implicit none
contains
  subroutine nested_grid_BC_apply_intT(var, istag, jstag, npx, npy, npz, bd, step, split, BC, bctype)
    type(fv_grid_bounds_type), intent(in) :: bd
    integer, intent(in) :: istag, jstag, npx, npy, npz
    real, intent(inout) :: var(bd%isd:bd%ied+istag, bd%jsd:bd%jed+jstag, npz)
    real, intent(in) :: step, split
    type(fv_nest_BC_type_3D), intent(in) :: BC
    integer, intent(in), optional :: bctype
    error stop 'nested grids are not supported in the single-face build'
  end subroutine nested_grid_BC_apply_intT
end module boundary_mod

module fv_grid_utils_mod
  use mpp_domains_mod, only: domain2d
  implicit none
  real, parameter :: big_number = 1.d8
  real, parameter :: ptop_min = 1.d-8
contains
  real function g_sum(domain, p, ifirst, ilast, jfirst, jlast, ngc, area, mode, reproduce)
    type(domain2d), intent(in) :: domain
    integer, intent(in) :: ifirst, ilast, jfirst, jlast, ngc, mode
    real, intent(in) :: p(ifirst:ilast, jfirst:jlast)
    real(kind=8), intent(in) :: area(ifirst-ngc:ilast+ngc, jfirst-ngc:jlast+ngc)
    logical, intent(in), optional :: reproduce
    g_sum = sum(p*area(ifirst:ilast, jfirst:jlast))
  end function g_sum
end module fv_grid_utils_mod

module fv_cmp_mod
  implicit none
contains
  subroutine qs_init(kmp)
    integer, intent(in) :: kmp
  end subroutine qs_init
  subroutine fv_sat_adj(mdt, zvir, is, ie, js, je, ng, hydrostatic, consv_te, te0, qv, ql, qi, qr, qs, qg, &
                        hs, dpln, delz, pt, dp, q_con, cappa, area, dtdt, out_dt, last_step, do_qa, qa)
    integer, intent(in) :: is, ie, js, je, ng
    real, intent(in) :: mdt, zvir
    logical, intent(in) :: hydrostatic, consv_te, out_dt, last_step, do_qa
    real :: te0(*), qv(*), ql(*), qi(*), qr(*), qs(*), qg(*), hs(*), dpln(*)
    real :: delz(*), pt(*), dp(*), q_con(*), cappa(*), dtdt(*), qa(*)
    real(kind=8) :: area(*)
    error stop 'fv_sat_adj is not supported in the single-face build'
  end subroutine fv_sat_adj
end module fv_cmp_mod
