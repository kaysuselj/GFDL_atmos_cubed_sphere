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

!>@brief The module 'adj_face_mod' holds the passive inputs of one advection step on one cube face
!! (grid metrics, C-grid winds with halo, pressure edges, ak/bk), read from the files written by
!! prep/make_face_inputs.py.
module adj_face_mod

  use fv_arrays_mod,   only: fv_grid_type, fv_flags_type, fv_grid_bounds_type
  use fv_mp_mod,       only: ng, group_halo_update_type
  use mpp_domains_mod, only: domain2d, adj_set_cgrid_halo

  implicit none
  public

  integer, parameter :: n = 24, npx = n + 1, npy = n + 1, npz = 72, nq = 1
  integer, parameter :: is = 1, ie = n, js = 1, je = n
  integer, parameter :: isd = is - ng, ied = ie + ng, jsd = js - ng, jed = je + ng
  real, parameter :: dt = 1800.

  type(fv_grid_bounds_type), save :: bd
  type(fv_grid_type), target, save :: gridstruct
  type(fv_flags_type), save :: flagstruct
  type(domain2d), save :: domain
  type(group_halo_update_type), save :: q_pack

  !> C-grid winds: x-faces (is:ie+1, jsd:jed), y-faces (isd:ied, js:je+1), with halo
  real, allocatable, save :: cx_e(:,:,:), cy_e(:,:,:), mfx_e(:,:,:), mfy_e(:,:,:)
  real, allocatable, save :: ple0(:,:,:), ple1(:,:,:), pleadv_gchp(:,:,:)
  real, allocatable, save :: q_in(:,:,:), q_out_gchp(:,:,:)
  real, save :: ak(npz+1), bk(npz+1), ptop

  !> level used by the single-level kernel test (fv_tp_2d)
  integer, save :: k_test = 40

contains

  subroutine read_r8(dir, name, a, nelem)
    character(len=*), intent(in) :: dir, name
    integer, intent(in) :: nelem
    real, intent(out) :: a(nelem)
    integer :: u, ios, fsize
    open(newunit=u, file=trim(dir)//'/'//trim(name)//'.r8', access='stream', form='unformatted', &
         status='old', action='read', iostat=ios)
    if (ios /= 0) then
       write(*,*) 'cannot open ', trim(dir)//'/'//trim(name)//'.r8'
       error stop
    endif
    inquire(unit=u, size=fsize)
    if (fsize /= 8*nelem) then
       write(*,*) trim(name), ': file has ', fsize/8, ' values, expected ', nelem
       error stop
    endif
    read(u) a
    close(u)
  end subroutine read_r8

  subroutine load_face(dir)
    character(len=*), intent(in) :: dir

    bd%is = is; bd%ie = ie; bd%js = js; bd%je = je
    bd%isd = isd; bd%ied = ied; bd%jsd = jsd; bd%jed = jed

    allocate(gridstruct%area(isd:ied,jsd:jed), gridstruct%rarea(isd:ied,jsd:jed))
    allocate(gridstruct%dxa(isd:ied,jsd:jed), gridstruct%dya(isd:ied,jsd:jed))
    allocate(gridstruct%dx(isd:ied,jsd:jed+1), gridstruct%dy(isd:ied+1,jsd:jed))
    allocate(gridstruct%rdxc(isd:ied+1,jsd:jed), gridstruct%rdyc(isd:ied,jsd:jed+1))
    allocate(gridstruct%del6_u(isd:ied,jsd:jed+1), gridstruct%del6_v(isd:ied+1,jsd:jed))
    allocate(gridstruct%cosa_s(isd:ied,jsd:jed), gridstruct%rsin2(isd:ied,jsd:jed))
    allocate(gridstruct%sin_sg(isd:ied,jsd:jed,9), gridstruct%area_64(isd:ied,jsd:jed))
    call read_r8(dir, 'area',   gridstruct%area,   size(gridstruct%area))
    call read_r8(dir, 'dxa',    gridstruct%dxa,    size(gridstruct%dxa))
    call read_r8(dir, 'dya',    gridstruct%dya,    size(gridstruct%dya))
    call read_r8(dir, 'dx',     gridstruct%dx,     size(gridstruct%dx))
    call read_r8(dir, 'dy',     gridstruct%dy,     size(gridstruct%dy))
    call read_r8(dir, 'sin_sg', gridstruct%sin_sg, size(gridstruct%sin_sg))
    gridstruct%rarea   = 1. / gridstruct%area
    gridstruct%area_64 = gridstruct%area
    gridstruct%rdxc = 0.; gridstruct%rdyc = 0.; gridstruct%del6_u = 0.; gridstruct%del6_v = 0.
    gridstruct%cosa_s = 0.; gridstruct%rsin2 = 1.
    gridstruct%da_min = minval(gridstruct%area(is:ie,js:je))

    allocate(cx_e(is:ie+1,jsd:jed,npz), mfx_e(is:ie+1,jsd:jed,npz))
    allocate(cy_e(isd:ied,js:je+1,npz), mfy_e(isd:ied,js:je+1,npz))
    call read_r8(dir, 'cx',  cx_e,  size(cx_e))
    call read_r8(dir, 'mfx', mfx_e, size(mfx_e))
    call read_r8(dir, 'cy',  cy_e,  size(cy_e))
    call read_r8(dir, 'mfy', mfy_e, size(mfy_e))
    call adj_set_cgrid_halo(is, ie, js, je, ng, npz, cx_e, cy_e, mfx_e, mfy_e)

    allocate(ple0(is:ie,js:je,npz+1), ple1(is:ie,js:je,npz+1), pleadv_gchp(is:ie,js:je,npz+1))
    call read_r8(dir, 'ple0',   ple0,        size(ple0))
    call read_r8(dir, 'ple1',   ple1,        size(ple1))
    call read_r8(dir, 'pleadv', pleadv_gchp, size(pleadv_gchp))
    allocate(q_in(isd:ied,jsd:jed,npz), q_out_gchp(is:ie,js:je,npz))
    call read_r8(dir, 'q_in',       q_in,       size(q_in))
    call read_r8(dir, 'q_out_gchp', q_out_gchp, size(q_out_gchp))
    call read_r8(dir, 'ak', ak, npz+1)
    call read_r8(dir, 'bk', bk, npz+1)
    ptop = ak(1)
  end subroutine load_face

end module adj_face_mod
