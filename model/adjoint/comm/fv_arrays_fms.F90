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

!>@brief The part of fv_arrays_mod (model/fv_arrays.F90) that fv_mp_mod.F90 uses (domain_decomp),
!! for the FMS communication back-end of the TL/AD tests. The components have the names and types
!! of fv_arrays.F90.
module fv_arrays_mod
  use mpp_domains_mod, only: domain2d
  implicit none
  private
  public :: fv_grid_bounds_type, fv_grid_type, fv_atmos_type

  type fv_grid_bounds_type
     integer :: is,  ie,  js,  je
     integer :: isd, ied, jsd, jed
     integer :: isc, iec, jsc, jec
     integer :: ng
  end type fv_grid_bounds_type

  type fv_grid_type
     logical :: square_domain = .false.
  end type fv_grid_type

  type fv_atmos_type
     integer :: grid_number = 1
     integer, allocatable, dimension(:) :: pelist
     type(fv_grid_bounds_type) :: bd
     type(domain2D) :: domain
     type(domain2D) :: domain_for_coupler
     integer :: num_contact, npes_per_tile, tile, npes_this_grid
     type(fv_grid_type) :: gridstruct
  end type fv_atmos_type

end module fv_arrays_mod
