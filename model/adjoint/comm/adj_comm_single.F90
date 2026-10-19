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

!>@brief Single-face communication back-end (one process, no MPI) for the TL/AD test build.
!>@details The q halo update is the identity: the q halo is an input of the test operator. With
!! qhalo_fill it is replaced by passive data instead (qhalo, the neighbouring faces' q from the
!! file), for the forward check against GCHP; then the halo carries no sensitivity. Global sums and
!! maxima are the identity. See stubs/stubs.F90 for the TAF directives of the external routines.
module adj_backend_mod
  implicit none
  private
  public :: multiface, qhalo_fill, qhalo, backend_init, backend_finalize, backend_rank, &
            backend_nranks, backend_bounds, backend_cgrid_halo, backend_allsum, backend_allmax

  logical, parameter :: multiface = .false.
  logical, save :: qhalo_fill = .false.
  real, allocatable, save :: qhalo(:,:,:)        !< (ni, nj, nk) as passed to adj_halo_q

contains

  subroutine backend_init(layout)
    integer, intent(in) :: layout
  end subroutine backend_init

  subroutine backend_finalize()
  end subroutine backend_finalize

  integer function backend_rank()
    backend_rank = 0
  end function backend_rank

  integer function backend_nranks()
    backend_nranks = 1
  end function backend_nranks

  subroutine backend_bounds(is, ie, js, je, tile)
    integer, intent(out) :: is, ie, js, je, tile
    error stop 'backend_bounds: single-face back-end'
  end subroutine backend_bounds

  subroutine backend_cgrid_halo(xi, yi, xe, ye)
    real, intent(in) :: xi(:,:,:), yi(:,:,:)
    real, intent(inout) :: xe(:,:,:), ye(:,:,:)
    error stop 'backend_cgrid_halo: single-face back-end'
  end subroutine backend_cgrid_halo

  real function backend_allsum(x)
    real, intent(in) :: x
    backend_allsum = x
  end function backend_allsum

  real function backend_allmax(x)
    real, intent(in) :: x
    backend_allmax = x
  end function backend_allmax

end module adj_backend_mod

!> q halo update (black box for TAF; linear, its own TL)
subroutine adj_halo_q(q, ni, nj, nk)
  use adj_backend_mod, only: qhalo_fill, qhalo
  implicit none
  integer, intent(in) :: ni, nj, nk
  real, intent(inout) :: q(ni, nj, nk)
  integer, parameter :: ng = 3
  if (qhalo_fill) then
     q(1:ng, :, :) = qhalo(1:ng, :, :)
     q(ni-ng+1:ni, :, :) = qhalo(ni-ng+1:ni, :, :)
     q(ng+1:ni-ng, 1:ng, :) = qhalo(ng+1:ni-ng, 1:ng, :)
     q(ng+1:ni-ng, nj-ng+1:nj, :) = qhalo(ng+1:ni-ng, nj-ng+1:nj, :)
  endif
end subroutine adj_halo_q

!> adjoint of adj_halo_q: identity, or the halo adjoint is dropped (halo from passive data);
!! argument list as TAF generates the call (q, q_ad, then ni, nj, nk; run_taf.sh checks it)
subroutine adj_halo_q_ad(q, q_ad, ni, nj, nk)
  use adj_backend_mod, only: qhalo_fill
  implicit none
  integer, intent(in) :: ni, nj, nk
  real, intent(inout) :: q(ni, nj, nk), q_ad(ni, nj, nk)
  integer, parameter :: ng = 3
  if (qhalo_fill) then
     q_ad(1:ng, :, :) = 0.
     q_ad(ni-ng+1:ni, :, :) = 0.
     q_ad(ng+1:ni-ng, 1:ng, :) = 0.
     q_ad(ng+1:ni-ng, nj-ng+1:nj, :) = 0.
  endif
end subroutine adj_halo_q_ad

!> global sum over ranks (black box; linear and self-adjoint): the identity on one process
subroutine adj_global_sum(a, n)
  implicit none
  integer, intent(in) :: n
  real, intent(inout) :: a(n)
end subroutine adj_global_sum

!> global maximum over ranks (black box, passive): the identity on one process
subroutine adj_reduce_max(a, n)
  implicit none
  integer, intent(in) :: n
  real, intent(inout) :: a(n)
end subroutine adj_reduce_max
