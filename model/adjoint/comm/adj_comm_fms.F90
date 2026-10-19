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

!>@brief Six-face communication back-end (FMS, MPI) for the TL/AD test build.
!>@details The cubed-sphere domain comes from the real fv_mp_mod (domain_decomp, as in GCHP). The q
!! halo update is GCHP's: fv_mp_mod's start_group_halo_update + complete_group_halo_update. Its
!! adjoint is FMS's mpp_update_domains_ad (checked by test_halo: <H x, y> = <x, H^T y>). The global
!! sum is mpp_sum (self-adjoint), the global max fv_mp_mod's mp_reduce_max, and the C-grid wind
!! halos come from the mpp_get_boundary / mpp_update_domains calls of offline_tracer_advection.
!! See stubs/stubs.F90 for the TAF directives of the external routines.
module adj_backend_mod
  use mpp_mod,         only: mpp_pe, mpp_npes, mpp_sum, mpp_max, mpp_set_current_pelist, mpp_get_current_pelist
  use mpp_domains_mod, only: mpp_get_boundary, mpp_update_domains, mpp_update_domains_ad, CGRID_NE
  use fms_mod,         only: fms_init
  use fv_mp_mod,       only: mp_start, mp_assign_gid, domain_decomp, group_halo_update_type, &
                             start_group_halo_update, complete_group_halo_update
  use fv_arrays_mod,   only: fv_atmos_type
  implicit none
  private
  public :: multiface, qhalo_fill, qhalo, backend_init, backend_finalize, backend_rank, &
            backend_nranks, backend_bounds, backend_cgrid_halo, backend_allsum, backend_allmax
  public :: halo_q, halo_q_ad, global_sum, reduce_max

  logical, parameter :: multiface = .true.
  logical, save :: qhalo_fill = .false.              !< not used on six faces (real exchange)
  real, allocatable, save :: qhalo(:,:,:)

  integer, parameter :: npx = 25, npy = 25
  type(fv_atmos_type), target, save :: Atm
  type(group_halo_update_type), save :: q_group
  integer, save :: q_group_nk = -1

contains

  subroutine backend_init(layout)
    integer, intent(in) :: layout
    integer :: lay(2), io_layout(2), nregions, commid, i
    include 'mpif.h'
    ! as fv_control (fv_init / run_setup) does it, with halo_update_type = 1 (the default)
    call fms_init(MPI_COMM_WORLD)
    call mp_assign_gid()
    allocate(Atm%pelist(mpp_npes()))
    Atm%pelist = (/ (i, i = 0, mpp_npes()-1) /)
    Atm%npes_this_grid = mpp_npes()
    call mpp_set_current_pelist(Atm%pelist)
    call mpp_get_current_pelist(Atm%pelist, commID=commid)
    call mp_start(commid, 1)
    lay = (/ layout, layout /); io_layout = (/ 1, 1 /); nregions = 6
    if (mpp_npes() /= 6*layout*layout) then
       if (mpp_pe() == 0) write(*,*) 'backend_init: need 6*layout**2 ranks, have ', mpp_npes()
       error stop
    endif
    call domain_decomp(npx, npy, nregions, 0, .false., Atm, lay, io_layout)
  end subroutine backend_init

  !> MPI_Finalize directly: fms_end hangs in this standalone setting (after the clock report)
  subroutine backend_finalize()
    integer :: ierr
    include 'mpif.h'
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    call MPI_FINALIZE(ierr)
  end subroutine backend_finalize

  integer function backend_rank()
    backend_rank = mpp_pe()
  end function backend_rank

  integer function backend_nranks()
    backend_nranks = mpp_npes()
  end function backend_nranks

  subroutine backend_bounds(is, ie, js, je, tile)
    integer, intent(out) :: is, ie, js, je, tile
    is = Atm%bd%is; ie = Atm%bd%ie; js = Atm%bd%js; je = Atm%bd%je; tile = Atm%tile
  end subroutine backend_bounds

  !> C-grid halo of an (x, y) pair from the interior, as offline_tracer_advection does it:
  !! xi, yi on (is:ie, js:je, :); xe on (is:ie+1, jsd:jed, :), ye on (isd:ied, js:je+1, :)
  subroutine backend_cgrid_halo(xi, yi, xe, ye)
    real, intent(in) :: xi(:,:,:), yi(:,:,:)
    real, intent(inout) :: xe(:,:,:), ye(:,:,:)
    integer :: is, ie, js, je, isd, ied, jsd, jed, nk
    real, allocatable :: xl(:,:,:), yl(:,:,:), wb(:,:), eb(:,:), sb(:,:), nb(:,:)
    is = Atm%bd%is; ie = Atm%bd%ie; js = Atm%bd%js; je = Atm%bd%je
    isd = Atm%bd%isd; ied = Atm%bd%ied; jsd = Atm%bd%jsd; jed = Atm%bd%jed
    nk = size(xi, 3)
    allocate(xl(isd:ied+1, jsd:jed, nk), yl(isd:ied, jsd:jed+1, nk))
    allocate(wb(js:je, nk), eb(js:je, nk), sb(is:ie, nk), nb(is:ie, nk))
    xl = 0.; yl = 0.
    xl(is:ie, js:je, :) = xi
    yl(is:ie, js:je, :) = yi
    call mpp_get_boundary(xl, yl, Atm%domain, wbufferx=wb, ebufferx=eb, sbuffery=sb, nbuffery=nb, &
                          gridtype=CGRID_NE)
    xl(ie+1, js:je, :) = eb
    yl(is:ie, je+1, :) = nb
    call mpp_update_domains(xl, yl, Atm%domain, gridtype=CGRID_NE, complete=.true.)
    xe = xl(is:ie+1, jsd:jed, :)
    ye = yl(isd:ied, js:je+1, :)
  end subroutine backend_cgrid_halo

  real function backend_allsum(x)
    real, intent(in) :: x
    real :: a(1)
    a(1) = x
    call mpp_sum(a, 1)
    backend_allsum = a(1)
  end function backend_allsum

  real function backend_allmax(x)
    real, intent(in) :: x
    real :: a
    a = x
    call mpp_max(a)
    backend_allmax = a
  end function backend_allmax

  !> GCHP's q halo update (data domain arrays, any lower bounds)
  subroutine halo_q(q)
    real, intent(inout) :: q(:,:,:)
    if (q_group_nk /= -1 .and. q_group_nk /= size(q, 3)) error stop 'halo_q: one group per field shape'
    q_group_nk = size(q, 3)
    call start_group_halo_update(q_group, q, Atm%domain)
    call complete_group_halo_update(q_group, Atm%domain)
  end subroutine halo_q

  !> adjoint of halo_q
  subroutine halo_q_ad(q_ad)
    real, intent(inout) :: q_ad(:,:,:)
    call mpp_update_domains_ad(q_ad, Atm%domain)
  end subroutine halo_q_ad

  subroutine global_sum(a, n)
    integer, intent(in) :: n
    real, intent(inout) :: a(n)
    call mpp_sum(a, n)
  end subroutine global_sum

  subroutine reduce_max(a, n)
    use fv_mp_mod, only: mp_reduce_max
    integer, intent(in) :: n
    real, intent(inout) :: a(n)
    call mp_reduce_max(a, n)
  end subroutine reduce_max

end module adj_backend_mod

!> q halo update (black box for TAF; linear, its own TL)
subroutine adj_halo_q(q, ni, nj, nk)
  use adj_backend_mod, only: halo_q
  implicit none
  integer, intent(in) :: ni, nj, nk
  real, intent(inout) :: q(ni, nj, nk)
  call halo_q(q)
end subroutine adj_halo_q

!> adjoint of adj_halo_q; argument list as TAF generates the call (q, q_ad, then the DEPEND
!! arguments ni, nj, nk; run_taf.sh checks it)
subroutine adj_halo_q_ad(q, q_ad, ni, nj, nk)
  use adj_backend_mod, only: halo_q_ad
  implicit none
  integer, intent(in) :: ni, nj, nk
  real, intent(in) :: q(ni, nj, nk)
  real, intent(inout) :: q_ad(ni, nj, nk)
  call halo_q_ad(q_ad)
end subroutine adj_halo_q_ad

!> global sum over ranks (black box; linear and self-adjoint)
subroutine adj_global_sum(a, n)
  use adj_backend_mod, only: global_sum
  implicit none
  integer, intent(in) :: n
  real, intent(inout) :: a(n)
  call global_sum(a, n)
end subroutine adj_global_sum

!> global maximum over ranks (black box, passive)
subroutine adj_reduce_max(a, n)
  use adj_backend_mod, only: reduce_max
  implicit none
  integer, intent(in) :: n
  real, intent(inout) :: a(n)
  call reduce_max(a, n)
end subroutine adj_reduce_max
