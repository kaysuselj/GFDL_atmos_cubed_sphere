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

!>@brief The module 'fv_tracer2d.F90' performs sub-cycled tracer advection.
!>@see \cite lin2004vertically

! Modules Included:
! <table>
! <tr>
!     <th>Module Name</th>
!     <th>Functions Included</th>
!   </tr>
! <table>
!   <tr>
!     <td>boundary_mod</td>
!     <td>nested_grid_BC_apply_intT</td>
!   </tr>
!   <tr>
!     <td>fv_arrays_mod</td>
!     <td>fv_grid_type, fv_nest_type, fv_atmos_type, fv_grid_bounds_type</td>
!   </tr>
!   <tr>
!   <tr>
!     <td>fv_mp_mod</td>
!     <td>mp_reduce_max, ng, mp_gather, is_master, group_halo_update_type, 
!         start_group_halo_update, complete_group_halo_update</td>
!   </tr>
!    <tr>
!     <td>fv_timing_mod</td>
!     <td>timing_on, timing_off</td>
!   </tr>
!  <tr>
!     <td>mpp_mod</td>
!     <td>mpp_error, FATAL, mpp_broadcast, mpp_send, mpp_recv, mpp_sum, mpp_max</td>
!   </tr>
!   <tr>
!     <td>mpp_domains_mod</td>
!     <td>mpp_update_domains, CGRID_NE, domain2d</td>
!   </tr>
!   <tr>
!     <td>tp_core_mod</td>
!     <td>fv_tp_2d, copy_corners</td>
!   </tr>
! </table>

module fv_tracer2d_mod
   use tp_core_mod,       only: fv_tp_2d, copy_corners
   use fv_mp_mod,         only: mp_reduce_max
   use fv_mp_mod,         only: ng, mp_gather, is_master
   use fv_mp_mod,         only: group_halo_update_type
   use fv_mp_mod,         only: start_group_halo_update, complete_group_halo_update
   use mpp_domains_mod,   only: mpp_update_domains, CGRID_NE, domain2d, mpp_get_boundary
   use mpp_domains_mod,   only: NON_BITWISE_EXACT_SUM, BITWISE_EXACT_SUM, BITWISE_EFP_SUM
   use fv_timing_mod,     only: timing_on, timing_off
   use boundary_mod,      only: nested_grid_BC_apply_intT
   use fv_arrays_mod,     only: fv_grid_type, fv_flags_type, fv_nest_type, fv_atmos_type, fv_grid_bounds_type, REAL8
   use mpp_mod,           only: mpp_error, FATAL, mpp_broadcast, mpp_send, mpp_recv, mpp_sum, mpp_max
   use fv_grid_utils_mod, only: g_sum_r8

implicit none
private

public :: tracer_2d, tracer_2d_nested, tracer_2d_1L, offline_tracer_advection

real, allocatable, dimension(:,:,:) :: nest_fx_west_accum, nest_fx_east_accum, nest_fx_south_accum, nest_fx_north_accum

! GCHP: parameters rescale_trop_only and show_scaling
logical, parameter :: rescale_trop_only=.True.
logical, parameter :: show_scaling=.False.

contains

!>@brief The subroutine 'tracer_2d_1L' performs 2-D horizontal-to-lagrangian transport.
!>@details This subroutine is called if 'z_tracer = .true.'
!! It modifies 'tracer_2d' so that each layer uses a different diagnosed number 
!! of split tracer timesteps. This potentially accelerates tracer advection when there
!! is a large difference in layer-maximum wind speeds (cf. polar night jet).
subroutine tracer_2d_1L(q, dp1, mfx, mfy, cx, cy, gridstruct, bd, domain, npx, npy, npz,   &
                        nq,  hord, dt, id_divg, q_pack, nord_tr, trdm, lim_fac, dpA)

      type(fv_grid_bounds_type), intent(IN) :: bd
      integer, intent(IN) :: npx
      integer, intent(IN) :: npy
      integer, intent(IN) :: npz
      integer, intent(IN) :: nq    !< number of tracers to be advected
      integer, intent(IN) :: hord, nord_tr
      integer, intent(IN) :: id_divg
      real   , intent(IN) :: dt, trdm
      real   , intent(IN) :: lim_fac
      type(group_halo_update_type), intent(inout) :: q_pack
      real   , intent(INOUT) :: q(bd%isd:bd%ied,bd%jsd:bd%jed,npz,nq)   !< Tracers
      real   , intent(INOUT) :: dp1(bd%isd:bd%ied,bd%jsd:bd%jed,npz)    !< DELP before dyn_core
      real   , intent(IN   ) :: mfx(bd%is:bd%ie+1,bd%js:bd%je,  npz)    !< Mass Flux X-Dir
      real   , intent(IN   ) :: mfy(bd%is:bd%ie  ,bd%js:bd%je+1,npz)    !< Mass Flux Y-Dir
      real   , intent(IN   ) ::  cx(bd%is:bd%ie+1,bd%jsd:bd%jed  ,npz)  !< Courant Number X-Dir
      real   , intent(IN   ) ::  cy(bd%isd:bd%ied,bd%js :bd%je +1,npz)  !< Courant Number Y-Dir
      real   , optional, intent(OUT) :: dpA(bd%is:bd%ie,bd%js:bd%je,npz)    ! DELP after advection
      type(fv_grid_type), intent(IN), target :: gridstruct
      type(domain2d), intent(INOUT) :: domain

! Local Arrays
      real :: qn2(bd%isd:bd%ied,bd%jsd:bd%jed,nq)   !< 3D tracers
      real :: dp2(bd%is:bd%ie,bd%js:bd%je)
      real :: fx(bd%is:bd%ie+1,bd%js:bd%je )
      real :: fy(bd%is:bd%ie , bd%js:bd%je+1)
      real :: ra_x(bd%is:bd%ie,bd%jsd:bd%jed)
      real :: ra_y(bd%isd:bd%ied,bd%js:bd%je)
      real :: xfx(bd%is:bd%ie+1,bd%jsd:bd%jed  ,npz)
      real :: yfx(bd%isd:bd%ied,bd%js: bd%je+1, npz)
      real :: mfx2(bd%is:bd%ie+1,bd%js:bd%je,npz)
      real :: mfy2(bd%is:bd%ie,bd%js:bd%je+1,npz)
      real ::  cx2(bd%is:bd%ie+1,bd%jsd:bd%jed, npz)
      real ::  cy2(bd%isd:bd%ied,bd%js :bd%je +1, npz)
      real :: c2d(bd%is:bd%ie,bd%js:bd%je)
      real :: cmax(npz), cmin(npz)
      real :: qmax(npz*(nq+1)), qmin(npz*(nq+1))
      integer :: icount(npz,nq)
      real :: frac
      integer :: nsplt
      integer :: i,j,k,n,it,iq

      real, pointer, dimension(:,:) :: area, rarea
      real, pointer, dimension(:,:,:) :: sin_sg
      real, pointer, dimension(:,:) :: dxa, dya, dx, dy

      integer :: is,  ie,  js,  je
      integer :: isd, ied, jsd, jed

      is  = bd%is
      ie  = bd%ie
      js  = bd%js
      je  = bd%je
      isd = bd%isd
      ied = bd%ied
      jsd = bd%jsd
      jed = bd%jed

       area => gridstruct%area
      rarea => gridstruct%rarea

      sin_sg => gridstruct%sin_sg
      dxa    => gridstruct%dxa 
      dya    => gridstruct%dya 
      dx     => gridstruct%dx  
      dy     => gridstruct%dy  

!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,jsd,jed,npz,cx,xfx,dxa,dy, &
!$OMP                                  sin_sg,cy,yfx,dya,dx) private(i,j,k)
  do k=1,npz
     do j=jsd,jed
        do i=is,ie+1
           if (cx(i,j,k) > 0.) then
              xfx(i,j,k) = cx(i,j,k)*dxa(i-1,j)*dy(i,j)*sin_sg(i-1,j,3)
           else
              xfx(i,j,k) = cx(i,j,k)*dxa(i,  j)*dy(i,j)*sin_sg(i,  j,1)
           endif
        enddo
     enddo
     do j=js,je+1
        do i=isd,ied
           if (cy(i,j,k) > 0.) then
              yfx(i,j,k) = cy(i,j,k)*dya(i,j-1)*dx(i,j)*sin_sg(i,j-1,4)
           else
              yfx(i,j,k) = cy(i,j,k)*dya(i,j  )*dx(i,j)*sin_sg(i,j,  2)
           endif
        enddo
     enddo
  enddo  ! k-loop

                               call timing_on('COMM_TOTAL')
                         call timing_on('COMM_TRACER')
  call complete_group_halo_update(q_pack, domain)
                        call timing_off('COMM_TRACER')
                              call timing_off('COMM_TOTAL')

! Check for levels where Q does not need to be advected
!$OMP parallel do default(none) shared(nq,npz,is,ie,js,je,q,qmax) private(iq,i,j,k,n)
   do iq=1,nq
      do k=1,npz
         n=(iq-1)*npz + k
         qmax(n) = 0.0
         do j=js,je
            do i=is,ie
               qmax(n) = max(qmax(n),q(i,j,k,iq))
            enddo
         enddo
      enddo
   enddo
!$OMP parallel do default(none) shared(is,ie,js,je,npz,nq,&
!$OMP                                  sin_sg,cx,cy,cmax,qmax) private(i,j,k,n)
   do k=1,npz
      cmax(k) = 0.
      if ( k < npz/6 ) then
           do j=js,je
              do i=is,ie
                cmax(k) = max( cmax(k), abs(cx(i,j,k)), abs(cy(i,j,k)) )
              enddo
           enddo
      else
           do j=js,je
              do i=is,ie
                cmax(k) = max( cmax(k), max(abs(cx(i,j,k)),abs(cy(i,j,k)))+1.-sin_sg(i,j,5) )
              enddo
           enddo
      endif
 !!!  if ( is_master() )  write(*,*) 'tracer_2d_1L: k, nsplt, cmax =', k, int(1. + cmax(k)), cmax(k)
    ! add to qmax for allreduce
     n=nq*npz + k
     qmax(n) = cmax(k)
   enddo  ! k-loop

                        call timing_on('COMM_TOTAL')
                            call timing_on('COMM_TRACER_MAX')
  call mp_reduce_max(qmax,(nq+1)*npz)
                           call timing_off('COMM_TRACER_MAX')
                       call timing_off('COMM_TOTAL')

  ! get cmax from allreduce array
  do k=1,npz
     n=nq*npz + k
     cmax(k) = qmax(n)
  enddo  ! k-loop

!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,jsd,jed,npz,cx,xfx, &
!$OMP                                  cy,yfx,mfx,mfy,cmax,mfx2,mfy2,cx2,cy2)   &
!$OMP                          private(nsplt, frac)
  do k=1,npz

     mfx2(:,:,k)=mfx(:,:,k)
     mfy2(:,:,k)=mfy(:,:,k)
     cx2(:,:,k)=cx(:,:,k)
     cy2(:,:,k)=cy(:,:,k)
     nsplt = int(1. + cmax(k))
     if ( nsplt > 1 ) then
        frac  = 1. / real(nsplt)
        do j=jsd,jed
           do i=is,ie+1
               cx2(i,j,k) =  cx(i,j,k) * frac
              xfx(i,j,k) = xfx(i,j,k) * frac
           enddo
        enddo
        do j=js,je
           do i=is,ie+1
              mfx2(i,j,k) = mfx(i,j,k) * frac
           enddo
        enddo
        do j=js,je+1
           do i=isd,ied
              cy2(i,j,k) =  cy(i,j,k) * frac
             yfx(i,j,k) = yfx(i,j,k) * frac
           enddo
        enddo
        do j=js,je+1
           do i=is,ie
              mfy2(i,j,k) = mfy(i,j,k) * frac
           enddo
        enddo
     endif

  enddo

! Begin k-independent tracer transport; can not be OpenMPed because the mpp_update call.
  do k=1,npz

!$OMP parallel do default(none) shared(k,is,ie,js,je,isd,ied,jsd,jed,xfx,area,yfx,ra_x,ra_y)
     do j=jsd,jed
        do i=is,ie
           ra_x(i,j) = area(i,j) + (xfx(i,j,k) - xfx(i+1,j,k))
        enddo
        if ( j>=js .and. j<=je ) then
           do i=isd,ied
              ra_y(i,j) = area(i,j) + (yfx(i,j,k) - yfx(i,j+1,k))
           enddo
        endif
     enddo

     nsplt = int(1. + cmax(k))

     do it=1,nsplt

!$OMP parallel do default(none) shared(k,is,ie,js,je,rarea,mfx2,mfy2,dp1,dp2)
        do j=js,je
           do i=is,ie
              dp2(i,j) = dp1(i,j,k) + ((mfx2(i,j,k)-mfx2(i+1,j,k))+(mfy2(i,j,k)-mfy2(i,j+1,k)))*rarea(i,j)
           enddo
        enddo

!$OMP parallel do default(none) shared(k,nsplt,it,is,ie,js,je,isd,ied,jsd,jed,npx,npy,cx2,xfx,hord,trdm, &
!$OMP                                  nord_tr,nq,gridstruct,bd,cy2,yfx,mfx2,mfy2,qn2,q,ra_x,ra_y,dp1,dp2,rarea,lim_fac) & 
!$OMP                          private(fx,fy)
        do iq=1,nq
        if ( nsplt /= 1 ) then
           if ( it==1 ) then
              do j=jsd,jed
                 do i=isd,ied
                    qn2(i,j,iq) = q(i,j,k,iq)
                 enddo
              enddo
           endif
           call fv_tp_2d(qn2(isd,jsd,iq), cx2(is,jsd,k), cy2(isd,js,k), &
                         npx, npy, hord, fx, fy, xfx(is,jsd,k), yfx(isd,js,k), &
                         gridstruct, bd, ra_x, ra_y, lim_fac, mfx=mfx2(is,js,k), mfy=mfy2(is,js,k))
           if ( it < nsplt ) then   ! not last call
              do j=js,je
              do i=is,ie
                 qn2(i,j,iq) = (qn2(i,j,iq)*dp1(i,j,k)+((fx(i,j)-fx(i+1,j))+(fy(i,j)-fy(i,j+1)))*rarea(i,j))/dp2(i,j)
              enddo
              enddo
           else
              do j=js,je
              do i=is,ie
                 q(i,j,k,iq) = (qn2(i,j,iq)*dp1(i,j,k)+((fx(i,j)-fx(i+1,j))+(fy(i,j)-fy(i,j+1)))*rarea(i,j))/dp2(i,j)
              enddo
              enddo
           endif
        else
           call fv_tp_2d(q(isd,jsd,k,iq), cx2(is,jsd,k), cy2(isd,js,k), &
                         npx, npy, hord, fx, fy, xfx(is,jsd,k), yfx(isd,js,k), &
                         gridstruct, bd, ra_x, ra_y, lim_fac, mfx=mfx2(is,js,k), mfy=mfy2(is,js,k))
           do j=js,je
              do i=is,ie
                 q(i,j,k,iq) = (q(i,j,k,iq)*dp1(i,j,k)+((fx(i,j)-fx(i+1,j))+(fy(i,j)-fy(i,j+1)))*rarea(i,j))/dp2(i,j)
              enddo
           enddo
        endif
        enddo   !  tracer-loop

        if ( it < nsplt ) then   ! not last call
             do j=js,je
                do i=is,ie
                   dp1(i,j,k) = dp2(i,j)
                enddo
             enddo
                               call timing_on('COMM_TOTAL')
                         call timing_on('COMM_TRACER')
             call mpp_update_domains(qn2, domain)
                        call timing_off('COMM_TRACER')
                              call timing_off('COMM_TOTAL')
        endif
     enddo  ! time-split loop

     if (present(dpA)) then
        dpA(:,:,k)=dp2
     endif

  enddo    ! k-loop

end subroutine tracer_2d_1L

!>@brief The subroutine 'tracer_2d' is the standard routine for sub-cycled tracer advection.
subroutine tracer_2d(q, dp1, mfx, mfy, cx, cy, gridstruct, bd, domain, npx, npy, npz,   &
                     nq,  hord, q_split, dt, id_divg, q_pack, nord_tr, trdm, lim_fac, dpA)

      type(fv_grid_bounds_type), intent(IN) :: bd
      integer, intent(IN) :: npx
      integer, intent(IN) :: npy
      integer, intent(IN) :: npz
      integer, intent(IN) :: nq    !< number of tracers to be advected
      integer, intent(IN) :: hord, nord_tr
      integer, intent(IN) :: q_split
      integer, intent(IN) :: id_divg
      real   , intent(IN) :: dt, trdm
      real   , intent(IN) :: lim_fac
      type(group_halo_update_type), intent(inout) :: q_pack
      real   , intent(INOUT) :: q(bd%isd:bd%ied,bd%jsd:bd%jed,npz,nq)   !< Tracers
      real   , intent(INOUT) :: dp1(bd%isd:bd%ied,bd%jsd:bd%jed,npz)    !< DELP before dyn_core
      real   , intent(IN   ) :: mfx(bd%is:bd%ie+1,bd%js:bd%je,  npz)    !< Mass Flux X-Dir
      real   , intent(IN   ) :: mfy(bd%is:bd%ie  ,bd%js:bd%je+1,npz)    !< Mass Flux Y-Dir
      real   , intent(IN   ) ::  cx(bd%is:bd%ie+1,bd%jsd:bd%jed  ,npz)  !< Courant Number X-Dir
      real   , intent(IN   ) ::  cy(bd%isd:bd%ied,bd%js :bd%je +1,npz)  !< Courant Number Y-Dir
      real   , optional, intent(OUT) :: dpA(bd%is:bd%ie,bd%js:bd%je,npz)! DELP after advection
      type(fv_grid_type), intent(IN), target :: gridstruct
      type(domain2d), intent(INOUT) :: domain

! Local Arrays
      real :: dp2(bd%is:bd%ie,bd%js:bd%je)
      real :: fx(bd%is:bd%ie+1,bd%js:bd%je )
      real :: fy(bd%is:bd%ie , bd%js:bd%je+1)
      real :: ra_x(bd%is:bd%ie,bd%jsd:bd%jed)
      real :: ra_y(bd%isd:bd%ied,bd%js:bd%je)
      real :: xfx(bd%is:bd%ie+1,bd%jsd:bd%jed  ,npz)
      real :: yfx(bd%isd:bd%ied,bd%js: bd%je+1, npz)
      real :: mfx2(bd%is:bd%ie+1,bd%js:bd%je,npz)
      real :: mfy2(bd%is:bd%ie,bd%js:bd%je+1,npz)
      real ::  cx2(bd%is:bd%ie+1,bd%jsd:bd%jed, npz)
      real ::  cy2(bd%isd:bd%ied,bd%js :bd%je +1, npz)
      real :: cmax(npz)
      real :: frac, rdt
      integer :: maxsplt, ksplt(npz)
      integer :: i,j,k,n,it,iq

      real, pointer, dimension(:,:) :: area, rarea
      real, pointer, dimension(:,:,:) :: sin_sg
      real, pointer, dimension(:,:) :: dxa, dya, dx, dy

      integer :: is,  ie,  js,  je
      integer :: isd, ied, jsd, jed

      is  = bd%is
      ie  = bd%ie
      js  = bd%js
      je  = bd%je
      isd = bd%isd
      ied = bd%ied
      jsd = bd%jsd
      jed = bd%jed

       area => gridstruct%area
      rarea => gridstruct%rarea

      sin_sg => gridstruct%sin_sg
      dxa    => gridstruct%dxa 
      dya    => gridstruct%dya 
      dx     => gridstruct%dx  
      dy     => gridstruct%dy  

!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,jsd,jed,npz,cx,xfx,dxa,dy, &  
!$OMP                                  sin_sg,cy,yfx,dya,dx,cmax,q_split)
    do k=1,npz
       do j=jsd,jed
          do i=is,ie+1
             if (cx(i,j,k) > 0.) then
                 xfx(i,j,k) = cx(i,j,k)*dxa(i-1,j)*dy(i,j)*sin_sg(i-1,j,3)
             else
                 xfx(i,j,k) = cx(i,j,k)*dxa(i,j)*dy(i,j)*sin_sg(i,j,1)
             endif
          enddo
       enddo
       do j=js,je+1
          do i=isd,ied
              if (cy(i,j,k) > 0.) then
                  yfx(i,j,k) = cy(i,j,k)*dya(i,j-1)*dx(i,j)*sin_sg(i,j-1,4)
              else
                  yfx(i,j,k) = cy(i,j,k)*dya(i,j)*dx(i,j)*sin_sg(i,j,2)
              endif
          enddo
       enddo
       if ( q_split == 0 ) then
         cmax(k) = 0.
         if ( k < npz/6 ) then
            do j=js,je
               do i=is,ie
                  cmax(k) = max( cmax(k), abs(cx(i,j,k)), abs(cy(i,j,k)) )
               enddo
            enddo
         else
            do j=js,je
               do i=is,ie
                  cmax(k) = max( cmax(k), max(abs(cx(i,j,k)),abs(cy(i,j,k)))+1.-sin_sg(i,j,5) )
               enddo
            enddo
         endif
       endif
    enddo

!--------------------------------------------------------------------------------

! Determine global cmax on levels:
    if ( q_split == 0 ) then
       call mp_reduce_max(cmax,npz)
    else
       cmax = q_split-1.
    endif

! Determine maxsplt for outer iteration loop. 
! Advection will be done base on levels using ksplt
    maxsplt = 0
    do k=1,npz
       ksplt(k) = int(1. + cmax(k))
       maxsplt = max(maxsplt,ksplt(k))
    enddo

!--------------------------------------------------------------------------------

!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,jsd,jed,npz,cx,cx2,xfx,mfx,mfx2,cy,cy2,yfx,mfy,mfy2,cmax,ksplt) &
!$OMP                          private( frac )
    do k=1,npz
        mfx2(:,:,k)=mfx(:,:,k)
        mfy2(:,:,k)=mfy(:,:,k)
        cx2(:,:,k)=cx(:,:,k)
        cy2(:,:,k)=cy(:,:,k)
        if( ksplt(k) /= 1 ) then
            frac  = 1. / real(ksplt(k))
            do j=jsd,jed
               do i=is,ie+1
                  cx2(i,j,k) =   cx(i,j,k) * frac
                  xfx(i,j,k) =  xfx(i,j,k) * frac
               enddo
            enddo
            do j=js,je
               do i=is,ie+1
                  mfx2(i,j,k) = mfx(i,j,k) * frac
               enddo
            enddo
            do j=js,je+1
               do i=isd,ied
                  cy2(i,j,k) =  cy(i,j,k) * frac
                  yfx(i,j,k) = yfx(i,j,k) * frac
               enddo
            enddo
            do j=js,je+1
               do i=is,ie
                  mfy2(i,j,k) = mfy(i,j,k) * frac
               enddo
            enddo
        endif
    enddo

    do it=1,maxsplt
                        call timing_on('COMM_TOTAL')
                            call timing_on('COMM_TRACER')
      call complete_group_halo_update(q_pack, domain)
                           call timing_off('COMM_TRACER')
                       call timing_off('COMM_TOTAL')

!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,jsd,jed,npz,dp1,mfx2,mfy2,rarea,nq,&
!$OMP                                  area,xfx,yfx,q,cx2,cy2,npx,npy,hord,gridstruct,bd,it,ksplt,nord_tr,trdm,lim_fac,dpA) &
!$OMP                          private(dp2, ra_x, ra_y, fx, fy, i, j, k, n, iq)
     do k=1,npz

       if ( it <= ksplt(k) ) then

         do j=js,je
            do i=is,ie
               dp2(i,j) = dp1(i,j,k) + ((mfx2(i,j,k)-mfx2(i+1,j,k))+(mfy2(i,j,k)-mfy2(i,j+1,k)))*rarea(i,j)
            enddo
         enddo

         do j=jsd,jed
            do i=is,ie
               ra_x(i,j) = area(i,j) + (xfx(i,j,k) - xfx(i+1,j,k))
               if (cx2(i,j,k) > 1.0)  write(*,*) 'cx2(i,j,k) > 1.0 : ', cx2(i,j,k), i, j, k
            enddo
         enddo
         do j=js,je
            do i=isd,ied
               ra_y(i,j) = area(i,j) + (yfx(i,j,k) - yfx(i,j+1,k))
               if (cy2(i,j,k) > 1.0)  write(*,*) 'cy2(i,j,k) > 1.0 : ', cy2(i,j,k), i, j, k
            enddo
         enddo

         do iq=1,nq
           if ( it==1 .and. trdm>1.e-4 ) then
            call fv_tp_2d(q(isd,jsd,k,iq), cx2(is,jsd,k), cy2(isd,js,k), &
                          npx, npy, hord, fx, fy, xfx(is,jsd,k), yfx(isd,js,k), &
                          gridstruct, bd, ra_x, ra_y, lim_fac, mfx=mfx2(is,js,k), mfy=mfy2(is,js,k),   &
                          mass=dp1(isd,jsd,k), nord=nord_tr, damp_c=trdm)
           else
            call fv_tp_2d(q(isd,jsd,k,iq), cx2(is,jsd,k), cy2(isd,js,k), &
                          npx, npy, hord, fx, fy, xfx(is,jsd,k), yfx(isd,js,k), &
                          gridstruct, bd, ra_x, ra_y, lim_fac, mfx=mfx2(is,js,k), mfy=mfy2(is,js,k))
           endif
           do j=js,je
              do i=is,ie
                 q(i,j,k,iq) = ( q(i,j,k,iq)*dp1(i,j,k) + &
                               ((fx(i,j)-fx(i+1,j))+(fy(i,j)-fy(i,j+1)))*rarea(i,j) )/dp2(i,j)
              enddo
           enddo
         enddo

         if ( it /= ksplt(k) ) then
            do j=js,je
               do i=is,ie
                  dp1(i,j,k) = dp2(i,j)
               enddo
            enddo
         elseif (it == ksplt(k)) then
            if (present(dpA)) then
               dpA(:,:,k)=dp2
            end if
         endif

       endif   ! ksplt(k)

     enddo ! npz

     if ( it /= maxsplt ) then
                      call timing_on('COMM_TOTAL')
                          call timing_on('COMM_TRACER')
           call start_group_halo_update(q_pack, q, domain)
                          call timing_off('COMM_TRACER')
                      call timing_off('COMM_TOTAL')
     endif

   enddo  ! maxsplt


end subroutine tracer_2d


subroutine tracer_2d_nested(q, dp1, mfx, mfy, cx, cy, gridstruct, bd, domain, npx, npy, npz,   &
                     nq,  hord, q_split, dt, id_divg, q_pack, nord_tr, trdm, &
                     k_split, neststruct, parent_grid, lim_fac)

      type(fv_grid_bounds_type), intent(IN) :: bd
      integer, intent(IN) :: npx
      integer, intent(IN) :: npy
      integer, intent(IN) :: npz
      integer, intent(IN) :: nq    !< number of tracers to be advected
      integer, intent(IN) :: hord, nord_tr
      integer, intent(IN) :: q_split, k_split
      integer, intent(IN) :: id_divg
      real   , intent(IN) :: dt, trdm
      real   , intent(IN) :: lim_fac
      type(group_halo_update_type), intent(inout) :: q_pack
      real   , intent(INOUT) :: q(bd%isd:bd%ied,bd%jsd:bd%jed,npz,nq)   !< Tracers
      real   , intent(INOUT) :: dp1(bd%isd:bd%ied,bd%jsd:bd%jed,npz)    !< DELP before dyn_core
      real   , intent(INOUT) :: mfx(bd%is:bd%ie+1,bd%js:bd%je,  npz)    !< Mass Flux X-Dir
      real   , intent(INOUT) :: mfy(bd%is:bd%ie  ,bd%js:bd%je+1,npz)    !< Mass Flux Y-Dir
      real   , intent(INOUT) ::  cx(bd%is:bd%ie+1,bd%jsd:bd%jed  ,npz)  !< Courant Number X-Dir
      real   , intent(INOUT) ::  cy(bd%isd:bd%ied,bd%js :bd%je +1,npz)  !< Courant Number Y-Dir
      type(fv_grid_type), intent(IN), target :: gridstruct
      type(fv_nest_type), intent(INOUT) :: neststruct
      type(fv_atmos_type), intent(INOUT) :: parent_grid
      type(domain2d), intent(INOUT) :: domain

! Local Arrays
      real :: dp2(bd%is:bd%ie,bd%js:bd%je)
      real :: fx(bd%is:bd%ie+1,bd%js:bd%je )
      real :: fy(bd%is:bd%ie , bd%js:bd%je+1)
      real :: ra_x(bd%is:bd%ie,bd%jsd:bd%jed)
      real :: ra_y(bd%isd:bd%ied,bd%js:bd%je)
      real :: xfx(bd%is:bd%ie+1,bd%jsd:bd%jed  ,npz)
      real :: yfx(bd%isd:bd%ied,bd%js: bd%je+1, npz)
      real :: cmax(npz)
      real :: cmax_t
      real :: c_global
      real :: frac, rdt
      integer :: nsplt, nsplt_parent, msg_split_steps = 1
      integer :: i,j,k,it,iq

      real, pointer, dimension(:,:) :: area, rarea
      real, pointer, dimension(:,:,:) :: sin_sg
      real, pointer, dimension(:,:) :: dxa, dya, dx, dy

      integer :: is,  ie,  js,  je
      integer :: isd, ied, jsd, jed

      is  = bd%is
      ie  = bd%ie
      js  = bd%js
      je  = bd%je
      isd = bd%isd
      ied = bd%ied
      jsd = bd%jsd
      jed = bd%jed

       area => gridstruct%area
      rarea => gridstruct%rarea

      sin_sg => gridstruct%sin_sg
      dxa    => gridstruct%dxa 
      dya    => gridstruct%dya 
      dx     => gridstruct%dx  
      dy     => gridstruct%dy  

!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,jsd,jed,npz,cx,xfx,dxa,dy, &
!$OMP                                  sin_sg,cy,yfx,dya,dx)
      do k=1,npz
         do j=jsd,jed
            do i=is,ie+1
               if (cx(i,j,k) > 0.) then
                  xfx(i,j,k) = cx(i,j,k)*dxa(i-1,j)*dy(i,j)*sin_sg(i-1,j,3)
               else
                  xfx(i,j,k) = cx(i,j,k)*dxa(i,j)*dy(i,j)*sin_sg(i,j,1)
               endif
            enddo
         enddo
         do j=js,je+1
            do i=isd,ied
               if (cy(i,j,k) > 0.) then
                  yfx(i,j,k) = cy(i,j,k)*dya(i,j-1)*dx(i,j)*sin_sg(i,j-1,4)
               else
                  yfx(i,j,k) = cy(i,j,k)*dya(i,j)*dx(i,j)*sin_sg(i,j,2)
               endif
            enddo
         enddo
      enddo

!--------------------------------------------------------------------------------
  if ( q_split == 0 ) then
! Determine nsplt

!$OMP parallel do default(none) shared(is,ie,js,je,npz,cmax,cx,cy,sin_sg) &
!$OMP                          private(cmax_t )
      do k=1,npz
         cmax(k) = 0.
         if ( k < 4 ) then
! Top layers: C < max( abs(c_x), abs(c_y) )
            do j=js,je
               do i=is,ie
                  cmax_t  = max( abs(cx(i,j,k)), abs(cy(i,j,k)) )
                  cmax(k) = max( cmax_t, cmax(k) )
               enddo
            enddo
         else
            do j=js,je
               do i=is,ie
                  cmax_t  = max(abs(cx(i,j,k)), abs(cy(i,j,k))) + 1.-sin_sg(i,j,5)
                  cmax(k) = max( cmax_t, cmax(k) )
               enddo
            enddo
         endif
      enddo
      call mp_reduce_max(cmax,npz)

! find global max courant number and define nsplt to scale cx,cy,mfx,mfy
      c_global = cmax(1)
      if ( npz /= 1 ) then                ! if NOT shallow water test case
         do k=2,npz
            c_global = max(cmax(k), c_global)
         enddo
      endif
      nsplt = int(1. + c_global)
!!    if ( is_master() .and. nsplt > 3 )  write(*,*) 'Tracer_2d_split=', nsplt, c_global
   else
      nsplt = q_split
      if (gridstruct%nested .and. neststruct%nestbctype > 1) msg_split_steps = max(q_split/parent_grid%flagstruct%q_split,1)
   endif

!--------------------------------------------------------------------------------

   frac  = 1. / real(nsplt)

      if( nsplt /= 1 ) then
!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,jsd,jed,npz,cx,frac,xfx,mfx,cy,yfx,mfy)
          do k=1,npz
             do j=jsd,jed
                do i=is,ie+1
                   cx(i,j,k) =  cx(i,j,k) * frac
                   xfx(i,j,k) = xfx(i,j,k) * frac
                enddo
             enddo
             do j=js,je
                do i=is,ie+1
                   mfx(i,j,k) = mfx(i,j,k) * frac
                enddo
             enddo

             do j=js,je+1
                do i=isd,ied
                   cy(i,j,k) =  cy(i,j,k) * frac
                  yfx(i,j,k) = yfx(i,j,k) * frac
                enddo
             enddo

             do j=js,je+1
                do i=is,ie
                  mfy(i,j,k) = mfy(i,j,k) * frac
                enddo
             enddo
          enddo
      endif


    do it=1,nsplt
       if ( gridstruct%nested ) then
          neststruct%tracer_nest_timestep = neststruct%tracer_nest_timestep + 1
       end if
                        call timing_on('COMM_TOTAL')
                            call timing_on('COMM_TRACER')
      call complete_group_halo_update(q_pack, domain)
                           call timing_off('COMM_TRACER')
                       call timing_off('COMM_TOTAL')
	    
      if (gridstruct%nested) then
            do iq=1,nq
                 call nested_grid_BC_apply_intT(q(isd:ied,jsd:jed,:,iq), &
                      0, 0, npx, npy, npz, bd, &
                      real(neststruct%tracer_nest_timestep)+real(nsplt*k_split), real(nsplt*k_split), &
                 neststruct%q_BC(iq), bctype=neststruct%nestbctype  )
           enddo
      endif


!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,jsd,jed,npz,dp1,mfx,mfy,rarea,nq, &
!$OMP                                  area,xfx,yfx,q,cx,cy,npx,npy,hord,gridstruct,bd,it,nsplt,nord_tr,trdm,lim_fac) &
!$OMP                          private(dp2, ra_x, ra_y, fx, fy)
      do k=1,npz

         do j=js,je
            do i=is,ie
               dp2(i,j) = dp1(i,j,k) + ((mfx(i,j,k)-mfx(i+1,j,k))+(mfy(i,j,k)-mfy(i,j+1,k)))*rarea(i,j)
            enddo
         enddo

         do j=jsd,jed
            do i=is,ie
               ra_x(i,j) = area(i,j) + (xfx(i,j,k) - xfx(i+1,j,k))
            enddo
         enddo
         do j=js,je
            do i=isd,ied
               ra_y(i,j) = area(i,j) + (yfx(i,j,k) - yfx(i,j+1,k))
            enddo
         enddo

         do iq=1,nq
         if ( it==1 .and. trdm>1.e-4 ) then
            call fv_tp_2d(q(isd,jsd,k,iq), cx(is,jsd,k), cy(isd,js,k), &
                          npx, npy, hord, fx, fy, xfx(is,jsd,k), yfx(isd,js,k), &
                          gridstruct, bd, ra_x, ra_y, lim_fac, mfx=mfx(is,js,k), mfy=mfy(is,js,k),   &
                          mass=dp1(isd,jsd,k), nord=nord_tr, damp_c=trdm)
         else
            call fv_tp_2d(q(isd,jsd,k,iq), cx(is,jsd,k), cy(isd,js,k), &
                          npx, npy, hord, fx, fy, xfx(is,jsd,k), yfx(isd,js,k), &
                          gridstruct, bd, ra_x, ra_y, lim_fac, mfx=mfx(is,js,k), mfy=mfy(is,js,k))
         endif
            do j=js,je
               do i=is,ie
                  q(i,j,k,iq) = ( q(i,j,k,iq)*dp1(i,j,k) + &
                                ((fx(i,j)-fx(i+1,j))+(fy(i,j)-fy(i,j+1)))*rarea(i,j) )/dp2(i,j)
               enddo
               enddo
          enddo
      enddo ! npz

      if ( it /= nsplt ) then
                      call timing_on('COMM_TOTAL')
                          call timing_on('COMM_TRACER')
           call start_group_halo_update(q_pack, q, domain)
                          call timing_off('COMM_TRACER')
                      call timing_off('COMM_TOTAL')
      endif
           !Apply nested-grid BCs
           if ( gridstruct%nested ) then
              do iq=1,nq


                 call nested_grid_BC_apply_intT(q(isd:ied,jsd:jed,:,iq), &
                      0, 0, npx, npy, npz, bd, &
                      real(neststruct%tracer_nest_timestep), real(nsplt*k_split), &
                 neststruct%q_BC(iq), bctype=neststruct%nestbctype  )

              end do
           end if


   enddo  ! nsplt

   if ( id_divg > 0 ) then
        rdt = 1./(frac*dt)

!$OMP parallel do default(none) shared(is,ie,js,je,npz,dp1,xfx,yfx,rarea,rdt)
        do k=1,npz
        do j=js,je
           do i=is,ie
              dp1(i,j,k) = ((xfx(i+1,j,k)-xfx(i,j,k)) + (yfx(i,j+1,k)-yfx(i,j,k)))*rarea(i,j)*rdt
           enddo
        enddo
        enddo
   endif

 end subroutine tracer_2d_nested

subroutine offline_tracer_advection(q, pleB, pleA, mfx, mfy, cx, cy, &
                                    gridstruct, flagstruct, bd, domain, &
                                    npx, npy, npz,   &
                                    nq, dt, q_split, pleadv) ! GCHP: added argument pleadv

      use fv_mapz_mod,        only: mapn_tracer
      use fv_fill_mod,        only: fillz

      integer, intent(IN) :: npx
      integer, intent(IN) :: npy
      integer, intent(IN) :: npz
      integer, intent(IN) :: nq    ! number of tracers to be advected
      integer, intent(IN) :: q_split
      type(fv_grid_bounds_type), intent(IN   ) :: bd
      type(fv_flags_type), intent(INOUT) :: flagstruct
      type(fv_grid_type), intent(IN), target :: gridstruct
      type(domain2D), intent(INOUT) :: domain

      real, intent(IN   ) :: dt
      real, intent(IN   ) ::pleB(bd%is:bd%ie,bd%js:bd%je,npz+1)      ! DELP before dyn_core
      real, intent(INOUT) ::pleA(bd%is:bd%ie,bd%js:bd%je,npz+1)      ! DELP after dyn_core
      real, intent(IN   ) ::  cx(bd%is:bd%ie,bd%js:bd%je,npz)        ! Courant Number X-Dir
      real, intent(IN   ) ::  cy(bd%is:bd%ie,bd%js:bd%je,npz)        ! Courant Number Y-Dir
      real, intent(IN   ) :: mfx(bd%is:bd%ie,bd%js:bd%je,npz)        ! Mass Flux X-Dir
      real, intent(IN   ) :: mfy(bd%is:bd%ie,bd%js:bd%je,npz)        ! Mass Flux Y-Dir
      real, intent(INOUT) ::   q(bd%is:bd%ie,bd%js:bd%je,npz,nq)     ! Tracers

      ! GCHP: add optional argument pleadv
      real, optional, intent(OUT) :: pleadv(bd%is:bd%ie,bd%js:bd%je,npz+1) ! DELP after adv_core

! Local Arrays
      real ::   xL(bd%isd:bd%ied+1,bd%jsd:bd%jed  ,npz)  ! X-Dir for MPP Updates
      real ::   yL(bd%isd:bd%ied  ,bd%jsd:bd%jed+1,npz)  ! Y-Dir for MPP Updates
      real ::  cxL(bd%is :bd%ie +1,bd%jsd:bd%jed  ,npz)  ! Courant Number X-Dir
      real ::  cyL(bd%isd:bd%ied  ,bd%js :bd%je +1,npz)  ! Courant Number Y-Dir
      real :: mfxL(bd%is :bd%ie +1,bd%js :bd%je   ,npz)  ! Mass Flux X-Dir
      real :: mfyL(bd%is :bd%ie   ,bd%js :bd%je +1,npz)  ! Mass Flux Y-Dir
      real ::  dpL(bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)  ! Pressure Thickness
      real ::  dpA(bd%is :bd%ie   ,bd%js :bd%je   ,npz)  ! Pressure Thickness
! Local Tracer Arrays
      real ::   q3(bd%isd:bd%ied,bd%jsd:bd%jed, npz,nq)! Ghosted 3D Tracers
      real ::   q2(bd%is :bd%ie ,               npz   )! 2D Tmp
      integer :: kord_tracers(nq)
! Local Buffer Arrarys
      real :: wbuffer(bd%js:bd%je,npz)
      real :: sbuffer(bd%is:bd%ie,npz)
      real :: ebuffer(bd%js:bd%je,npz)
      real :: nbuffer(bd%is:bd%ie,npz)
! Local Remap Arrays
      real  pe1(bd%is:bd%ie,npz+1)
      real  pe2(bd%is:bd%ie,npz+1)
      real  dp1(bd%is:bd%ie,npz)
      real  dp2(bd%is:bd%ie,npz)

! Local indices
      integer     :: i,j,k,n,iq

      real :: scalingFactor

      ! GCHP: define kStart
      integer :: kStart

      type(group_halo_update_type), save :: i_pack

      integer :: is,  ie,  js,  je
      integer :: isd, ied, jsd, jed

      is  = bd%is
      ie  = bd%ie
      js  = bd%js
      je  = bd%je
      isd = bd%isd
      ied = bd%ied
      jsd = bd%jsd
      jed = bd%jed

! Time-step
! Fill CX/CY C-Grid boundaries and update ghost regions
    xL(is:ie,js:je,:) = cx(:,:,:)
    yL(is:ie,js:je,:) = cy(:,:,:)
    call mpp_get_boundary(xL, yL, domain, &
                          wbufferx=wbuffer, ebufferx=ebuffer, &
                          sbuffery=sbuffer, nbuffery=nbuffer, &
                          gridtype=CGRID_NE )
    xL(ie+1,js:je,:) = ebuffer
    yL(is:ie,je+1,:) = nbuffer
    call mpp_update_domains( xL, yL, domain, gridtype=CGRID_NE, complete=.true.)
    cxL(is:ie+1,jsd:jed,:) = xL(is:ie+1,jsd:jed,:)
    cyL(isd:ied,js:je+1,:) = yL(isd:ied,js:je+1,:)

! Fill MFX/MFY C-Grid boundaries
    xL(is:ie,js:je,:) = mfx(:,:,:)
    yL(is:ie,js:je,:) = mfy(:,:,:)
    call mpp_get_boundary(xL, yL, domain, &
                          wbufferx=wbuffer, ebufferx=ebuffer, &
                          sbuffery=sbuffer, nbuffery=nbuffer, &
                          gridtype=CGRID_NE )
    xL(ie+1,js:je,:) = ebuffer
    yL(is:ie,je+1,:) = nbuffer
    mfxL(is:ie+1,js:je,:) = xL(is:ie+1,js:je,:)
    mfyL(is:ie,js:je+1,:) = yL(is:ie,js:je+1,:)

! Fill local tracers and pressure thickness
    dpL(bd%is:bd%ie,bd%js:bd%je,:) = pleB(:,:,2:npz+1) - pleB(:,:,1:npz)
    q3(is:ie,js:je,:,:) = q(is:ie,js:je,:,:)
    call start_group_halo_update(i_pack, q3, domain)

    if ( flagstruct%z_tracer .and. q_split==0 ) then
         call tracer_2d_1L(q3, dpL, mfxL, mfyL, cxL, cyL, &
                         gridstruct, bd, domain, npx, npy, npz, nq,    &
                         flagstruct%hord_tr, dt, 0, i_pack, &
                         flagstruct%nord_tr, flagstruct%trdm2, flagstruct%lim_fac, dpA=dpA)
    else
         call tracer_2d(q3, dpL, mfxL, mfyL, cxL, cyL, gridstruct, bd, domain, npx, npy, npz, nq,    &
                        flagstruct%hord_tr, q_split, dt, 0, i_pack, &
                        flagstruct%nord_tr, flagstruct%trdm2, flagstruct%lim_fac, dpA=dpA)
    endif

    ! ewl todo: ptop is no longer passed to offline_tracer_advection
    ! Need to add it back in, and define within AdvCore_GridCompMod
    ! For now, hard-code here (temporary!) to what is in DynCore_GridCompMod.F90
    ! NOTE: ptop might now be the same as pleB(:,:,1), or would it be PLEA(:,:,1)?
    ! GCHP: Build post-advection pressure edges if requested
    if (present(pleadv)) then
       pleadv(:,:,1) = 25499.234876157 ! Model top pressure in Pa (ptop)
       do k=2,npz+1
          pleadv(:,:,k) = pleadv(:,:,k-1) + dpA(:,:,k-1)
       end do
    end if

!------------------------------------------------------------------
! Re-Map constituents
!------------------------------------------------------------------

    ! ewl: bk is no longer passed in. Will need an alternate way to determine min k
    ! in troposphere. Or add it back as an argument. If kstart is always the same,
    ! then why not pass that as an argument instead?
    ! For now, hard-code to 41 based on GEOS 72-level grid BP values in GEOS-Chem
    ! Note that k is positive down, so bk is reverse of GEOS-Chem BP
    ! GCHP: rescale trop only
    kStart = 41
    !kStart = 0
    !if (rescale_trop_only) then
    !   do k=1,npz
    !      if ((kStart.le.0).and.(bk(k+1).gt.0.d0)) then
    !         kStart=k
    !         exit
    !      end if
    !   end do
    !end if

       kord_tracers = flagstruct%kord_tr
!$OMP parallel do default(none) shared(is,ie,isd,ied,js,je,jsd,jed,npz,nq, &
!$OMP                                  pleB,dpA,pleA,q3,flagstruct,kord_tracers) &
!$OMP                          private(i,j,k,pe1,dp1,pe2,dp2)
       do j=js,je
        ! pressures mapping from (dpA is new delp after tracer_2d)
          pe1(:,1) = pleB(:,j,1)
          do k=2,npz+1
            pe1(:,k) = pe1(:,k-1) + dpA(:,j,k-1)
          enddo
          do k=1,npz
             dp1(:,k) = pe1(:,k+1) - pe1(:,k)
          enddo
        ! pressures mapping to
          pe2 = pleA(:,j,:)
          do k=1,npz
             dp2(:,k) = pe2(:,k+1) - pe2(:,k)
          enddo
          call mapn_tracer(nq, npz, pe1, pe2, q3, dp2, kord_tracers, j,     &
                           is, ie, isd, ied, jsd, jed, 0., flagstruct%fill)
       enddo

       ! Rescale tracers based on pleA at destination timestep
       !------------------------------------------------------
       ! GCHP: rescale trop only
       if (rescale_trop_only) then
          do iq=1,nq
               ! Scale tracers
               !---------------
               scalingFactor = calcScalingFactorTrop(q(is:ie,js:je,1:npz,iq), q3(is:ie,js:je,1:npz,iq), pleB, pleA, &
                                 npz, domain, gridstruct, flagstruct, kstart, bd)
               q(is:ie,js:je,1:(kStart-1),iq) = q3(is:ie,js:je,1:(kStart-1),iq)
               q(is:ie,js:je,kStart:npz,iq) = q3(is:ie,js:je,kStart:npz,iq) * scalingFactor
               if (show_scaling.and.is_master()) then
                  write(*,'(a,I3,a,E16.5E4)') ' FV3 scaling error for tracer ', iq, ': ', &
                      (scalingFactor - 1.0)
               end if
          enddo
       else
          do iq=1,nq
               ! Scale tracers
               !---------------
               scalingFactor = calcScalingFactor(q(is:ie,js:je,1:npz,iq), q3(is:ie,js:je,1:npz,iq), pleB, pleA, &
                                 npz, domain, gridstruct, flagstruct, bd)
               q(is:ie,js:je,1:npz,iq) = q3(is:ie,js:je,1:npz,iq) * scalingFactor
               if (show_scaling.and.is_master()) then
                  write(*,'(a,I3,a,E16.5E4)') ' FV3 scaling error for tracer ', iq, ': ', &
                      (scalingFactor - 1.0)
               end if
          enddo
       endif

end subroutine offline_tracer_advection

!------------------------------------------------------------------------------------

         ! Added for GCHP to limit scaling to troposphere
         ! This is a copy of function calcScalingFactor except where differences are noted
         function calcScalingFactorTrop(q1, q2, ple1, ple2, npz, domain, gridstruct, flagstruct, kstart, bd) result(scaling)
         integer, intent(in) :: npz
         integer, intent(in) :: kstart
         type(fv_grid_bounds_type), intent(IN   ) :: bd
         real, intent(in) :: q1(bd%is:bd%ie,bd%js:bd%je,npz)
         real, intent(in) :: q2(bd%is:bd%ie,bd%js:bd%je,npz)
         real, intent(in) :: ple1(bd%is:bd%ie,bd%js:bd%je,npz+1)
         real, intent(in) :: ple2(bd%is:bd%ie,bd%js:bd%je,npz+1)
         type(domain2D), intent(INOUT) :: domain
         type(fv_grid_type), intent(IN   ) :: gridstruct
         type(fv_flags_type), intent(INOUT) :: flagstruct
         integer :: sflag
         real :: scaling

         integer :: k
         real(REAL8)     :: qsum1(bd%is:bd%ie,bd%js:bd%je)
         real(REAL8)     :: qsum2(bd%is:bd%ie,bd%js:bd%je)
         real(REAL8)     :: qsum3(bd%is:bd%ie,bd%js:bd%je)
         real(REAL8)     :: tempsum(bd%is:bd%ie,bd%js:bd%je)
         real(REAL8)     :: globalSums(3)
         real(REAL8), parameter :: TINY_DENOMINATOR = tiny(1.d0)
         real(REAL8)     :: scalingR8

         !-------
         ! Compute partial sum on local array first to minimize communication.
         ! This algorithm will not be strongly repdroducible under changes do domain
         ! decomposition, but uses far less communication bandwidth (and memory BW)
         ! then the preceding implementation.
         !-------
         qsum1(:,:) = 0.d0
         qsum2(:,:) = 0.d0
         qsum3(:,:) = 0.d0
         tempSum(:,:) = 0.d0
         do k=1,npz
            qsum1(:,:) = qsum1(:,:) + q1(:,:,k) * (ple1(:,:,k+1)-ple1(:,:,k))
            tempSum(:,:) = q2(:,:,k) * (ple2(:,:,k+1)-ple2(:,:,k))
            qsum2(:,:) = qsum2(:,:) + tempSum
            if (k.lt.kStart) then
               ! In stratosphere
               qsum3(:,:) = qsum3(:,:) + tempSum
            endif
         enddo

         ! Get total pre-advection mass, total post-advection mass, and total
         ! stratospheric mass after advection

         ! numerator: total mass based on advected air
         globalSums(1) = g_sum_r8(domain, qsum1, bd%is,bd%ie, bd%js,bd%je, 0, &
                                  gridstruct%area_64(bd%is:bd%ie,bd%js:bd%je), 1)

         ! denominator: total mass based on met field pressures after advection
         globalSums(2) = g_sum_r8(domain, qsum2, bd%is,bd%ie, bd%js,bd%je, 0, &
                                  gridstruct%area_64(bd%is:bd%ie,bd%js:bd%je), 1)

         ! Total stratospheric mass after advection
         globalSums(3) = g_sum_r8(domain, qsum3, bd%is,bd%ie, bd%js,bd%je, 0, &
                                  gridstruct%area_64(bd%is:bd%ie,bd%js:bd%je), 1)

         ! Subtract stratospheric mass
         globalSums(1) = globalSums(1) - globalSums(3)
         globalSums(2) = globalSums(2) - globalSums(3)

         if (globalSums(2) > TINY_DENOMINATOR) then
            scalingR8 =  globalSums(1) / globalSums(2)
            !#################################################################
            ! This line was added to ensure strong reproducibility of the code
            ! SDE 2017-02-17: Disabled as it compromises mass conservation
            ! EWL TODO: Is this still true?
            !#################################################################
            !scaling = REAL(scalingR8, KIND=kind(1.00))
         else
            scaling = 1.0
         end if

         end function calcScalingFactorTrop

!------------------------------------------------------------------------------------

         function calcScalingFactor(q1, q2, ple1, ple2, npz, domain, gridstruct, flagstruct, bd) result(scaling)
         integer, intent(in) :: npz
         type(fv_grid_bounds_type), intent(IN   ) :: bd
         real, intent(in) :: q1(bd%is:bd%ie,bd%js:bd%je,npz)
         real, intent(in) :: q2(bd%is:bd%ie,bd%js:bd%je,npz)
         real, intent(in) :: ple1(bd%is:bd%ie,bd%js:bd%je,npz+1)
         real, intent(in) :: ple2(bd%is:bd%ie,bd%js:bd%je,npz+1)
         type(domain2D), intent(INOUT) :: domain
         type(fv_grid_type), intent(IN   ) :: gridstruct
         type(fv_flags_type), intent(INOUT) :: flagstruct
         integer :: sflag
         real :: scaling

         integer :: k
         real(REAL8)     :: qsum1(bd%is:bd%ie,bd%js:bd%je)
         real(REAL8)     :: qsum2(bd%is:bd%ie,bd%js:bd%je)
         real(REAL8)     :: globalSums(2)
         real(REAL8), parameter :: TINY_DENOMINATOR = tiny(1.d0)
         real(REAL8)     :: scalingR8
         !-------
         ! Compute partial sum on local array first to minimize communication.
         ! This algorithm will not be strongly repdroducible under changes do domain
         ! decomposition, but uses far less communication bandwidth (and memory BW)
         ! then the preceding implementation.
         !-------
         qsum1(:,:) = 0.d0
         qsum2(:,:) = 0.d0
         do k=1,npz
            qsum1(:,:) = qsum1(:,:) + q1(:,:,k) * (ple1(:,:,k+1)-ple1(:,:,k))
            qsum2(:,:) = qsum2(:,:) + q2(:,:,k) * (ple2(:,:,k+1)-ple2(:,:,k))
         enddo

         ! numerator
         globalSums(1) = g_sum_r8(domain, qsum1, bd%is,bd%ie, bd%js,bd%je, 0, &
                                  gridstruct%area_64(bd%is:bd%ie,bd%js:bd%je), 1)
         ! denominator
         globalSums(2) = g_sum_r8(domain, qsum2, bd%is,bd%ie, bd%js,bd%je, 0, &
                                  gridstruct%area_64(bd%is:bd%ie,bd%js:bd%je), 1)

         if (globalSums(2) > TINY_DENOMINATOR) then
            scalingR8 =  globalSums(1) / globalSums(2)
            !#################################################################
            ! This line was added to ensure strong reproducibility of the code
            !#################################################################
            scaling = REAL(scalingR8, KIND=kind(1.00))
         else
            scaling = 1.0
         end if

         end function calcScalingFactor

end module fv_tracer2d_mod
