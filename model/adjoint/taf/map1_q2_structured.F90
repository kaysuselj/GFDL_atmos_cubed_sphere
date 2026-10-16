!>@brief GOTO-free form of map1_q2 (fv_mapz.F90), for the TAF input only.
!>@details The original jumps out of the nested search loops with GOTO 555 / GOTO 123, which makes
!! its control-flow graph irreducible, and TAF's data-dependence analysis then loses the q1 -> q2
!! dependence. This version uses flags instead, with the same arithmetic in the same order, so the
!! results are bit-identical (checked by build/test_map1_q2 on the trajectory data).
!! run_taf.sh puts it in place of map1_q2 in the preprocessed copy of fv_mapz.F90 that is sent to
!! TAF; fv_mapz.F90 itself is not changed.
 subroutine map1_q2(km,   pe1,   q1,            &
                    kn,   pe2,   q2,   dp2,     &
                    i1,   i2,    iv,   kord, j, &
                    ibeg, iend, jbeg, jend, q_min )


! INPUT PARAMETERS:
      integer, intent(in) :: j
      integer, intent(in) :: i1, i2
      integer, intent(in) :: ibeg, iend, jbeg, jend
      integer, intent(in) :: iv                !< Mode: 0 ==  constituents 1 == ???
      integer, intent(in) :: kord
      integer, intent(in) :: km                !< Original vertical dimension
      integer, intent(in) :: kn                !< Target vertical dimension

      real, intent(in) ::  pe1(i1:i2,km+1)     !< pressure at layer edges from model top to bottom surface in the original vertical coordinate
      real, intent(in) ::  pe2(i1:i2,kn+1)     !< pressure at layer edges from model top to bottom surface in the new vertical coordinate
      real, intent(in) ::  q1(ibeg:iend,jbeg:jend,km) !< Field input
      real, intent(in) ::  dp2(i1:i2,kn)
      real, intent(in) ::  q_min
! INPUT/OUTPUT PARAMETERS:
      real, intent(inout):: q2(i1:i2,kn) !< Field output
! LOCAL VARIABLES:
      real   qs(i1:i2)
      real   dp1(i1:i2,km)
      real   q4(4,i1:i2,km)
      real   pl, pr, qsum, dp, esl

      integer i, k, l, m, k0
      logical found, within, mdone

      do k=1,km
         do i=i1,i2
             dp1(i,k) = pe1(i,k+1) - pe1(i,k)
            q4(1,i,k) = q1(i,j,k)
         enddo
      enddo

! Compute vertical subgrid distribution
   if ( kord >7 ) then
        call  scalar_profile( qs, q4, dp1, km, i1, i2, iv, kord, q_min )
   else
        call ppm_profile( q4, dp1, km, i1, i2, iv, kord )
   endif

! Mapping
      do i=i1,i2
         k0 = 1
         do k=1,kn
            found  = .false.
            within = .false.
            do l=k0,km
! locate the top edge: pe2(i,k)
               if ( .not. found ) then
               if(pe2(i,k) >= pe1(i,l) .and. pe2(i,k) <= pe1(i,l+1)) then
                  found = .true.
                  pl = (pe2(i,k)-pe1(i,l)) / dp1(i,l)
                  if(pe2(i,k+1) <= pe1(i,l+1)) then
! entire new grid is within the original grid
                     pr = (pe2(i,k+1)-pe1(i,l)) / dp1(i,l)
                     q2(i,k) = q4(2,i,l) + 0.5*(q4(4,i,l)+q4(3,i,l)-q4(2,i,l))  &
                                *(pr+pl)-q4(4,i,l)*r3*(pr*(pr+pl)+pl**2)
                     within = .true.
                  else
! Fractional area...
                     qsum = (pe1(i,l+1)-pe2(i,k))*(q4(2,i,l)+0.5*(q4(4,i,l)+   &
                             q4(3,i,l)-q4(2,i,l))*(1.+pl)-q4(4,i,l)*           &
                              (r3*(1.+pl*(1.+pl))))
                     mdone = .false.
                     do m=l+1,km
                        if ( .not. mdone ) then
! locate the bottom edge: pe2(i,k+1)
                        if(pe2(i,k+1) > pe1(i,m+1) ) then
                                                   ! Whole layer..
                           qsum = qsum + dp1(i,m)*q4(1,i,m)
                        else
                           dp = pe2(i,k+1)-pe1(i,m)
                           esl = dp / dp1(i,m)
                           qsum = qsum + dp*(q4(2,i,m)+0.5*esl*               &
                               (q4(3,i,m)-q4(2,i,m)+q4(4,i,m)*(1.-r23*esl)))
                           k0 = m
                           mdone = .true.
                        endif
                        endif
                     enddo
                  endif
                  if ( within ) k0 = l
               endif
               endif
            enddo
            if ( .not. within ) q2(i,k) = qsum / dp2(i,k)
         enddo
      enddo

 end subroutine map1_q2
