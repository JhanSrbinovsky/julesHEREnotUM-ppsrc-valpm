! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  SUBROUTINE SIEVE--------------------------------------------------

!  PURPOSE : TO CALCULATE THE THROUGHFALL OF WATER FALLING
!            THROUGH THE SURFACE CANOPY

!  SUITABLE FOR SINGLE COLUMN MODEL USE

!  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER NO 25
!                  SECTION (3B(II)), EQN(P252.9)

!--------------------------------------------------------------------

!    ARGUMENTS---------------------------------------------------------

SUBROUTINE sieve (                                                &
 npnts,tile_pts,tile_index,area,can_cpy,r,frac,timestep,          &
 can_wcnt,tot_tfall                                               &
 )


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

INTEGER                                                           &
 npnts                                                            &
                      ! IN Total number of land points.
,tile_pts                                                         &
                      ! IN Number of tile points.
,tile_index(npnts)    ! IN Index of tile points.

REAL                                                              &
 area                                                             &
                      ! IN Fractional area of gridbox over which
                      !    water falls (%).
,can_cpy(npnts)                                                   &
                      ! IN Canopy capacity (kg/m2).
,r(npnts)                                                         &
                      ! IN Water fall rate (kg/m2/s).
,frac(npnts)                                                      &
                      ! IN Tile fraction.
,timestep             ! IN Timestep (s).

REAL                                                              &
 can_wcnt(npnts)                                                  &
                      ! INOUT Canopy water content (kg/m2).
,tot_tfall(npnts)     ! INOUT Cummulative canopy throughfall
!                           !       (kg/m2/s).

!  Workspace --------------------------------------------------------
REAL                                                              &
 aexp                                                             &
                      ! Used in calculation of exponential
                      ! in throughfall formula.
,can_ratio                                                        &
                      ! CAN_WCNT / CAN_CPY
,tfall(npnts)                                                     &
                      ! Local throughfall (kg/m2/s).
,smallp                                                           &
                      ! Small positive number << 1
,smallestp            ! Smallest +ve real which can be represented

INTEGER                                                           &
 i                                                                &
                      ! Land point index.
,j                    ! Counter for loop over tile points.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('SIEVE',zhook_in,zhook_handle)

smallestp = TINY(1.0)
smallp = EPSILON( r )

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j, i, aexp, can_ratio)
DO j=1,tile_pts
  i = tile_index(j)
  IF (can_cpy(i)  >   0.0 .AND. r(i)  >  smallp ) THEN
     aexp = area*can_cpy(i)/(r(i)*timestep)
     ! Only calculate if AEXP is small enough to avoid underflow
     IF (aexp < -LOG(smallestp)) THEN
        aexp = EXP(-aexp)
     ELSE
        aexp = 0.0
     END IF

     can_ratio = can_wcnt(i) / can_cpy(i)

     can_ratio = MIN(can_ratio,1.0)
     tfall(i) = r(i) * ((1.0-can_ratio)*aexp + can_ratio)
  ELSE
     tfall(i) = r(i)
  END IF

  can_wcnt(i) = can_wcnt(i) + (r(i) - tfall(i))*timestep
  tot_tfall(i) = tot_tfall(i) + frac(i)*tfall(i)
END DO
!$OMP END PARALLEL DO

IF (lhook) CALL dr_hook('SIEVE',zhook_out,zhook_handle)
RETURN
END SUBROUTINE sieve
