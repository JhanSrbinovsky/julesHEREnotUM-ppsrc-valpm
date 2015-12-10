! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  SUBROUTINE FRUNOFF------------------------------------------------

!  PURPOSE : TO CALCULATE SURFACE RUNOFF

!  SUITABLE FOR SINGLE COLUMN MODEL USE

!  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER NO 25
!                  SECTION (3B(II)), EQN(P252.14)
!--------------------------------------------------------------------

!    ARGUMENTS---------------------------------------------------------

SUBROUTINE frunoff (                                              &
 npnts,tile_pts,tile_index,area,                                  &
 can_cpy,can_wcnt,infil,r,frac,timestep,                          &
 surf_roff                                                        &
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
,can_wcnt(npnts)                                                  &
                      ! IN Canopy water content (kg/m2).
,infil(npnts)                                                     &
                      ! IN Infiltration rate (kg/m2/s).
,r(npnts)                                                         &
                      ! IN Water fall rate (kg/m2/s).
,frac(npnts)                                                      &
                      ! IN Tile fraction.
,timestep             ! IN Timestep (s).

REAL                                                              &
 surf_roff(npnts)     ! OUT Cummulative surface runoff (kg/m2/s).

!  Workspace --------------------------------------------------------
REAL                                                              &
 aexp                                                             &
                      ! Used in the calculation of exponential
,aexp1                                                            &
                      ! terms in the surface runoff formula.
,aexp2                                                            &
                      !
,cm                                                               &
                      ! (CAN_CPY - CAN_WCNT)/TIMESTEP
,can_ratio                                                        &
                      ! CAN_WCNT / CAN_CPY
,runoff                                                           &
                      ! Local runoff.
,smallestp            ! Smallest +ve number that can be represented

INTEGER                                                           &
 i                                                                &
                      ! Land point index.
,j                    ! Counter for loop over tile points.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('FRUNOFF',zhook_in,zhook_handle)

smallestp = TINY(1.0)



!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, runoff, aexp, aexp1, aexp2,&
!$OMP& can_ratio, cm)
!CDIR NODEP
DO j=1,tile_pts
  i = tile_index(j)
  runoff = 0.
  IF ( r(i) > EPSILON(r(i)) ) THEN
    IF ( infil(i)*timestep <= can_wcnt(i)                         &
                                   .AND. can_cpy(i) >  0.0 ) THEN
! Infiltration in timestep < or = canopy water content
      aexp = area*can_cpy(i)/r(i)
      IF ( can_wcnt(i) > EPSILON(can_wcnt(i)) ) THEN
        aexp1 = EXP( -aexp*infil(i)/can_wcnt(i))
      ELSE
        aexp1 = 0.0
      END IF
      aexp2 = EXP( -aexp/timestep)
      can_ratio = can_wcnt(i)/can_cpy(i)
      can_ratio = MIN(can_ratio,1.0)
      runoff = r(i) * ( can_ratio*aexp1 +                        &
                                       (1. - can_ratio)*aexp2 )
!                                                        ... P252.14A
    ELSE
! Infiltration in timestep > canopy water content
      cm = (can_cpy(i)-can_wcnt(i))/timestep
      cm = MAX(cm,0.0)
      ! Only compute AEXP if will not generate an underflow error
      IF ( area*(infil(i)+cm)/r(i) < - LOG(smallestp) ) THEN
         aexp = EXP( -area*(infil(i)+cm)/r(i))
      ELSE
         aexp = 0.0
      END IF
      runoff = r(i)*aexp                    !     ... P252.14B
    END IF
  END IF
  surf_roff(i) = surf_roff(i) + frac(i)*runoff
END DO
!$OMP END PARALLEL DO

IF (lhook) CALL dr_hook('FRUNOFF',zhook_out,zhook_handle)
RETURN
END SUBROUTINE frunoff
