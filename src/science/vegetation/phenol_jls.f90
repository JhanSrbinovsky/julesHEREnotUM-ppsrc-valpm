
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine PHENOL -------------------------------------------------
!
!
! Purpose :  Parametrizes leaf phenological changes and updates the
!            leaf area index and the leaf turnover rate.
!
! -------------------------------------------------------------------
SUBROUTINE phenol (land_pts,veg_pts,veg_index,n,g_leaf,ht         &
,                  dtime_phen,g_leaf_phen,lai)

USE pftparm
USE trif

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

INTEGER                                                           &
 land_pts                                                         &
                            ! IN Total number of land points.
,veg_pts                                                          &
                            ! IN Number of vegetated points.
,veg_index(land_pts)                                              &
                            ! IN Index of vegetated points
!                                 !    on the land grid.
,n                          ! IN Plant functional type.

REAL                                                              &
 g_leaf(land_pts)                                                 &
                            ! IN Rate of leaf turnover (/360days).
,ht(land_pts)                                                     &
                            ! IN Canopy height (m).
,dtime_phen                                                       &
                            ! IN Timestep (years).
,g_leaf_phen(land_pts)                                            &
                            ! OUT Rate of leaf turnover
!                                 !     including leaf phenology
!                                 !     (/360days).
,lai(land_pts)                                                    &
                            ! INOUT Leaf area index.
,dphen                                                            &
                            ! WORK Increment to phenological
!                                 !      state.
,lai_bal(land_pts)                                                &
                            ! WORK Balanced growth LAI.
,phen(land_pts)             ! WORK Phenological state.

INTEGER                                                           &
 j,l                        ! Loop counters

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


!-----------------------------------------------------------------------
! Diagnose the phenological state
!-----------------------------------------------------------------------
IF (lhook) CALL dr_hook('PHENOL',zhook_in,zhook_handle)
DO j=1,veg_pts
  l = veg_index(j)
  lai_bal(l) = (a_ws(n)*eta_sl(n)*ht(l)                           &
               /a_wl(n))**(1.0/(b_wl(n)-1))
  phen(l) = lai(l)/lai_bal(l)
END DO

!-----------------------------------------------------------------------
! Update the phenological state and output the leaf turnover rate in
! terms of the balanced growth LAI
!-----------------------------------------------------------------------
DO j=1,veg_pts
  l = veg_index(j)

  IF (g_leaf(l) >  2*g_leaf_0(n)) THEN
    dphen = -dtime_phen*g_grow(n)
    dphen = MAX(dphen,(0.01-phen(l)))
    g_leaf_phen(l) = -dphen/dtime_phen
  ELSE
    dphen = dtime_phen*g_grow(n)*(1.0-phen(l))
    dphen = MIN(dphen,(1.0-phen(l)))
    g_leaf_phen(l) = phen(l)*g_leaf(l)
  END IF

!-----------------------------------------------------------------------
! Update the leaf area index
!-----------------------------------------------------------------------
  phen(l) = phen(l) + dphen
  lai(l) = phen(l)*lai_bal(l)

END DO

IF (lhook) CALL dr_hook('PHENOL',zhook_out,zhook_handle)
RETURN

END SUBROUTINE phenol
