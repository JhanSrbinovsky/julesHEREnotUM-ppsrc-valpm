! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Purpose:
! Calculates the leaf turnover rate as a function of temperature and
! soil water availability

!***********************************************************************
SUBROUTINE leaf_lit (land_pts,veg_pts,veg_index,n,fsmc,tstar      &
,                    g_leaf)

USE pftparm

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
 fsmc(land_pts)                                                   &
                            ! IN Soil moisture availability
!                                 !    factor.
,tstar(land_pts)                                                  &
                            ! IN Surface temperature (K).
,g_leaf(land_pts)                                                 &
                            ! OUT Rate of leaf turnover
!                                 !     (/360days).
,fm,ft                      ! WORK Soil moisture and leaf
!                                        temperature amplifiers of
!                                        leaf turnover.

INTEGER                                                           &
 j,l                        ! Loop counters

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------
! Calculate the leaf turnover rate
!-----------------------------------------------------------------------
IF (lhook) CALL dr_hook('LEAF_LIT',zhook_in,zhook_handle)
DO j=1,veg_pts
  l = veg_index(j)

  ft = 1.0
  fm = 1.0
  IF (tstar(l)  <   tleaf_of(n)) THEN
    ft = 1.0 + dgl_dt(n)*(tleaf_of(n)-tstar(l))
  ELSE IF (fsmc(l)  <   fsmc_of(n)) THEN
    fm = 1.0 + dgl_dm(n)*(fsmc_of(n)-fsmc(l))
  END IF

  g_leaf(l) = g_leaf_0(n)*ft*fm

END DO

IF (lhook) CALL dr_hook('LEAF_LIT',zhook_out,zhook_handle)
RETURN

END SUBROUTINE leaf_lit
