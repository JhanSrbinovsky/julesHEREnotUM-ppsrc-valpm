! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine CANCAP ------------------------------------------------

! Purpose : Calculate the heat capacity of a given PFT from its LAI

! -----------------------------------------------------------------
SUBROUTINE cancap (land_pts,veg_pts,veg_index,can_model,ft        &
,                  ht,lai,canhc,vfrac)

USE trif
USE pftparm
USE surf_param, ONLY : hleaf,hwood

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

 INTEGER                                                          &
 land_pts                                                         &
                            ! IN Total number of land points.
,veg_pts                                                          &
                            ! IN Number of vegetated points.
,veg_index(land_pts)                                              &
                            ! IN Index of vegetated points.
,can_model                  ! IN Swith for thermal vegetation
!                                 !    canopy

INTEGER                                                           &
 ft                         ! IN Plant functional type.

REAL                                                              &
 ht(land_pts)                                                     &
                            ! IN Vegetation height (m).
,lai(land_pts)                                                    &
                            ! IN Leaf area index.
,canhc(land_pts)                                                  &
                            ! OUT Areal heat capacity of
!                                 !     vegetation canopy (J/K/m2).
,vfrac(land_pts)            ! OUT Fractional canopy coverage.

REAL                                                              &
 lai_bal(land_pts)                                                &
                            ! WORK Leaf area index in balanced
!                                 !      growth state.
,leaf(land_pts)                                                   &
                            ! WORK Leaf biomass (kg C/m2).
,wood(land_pts)             ! WORK Woody biomass (kg C/m2).

INTEGER                                                           &
 j,l                        ! WORK Loop counters.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


IF (lhook) CALL dr_hook('CANCAP',zhook_in,zhook_handle)
DO j=1,veg_pts
  l = veg_index(j)
  canhc(l) = 0.
  vfrac(l) = 0.
END DO

IF (can_model  ==  2) THEN
!     Radiative canopy without heat capacity
DO j=1,veg_pts
  l = veg_index(j)
  canhc(l) = 0.
  vfrac(l) = 1. - EXP(-kext(ft)*lai(l))
END DO

ELSE IF (can_model == 3 .OR. can_model == 4) THEN
!     Radiative canopy with heat capacity
  DO j=1,veg_pts
    l = veg_index(j)
    lai_bal(l) = ( a_ws(ft)*eta_sl(ft)*ht(l) /                    &
                   a_wl(ft) )**(1.0/(b_wl(ft)-1))
    leaf(l) = sigl(ft)*lai_bal(l)
    wood(l) = a_wl(ft)*(lai_bal(l)**b_wl(ft))
    canhc(l) = hleaf*leaf(l) + hwood*wood(l)
    vfrac(l) = 1. - EXP(-kext(ft)*lai(l))
  END DO

END IF

IF (lhook) CALL dr_hook('CANCAP',zhook_out,zhook_handle)
RETURN
END SUBROUTINE cancap
