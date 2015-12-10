! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Description:
! Routine to calculate the aerodynamic resistance

!**********************************************************************
SUBROUTINE raero (land_pts,land_index,veg_pts,veg_index           &
,                 rib,wind,z0h,z0m,z1,ra)

USE atm_fields_bounds_mod
USE theta_field_sizes, ONLY : t_i_length

USE c_vkman, ONLY : karman=>vkman
USE surf_param, ONLY : ah,cz

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

INTEGER                                                           &
 land_pts                                                         &
                            ! IN Total number of land points.
,land_index(land_pts)                                             &
                            ! IN Index of land points on the
!                                 !    P-grid.
,veg_pts                                                          &
                            ! IN Number of vegetated points.
,veg_index(land_pts)                                              &
                            ! IN Index of vegetated points
!                                 !    on the land grid.
,i,j,k,l                    ! WORK Loop counters.

REAL                                                              &
 rib(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)         &
                            ! IN Bulk Richardson number.
,wind(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)        &
                            ! IN Windspeed (m/s).
,z0h(land_pts)                                                    &
                            ! IN Roughness length for heat (m).
,z0m(land_pts)                                                    &
                            ! IN Roughness length for momentum (m)
,z1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)          &
                            ! IN Reference level (m).
,ra(land_pts)                                                     &
                            ! OUT Aerodynamic resistance (s/m).
,bh                                                               &
                            ! WORK Stability coefficient.
,chn(land_pts)                                                    &
                            ! WORK Neutral drag coefficient.
,fh(land_pts)                                                     &
                            ! WORK Stability factor.
,zetah,zetam                ! WORK Tempories in calculation of
!                                 !      CHN.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------
! Calculate the neutral bulk tranfer coefficient.
!-----------------------------------------------------------------------
IF (lhook) CALL dr_hook('RAERO',zhook_in,zhook_handle)

DO k=1,veg_pts
  l = veg_index(k)
  j=(land_index(l)-1)/t_i_length + 1
  i = land_index(l) - (j-1)*t_i_length
  zetam = LOG((z1(i,j) + z0m(l)) / z0m(l))
  zetah = LOG((z1(i,j) + z0m(l)) / z0h(l))
  chn(l) = (karman * karman) / (zetah * zetam)
END DO

!-----------------------------------------------------------------------
! Calculate the stability factor.
!-----------------------------------------------------------------------
DO k=1,veg_pts
  l = veg_index(k)
  j=(land_index(l)-1)/t_i_length + 1
  i = land_index(l) - (j-1)*t_i_length
  bh = ah * chn(l) * cz * SQRT (z1(i,j) / z0h(l))
  IF (rib(i,j)  >=  0.0) THEN
    fh(l) = 1.0 / (1 + ah * rib(i,j))
  ELSE
    fh(l) = 1 - ah * rib(i,j) / (1 + bh * SQRT(-rib(i,j)))
  END IF
END DO

!-----------------------------------------------------------------------
! Calculate the aerodynamic resistance.
!-----------------------------------------------------------------------
DO k=1,veg_pts
  l = veg_index(k)
  j=(land_index(l)-1)/t_i_length + 1
  i = land_index(l) - (j-1)*t_i_length
  ra(l) = 1.0 / (fh(l) * chn(l) * wind(i,j))
END DO

IF (lhook) CALL dr_hook('RAERO',zhook_out,zhook_handle)
RETURN
END SUBROUTINE raero
