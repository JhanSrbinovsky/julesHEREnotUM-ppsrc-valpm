
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! SUBROUTINE pft_sparm
!
! Purpose:
! Routine to calculate the land surface parameters of a given PFT from
! its areal fraction and structural properties.

!**********************************************************************
SUBROUTINE pft_sparm  (land_pts,n,tile_index,tile_pts             &
,                      ht,lai,satcon,catch_t,infil_t,z0_t)


USE pftparm

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

INTEGER                                                           &
 land_pts                                                         &
                            ! IN Number of land points.
,n                                                                &
                            ! IN Plant functional type.
,tile_pts                                                         &
                            ! IN Number of land points which
!                                 !    include the surface type.
,tile_index(land_pts)                                             &
                            ! IN Indices of land points which
!                                 !    include the surface type.
,j,l                        ! WORK Loop counters.

REAL                                                              &
 ht(land_pts)                                                     &
                            ! IN Vegetation height (m).
,lai(land_pts)                                                    &
                            ! IN Leaf area index.
,satcon(land_pts)                                                 &
                            ! IN Saturated hydraulic conductivity
!                                 !    (kg/m2/s).
,catch_t(land_pts)                                                &
                            ! OUT Canopy capacity (kg/m2).
,infil_t(land_pts)                                                &
                            ! OUT Maximum surface infiltration
!                                 !     rate (kg/m2/s).
,z0_t(land_pts)             ! OUT Roughness length (m).

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('PFT_SPARM',zhook_in,zhook_handle)

DO j=1,tile_pts
  l = tile_index(j)
  z0_t(l) = dz0v_dh(n) * ht(l)
  catch_t(l) = catch0(n) + dcatch_dlai(n) * lai(l)
  infil_t(l) = infil_f(n) * satcon(l)
END DO


IF (lhook) CALL dr_hook('PFT_SPARM',zhook_out,zhook_handle)
RETURN
END SUBROUTINE pft_sparm
