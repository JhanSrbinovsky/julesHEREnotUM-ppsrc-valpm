
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine FLAKE_INIT ----------------------------------
!
! Description: Dummy routine for UM without FLake.
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.
!
!   Code Owner: See Unified Model Code Owner's HTML page
!   This file belongs in section: Land

SUBROUTINE flake_interface ( land_points, tile_pts, tile_index                 &
                           , u_s_lake, surf_ht_flux_lk, sw_down                &
                           , lake_depth, lake_fetch, coriolis_param, timestep  &
                           , lake_albedo                                       &
                           , lake_t_snow, lake_t_ice, lake_t_mean, lake_t_mxl  &
                           , lake_shape_factor                                 &
                           , lake_h_snow, lake_h_ice, lake_h_mxl               &
                           , lake_t_sfc, ts1_lake, g_dt                        &
                           , trap_frozen, trap_unfrozen )

USE ereport_mod, ONLY : ereport
IMPLICIT NONE
INTEGER, INTENT(IN) ::   land_points     &
                       , tile_pts        &
                       , trap_frozen     &
                       , trap_unfrozen

INTEGER, INTENT(IN), DIMENSION (land_points) :: tile_index

REAL, INTENT(IN) ::  timestep

REAL, INTENT(IN), DIMENSION (land_points) ::  &
                          surf_ht_flux_lk     &
                         , u_s_lake           &
                         , sw_down            &
                         , lake_depth         &
                         , lake_fetch         &
                         , coriolis_param     &
                         , lake_albedo        &
                         , lake_t_snow        &
                         , lake_t_ice         &
                         , lake_t_mean        &
                         , lake_t_mxl         &
                         , lake_shape_factor  &
                         , lake_h_snow        &
                         , lake_h_ice         &
                         , lake_h_mxl         &
                         , lake_t_sfc         &
                         , ts1_lake           &
                         , g_dt

! Error reporting

  CHARACTER*(*), PARAMETER :: RoutineName = 'flake_interface_0a'
  CHARACTER(LEN=256)       :: message
  INTEGER                  :: errorstatus

!---------------------------------------------------
! This routine is a dummy and should not be called.
! If it has been called, it indicates that
! L_FLAKE_MODEL is TRUE without the UM_FLAKE if-def
! being set.
!---------------------------------------------------
  errorstatus = 42
  message = "The FLake switch should only be TRUE with FLake-compatible code."
  CALL ereport(RoutineName, errorstatus, message)


RETURN

END SUBROUTINE flake_interface
