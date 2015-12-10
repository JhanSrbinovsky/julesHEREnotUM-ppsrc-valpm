
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


SUBROUTINE flake_init(  land_field,ntiles,sm_levels,land_index  &
                       ,frac_land                               &
                       ,tstar_tile,snow_tile,deep_soil_temp )


  USE nstypes, ONLY : ntype

IMPLICIT NONE

! Subroutine arguments

!   Scalar arguments with intent(in):
  INTEGER, INTENT(IN) ::                                &
  land_field,                                           &
  ntiles,                                               &
  sm_levels

  INTEGER, INTENT(IN) ::                                &
  land_index(land_field)

  REAL, INTENT(IN) ::                                   &
  frac_land(land_field,ntype)

  REAL, INTENT(IN) ::                                   &
  tstar_tile(land_field,ntiles)                         &
 ,snow_tile( land_field,ntiles)                         &
 ,deep_soil_temp(land_field,sm_levels)

! Error reporting

  CHARACTER*(*), PARAMETER :: RoutineName = 'flake_init_0a'
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
! DEPENDS ON: EREPORT
  CALL ereport(RoutineName, errorstatus, message)


RETURN

END SUBROUTINE flake_init
