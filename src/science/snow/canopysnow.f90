! *****************************COPYRIGHT*******************************

! (c) [University of Edinburgh] [2009]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the JULES collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC237]

! *****************************COPYRIGHT*******************************
!  SUBROUTINE CANOPYSNOW-----------------------------------------------

! Description:
!     Partition snowfall into canopy interception, throughfall and unloading.

! Subroutine Interface:
SUBROUTINE canopysnow ( land_pts,tile_pts,timestep,cansnowtile    &
,                       tile_index,catch_snow,con_snow,ls_snow    &
,                       melt_tile,snow_can,snowfall )

USE snow_param, ONLY :                                            &
!  imported scalars with intent(in)
 snowinterceptfact                                                &
                    !   constant in relationship between mass of
!                             intercepted snow and snowfall rate
,snowunloadfact      !   constant in relationship between canopy snow
!                             unloading and canopy snow melt rate

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Scalar arguments with intent(in)
INTEGER,INTENT(IN) ::                                             &
 land_pts                                                         &
                           !  Number of land points
,tile_pts                  !  Number of tile points

REAL, INTENT(IN) ::                                               &
 timestep                  !  Timestep (s)

LOGICAL, INTENT(IN) ::                                            &
 cansnowtile               !  Switch for canopy snow model

! Array arguments with intent(in)
INTEGER, INTENT(IN) ::                                            &
 tile_index(land_pts)      !  Index of tile points

REAL, INTENT(IN) ::                                               &
 catch_snow(land_pts)                                             &
                           !  Canopy snow capacity (kg/m2)
,con_snow(land_pts)                                               &
                           !  Convective snowfall rate (kg/m2/s)
,ls_snow(land_pts)                                                &
                           !  Large-scale snowfall rate (kg/m2/s)
,melt_tile(land_pts)       !  Canopy snow melt rate (kg/m2/s)

! Array arguments with intent(inout)
REAL, INTENT(INOUT) ::                                            &
 snow_can(land_pts)      !  Canopy snow load, if the
!                              !  canopy snow model is selected (kg/m2)

! Array arguments with intent(out)
REAL, INTENT(OUT) ::                                              &
 snowfall(land_pts)      !  Snow reaching the ground in timestep
                         !  (kg/m2)

! Local scalars
INTEGER ::                                                        &
 i                                                                &
                         !  Land point index
,k                       !  Tile point index

REAL ::                                                           &
 intercept                                                        &
                       !  Snow intercepted by canopy in timestep
                       !  (kg/m2)
,unload                !  Canopy snow unloaded in timestep (kg/m2)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!-------------------------------------------------------------------------------

IF (lhook) CALL dr_hook('CANOPYSNOW',zhook_in,zhook_handle)

IF ( cansnowtile ) THEN

  DO k=1,tile_pts
    i = tile_index(k)
    snowfall(i) = ( ls_snow(i) + con_snow(i) ) * timestep
    intercept = snowinterceptfact * (catch_snow(i) - snow_can(i)) &
                  * (1.0 - EXP(-snowfall(i)/catch_snow(i)))
    unload = snowunloadfact * melt_tile(i) * timestep
!-----------------------------------------------------------------------
!         At this point, the value of unload can be larger than the
!         amount of snow on the canopy (which has already had melt
!         and sublimation removed), so we need to limit to the amount
!         of snow available. However, snow_can can be <0 (with small
!         absolute value) because of "issues" in the surface flux code,
!         so we also restrict unload to be >=0.
!-----------------------------------------------------------------------
    unload = MAX( MIN( unload, snow_can(i) ), 0.0 )
    snow_can(i) = snow_can(i) + intercept - unload
    snowfall(i) = snowfall(i) - intercept + unload
  END DO

ELSE

  DO k=1,tile_pts
    i = tile_index(k)
    snowfall(i) = ( ls_snow(i) + con_snow(i) ) * timestep
  END DO

END IF

IF (lhook) CALL dr_hook('CANOPYSNOW',zhook_out,zhook_handle)
RETURN

END SUBROUTINE canopysnow
