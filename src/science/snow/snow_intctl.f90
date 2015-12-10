
! *****************************COPYRIGHT*******************************

! (c) [University of Edinburgh] [2009]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the JULES collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC237]

! *****************************COPYRIGHT*******************************
!  SUBROUTINE SNOW_INTCTL-----------------------------------------------

! Description:
!     Calling routine for snow module

! Subroutine Interface:
SUBROUTINE snow_intctl ( land_pts,timestep,stf_hf_snow_melt,      &
                         ntiles,tile_pts,tile_index,catch_snow,   &
                         con_snow,tile_frac,ls_snow,ei_tile,hcaps,&
                         hcons,melt_tile,smcl1,sthf1,             &
                         surf_htf_tile,tsoil1,tstar_tile,v_sat1,  &
                         rgrain,snow_grnd,snow_melt,snow_tile,    &
                         hf_snow_melt,lying_snow,                 &
                         snomlt_sub_htf,surf_ht_flux_ld )

USE prognostics, ONLY: rgrainl,rho_snow_grnd,sice,sliq,          &
                       snowdepth,tsnow,nsnow,rho_snow


USE snow_param, ONLY : ds

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Scalar arguments with intent(in)
INTEGER, INTENT(IN) ::                                            &
 land_pts              ! Total number of land points

REAL, INTENT(IN) ::                                               &
 timestep              ! Timestep length (s)

LOGICAL, INTENT(IN) ::                                            &
 stf_hf_snow_melt      ! Stash flag for snowmelt heat flux

! Array arguments with intent(in)
INTEGER, INTENT(IN) ::                                            &
 ntiles                                                           &
                       !  Number of land tiles
,tile_pts(ntiles)                                                 &
                       ! Number of tile points
,tile_index(land_pts,ntiles)
                       ! Index of tile points

REAL, INTENT(IN) ::                                               &
 catch_snow(land_pts,ntiles)                                      &
                       ! canopy snow capacity (kg/m2)
,con_snow(land_pts)                                               &
                       ! Convective snowfall rate (kg/m2/s)
,tile_frac(land_pts,ntiles)                                       &
                       ! Tile fractions
,ls_snow(land_pts)                                                &
                       ! Large-scale snowfall rate (kg/m2/s)
,ei_tile(land_pts,ntiles)                                         &
                       ! Sublimation of snow (kg/m2/s)
,hcaps(land_pts)                                                  &
                       ! Soil heat capacity of top layer(J/K/m3).
,hcons(land_pts)                                                  &
                       ! Thermal conductivity of top soil layer,
!                            ! including water and ice (W/m/K)
,smcl1(land_pts)                                                  &
                       ! Moisture content of surface soil
!                            ! layer (kg/m2)
,sthf1(land_pts)                                                  &
                       ! Frozen soil moisture content of
!                            ! surface layer as a fraction of saturation.
,surf_htf_tile(land_pts,ntiles)                                   &
                       ! Surface heat flux (W/m2)
,tsoil1(land_pts)                                                 &
                       ! Soil surface layer temperature (K)
,tstar_tile(land_pts,ntiles)                                      &
                       ! Tile surface temperature (K)
,v_sat1(land_pts)      ! Surface soil layer volumetric
!                            ! moisture concentration at saturation

! Array arguments with intent(inout)
REAL, INTENT(INOUT) ::                                            &
 melt_tile(land_pts,ntiles)                                       &
                       ! Snowmelt rate (kg/m2/s)
,rgrain(land_pts,ntiles)                                          &
                       ! Snow surface grain size (microns)
,snow_grnd(land_pts,ntiles)                                       &
                       ! Snow mass on ground (kg/m2)
,snow_tile(land_pts,ntiles)
                       ! Lying snow on tiles (kg/m2)

! Array arguments with intent(out)
REAL, INTENT(OUT) ::                                              &
 hf_snow_melt(land_pts)                                           &
                          ! Gridbox snowmelt heat flux (W/m2)
,lying_snow(land_pts)                                             &
                          ! Gridbox snow mass (kg m-2)
,snomlt_sub_htf(land_pts)                                         &
                          ! Sub-canopy snowmelt heat flux (W/m2)
,snow_melt(land_pts)                                              &
                          ! Gridbox snowmelt (kg/m2/s)
,surf_ht_flux_ld(land_pts)   ! Surface heat flux on land (W/m2).

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------
! External routines called:
!-----------------------------------------------------------------------
EXTERNAL snow

IF (lhook) CALL dr_hook('SNOW_INTCTL',zhook_in,zhook_handle)

!! DEPENDS ON: snow
!CALL snow ( land_pts,timestep,stf_hf_snow_melt,ntiles,tile_pts,   &
!            tile_index,catch_snow,con_snow,tile_frac,ls_snow,     &
!            ei_tile,hcaps,hcons,melt_tile,smcl1,sthf1,            &
!            surf_htf_tile,tsoil1,tstar_tile,v_sat1,               &
!            rgrain,rgrainl,rho_snow_grnd,sice,sliq,               &
!            snow_grnd,snow_tile,snowdepth,                        &
!            tsnow,nsnow,ds,hf_snow_melt,lying_snow,rho_snow,      &
!            snomlt_sub_htf,snow_melt,                             &
!            surf_ht_flux_ld )

IF (lhook) CALL dr_hook('SNOW_INTCTL',zhook_out,zhook_handle)
RETURN
END SUBROUTINE snow_intctl
