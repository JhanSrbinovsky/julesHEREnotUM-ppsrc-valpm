! *****************************COPYRIGHT*******************************

! (c) [University of Edinburgh] [2009]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the JULES collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC237]

! *****************************COPYRIGHT*******************************
!  SUBROUTINE SNOW-----------------------------------------------------

! Description:
!     Calling routine for snow module

! Subroutine Interface:
SUBROUTINE snow ( land_pts,timestep,stf_hf_snow_melt,ntiles,      &
                  tile_pts,tile_index,catch_snow,con_snow,        &
                  tile_frac,ls_snow,ei_tile,hcaps,hcons,melt_tile,&
                  smcl1,sthf1,surf_htf_tile,tsoil1,tstar_tile,    &
                  v_sat1,rgrain,rgrainl,rho_snow_grnd,sice,       &
                  sliq,snow_grnd,snow_tile,snowdepth,             &
                  tsnow,nsnow,ds,hf_snow_melt,lying_snow,         &
                  rho_snow,snomlt_sub_htf,snow_melt,              &
                  surf_ht_flux_ld )

USE c_lheat, ONLY :                                               &
!  imported scalar parameters
 lf                 !  latent heat of fusion of water at 0degc

USE ancil_info, ONLY :                                            &
!  imported arrays with intent(in)
 nsmax              !  Maximum possible number of snow layers

USE snow_param, ONLY :                                            &
 cansnowtile        !  switch for canopy snow model

USE switches, ONLY :                                              &
!  imported scalars with intent(in)
 l_snow_albedo        ! Switch for prognostic snow albedo


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
,tstar_tile(land_pts,ntiles)                                      &
                       ! Tile surface temperature (K)
,v_sat1(land_pts)      ! Surface soil layer volumetric
!                            ! moisture concentration at saturation

! Array arguments with intent(inout)
REAL, INTENT(INOUT) ::                                            &
 melt_tile(land_pts,ntiles)                                       &
                       ! Surface or canopy snowmelt rate (kg/m2/s)
                       ! On output, this is the total melt rate
                       ! for the tile (i.e. sum of  melt on canopy
                       ! and ground).
,tsoil1(land_pts)                                                 &
                       ! Soil surface layer temperature (K)
,rgrain(land_pts,ntiles)                                          &
                       ! Snow surface grain size (microns)
,rgrainl(land_pts,ntiles,nsmax)                                   &
                       ! Snow layer grain size (microns)
,rho_snow_grnd(land_pts,ntiles)                                   &
                       ! Snowpack bulk density (kg/m3)
,sice(land_pts,ntiles,nsmax)                                      &
                       ! Ice content of snow layers (kg/m2)
,sliq(land_pts,ntiles,nsmax)                                      &
                       ! Liquid content of snow layers (kg/m2)
,snow_grnd(land_pts,ntiles)                                       &
                       ! Snow beneath canopy (kg/m2)
,snow_tile(land_pts,ntiles)                                       &
                       ! Snow mass (Kg/m2)
,snowdepth(land_pts,ntiles)                                       &
                       ! Snow depth (m)
,tsnow(land_pts,ntiles,nsmax)
                       ! Snow layer temperatures (K)

! Array arguments with intent(out)
INTEGER, INTENT(OUT) ::                                           &
 nsnow(land_pts,ntiles)   ! Number of snow layers

REAL, INTENT(OUT) ::                                              &
 ds(land_pts,ntiles,nsmax)                                        &
                          ! Snow layer thicknesses (m)
,hf_snow_melt(land_pts)                                           &
                          ! Gridbox snowmelt heat flux (W/m2)
,lying_snow(land_pts)                                             &
                          ! Gridbox snow mass (kg m-2)
,rho_snow(land_pts,ntiles,nsmax)                                  &
                          ! Snow layer densities (kg/m3)
,snomlt_sub_htf(land_pts)                                         &
                          ! Sub-canopy snowmelt heat flux (W/m2)
,snow_melt(land_pts)                                              &
                          ! Gridbox snowmelt (Kg/m2/s)
,surf_ht_flux_ld(land_pts)   ! Surface heat flux on land (W/m2).

! Local scalars
INTEGER ::                                                        &
 i                                                                &
                      ! land point index and loop counter
,j                                                                &
                      ! tile pts loop counter
,k                                                                &
                      ! tile number
,n                    ! tile loop counter

! Local arrays
REAL ::                                                           &
 csnow(land_pts,nsmax)                                            &
                      ! Areal heat capacity of layers (J/K/m2)
,ksnow(land_pts,nsmax)                                            &
                      ! Thermal conductivity of layers (W/m/K)
,rho0(land_pts)                                                   &
                      ! Density of fresh snow (kg/m3)
!                         ! Where NSNOW=0, rho0 is the density
!                         ! of the snowpack.
,snowfall(land_pts)                                               &
                      ! Snowfall reaching the ground in timestep
                      ! (kg/m2) - includes any canopy unloading
,snowmass(land_pts)                                               &
                      ! Snow mass on the ground (Kg/m2)
,rgrain0(land_pts)                                                &
                      ! Fresh snow grain size (microns)
,sice0(land_pts)                                                  &
                      ! Ice content of fresh snow (kg/m2)
!                           ! Where NSNOW=0, SICE0 is the mass of
!                           ! the snowpack.
,snow_can(land_pts,ntiles)                                        &
                      ! Canopy snow load (Kg/m2)
,snow_soil_htf(land_pts,ntiles)                                   &
                      ! Heat flux into the uppermost
!                           ! subsurface layer (W/m2)
!                           ! i.e. snow to ground, or into
!                           ! snow/soil composite layer
,tsnow0(land_pts)     ! Temperature of fresh snow (K)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------
! External routines called:
!-----------------------------------------------------------------------
EXTERNAL canopysnow,layersnow,relayersnow,snowtherm,snowpack,     &
         snowgrain,compactsnow

!-----------------------------------------------------------------------
! Initialise gridbox variables.
!-----------------------------------------------------------------------
IF (lhook) CALL dr_hook('SNOW',zhook_in,zhook_handle)

lying_snow(:) = 0.0
snomlt_sub_htf(:) = 0.0
snow_melt(:) = 0.0

DO n=1,ntiles
!-----------------------------------------------------------------------
! Set snow mass variables
!-----------------------------------------------------------------------

  IF ( cansnowtile(n) ) THEN
!-----------------------------------------------------------------------
!         With the canopy snow model, snow_tile is held on the canopy,
!         while the mass of the snowpack (on the ground) is snow_grnd.
!-----------------------------------------------------------------------
    DO k=1,tile_pts(n)
      i = tile_index(k,n)
      snow_can(i,n) = snow_tile(i,n)
      snowmass(i) = snow_grnd(i,n)

!-----------------------------------------------------------------------
! Subtract sublimation and melt of canopy snow.
!-----------------------------------------------------------------------
      snow_can(i,n) = snow_can(i,n) -                             &
                     ( ei_tile(i,n) + melt_tile(i,n) ) * timestep
    END DO

  ELSE

!-----------------------------------------------------------------------
! Without the snow canopy model, all the snow is in the snowpack
!-----------------------------------------------------------------------
    DO k=1,tile_pts(n)
      i = tile_index(k,n)
      snowmass(i) = snow_tile(i,n)
    END DO

  END IF

!-----------------------------------------------------------------------
! Canopy interception, throughfall and unloading of snow
!-----------------------------------------------------------------------
! DEPENDS ON: canopysnow
  CALL canopysnow ( land_pts,tile_pts(n),timestep,cansnowtile(n), &
                    tile_index(:,n),catch_snow(:,n),con_snow,     &
                    ls_snow,melt_tile(:,n),snow_can(:,n),         &
                    snowfall )

!-----------------------------------------------------------------------
! Divide snow pack into layers
!-----------------------------------------------------------------------
! DEPENDS ON: layersnow
  CALL layersnow ( land_pts,tile_pts(n),tile_index(:,n),          &
                   snowdepth(:,n),nsnow(:,n),ds(:,n,:) )

!-----------------------------------------------------------------------
! Thermal properties of snow layers
!-----------------------------------------------------------------------
  IF ( nsmax > 0 )                                                &
! DEPENDS ON: snowtherm
    CALL snowtherm ( land_pts,tile_pts(n),nsnow(:,n),             &
                     tile_index(:,n),ds(:,n,:),sice(:,n,:),       &
                     sliq(:,n,:),csnow,ksnow )

!-----------------------------------------------------------------------
! Snow thermodynamics and hydrology
!-----------------------------------------------------------------------
! DEPENDS ON: snowpack
  CALL snowpack ( land_pts,tile_pts(n),timestep,cansnowtile(n),   &
                  nsnow(:,n),tile_index(:,n),csnow,ei_tile(:,n),  &
                  hcaps,hcons,ksnow,rho_snow_grnd(:,n),smcl1,     &
                  snowfall,sthf1,surf_htf_tile(:,n),              &
                  tile_frac(:,n),v_sat1,ds(:,n,:),melt_tile(:,n), &
                  sice(:,n,:),sliq(:,n,:),snomlt_sub_htf,         &
                  snowdepth(:,n),snowmass,tsnow(:,n,:),           &
                  tsoil1,snow_soil_htf(:,n),rho_snow(:,n,:),      &
                  rho0,sice0,tsnow0 )

!-----------------------------------------------------------------------
! Growth of snow grains
!-----------------------------------------------------------------------
  IF ( l_snow_albedo )                                            &
! DEPENDS ON: snowgrain
    CALL snowgrain ( land_pts,tile_pts(n),timestep,nsnow(:,n),    &
                     tile_index(:,n),sice(:,n,:),snowfall,        &
                     snowmass,tsnow(:,n,:),tstar_tile(:,n),       &
                     rgrain(:,n),rgrainl(:,n,:),rgrain0 )

  IF ( nsmax > 0 ) THEN
!-----------------------------------------------------------------------
! Mechanical compaction of snow
!-----------------------------------------------------------------------
! DEPENDS ON: compactsnow
    CALL compactsnow ( land_pts,tile_pts(n),timestep,nsnow(:,n),  &
                       tile_index(:,n),sice(:,n,:),sliq(:,n,:),   &
                       tsnow(:,n,:),rho_snow(:,n,:),ds(:,n,:) )

!-----------------------------------------------------------------------
! Redivide snowpack after changes in depth, conserving mass and energy
!-----------------------------------------------------------------------
! DEPENDS ON: relayersnow
    CALL relayersnow ( land_pts,tile_pts(n),tile_index(:,n),      &
                       rgrain0,rho0,sice0,snowfall,               &
                       snowmass,tsnow0,nsnow(:,n),ds(:,n,:),      &
                       rgrain(:,n),rgrainl(:,n,:),sice(:,n,:),    &
                       rho_snow_grnd(:,n),sliq(:,n,:),            &
                       tsnow(:,n,:),rho_snow(:,n,:),              &
                       snowdepth(:,n) )
  END IF  !  NSMAX>0

!-----------------------------------------------------------------------
! Copy into final snow mass variables
!-----------------------------------------------------------------------
  IF ( cansnowtile(n) ) THEN
    DO k=1,tile_pts(n)
      i = tile_index(k,n)
      snow_grnd(i,n) = snowmass(i)
      snow_tile(i,n) = snow_can(i,n)
    END DO
  ELSE
    DO k=1,tile_pts(n)
      i = tile_index(k,n)
      snow_tile(i,n) = snowmass(i)
    END DO
  END IF

!-----------------------------------------------------------------------
! Increment gridbox lying snow and snow melt.
!-----------------------------------------------------------------------
  DO k=1,tile_pts(n)
    i = tile_index(k,n)
    lying_snow(i) = lying_snow(i) +                               &
                       tile_frac(i,n) * snow_tile(i,n)
!     Add snow beneath canopy.
    IF ( cansnowtile(n) )                                         &
      lying_snow(i) = lying_snow(i) +                             &
                         tile_frac(i,n) * snow_grnd(i,n)

!     Snow melt.
    snow_melt(i) = snow_melt(i) + tile_frac(i,n) * melt_tile(i,n)
  END DO

END DO  !  tiles

!-----------------------------------------------------------------------
! Calculate the total snowmelt heat flux.
!-----------------------------------------------------------------------
IF ( stf_hf_snow_melt ) THEN
  DO i=1,land_pts
    hf_snow_melt(i) = lf * snow_melt(i)
  END DO
END IF

!-----------------------------------------------------------------------
! Calculate gridbox surface heat flux over land.
!-----------------------------------------------------------------------
surf_ht_flux_ld(:) = 0.0
DO n=1,ntiles
  DO j=1,tile_pts(n)
    i = tile_index(j,n)
    surf_ht_flux_ld(i) = surf_ht_flux_ld(i) +                     &
                          tile_frac(i,n) * snow_soil_htf(i,n)
  END DO
END DO
IF (lhook) CALL dr_hook('SNOW',zhook_out,zhook_handle)
RETURN

END SUBROUTINE snow
