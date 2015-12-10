! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  SUBROUTINE SF_EVAP------------------------------------------------

!  Purpose: Calculate surface evaporation and sublimation amounts
!           (without applying them to the surface stores).

!  Suitable for single column usage.

!  Documentation: UMDP 24

!--------------------------------------------------------------------
SUBROUTINE sf_evap (                                              &
 land_pts,ntiles,                                                 &
 land_index,tile_index,tile_pts,nshyd,fland,                      &
 ashtf_prime_tile,canopy,dtrdz_1,flake,fraca,snow_tile,resfs,     &
 resft,rhokh_1,tile_frac,smc,wt_ext_tile,timestep,GAMMA,          &
 fqw_1,fqw_tile,ftl_1,ftl_tile,tstar_tile,                        &
 ecan,ecan_tile,elake_tile,esoil,esoil_tile,ei_tile,ext           &
 )

USE atm_fields_bounds_mod
USE theta_field_sizes, ONLY : t_i_length

USE c_r_cp
USE c_lheat
USE c_0_dg_c

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

INTEGER                                                           &
 land_pts                                                         &
                       ! IN Number of land points to be processed.
,ntiles                                                           &
                       ! IN Number of tiles per land point.
,land_index(land_pts)                                             &
                       ! IN Index of land points.
,tile_index(land_pts,ntiles)                                      &
!                            ! IN Index of tile points.
,tile_pts(ntiles)                                                 &
                       ! IN Number of tile points.
,nshyd                 ! IN Number of soil moisture levels.

REAL                                                              &
 fland(land_pts)                                                  &
                       ! IN Fraction of gridbox which is land.
,ashtf_prime_tile(land_pts,ntiles)                                &
!                            ! IN Adjusted SEB coefficient
,canopy(land_pts,ntiles)                                          &
!                            ! IN Surface/canopy water on land
!                            !    tiles (kg/m2).
,dtrdz_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)     &
!                            ! IN -g.dt/dp for surface layer
,flake(land_pts,ntiles)                                           &
                       ! IN Lake fraction.
,fraca(land_pts,ntiles)                                           &
                       ! IN Fraction of surface moisture flux
!                            !    with only aerodynamic resistance
!                            !    for land tiles.
,snow_tile(land_pts,ntiles)                                       &
!                            ! IN Lying snow amount on tiles (kg/m2).
,resfs(land_pts,ntiles)                                           &
                       ! IN Combined soil, stomatal and aerodynam.
!                            !    resistance factor for fraction 1-FRACA
!                            !    of land tiles.
,resft(land_pts,ntiles)                                           &
                       ! IN Total resistance factor
!                            !    FRACA+(1-FRACA)*RESFS.
,rhokh_1(land_pts,ntiles)                                         &
!                            ! IN Surface exchange coefficients.
,tile_frac(land_pts,ntiles)                                       &
!                            ! IN Tile fractions.
,smc(land_pts)                                                    &
                       ! IN Available soil moisture (kg/m2).
,wt_ext_tile(land_pts,nshyd,ntiles)                               &
!                            ! IN Fraction of transpiration
!                            !    extracted from each soil layer
!                            !    by each tile.
,timestep                                                         &
                       ! IN Timestep in seconds.
,GAMMA                 ! IN implicit weight in level 1

REAL                                                              &
 fqw_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)       &
                       ! INOUT Surface moisture flux (kg/m2/s).
,fqw_tile(land_pts,ntiles)                                        &
!                            ! INOUT Local FQW_1 for tiles.
,ftl_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)       &
                       ! INOUT Surface sensible heat flux (W/m2).
,ftl_tile(land_pts,ntiles)                                        &
!                            ! INOUT Local FTL_1 for tiles.
,tstar_tile(land_pts,ntiles)
!                            ! INOUT Tile surface temperatures (K).

REAL                                                              &
 ecan(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)        &
                       ! OUT Gridbox mean evaporation from canopy/
!                            !     surface store (kg per sq m per s).
!                            !     Zero over sea and sea-ice.
,ecan_tile(land_pts,ntiles)                                       &
!                            ! OUT ECAN for land tiles.
,elake_tile(land_pts,ntiles)                                      &
!                            ! OUT Lake evaporation.
,esoil(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)       &
                       ! OUT Gridbox mean evapotranspiration from
!                            !     soil moisture (kg per sq m per s).
!                            !     Zero over sea and sea-ice.
,esoil_tile(land_pts,ntiles)                                      &
!                            ! OUT ESOIL for land tiles.
,ei_tile(land_pts,ntiles)                                         &
!                            ! OUT Sublimation from snow or land-ice
!                            !     (kg per sq m per s).
,ext(land_pts,nshyd)   ! OUT Extraction of water from each
!                            !     soil layer (kg/m2/s).

REAL                                                              &
 dfqw(land_pts)                                                   &
                       ! Increment in GBM moisture flux.
,dftl(land_pts)                                                   &
                       ! Increment in GBM sensible heat flux.
,e_tile_old(land_pts,ntiles)                                      &
!                            ! Surface moisture flux before adjustment.
,le_tile_old(land_pts,ntiles)
!                            ! Surf latent heat flux before adjustment.

REAL                                                              &
 diff_lat_htf                                                     &
                       ! Increment in local latent heat flux.
,diff_sens_htf                                                    &
                       ! Increment in local sensible heat flux.
,dtstar                                                           &
                       ! Increment in local surface temperature.
,edt                                                              &
                       ! Moisture flux x timestep
,rhokh1_prime          ! Modified forward time-weighted
                       ! transfer coefficient.

INTEGER                                                           &
 i,j                                                              &
             ! Loop counter (horizontal field index).
,k                                                                &
             ! Loop counter (land, snow or land-ice field index).
,m                                                                &
             ! Loop counter (soil level index).
,l                                                                &
             ! Loop counter (land point field index).
,n           ! Loop counter (tile index).

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('SF_EVAP',zhook_in,zhook_handle)

DO n=1,ntiles
  DO k=1,tile_pts(n)
    l = tile_index(k,n)
    e_tile_old(l,n) = fqw_tile(l,n)
    IF (snow_tile(l,n)  >   0.) THEN
      le_tile_old(l,n) = (lc + lf)*fqw_tile(l,n)
    ELSE
      le_tile_old(l,n) = lc*fqw_tile(l,n)
    END IF
  END DO
END DO

DO n=1,ntiles
  DO l=1,land_pts
    ecan_tile(l,n) = 0.
    esoil_tile(l,n) = 0.
    elake_tile(l,n) = 0.
    ei_tile(l,n) = 0.
  END DO
END DO

!-----------------------------------------------------------------------
! Sublimation from snow-covered land tiles
!-----------------------------------------------------------------------
DO n=1,ntiles
  DO k=1,tile_pts(n)
    l = tile_index(k,n)
    IF (snow_tile(l,n)  >   0.) THEN
      ei_tile(l,n) =  fqw_tile(l,n)
      edt = ei_tile(l,n)*timestep
      IF ( edt  >   snow_tile(l,n) )                              &
        ei_tile(l,n) = snow_tile(l,n) / timestep
        fqw_tile(l,n) = fqw_tile(l,n) -  ei_tile(l,n)
    END IF
  END DO
END DO

!-----------------------------------------------------------------------
! Surface evaporation from and condensation onto snow-free land
!-----------------------------------------------------------------------
DO j=tdims%j_start,tdims%j_end
  DO i=tdims%i_start,tdims%i_end
    ecan(i,j) = 0.
    esoil(i,j) = 0.
  END DO
END DO

DO n=1,ntiles
  DO k=1,tile_pts(n)
    l = tile_index(k,n)
    j=(land_index(l)-1)/t_i_length + 1
    i = land_index(l) - (j-1)*t_i_length
    IF ( fqw_tile(l,n)  >   0.0 ) THEN
      ecan_tile(l,n) = (1. - flake(l,n)) *                        &
                       fraca(l,n) * fqw_tile(l,n) / resft(l,n)
      esoil_tile(l,n) = (1. - flake(l,n)) *                       &
                        (1. - fraca(l,n))*resfs(l,n)*fqw_tile(l,n)&
                                                      / resft(l,n)
      elake_tile(l,n) = flake(l,n)*fqw_tile(l,n) / resft(l,n)
      edt = ecan_tile(l,n)*timestep
      IF ( edt  >   canopy(l,n) ) THEN
        esoil_tile(l,n) =  (1. - flake(l,n)) *                    &
                           (1. - fraca(l,n)*canopy(l,n)/edt) *    &
                               resfs(l,n)*fqw_tile(l,n)/resft(l,n)
        ecan_tile(l,n) = canopy(l,n) / timestep
      END IF
    ELSE IF (snow_tile(l,n) <= 0.) THEN
      IF (tstar_tile(l,n) >= tm) THEN
        ecan_tile(l,n) = (1. - flake(l,n))*fqw_tile(l,n)
        elake_tile(l,n) = flake(l,n)*fqw_tile(l,n)
      ELSE
        ei_tile(l,n) =  fqw_tile(l,n)
      END IF
    END IF
    ecan(i,j) = ecan(i,j) + tile_frac(l,n)*ecan_tile(l,n)
    esoil(i,j) = esoil(i,j) + tile_frac(l,n)*esoil_tile(l,n)
  END DO
END DO

!-----------------------------------------------------------------------
! Soil evapotranspiration
!-----------------------------------------------------------------------
DO l=1,land_pts
  j=(land_index(l)-1)/t_i_length + 1
  i = land_index(l) - (j-1)*t_i_length
  edt = esoil(i,j)*timestep
  IF ( edt  >   smc(l) ) THEN
    DO n=1,ntiles
      esoil_tile(l,n) = smc(l)*esoil_tile(l,n) / edt
    END DO
    esoil(i,j) = smc(l) / timestep
  END IF
END DO

DO m=1,nshyd
  DO l=1,land_pts
    ext(l,m) = 0.
  END DO
END DO

DO m=1,nshyd
  DO n=1,ntiles
    DO k=1,tile_pts(n)
      l = tile_index(k,n)
      ext(l,m) = ext(l,m) + tile_frac(l,n)*wt_ext_tile(l,m,n)     &
                                          *esoil_tile(l,n)
    END DO
  END DO
END DO

!-----------------------------------------------------------------------
! Calculate increments to surface heat fluxes, moisture fluxes and
! temperatures
!-----------------------------------------------------------------------
DO l=1,land_pts
  dftl(l) = 0.
  dfqw(l) = 0.
END DO

DO n=1,ntiles
  DO k=1,tile_pts(n)
    l = tile_index(k,n)
    j=(land_index(l)-1)/t_i_length + 1
    i = land_index(l) - (j-1)*t_i_length
    rhokh1_prime = 1. / ( 1. / rhokh_1(l,n)                       &
                       + GAMMA*dtrdz_1(i,j) )
    diff_lat_htf = (lc + lf)*ei_tile(l,n) + lc*ecan_tile(l,n)     &
                    + lc*esoil_tile(l,n) + lc*elake_tile(l,n)     &
                    - le_tile_old(l,n)
    dtstar = - diff_lat_htf /                                     &
                ( cp*rhokh1_prime + ashtf_prime_tile(l,n) )
    diff_sens_htf = cp * rhokh1_prime * dtstar
    ftl_tile(l,n) = ftl_tile(l,n) + diff_sens_htf
    tstar_tile(l,n) = tstar_tile(l,n) + dtstar
    dftl(l) = dftl(l) + tile_frac(l,n)*diff_sens_htf
    dfqw(l) = dfqw(l) + tile_frac(l,n)*( ecan_tile(l,n) +         &
                  esoil_tile(l,n) + ei_tile(l,n) + elake_tile(l,n)&
                  - e_tile_old(l,n) )
  END DO
END DO

!-----------------------------------------------------------------------
! Update level 1 temperature and humidity and GBM heat and moisture
! fluxes due to limited moisture availability
!-----------------------------------------------------------------------
DO l=1,land_pts
  j=(land_index(l)-1)/t_i_length + 1
  i = land_index(l) - (j-1)*t_i_length
  ftl_1(i,j) = ftl_1(i,j) + fland(l)*dftl(l)
  fqw_1(i,j) = fqw_1(i,j) + fland(l)*dfqw(l)
END DO

IF (lhook) CALL dr_hook('SF_EVAP',zhook_out,zhook_handle)
RETURN
END SUBROUTINE sf_evap
