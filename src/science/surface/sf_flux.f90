! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!     SUBROUTINE SF_FLUX -----------------------------------------------
! Description:
! Subroutines SF_FLUX to calculate explicit surface fluxes of
! heat and moisture
!-----------------------------------------------------------------------

SUBROUTINE sf_flux (                                              &
 points,tile_pts,fld_sea,pts_index,tile_index,                    &
 nsnow,n,canhc,dzsurf,hcons,qs1_elev,qstar,q_elev,radnet,resft,   &
 rhokh_1,smvcst,snowdepth,tile_frac,timestep,                     &
 t_elev,ts1_elev,tstar,vfrac,rhokh_can,                           &
 z0h,z0m_eff,z1_tq,lh0,emis_tile,emis_soil,                       &
 salinityfactor,anthrop_heat,scaling_urban,                       &
 fqw_1_gb,ftl_1_gb,                                               &
 alpha1,ashtf_prime,fqw_1,epot,ftl_1,dtstar                       &
 ,sea_point                                                       &
 )

USE atm_fields_bounds_mod
USE theta_field_sizes, ONLY : t_i_length

USE csigma
USE c_r_cp
USE c_lheat
USE c_g
USE c_0_dg_c
USE c_epslon, ONLY : c_virtual,epsil=>repsilon
USE snow_param, ONLY : snow_hcon
USE surf_param, ONLY : grcp,ls
USE nstypes, ONLY : urban_canyon, urban_roof
USE switches_urban, ONLY : l_moruses_storage
USE switches, ONLY : l_aggregate                                  &
                    ,l_epot_corr

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

INTEGER                                                           &
 points                                                           &
                      ! IN Total number of points.
,tile_pts                                                         &
                      ! IN Number of tile points.
,pts_index(points)                                                &
                      ! IN Index of points.
,tile_index(points)                                               &
                      ! IN Index of tile points.
,nsnow(points)                                                    &
                      ! IN Number of snow layers
,n                    ! IN Tile number.
                      ! For sea and sea-ice this = 0

REAL                                                              &
 fld_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)     &
!                          ! IN Fraction of land or sea
,canhc(points)                                                    &
!                          ! IN Areal heat capacity of canopy (J/K/m2).
,dzsurf(points)                                                   &
!                          ! IN Surface layer thickness (m).
,hcons(points)                                                    &
!                          ! IN Soil thermal conductivity (W/m/K).
,qs1_elev(points)                                                 &
!                          ! IN Sat. specific humidity qsat(T_ELEV,PSTAR
,qstar(points)                                                    &
!                          ! IN Surface qsat.
,q_elev(points)                                                   &
!                          ! IN Total water content of lowest
!                          !    atmospheric layer (kg per kg air).
,radnet(points)                                                   &
!                          ! IN Net surface radiation (W/m2) positive
!                          !    downwards
,resft(points)                                                    &
!                          ! IN Total resistance factor.
,rhokh_1(points)                                                  &
!                          ! IN Surface exchange coefficient.
,smvcst(points)                                                   &
!                          ! IN Volumetric saturation point
!                          !    - zero at land-ice points.
,snowdepth(points)                                                &
!                          ! IN Snow depth (on ground) (m)
,tile_frac(points)                                                &
!                          ! IN Tile fraction.
,timestep                                                         &
!                          ! IN Timestep (s).
,t_elev(points)                                                   &
!                          ! IN Liquid/frozen water temperature for
!                          !     lowest atmospheric layer (K).
,ts1_elev(points)                                                 &
!                          ! IN Temperature of surface layer (K).
,tstar(points)                                                    &
!                          ! IN Surface temperature (K).
,vfrac(points)                                                    &
!                          ! IN Fractional canopy coverage.
,rhokh_can(points)                                                &
!                          ! IN Exchange coefficient for canopy air
!                          !     to surface
,z0h(points)                                                      &
!                          ! IN Roughness length for heat and moisture
,z0m_eff(points)                                                  &
!                          ! IN Effective roughness length for momentum
,z1_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)       &
!                          ! IN Height of lowest atmospheric level (m).
,emis_tile(points)                                                &
!                          ! IN Emissivity for land tiles
,emis_soil(points)                                                &
!                          ! IN Emissivity of underlying soil
,lh0                                                              &
!                          ! IN Latent heat for snow free surface
!                          !    =LS for sea-ice, =LC otherwise
,salinityfactor                                                   &
!                          ! IN Factor allowing for the effect of the
!                          !    salinity of sea water on the
!                          !    evaporative flux.
,anthrop_heat(points)                                             &
!                          ! IN Anthropogenic contribution to surface
!                          !    heat flux (W/m2). Zero except for
!                          !    urban and L_ANTHROP_HEAT=.true.
!                          !    or for urban_canyon & urban_roof when
!                          !    l_urban2T is also .TRUE.
,scaling_urban(points)
!                          ! IN MORUSES: ground heat flux scaling;
!                          ! canyon tile only coupled to soil.
!                          ! This equals 1.0 except for urban tiles when
!                          ! MORUSES is used.


REAL                                                              &
 fqw_1_gb(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)    &
!                          ! INOUT GBM surface flux of
!                          !       QW (kg/m2/s).
,ftl_1_gb(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                           ! INOUT GBM surface flux of TL.

REAL                                                              &
 ashtf_prime(points)                                              &
!                          ! OUT Adjusted SEB coefficient
,alpha1(points)                                                   &
!                          ! OUT Gradient of saturated specific humidity
!                          !     with respect to temperature between the
!                          !     bottom model layer and the surface.
,fqw_1(points)                                                    &
!                          ! OUT Local surface flux of QW (kg/m2/s).
,epot(points)                                                     &
                     ! OUT
,ftl_1(points)                                                    &
!                          ! OUT Local surface flux of TL.
,dtstar(points)      ! OUT Change in TSTAR over timestep

REAL                                                              &
 sea_point           ! =1.0 IF SEA POINT, =0.0 OTHERWISE

! Workspace
REAL                                                              &
 ashtf(points)                                                    &
!                          ! Coefficient to calculate surface
!                          ! heat flux into soil (W/m2/K).
,dtstar_pot(points)                                               &
                           ! Change in TSTAR over timestep that is
!                          ! appropriate for the potential evaporation
,surf_ht_flux        ! Flux of heat from surface to sub-surface

! Scalars
INTEGER                                                           &
 i,j                                                              &
!                          ! Horizontal field index.
,k                                                                &
!                          ! Tile field index.
,l                   ! Points field index.

REAL                                                              &
 dq1                                                              &
!                          ! (qsat(TL_1,PSTAR)-QW_1) + g/cp*alpha1*Z1
,ds_ratio                                                         &
!                          ! 2 * snowdepth / depth of top soil layer.
,d_t                                                              &
!                          ! Temporary in calculation of alpha1.
,lh                  ! Latent heat (J/K/kg).

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('SF_FLUX',zhook_in,zhook_handle)

!-----------------------------------------------------------------------
!!  0 initialise
!-----------------------------------------------------------------------
ashtf( :) = 0.0
alpha1(:) = 0.0
ftl_1( :) = 0.0
epot(  :) = 0.0

!-----------------------------------------------------------------------
!!  1 Calculate gradient of saturated specific humidity for use in
!!    calculation of surface fluxes
!-----------------------------------------------------------------------
DO k=1,tile_pts
  l = tile_index(k)
  d_t = tstar(l) - t_elev(l)
  IF (d_t > 0.05 .OR. d_t < -0.05) THEN
    alpha1(l) = (qstar(l) - qs1_elev(l)) / d_t
  ELSE IF (t_elev(l) > tm) THEN
    alpha1(l) = epsil*lc*qs1_elev(l)*                             &
                (1. + c_virtual*qs1_elev(l)) /                    &
             ( r*t_elev(l)*t_elev(l))
  ELSE
    alpha1(l) = epsil*ls*qs1_elev(l)*                             &
                (1. + c_virtual*qs1_elev(l)) /                    &
             ( r*t_elev(l)*t_elev(l))
  END IF
END DO

! MORUSES: This only affects land tiles and not land-ice, sea or sea-ice tiles
! as the test on smvcst is used. This slightly revised version of the original
! coding.
! Note: This modification also only works when vf_nvg = 0. This is set in
! physiol for urban tiles.

DO k=1,tile_pts
  l = tile_index(k)
  IF ( .NOT. l_aggregate .AND. l_moruses_storage                  &
       .AND. smvcst(l) > EPSILON(smvcst)                          &
       .AND. n == urban_roof ) THEN
    hcons(l) = 0.0
  END IF

! Except when n == urban_canyon when MORUSES is used
! scaling_urban(l) = 1.0
  ashtf(l) = 2.0 * hcons(l) / dzsurf(l)
  ashtf(l) = ashtf(l) * scaling_urban(l)

  IF (snowdepth(l) > 0.0 .AND. smvcst(l) > EPSILON(smvcst(l))     &
     .AND. nsnow(l) == 0 ) THEN
    IF ( l_moruses_storage .AND. n == urban_roof ) THEN
! This required as HCONS(L) = 0 in this case. No clash with ice-tile here.
      ashtf(l) =  0.0
    ELSE
      ds_ratio = 2.0 * snowdepth(l) / dzsurf(l)
      IF (ds_ratio <= 1.0) THEN
        ashtf(l) =  ashtf(l) /                                    &
           (1. + ds_ratio*(hcons(l)/snow_hcon - 1.))
      ELSE
        ashtf(l) =  ashtf(l)*snow_hcon / hcons(l)
      END IF
    END IF
  END IF
END DO

DO k=1,tile_pts
  l = tile_index(k)
  j=(pts_index(l)-1)/t_i_length + 1
  i = pts_index(l) - (j-1)*t_i_length

  lh = lh0
  IF (snowdepth(l) > 0.) lh = ls

  ftl_1(l) = rhokh_1(l) * (tstar(l) - t_elev(l) -                 &
                  grcp*(z1_tq(i,j) + z0m_eff(l) - z0h(l)))
  epot(l) = rhokh_1(l) * (salinityfactor*qstar(l) - q_elev(l))
  fqw_1(l) = resft(l) * epot(l)

  surf_ht_flux = ((1.0 - vfrac(l)) * ashtf(l) +                   &
                    vfrac(l) * rhokh_can(l)) *                    &
                                   (tstar(l) - ts1_elev(l)) +     &
                 vfrac(l)*emis_soil(l)*emis_tile(l)*sbcon*        &
                  (tstar(l)**4.0 - ts1_elev(l)**4.0)


  ashtf_prime(l) = 4.0*(1. + emis_soil(l)*vfrac(l))*              &
                             emis_tile(l)*sbcon*tstar(l)**3.0 +   &
                             vfrac(l)*rhokh_can(l) +              &
                (1.0 - vfrac(l))*ashtf(l) + canhc(l)/timestep

  dtstar(l) = (radnet(l) + anthrop_heat(l) - cp*ftl_1(l) -        &
                         lh*fqw_1(l) - surf_ht_flux)  /           &
               ( rhokh_1(l)*(cp + lh*alpha1(l)*resft(l)) +        &
                    ashtf_prime(l) )

  ftl_1(l) = ftl_1(l) + rhokh_1(l) * dtstar(l) * (1.0 - sea_point)
  fqw_1(l) = fqw_1(l) + resft(l) * rhokh_1(l) * alpha1(l) *       &
                        dtstar(l) * (1.0 - sea_point)

  IF (l_epot_corr) THEN
    dtstar_pot(l) = (radnet(l) + anthrop_heat(l) - cp*ftl_1(l) -  &
                         lh*epot(l) - surf_ht_flux)  /            &
               ( rhokh_1(l)*(cp + lh*alpha1(l)) + ashtf_prime(l) )
    epot(l) = epot(l) + rhokh_1(l) * alpha1(l) * dtstar_pot(l)
  ELSE
    epot(l) = epot(l) + rhokh_1(l) * alpha1(l) * dtstar(l)
  END IF

  ftl_1_gb(i,j) = ftl_1_gb(i,j) +                                 &
                  fld_sea(i,j)*tile_frac(l)*ftl_1(l)
  fqw_1_gb(i,j) = fqw_1_gb(i,j) +                                 &
                  fld_sea(i,j)*tile_frac(l)*fqw_1(l)

END DO

IF (lhook) CALL dr_hook('SF_FLUX',zhook_out,zhook_handle)
RETURN
END SUBROUTINE sf_flux
