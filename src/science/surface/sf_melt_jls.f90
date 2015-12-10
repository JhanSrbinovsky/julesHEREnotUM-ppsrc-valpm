! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! SUBROUTINE SF_MELT----------------------------------------------------
!
! Purpose : Calculates surface melting (snow and sea-ice) and increments
!           surface fluxes to satisfy energy balance.
!           Sub-surface snowmelt is calculated and snowdepth incremented
!           by melt and sublimation in P251.
!-----------------------------------------------------------------------
SUBROUTINE sf_melt (                                              &
 points,pts_index                                                 &
,tile_index,tile_pts,fld_sea                                      &
,alpha1,ashtf_prime,dtrdz_1                                       &
,resft,rhokh_1,tile_frac,timestep,GAMMA                           &
,ei_tile,fqw_1,ftl_1,fqw_tile,ftl_tile                            &
,tstar_tile,snow_tile,snowdepth                                   &
,melt_tile                                                        &
 )

USE atm_fields_bounds_mod
USE theta_field_sizes, ONLY : t_i_length

USE c_r_cp
USE c_lheat
USE c_0_dg_c

USE switches, ONLY :                                              &
 frac_snow_subl_melt                                              &
,l_snowdep_surf

USE rad_param, ONLY :                                             &
 maskd

USE snow_param, ONLY :                                            &
 rho_snow_const

USE water_constants_mod, ONLY:                                    &
 rho_water

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

INTEGER                                                           &
 points                                                           &
                      ! IN Total number of points.
,pts_index(points)                                                &
                      ! IN Index of points.
,tile_index(points)                                               &
!                           ! IN Index of tile points.
,tile_pts             ! IN Number of tile points.

 REAL                                                             &
 fld_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)     &
!                           ! IN Fraction of land or sea.
,alpha1(points)                                                   &
!                           ! IN Gradients of saturated specific
!                           !    humidity with respect to temp.
!                           !    between the bottom model layer
!                           !    and surface.
,ashtf_prime(points)                                              &
!                           ! IN Adjusted SEB coefficient
,dtrdz_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)     &
!                           ! IN -g.dt/dp for surface layer
,resft(points)                                                    &
                       !IN Resistance factor.
,rhokh_1(points)                                                  &
!                           ! IN Surface exchange coefficient.
,tile_frac(points)                                                &
!                           ! IN Tile fractions.
,timestep                                                         &
                      ! IN Timestep (sec).
,GAMMA                ! IN implicit weight in level 1

REAL                                                              &
 ei_tile(points)                                                  &
!                           ! INOUT Sublimation for tile (kg/m2/s)
,fqw_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)       &
!                           ! INOUT GBM surface moisture flux (kg/m2/s).
,ftl_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)       &
!                           ! INOUT GBM surface sens. heat flux (W/m2).
,fqw_tile(points)                                                 &
!                           ! INOUT FQW for tile.
,ftl_tile(points)                                                 &
!                           ! INOUT FTL for tile.
,tstar_tile(points)                                               &
!                           ! INOUT Tile surface temperatures (K).
,snow_tile(points)                                                &
!                           ! INOUT Lying snow on tile (kg/m2).
,snowdepth(points)
!                           ! INOUT Depth of snow on tile (m).

REAL                                                              &
 melt_tile(points)
!                           ! OUT Surface snowmelt on tiles (kg/m2/s).

REAL                                                              &
 dfqw                                                             &
                      ! Moisture flux increment.
,dftl                                                             &
                      ! Sensible heat flux increment.
,dtstar                                                           &
                      ! Surface temperature increment.
,lcmelt                                                           &
                      ! Temporary in melt calculations.
,lsmelt                                                           &
                      ! Temporary in melt calculations.
,rhokh1_prime                                                     &
                      ! Modified forward time-weighted
!                           ! transfer coefficient.
,snow_max                                                         &
                      ! Snow available for melting.
,snow_density                                                     &
                      ! Density of snow on input
,tstarmax             ! Maximum gridbox mean surface temperature
!                           ! at sea points with ice.
INTEGER                                                           &
 i,j                                                              &
                      ! Loop counter - full horizontal field.
,k                                                                &
                      ! Loop counter - tile field
,l                                                                &
                      ! Loop counter - land field.
,n                    ! Loop counter - tile index.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('SF_MELT',zhook_in,zhook_handle)

melt_tile(:) = 0.

!-----------------------------------------------------------------------
!  Melt snow on tile if TSTAR_TILE is greater than TM.
!-----------------------------------------------------------------------
DO k=1,tile_pts
  l = tile_index(k)
  j=(pts_index(l)-1)/t_i_length + 1
  i = pts_index(l) - (j-1)*t_i_length
!
  IF (snowdepth(l) > SQRT(TINY(1.0))) THEN
    snow_density = snow_tile(l) / snowdepth(l)
  ELSE
    snow_density = rho_snow_const
  ENDIF
!
! The l_snowdep_surf if-test below is temporary to preserve
! bit comparison. It should be removed when this switch becomes
! default .TRUE.
!
  IF (l_snowdep_surf) THEN
    snow_density = MAX(rho_snow_const,snow_density)
    snow_density = MIN(rho_water     ,snow_density)
  ELSE
    snow_density = MAX(   1.0,snow_density)
    snow_density = MIN(1000.0,snow_density)
  ENDIF
!
  snow_max = MAX( 0.0, snow_tile(l) - ei_tile(l)*timestep )
  IF ( snow_max.gt.0.0 .AND. tstar_tile(l).gt.tm ) THEN
    rhokh1_prime = 1. / ( 1. / rhokh_1(l)                         &
                       + GAMMA*dtrdz_1(i,j) )
    lcmelt = (cp + lc*alpha1(l)*resft(l))*rhokh1_prime            &
             + ashtf_prime(l)
    lsmelt = lcmelt + lf*alpha1(l)*rhokh1_prime
    IF (frac_snow_subl_melt == 1) THEN
      dtstar = - MIN( (tstar_tile(l) - tm) *                      &
               (1.0 - EXP(-maskd*snow_max/snow_density)),         &
               lf*snow_max / (lcmelt*timestep) )
    ELSE
      dtstar = - MIN( tstar_tile(l) - tm ,                        &
                      lf*snow_max / (lcmelt*timestep) )
    END IF
    tstar_tile(l) = tstar_tile(l) + dtstar
    melt_tile(l) = - lsmelt*dtstar / lf
    dftl = cp*rhokh1_prime*dtstar
    dfqw = alpha1(l)*resft(l)*rhokh1_prime*dtstar
    ftl_tile(l) = ftl_tile(l) + dftl
    fqw_tile(l) = fqw_tile(l) + dfqw
    ei_tile(l) = ei_tile(l) + dfqw
!-----------------------------------------------------------------------
!  Update gridbox-mean quantities
!-----------------------------------------------------------------------
    dftl = tile_frac(l)*dftl
    dfqw = tile_frac(l)*dfqw
    ftl_1(i,j) = ftl_1(i,j) + fld_sea(i,j)*dftl
    fqw_1(i,j) = fqw_1(i,j) + fld_sea(i,j)*dfqw
  END IF
END DO

IF (lhook) CALL dr_hook('SF_MELT',zhook_out,zhook_handle)
RETURN
END SUBROUTINE sf_melt
