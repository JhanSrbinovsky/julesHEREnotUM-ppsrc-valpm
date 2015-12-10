! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  SUBROUTINE IM_SF_PT2 ---------------------------------------------

!  Purpose: Calculate implicit increments to surface variables
!           for the unconditionally stable and non-oscillatory
!           BL numerical solver.

!---------------------------------------------------------------------
!  Arguments :-
SUBROUTINE im_sf_pt2 (                                            &
 land_pts,land_index,ntiles,tile_index,tile_pts                   &
,flandg,tile_frac,snow_tile,nice_use,ice_fract,ice_fract_cat      &
,gamma_in,gamma1_in,gamma2_in,alpha1,alpha1_sice                  &
,ashtf_prime,ashtf_prime_tile                                     &
,resft,dtstar_tile,dtstar                                         &
,rhokm_u_1,rhokm_v_1,rhokh_1,rhokh1_sice                          &
,ctctq1,dqw1_1,dtl1_1,cq_cm_u_1                                   &
,cq_cm_v_1,du_1,dv_1,du_star1,dv_star1,flandg_u,flandg_v          &
,fqw_gb,ftl_gb                                                    &
,taux_1,taux_land,taux_land_star,taux_ssi,taux_ssi_star,tauy_1    &
,tauy_land,tauy_land_star,tauy_ssi,tauy_ssi_star                  &
,fqw_tile,epot_tile,ftl_tile,fqw_ice,ftl_ice,e_sea,h_sea          &
,l_correct,l_flux_bc                                              &
)

USE dynamics_grid_mod, ONLY: l_vatpoles

USE atm_fields_bounds_mod
USE theta_field_sizes, ONLY : t_i_length

USE c_r_cp
USE c_lheat
USE surf_param, ONLY : ls
USE switches  , ONLY : l_epot_corr

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

LOGICAL                                                           &
 l_correct                                                        &
                             ! Flag for BL solver
,l_flux_bc                   ! Logical for prescribed
                             ! surface fluxes (SCM)

INTEGER                                                           &
 land_pts                                                         &
                             ! IN Total number of land points.
,land_index(land_pts)                                             &
                             ! IN Index of land points.
,ntiles                                                           &
                             ! IN Number of land surface tiles.
,tile_index(land_pts,ntiles)                                      &
                             ! IN Index of tile points.
,tile_pts(ntiles)                                                 &
                             ! IN Number of tiles.
,nice_use                    ! IN Number of sea ice categories fully
                             !    used in surface exchange

REAL                                                              &
 flandg(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)      &
                             ! IN Land fraction
,tile_frac(land_pts,ntiles)                                       &
                             ! IN Tile fraction
,snow_tile(land_pts,ntiles)                                       &
                             ! IN Lying snow on land tiles (kg/m2)
,ice_fract(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)   &
                             ! IN Fraction of grid-box which is
!                                  !    sea-ice (decimal fraction).
,ice_fract_cat(tdims%i_start:tdims%i_end,                         &
               tdims%j_start:tdims%j_end,nice_use) &
                             ! IN Fraction of grid-box which is
!                            !    sea-ice (decimal fraction) per 
                             !    category if nice_use>1
,gamma_in                                                         &
                             ! IN Implicit weighting.
,gamma1_in(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)   &
,gamma2_in(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)   &
                             ! IN Implicit weights for uncond.
                             !    stable non-oscillatory BL solver
,alpha1(land_pts,ntiles)                                          &
                             ! IN Gradient of saturated specific
!                                  !    humidity with respect to
!                                  !    temperature between the bottom
!                                  !    model layer and the surface.
,alpha1_sice(tdims%i_start:tdims%i_end,                           &
             tdims%j_start:tdims%j_end,nice_use)                  &
                             ! IN ALPHA1 for sea-ice
,ashtf_prime(tdims%i_start:tdims%i_end,                           & 
             tdims%j_start:tdims%j_end,nice_use)                  &
                             ! IN Coefficient to calculate surface
!                                  !    heat flux into soil or sea-ice
!                                  !    (W/m2/K).

,ashtf_prime_tile(land_pts,ntiles)                                &
                             ! IN Coefficient to calculate heat
!                                  !    flux into land tiles (W/m2/K).
,dtstar_tile(land_pts,ntiles)                                     &
                             ! INOUT Change in TSTAR over timestep
!                                  !    for land tiles
,dtstar(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use)&
                             ! INOUT Change in TSTAR over timestep
!                                  !    for sea-ice
,resft(land_pts,ntiles)                                           &
                             ! IN Total resistance factor
,epot_tile(land_pts,ntiles)                                       &
                             ! IN surface tile potential
!                                  !    evaporation
,e_epot_tile(land_pts,ntiles)                                     &
                             ! Work ratio of explicit
!                                  !      EPOT/E
!                                  !    evaporation

,rhokm_u_1(udims%i_start:udims%i_end,udims%j_start:udims%j_end)   &
                             ! IN Level 1 exchange coefficient for
!                                  !    momentum
,rhokm_v_1(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end) &
                             ! IN Level 1 exchange coefficient for
!                                  !    momentum
,rhokh_1(land_pts,ntiles)                                         &
                             ! IN Surface exchange coeffs for FTL

,rhokh1_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end) &
                             ! IN Sea and sea-ice surface exchange
,ctctq1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)      &
,cq_cm_u_1(udims%i_start:udims%i_end,udims%j_start:udims%j_end)   &
                             ! IN Coefficient in U and V
!                                  !     tri-diagonal implicit matrix
,cq_cm_v_1(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end) &
                             ! IN Coefficient in U and V
!                                  !     tri-diagonal implicit matrix
,dqw1_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)      &
                             ! IN Incr obtained by semi-implicit inte
,dtl1_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)      &
                             ! IN Incr obtained by semi-implicit inte
,du_1(udims_s%i_start:udims_s%i_end,                              &
      udims_s%j_start:udims_s%j_end)                              &
                             ! IN Level 1 increment to u wind
!                                  !    field
,dv_1(vdims_s%i_start:vdims_s%i_end,                              &
      vdims_s%j_start:vdims_s%j_end)                              &
                             ! IN Level 1 increment to v wind
!                                  !    field
,du_star1(udims_s%i_start:udims_s%i_end,                          &
          udims_s%j_start:udims_s%j_end)                          &
,dv_star1(vdims_s%i_start:vdims_s%i_end,                          &
          vdims_s%j_start:vdims_s%j_end)                          &
,flandg_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end)    &
                             ! IN Land fraction on U grid.
,flandg_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)
                             ! IN Land fraction on V grid.


REAL                                                              &
 fqw_gb(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)      &
                             ! INOUT Grid-box value of QW flux at
!                                  !       Kg/sq m/s
,ftl_gb(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)      &
                             ! INOUT Grid-box value of TL flux at
!                                  !       i.e. H/Cp where H is sensible
!                                  !       in W per sq m).
,taux_1(udims%i_start:udims%i_end,udims%j_start:udims%j_end)      &
                             ! OUT   x-component of turbulent
!                                  !       stress at surface.
,taux_land(udims%i_start:udims%i_end,udims%j_start:udims%j_end)   &
                             ! INOUT x-component of turbulent
!                                  !       stress at land surface.
,taux_ssi(udims%i_start:udims%i_end,udims%j_start:udims%j_end)    &
                             ! INOUT x-component of turbulent
!                                  !       stress at sea surface.
,tauy_1(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)      &
                             ! OUT   y-component of turbulent
!                                  !       stress at surface.
,tauy_land(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)   &
                             ! INOUT y-component of turbulent
!                                  !       stress at land surface.
,tauy_ssi(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)    &
                             ! INOUT y-component of turbulent
!                                  !       stress at sea surface.
,taux_land_star(udims%i_start:udims%i_end,                        &
                udims%j_start:udims%j_end)                        &
,taux_ssi_star(udims%i_start:udims%i_end,                         &
               udims%j_start:udims%j_end)                         &
,tauy_land_star(vdims%i_start:vdims%i_end,                        &
                vdims%j_start:vdims%j_end)                        &
,tauy_ssi_star(vdims%i_start:vdims%i_end,                         &
               vdims%j_start:vdims%j_end)                         &
                             ! INOUT as above, temporarily needed
                             !       by predictor stage of the
                             !       new uncond stable BL solver
,fqw_tile(land_pts,ntiles)                                        &
                             ! INOUT Tile flux of QW. Kg/sq m/s
,ftl_tile(land_pts,ntiles)                                        &
                             ! INOUT Tile flux of TL
,e_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)       &
                             ! OUT Evaporation from sea times
!                                  !       leads fraction (kg/m2/s).
!                                  !       Zero over land.
,h_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)       &
                                   ! OUT Surface sensible heat flux ov
!                                  !       sea times leads fraction (W/m
!                                  !       Zero over land.
,fqw_ice(tdims%i_start:tdims%i_end,                               &
         tdims%j_start:tdims%j_end,nice_use)                      &
                             ! INOUT  surface flux of QW for
!                                  !  sea-ice fraction of gridsquare.
,ftl_ice(tdims%i_start:tdims%i_end,                               &
         tdims%j_start:tdims%j_end,nice_use)
                             ! INOUT surface flux of TL for
!                                  !  sea-ice fraction of gridsquare.

!  Local and other symbolic constants :-

! Workspace :-
REAL                                                              &
 rhokpm(land_pts,ntiles)                                          &
                             !  Surface exchange coeff for tiles
,rhokpm_sice(tdims%i_start:tdims%i_end,                           &
             tdims%j_start:tdims%j_end,nice_use)                  &
                             !  Sea-ice surface exchange coeff
,lat_ht                                                           &
          ! Latent heat of evaporation for snow-free land
!               ! or sublimation for snow-covered land and ice.
,apart(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,2)     &
                             ! Temporary array
,bpart(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,2)     &
                             ! Temporary array
,recip(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)       &
                             ! Temporary array
,ftl_land(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)    &
                             ! Temporary array
,fqw_land(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)    &
                             ! Temporary array
,ftl_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)    &
                             ! Temporary array (sum over all cats)
,fqw_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                             ! Temporary array (sum over all cats)

LOGICAL                                                           &
 epot_calc(land_pts,ntiles)  ! flag to connect first and second
                             ! epot calculations

!  Local scalars :-
INTEGER                                                           &
 i,j                                                              &
          ! Loop counter (horizontal field index).
,k                                                                &
          ! Loop counter (tile index).
,l                                                                &
          ! Loop counter (horizontal land index).
,n                                                                &
          ! Loop counter (tile counter).
,offset   ! offset

REAL                                                              &
 ftl_old                                                          &
          ! Used to hold current value of FTL_GB before updating
,gamma_1                                                          &
,gamma1                                                           &
,gamma2

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('IM_SF_PT2',zhook_in,zhook_handle)

IF (l_vatpoles) THEN
  offset=1
ELSE
  offset=0
END IF

!-------------------------------------------------------------------------
! initialise epot_calc
!-------------------------------------------------------------------------
epot_calc(:,:) = .FALSE.

!-------------------------------------------------------------------------

! Now compute surface stresses for the 2nd stage (predictor) of
! the new scheme (BL solver) using its discretization.

!-------------------------------------------------------------------------

IF ( .NOT. l_correct ) THEN
!-------------------------------------------------------------------------

! Compute surface stresses for the 1st stage (predictor) of
! the new scheme (BL solver) using its discretization.

!-------------------------------------------------------------------------
! Copy gamma into an array defined on u-points
! In principle, gamma should be interpolated, but sensitivity is expected
! to be small so the adjacent p-point is used to avoid message passing
! This makes v-at-poles/ENDGAME formulation consistent with the current
!  formulation
  DO j=udims%j_start,udims%j_end
    DO i=udims%i_start,udims%i_end
        gamma1 = gamma1_in(i+offset,j)
        gamma2 = gamma2_in(i+offset,j)
      IF ( flandg_u(i,j) > 0.0 ) THEN
        taux_land_star(i,j) = ( gamma2*taux_land(i,j) +           &
                 gamma1*rhokm_u_1(i,j)*du_1(i,j) ) /              &
                ( 1.0 + gamma1*rhokm_u_1(i,j)*cq_cm_u_1(i,j) )
      ELSE
        taux_land_star(i,j) = 0.0
      END IF
      IF ( flandg_u(i,j) < 1.0 ) THEN
        taux_ssi_star(i,j) = ( gamma2*taux_ssi(i,j) +             &
                 gamma1*rhokm_u_1(i,j)*du_1(i,j) ) /              &
                ( 1.0 + gamma1*rhokm_u_1(i,j)*cq_cm_u_1(i,j) )
      ELSE
        taux_ssi_star(i,j) = 0.0
      END IF
      taux_1(i,j) = flandg_u(i,j)*taux_land_star(i,j)             &
                  + ( 1.0-flandg_u(i,j))*taux_ssi_star(i,j)
    END DO
  END DO

! Copy gamma into an array defined on u-points
! In principle, gamma should be interpolated, but sensitivity is expected
! to be small so the adjacent p-point is used to avoid message passing
! This makes v-at-poles/ENDGAME formulation consistent with the current
! formulation
  IF (l_vatpoles) THEN
  DO j=vdims%j_start,vdims%j_end
    DO i=vdims%i_start,vdims%i_end
        IF ( j .EQ. vdims%j_end ) THEN
          gamma1 = gamma1_in(i,udims%j_end)
          gamma2 = gamma2_in(i,udims%j_end)
        ELSE
          gamma1 = gamma1_in(i,j+1)
          gamma2 = gamma2_in(i,j+1)
        END IF
      IF ( flandg_v(i,j) > 0.0 ) THEN
        tauy_land_star(i,j) = ( gamma2*tauy_land(i,j) +           &
                 gamma1*rhokm_v_1(i,j)*dv_1(i,j) ) /              &
                ( 1.0 + gamma1*rhokm_v_1(i,j)*cq_cm_v_1(i,j) )
      ELSE
        tauy_land_star(i,j) = 0.0
      END IF
      IF ( flandg_v(i,j) < 1.0 ) THEN
        tauy_ssi_star(i,j) = ( gamma2*tauy_ssi(i,j) +             &
                 gamma1*rhokm_v_1(i,j)*dv_1(i,j) ) /              &
                ( 1.0 + gamma1*rhokm_v_1(i,j)*cq_cm_v_1(i,j) )
      ELSE
        tauy_ssi_star(i,j) = 0.0
      END IF
      tauy_1(i,j) = flandg_v(i,j)*tauy_land_star(i,j)             &
                + ( 1.0-flandg_v(i,j))*tauy_ssi_star(i,j)
    END DO
  END DO
  ELSE
  DO j=vdims%j_start,vdims%j_end
    DO i=vdims%i_start,vdims%i_end
        gamma1 = gamma1_in(i,j)
        gamma2 = gamma2_in(i,j)
      IF ( flandg_v(i,j) > 0.0 ) THEN
        tauy_land_star(i,j) = ( gamma2*tauy_land(i,j) +           &
                 gamma1*rhokm_v_1(i,j)*dv_1(i,j) ) /              &
                ( 1.0 + gamma1*rhokm_v_1(i,j)*cq_cm_v_1(i,j) )
      ELSE
        tauy_land_star(i,j) = 0.0
      END IF
      IF ( flandg_v(i,j) < 1.0 ) THEN
        tauy_ssi_star(i,j) = ( gamma2*tauy_ssi(i,j) +             &
                 gamma1*rhokm_v_1(i,j)*dv_1(i,j) ) /              &
                ( 1.0 + gamma1*rhokm_v_1(i,j)*cq_cm_v_1(i,j) )
      ELSE
        tauy_ssi_star(i,j) = 0.0
      END IF
      tauy_1(i,j) = flandg_v(i,j)*tauy_land_star(i,j)             &
                + ( 1.0-flandg_v(i,j))*tauy_ssi_star(i,j)
    END DO
  END DO
  END IF

  gamma_1 = gamma_in

! time weights for specified scalar fluxes

  IF ( l_flux_bc ) gamma_1 = 0.0

!-------------------------------------------------------------------------

! Compute scalar surface fluxes using the standard scheme described
! in MOSES 2.2 technical documentation (hadley centre tecnical note 30).
! This needs to be done only in the first stage of the scheme (predictor).
! The same fluxes will be used for the 2nd stage (corrector).
! NOTE: scalar surface fluxes could be computed using the new scheme
!       but with a penalty in code complexity.

!------------------------------------------------------------------------

! Initialise APART and BPART to zero
  DO j=tdims%j_start,tdims%j_end
    DO i=tdims%i_start,tdims%i_end
      apart(i,j,1)=0.0
      apart(i,j,2)=0.0
      bpart(i,j,1)=0.0
      bpart(i,j,2)=0.0
      ftl_land(i,j)=0.0
      fqw_land(i,j)=0.0
      ftl_sice(i,j)=0.0
      fqw_sice(i,j)=0.0
    END DO
  END DO

! Land tiles
  DO n=1,ntiles
    DO k=1,tile_pts(n)
      l = tile_index(k,n)
      j=(land_index(l)-1)/t_i_length + 1
      i = land_index(l) - (j-1)*t_i_length
      lat_ht = lc
      IF (snow_tile(l,n) >  0.) lat_ht = ls

      rhokpm(l,n) = rhokh_1(l,n) / ( ashtf_prime_tile(l,n) +      &
               rhokh_1(l,n)*(lat_ht*alpha1(l,n)*resft(l,n) + cp) )

      apart(i,j,1)=apart(i,j,1) - tile_frac(l,n) *                &
               gamma_1 * rhokpm(l,n) *                            &
            ( lat_ht*resft(l,n)*rhokh_1(l,n)*alpha1(l,n) +        &
                         ashtf_prime_tile(l,n) )
      apart(i,j,2)=apart(i,j,2) + tile_frac(l,n) *                &
               gamma_1 * rhokpm(l,n) *                            &
               lat_ht*resft(l,n)*rhokh_1(l,n)
      bpart(i,j,1)=bpart(i,j,1) + tile_frac(l,n) *                &
               gamma_1 * resft(l,n)*rhokpm(l,n) *                 &
               cp*rhokh_1(l,n)*alpha1(l,n)
      bpart(i,j,2)=bpart(i,j,2) - tile_frac(l,n) *                &
               gamma_1 * resft(l,n)*rhokpm(l,n) *                 &
               ( cp*rhokh_1(l,n) + ashtf_prime_tile(l,n) )
    END DO
  END DO

! Sea points
  DO j=tdims%j_start,tdims%j_end
    DO i=tdims%i_start,tdims%i_end

      IF ( flandg(i,j) < 1.0 ) THEN
! Note, for icy points, this is the leads region in which case ice_fract=0 
        apart(i,j,1)= flandg(i,j)*apart(i,j,1)                    &
       - gamma_1 * ( 1.0 - flandg(i,j) ) * rhokh1_sice(i,j)       &
       * ( 1.0 - ice_fract(i,j) ) 
        apart(i,j,2)= flandg(i,j)*apart(i,j,2)
        bpart(i,j,1)= flandg(i,j)*bpart(i,j,1)
        bpart(i,j,2)= flandg(i,j)*bpart(i,j,2)                    &
       - gamma_1*( 1.0 - flandg(i,j) )*rhokh1_sice(i,j)           &
       * ( 1.0 - ice_fract(i,j) )    
 
      END IF
    END DO
  END DO

! Seaice points
  DO n=1,nice_use
    DO j=tdims%j_start,tdims%j_end
      DO i=tdims%i_start,tdims%i_end

        IF ( flandg(i,j) < 1.0 .AND. ice_fract_cat(i,j,n) >  0.0 ) THEN

          rhokpm_sice(i,j,n) = rhokh1_sice(i,j) / ( ashtf_prime(i,j,n) &
                 + rhokh1_sice(i,j)*(ls*alpha1_sice(i,j,n) + cp) )

          apart(i,j,1)=apart(i,j,1)                                    &
          - gamma_1 * (1.0-flandg(i,j)) * ice_fract_cat(i,j,n)         &
          * rhokpm_sice(i,j,n) *                                       &
          ( ls*rhokh1_sice(i,j)*alpha1_sice(i,j,n) + ashtf_prime(i,j,n) )

          apart(i,j,2)=apart(i,j,2)                                    &
          + gamma_1 * (1.0-flandg(i,j)) *ice_fract_cat(i,j,n)          &
          * rhokpm_sice(i,j,n) * ls*rhokh1_sice(i,j)

          bpart(i,j,1)=bpart(i,j,1)                                    &
          + gamma_1 * ice_fract_cat(i,j,n) * ( 1.0 - flandg(i,j) )     &
          * rhokpm_sice(i,j,n) *cp*rhokh1_sice(i,j)*alpha1_sice(i,j,n)

          bpart(i,j,2)=bpart(i,j,2)                                    &
         - gamma_1 * ice_fract_cat(i,j,n) * ( 1.0 - flandg(i,j) )      &
         * rhokpm_sice(i,j,n) *                                        &
         ( cp*rhokh1_sice(i,j) + ashtf_prime(i,j,n) ) 
 
       END IF
     END DO
   END DO
 END DO

! Land tiles
DO n=1,ntiles
  DO k=1,tile_pts(n)
    l = tile_index(k,n)
    j=(land_index(l)-1)/t_i_length + 1
    i = land_index(l) - (j-1)*t_i_length
    e_epot_tile(l,n)=1.0
    IF(      (epot_tile(l,n) >  0.0)                              &
        .AND.(fqw_tile(l,n)  > SQRT(TINY(0.0))) )THEN
      e_epot_tile(l,n)=epot_tile(l,n)/fqw_tile(l,n)
      epot_calc(l,n)=.TRUE.
    ENDIF
  END DO
END DO

! Calculate grid-box fluxes of heat and moisture
  DO j=tdims%j_start,tdims%j_end
    DO i=tdims%i_start,tdims%i_end

      recip(i,j)=( 1.0 + ctctq1(i,j)*apart(i,j,1) ) *             &
           ( 1.0 + ctctq1(i,j)*bpart(i,j,2) ) -                   &
           ctctq1(i,j)*apart(i,j,2)*ctctq1(i,j)*bpart(i,j,1)

      ftl_old=ftl_gb(i,j)

      ftl_gb(i,j) = ( ( 1.0 + ctctq1(i,j)*bpart(i,j,2) ) *        &
               ( ftl_old + apart(i,j,1)*dtl1_1(i,j) +             &
                 apart(i,j,2)*dqw1_1(i,j)) -                      &
                 ctctq1(i,j)*apart(i,j,2) * ( fqw_gb(i,j) +       &
                 bpart(i,j,1)*dtl1_1(i,j) +                       &
                 bpart(i,j,2)*dqw1_1(i,j)) ) / recip(i,j)

      fqw_gb(i,j) = ( ( 1.0 + ctctq1(i,j)*apart(i,j,1) ) *        &
                ( fqw_gb(i,j) + bpart(i,j,1)*dtl1_1(i,j) +        &
                  bpart(i,j,2)*dqw1_1(i,j)) -                     &
                  ctctq1(i,j)*bpart(i,j,1) * ( ftl_old +          &
                  apart(i,j,1)*dtl1_1(i,j) +                      &
                  apart(i,j,2)*dqw1_1(i,j)) ) / recip(i,j)

    END DO
  END DO

! Make implicit correction to tile fluxes

! Land tiles
  DO n=1,ntiles
    DO k=1,tile_pts(n)
      l = tile_index(k,n)
      j=(land_index(l)-1)/t_i_length + 1
      i = land_index(l) - (j-1)*t_i_length
      lat_ht = lc
      IF (snow_tile(l,n) >  0.) lat_ht = ls

      ftl_tile(l,n)=ftl_tile(l,n) -                               &
               gamma_1 * rhokpm(l,n) *                            &
            ( lat_ht*resft(l,n)*rhokh_1(l,n)*alpha1(l,n) +        &
                       ashtf_prime_tile(l,n) ) *                  &
         ( dtl1_1(i,j) - ctctq1(i,j)*ftl_gb(i,j) ) +              &
               gamma_1 * rhokpm(l,n) *                            &
               lat_ht*resft(l,n)*rhokh_1(l,n) *                   &
         ( dqw1_1(i,j) - ctctq1(i,j)*fqw_gb(i,j) )

      fqw_tile(l,n)=fqw_tile(l,n) +                               &
               gamma_1 * resft(l,n)*rhokpm(l,n) *                 &
               cp*rhokh_1(l,n)*alpha1(l,n) *                      &
         ( dtl1_1(i,j) - ctctq1(i,j)*ftl_gb(i,j) ) -              &
               gamma_1 * resft(l,n)*rhokpm(l,n) *                 &
               ( cp*rhokh_1(l,n) + ashtf_prime_tile(l,n) ) *      &
         ( dqw1_1(i,j) - ctctq1(i,j)*fqw_gb(i,j) )

      IF (l_epot_corr) THEN
        IF (epot_calc(l,n)) THEN
          epot_tile(l,n)=fqw_tile(l,n)*e_epot_tile(l,n)
        END IF
      ELSE
        epot_tile(l,n)=fqw_tile(l,n)*e_epot_tile(l,n)
      END IF

      fqw_land(i,j)=fqw_land(i,j)+fqw_tile(l,n)*tile_frac(l,n)
      ftl_land(i,j)=ftl_land(i,j)+ftl_tile(l,n)*tile_frac(l,n)

      dtstar_tile(l,n) = dtstar_tile(l,n) + gamma_1 *             &
             ( cp * rhokh_1(l,n) *                                &
               ( dtl1_1(i,j) - ctctq1(i,j) * ftl_gb(i,j) ) +      &
               lat_ht * resft(l,n) * rhokh_1(l,n) *               &
               ( dqw1_1(i,j) - ctctq1(i,j) * fqw_gb(i,j) ) ) /    &
             ( cp * rhokh_1(l,n) *                                &
               ( cp + lat_ht * resft(l,n) * alpha1(l,n) ) +       &
               ashtf_prime_tile(l,n) )

    END DO
  END DO

! Make implicit correction for sea ice category fluxes
  DO n=1,nice_use
    DO j=tdims%j_start,tdims%j_end
      DO i=tdims%i_start,tdims%i_end

       IF(flandg(i,j) <  1.0 .AND. ice_fract_cat(i,j,n) >  0.0) THEN

        ftl_ice(i,j,n)=ftl_ice(i,j,n) -                           &
             gamma_1 * rhokpm_sice(i,j,n) *                       &
            ( ls*rhokh1_sice(i,j)*alpha1_sice(i,j,n) +            &
                     ashtf_prime(i,j,n) ) *                       &
            ( dtl1_1(i,j) - ctctq1(i,j)*ftl_gb(i,j) ) +           &
                     gamma_1 * rhokpm_sice(i,j,n) *               &
                     ls*rhokh1_sice(i,j) *                        &
            ( dqw1_1(i,j) - ctctq1(i,j)*fqw_gb(i,j) )

        fqw_ice(i,j,n)=fqw_ice(i,j,n) +                           &
             gamma_1 * rhokpm_sice(i,j,n) *                       &
             cp*rhokh1_sice(i,j)*alpha1_sice(i,j,n) *             &
            ( dtl1_1(i,j) - ctctq1(i,j)*ftl_gb(i,j) ) -           &
                 gamma_1 *rhokpm_sice(i,j,n) *                    &
            ( cp*rhokh1_sice(i,j) + ashtf_prime(i,j,n) ) *        &
            ( dqw1_1(i,j) - ctctq1(i,j)*fqw_gb(i,j) )

        ! Calculate total sea ice fluxes
        fqw_sice(i,j)=fqw_sice(i,j) + fqw_ice(i,j,n)              &
                            * ice_fract_cat(i,j,n)/ice_fract(i,j)
        ftl_sice(i,j)=ftl_sice(i,j) + ftl_ice(i,j,n)              &
                            * ice_fract_cat(i,j,n)/ice_fract(i,j)

        dtstar(i,j,n) = dtstar(i,j,n) + gamma_1 *                 &
                ( cp * rhokh1_sice(i,j) *                         &
                ( dtl1_1(i,j) - ctctq1(i,j) * ftl_gb(i,j) ) +     &
                   ls * rhokh1_sice(i,j) *                        &
                ( dqw1_1(i,j) - ctctq1(i,j) * fqw_gb(i,j) ) ) /   &
                ( cp * rhokh1_sice(i,j) *                         &
                ( cp + ls * alpha1_sice(i,j,n) ) +                &
                   ashtf_prime(i,j,n) )

       ENDIF

      END DO
    END DO
  END DO

! Sea points  (= GBM - land fluxes - sea ice fluxes)
  DO j=tdims%j_start,tdims%j_end
    DO i=tdims%i_start,tdims%i_end

      h_sea(i,j) = 0.0
      e_sea(i,j) = 0.0

      IF(flandg(i,j) <  1.0 .and. (1.0-ice_fract(i,j)) /= 0.0 ) THEN
          h_sea(i,j)=( (ftl_gb(i,j)                           &
             - ftl_land(i,j)*flandg(i,j))/(1.-flandg(i,j))            &
             - ftl_sice(i,j) * ice_fract(i,j) ) /           &
             (1.0 - ice_fract(i,j))

          e_sea(i,j)=( (fqw_gb(i,j)                             &
             - fqw_land(i,j)*flandg(i,j))/(1.-flandg(i,j))             &
             - fqw_sice(i,j) * ice_fract(i,j) ) /          &
             (1.0 - ice_fract(i,j))

      END IF
    END DO
  END DO

ELSE ! IF L_CORRECT==TRUE THEN:

! Copy gamma into an array defined on u-points
! In principle, gamma should be interpolated, but sensitivity is expected
! to be small so the adjacent p-point is used to avoid message passing
! This makes v-at-poles/ENDGAME formulation consistent with the current
! formulation
  DO j=udims%j_start,udims%j_end
    DO i=udims%i_start,udims%i_end
        gamma1 = gamma1_in(i+offset,j)
        gamma2 = gamma2_in(i+offset,j)
      IF ( flandg_u(i,j) >  0.0 ) THEN
        taux_land(i,j) = ( gamma2*(taux_land(i,j)+                &
                               rhokm_u_1(i,j)*du_star1(i,j)) +    &
                 gamma1*rhokm_u_1(i,j)*du_1(i,j) ) /              &
                ( 1.0 + gamma1*rhokm_u_1(i,j)*cq_cm_u_1(i,j) )
      ELSE
        taux_land(i,j) = 0.0
      END IF
      IF ( flandg_u(i,j) <  1.0 ) THEN
        taux_ssi(i,j) = ( gamma2*(taux_ssi(i,j)+                  &
                              rhokm_u_1(i,j)*du_star1(i,j)) +     &
                 gamma1*rhokm_u_1(i,j)*du_1(i,j) ) /              &
                ( 1.0 + gamma1*rhokm_u_1(i,j)*cq_cm_u_1(i,j) )
      ELSE
        taux_ssi(i,j) = 0.0
      END IF
      taux_1(i,j) = flandg_u(i,j)*taux_land(i,j)                  &
                  + ( 1.0-flandg_u(i,j))*taux_ssi(i,j)
    END DO
  END DO

! Copy gamma into an array defined on u-points
! In principle, gamma should be interpolated, but sensitivity is expected
! to be small so the adjacent p-point is used to avoid message passing
! This makes v-at-poles/ENDGAME formulation consistent with the current
! formulation
  IF (l_vatpoles) THEN
  DO j=vdims%j_start,vdims%j_end
    DO i=vdims%i_start,vdims%i_end
        IF ( j .EQ. vdims%j_end ) THEN
          gamma1 = gamma1_in(i,udims%j_end)
          gamma2 = gamma2_in(i,udims%j_end)
        ELSE
          gamma1 = gamma1_in(i,j+1)
          gamma2 = gamma2_in(i,j+1)
        END IF
      IF ( flandg_v(i,j) >  0.0 ) THEN
        tauy_land(i,j) = ( gamma2*(tauy_land(i,j)+                &
                               rhokm_v_1(i,j)*dv_star1(i,j)) +    &
                 gamma1*rhokm_v_1(i,j)*dv_1(i,j) ) /              &
                ( 1.0 + gamma1*rhokm_v_1(i,j)*cq_cm_v_1(i,j) )
      ELSE
        tauy_land(i,j) = 0.0
      END IF
      IF ( flandg_v(i,j) <  1.0 ) THEN
        tauy_ssi(i,j) = ( gamma2*(tauy_ssi(i,j)+                  &
                               rhokm_v_1(i,j)*dv_star1(i,j)) +    &
                 gamma1*rhokm_v_1(i,j)*dv_1(i,j) ) /              &
                ( 1.0 + gamma1*rhokm_v_1(i,j)*cq_cm_v_1(i,j) )
      ELSE
        tauy_ssi(i,j) = 0.0
      END IF
      tauy_1(i,j) = flandg_v(i,j)*tauy_land(i,j)                  &
                  + ( 1.0-flandg_v(i,j))*tauy_ssi(i,j)
    END DO
  END DO
  ELSE
  DO j=vdims%j_start,vdims%j_end
    DO i=vdims%i_start,vdims%i_end
        gamma1 = gamma1_in(i,j)
        gamma2 = gamma2_in(i,j)
      IF ( flandg_v(i,j) >  0.0 ) THEN
        tauy_land(i,j) = ( gamma2*(tauy_land(i,j)+                &
                               rhokm_v_1(i,j)*dv_star1(i,j)) +    &
                 gamma1*rhokm_v_1(i,j)*dv_1(i,j) ) /              &
                ( 1.0 + gamma1*rhokm_v_1(i,j)*cq_cm_v_1(i,j) )
      ELSE
        tauy_land(i,j) = 0.0
      END IF
      IF ( flandg_v(i,j) <  1.0 ) THEN
        tauy_ssi(i,j) = ( gamma2*(tauy_ssi(i,j)+                  &
                               rhokm_v_1(i,j)*dv_star1(i,j)) +    &
                 gamma1*rhokm_v_1(i,j)*dv_1(i,j) ) /              &
                ( 1.0 + gamma1*rhokm_v_1(i,j)*cq_cm_v_1(i,j) )
      ELSE
        tauy_ssi(i,j) = 0.0
      END IF
      tauy_1(i,j) = flandg_v(i,j)*tauy_land(i,j)                  &
                  + ( 1.0-flandg_v(i,j))*tauy_ssi(i,j)
    END DO
  END DO
  END IF

END IF ! L_CORRECT

IF (lhook) CALL dr_hook('IM_SF_PT2',zhook_out,zhook_handle)
RETURN
END SUBROUTINE im_sf_pt2
