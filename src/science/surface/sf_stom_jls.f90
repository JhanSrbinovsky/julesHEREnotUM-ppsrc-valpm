! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!**********************************************************************
! Routine to calculate the bulk stomatal resistance and the canopy
! CO2 fluxes

!**********************************************************************
SUBROUTINE sf_stom  (land_pts,land_index                          &
,                    veg_pts,veg_index                            &
,                    ft,co2,co2_3d,co2_dim_len                    &
,                    co2_dim_row,l_co2_interactive                &
,                    fsmc,ht,ipar,lai,pstar                       &
,                    q1,ra,tstar,o3                               &
,                    can_rad_mod,ilayers,faparv                   &
,                    gpp,npp,resp_p,resp_w,gc                     &
,                    fapar_sun,fapar_shd,fsun                     &
,                    flux_o3,fo3                                  &
,                    isoprene,terpene,methanol,acetone)

USE atm_fields_bounds_mod
USE theta_field_sizes, ONLY : t_i_length

USE switches, ONLY : l_bvoc_emis

USE pftparm, ONLY :                                               &
!      imported arrays with intent(in)
    a_wl,a_ws,b_wl,eta_sl,kpar,nl0,nr_nl,ns_nl,omega,r_grow,      &
    sigl,glmin

USE ccarbon, ONLY :                                               &
!      imported scalar parameters
   epco2,epo2
   
USE c_rmol, ONLY : rmol

USE surf_param, ONLY :                                            &
!      imported scalar parameters
   iter,kn,o2

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE


INTEGER                                                           &
 land_pts                                                         &
                            ! IN Number of land points to be
!                                 !    processed.
,land_index(land_pts)                                             &
                            ! IN Index of land points on the
!                                 !    P-grid.
,veg_pts                                                          &
                            ! IN Number of vegetated points.
,veg_index(land_pts)                                              &
                            ! IN Index of vegetated points
!                                 !    on the land grid.
,co2_dim_len                                                      &
                            ! IN Length of a CO2 field row.
,co2_dim_row                ! IN Number of CO2 field rows.

INTEGER                                                           &
 ft                         ! IN Plant functional type.

LOGICAL l_co2_interactive   ! switch for 3D CO2 field

INTEGER                                                           &
  can_rad_mod                                                     &
!                                  !Switch for canopy radiation model
 ,ilayers
!                                  !No of layers in canopy radiation model

REAL                                                              &
 co2                                                              &
                            ! IN Atmospheric CO2 concentration
,co2_3d(co2_dim_len,co2_dim_row)                                  &
!                                 ! IN 3D atmos CO2 concentration
!                                 !    (kg CO2/kg air).
,fsmc(land_pts)                                                   &
                            ! IN Soil water factor.
,ht(land_pts)                                                     &
                            ! IN Canopy height (m).
,ipar(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)        &
                            ! IN Incident PAR (W/m2).
,lai(land_pts)                                                    &
                            ! IN Leaf area index.
,pstar(land_pts)                                                  &
                            ! IN Surface pressure (Pa).
,faparv(land_pts,ilayers)                                         &
                            ! IN Profile of absorbed PAR.
,fapar_shd(land_pts,ilayers)                                      &
                            ! IN Profile of absorbed DIFF_PAR.
,fapar_sun(land_pts,ilayers)                                      &
                            ! IN Profile of absorbed DIR_PAR.
,fsun(land_pts,ilayers)                                           &
                            ! IN fraction of sunlit leaves
,q1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)          &
                            ! IN Specific humidity at level 1
,ra(land_pts)                                                     &
                            ! IN Aerodynamic resistance (s/m).
,tstar(land_pts)                                                  &
                            ! IN Surface temperature (K).
,o3(land_pts)                                                     &
                            ! IN Surface ozone concentration (ppb).
,gpp(land_pts)                                                    &
                            ! OUT Gross Primary Productivity
!                                 !     (kg C/m2/s).
,npp(land_pts)                                                    &
                            ! OUT Net Primary Productivity
!                                 !     (kg C/m2/s).
,resp_p(land_pts)                                                 &
                            ! OUT Plant respiration rate
!                                 !     (kg C/m2/sec).
,resp_w(land_pts)                                                 &
                            ! OUT Wood respiration rate
!                                 !     (kg C/m2/sec).
,flux_o3(land_pts)                                                &
                            ! OUT Flux of O3 to stomata (nmol O3/m2/s).
,fo3(land_pts)                                                    &
                            ! OUT Ozone exposure factor.
,gc(land_pts)
                            ! INOUT Canopy resistance to H2O
!                                 !       (m/s).

! BVOC variables
REAL, INTENT(OUT) ::                                              &
 isoprene(land_pts)                                               &
                   ! OUT Isoprene Emission Flux (kgC/m2/s)
,terpene(land_pts)                                                &
                   ! OUT (Mono-)Terpene Emission Flux (kgC/m2/s)
,methanol(land_pts)                                               &
                   ! OUT Methanol Emission Flux (kgC/m2/s)
,acetone(land_pts)
                   ! OUT Acetone Emission Flux (kgC/m2/s)


!  External routines called :-
EXTERNAL qsat,leaf,leaf_limits



REAL                                                              &
 anetc(land_pts)                                                  &
                            ! WORK Net canopy photosynthesis
!                                 !     (mol CO2/m2/s).
,apar_crit(land_pts)                                              &
                            ! WORK Critical APAR below which
!                                 !      light is limiting (W/m2)
,co2c(land_pts)                                                   &
                            ! WORK Canopy level CO2 concentration
!                                 !      (kg CO2/kg air).
,ci(land_pts)                                                     &
                            ! WORK Internal CO2 pressure (Pa).
,dq(land_pts)                                                     &
                            ! WORK Specific humidity deficit
!                                 !      (kg H2O/kg air).
,dqc(land_pts)                                                    &
                            ! WORK Canopy level specific humidity
!                                 !      deficit (kg H2O/kg air).
,fpar(land_pts)                                                   &
                            ! WORK PAR absorption factor.
,lai_bal(land_pts)                                                &
                            ! WORK Leaf area index in balanced
!                                 !      growth state.
,nl(land_pts)                                                     &
                            ! WORK Mean leaf nitrogen
!                                 !      concentration (kg N/kg C).
,nl_bal(land_pts)                                                 &
                            ! WORK Mean leaf nitrogen
!                                 !      concentration in balanced
!                                 !      growth state (kg N/kg C).
,n_leaf(land_pts)                                                 &
                            ! WORK Nitrogen contents of the leaf,
,n_root(land_pts)                                                 &
                            !      root,
,n_stem(land_pts)                                                 &
                            !      and stem (kg N/m2).
,qs(land_pts)                                                     &
                            ! WORK Saturated specific humidity
!                                 !      (kg H2O/kg air).
,ra_rc(land_pts)                                                  &
                            ! WORK Ratio of aerodynamic resistance
!                                 !      to canopy resistance.
,rdc(land_pts)                                                    &
                            ! WORK Canopy dark respiration,
!                                 !      without soil water dependence
!                                 !      (mol CO2/m2/s).
,resp_p_g(land_pts)                                               &
                            ! WORK Plant growth respiration rate
!                                 !      (kg C/m2/sec).
,resp_p_m(land_pts)                                               &
                            ! WORK Plant maintenance respiration
!                                 !      rate (kg C/m2/sec).
,root(land_pts)                                                   &
                            ! WORK Root carbon (kg C/m2).
,faparv_layer(land_pts,ilayers)                                   &
                            ! WORK absorbed par(layers)
,flux_o3_l(land_pts)                                              &
                            ! WORK Flux of O3 to stomata (nmol O3/m2/s).
,flux_o3_l_sun(land_pts)                                          &
                            ! WORK Flux of O3 to stomata
                            !      for sunlit leaves
                            !      (for can_rad_mod=5)
                            !      (nmol O3/m2/s).
,flux_o3_l_shd(land_pts)                                          &
                            ! WORK Flux of O3 to stomata
                            !      for shaded leaves
                            !      (for can_rad_mod=5)
                            !      (nmol O3/m2/s).
,fo3_l(land_pts)                                                  &
                            ! WORK Ozone exposure factor.
,fo3_l_sun(land_pts)                                              &
                            ! WORK Ozone exposure factor
                            !      for sunlit leaves
                            !      (for can_rad_mod=5)
,fo3_l_shd(land_pts)                                              &
                            ! WORK Ozone exposure factor
                            !      for shaded leaves
                            !      (for can_rad_mod=5)
,o3mol(land_pts)            ! WORK Surface ozone concentration (moles).


INTEGER                                                           &
 i,j,k,l,m,n                  ! WORK Loop counters.


!-----------------------------------------------------------------------
! Parameters
!-----------------------------------------------------------------------
REAL, PARAMETER      :: cconu = 12.0e-3
                            ! kg C in 1 mol CO2
!  (mol/sec) / (watts) conversion for PAR:
REAL, PARAMETER      ::  conpar = 2.19e5


REAL                                                              &
 anetl(land_pts)                                                  &
                            ! WORK Net leaf photosynthesis
!                                 !      (mol CO2/m2/s/LAI).
,anetl_sun(land_pts)                                              &
!                                 ! WORK Net leaf photosynthesis of
!                                 !      sunlit leaves
!                                 !      (mol CO2/m2/s/LAI)
,anetl_shd(land_pts)                                              &
!                                 ! WORK Net leaf photosynthesis of
!                                 !      shaded leaves
!                                 !      (mol CO2/m2/s/LAI
,apar(land_pts)                                                   &
                            ! WORK PAR absorbed by the top leaf
!                                 !      (W/m2).
,acr(land_pts)                                                    &
                            ! WORK Absorbed PAR
!                                 !      (mol photons/m2/s).
,ca(land_pts)                                                     &
                            ! WORK Canopy level CO2 pressure
!                                 !      (Pa).
,gl(land_pts)                                                     &
                            ! WORK Leaf conductance for H2O
!                                 !      (m/s).
,gl_sun(land_pts)                                                 &
                            ! WORK Leaf conductance for H2O of
!                                 !      sunlit leaves (m/s).
,gl_shd(land_pts)                                                 &
                            ! WORK Leaf conductance for H2O of
!                                 !      shaded leaves (m/s).
,oa(land_pts)                                                     &
                            ! WORK Atmospheric O2 pressure
!                                 !      (Pa).
,rd(land_pts)                                                     &
                            ! WORK Dark respiration of top leaf
!                                 !      (mol CO2/m2/s).
,rd_sun(land_pts)                                                 &
                            ! WORK Dark respiration of sunlit leaves
!                                 !      (mol CO2/m2/s).
,rd_shd(land_pts)                                                 &
                            ! WORK Dark respiration of shaded leaves
!                                 !      (mol CO2/m2/s).

,wcarb(land_pts)                                                  &
                            ! WORK Carboxylation, ...
,wlite(land_pts)                                                  &
                            !      ... Light, and ...
,wexpt(land_pts)                                                  &
                            !      ... export limited gross ...
!                                 !      ... photosynthetic rates ...
!                                 !      ... (mol CO2/m2/s).
,wlitev(land_pts)                                                 &
!                                 ! WORK Light limited gross
!                                 !      photosynthetic rates
!                                 !      for each layer
!                                 !      (mol CO2/m2/s).
,wlitev_sun(land_pts)                                             &
!                                 ! WORK Light limited gross
!                                 !      photosynthetic rates
!                                 !      for sunlit leaves
!                                 !      (mol CO2/m2/s).
,wlitev_shd(land_pts)
!                                 ! WORK Light limited gross
!                                 !      photosynthetic rates
!                                 !      for shaded leaves
!                                 !      (mol CO2/m2/s).

REAL                                                              &
 aparv(land_pts)                                                  &
                            ! WORK APAR for each leaf layer
!                                 !      (W/m2).
,apar_lit(land_pts)                                               &
                            ! WORK Mean APAR for non-light
!                                 !      limited leaves (W/m2/LAI).
,apar_unlit(land_pts)                                             &
                            ! WORK Mean APAR for light
!                                 !      limited leaves (W/m2/LAI).
,dlai(land_pts)                                                   &
                            ! WORK LAI Increment.
,lai_lit(land_pts)                                                &
                            ! WORK Total LAI of non-light
!                                 !      limited leaves.
,lai_unlit(land_pts)        ! WORK Total LAI of light
!                                 !      limited leaves.

INTEGER                                                           &
 clos_index(land_pts)                                             &
                            ! WORK Index of land points
!                                 !      with closed stomata.
,clos_pts                                                         &
                            ! WORK Number of land points
!                                 !      with closed stomata.
,open_index(land_pts)                                             &
                            ! WORK Index of land points
!                                 !      with open stomata.
,open_pts                   ! WORK Number of land points
!                                 !      with open stomata.

REAL                                                              &
 nl_tmp(land_pts)           ! WORK Temporary to pass down to leaf limits

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('SF_STOM',zhook_in,zhook_handle)

!-----------------------------------------------------------------------
! Initialisation.
!-----------------------------------------------------------------------
nl(:)        = 0.0
anetl(:)     = 0.0
anetc(:)     = 0.0
gl(:)        = 0.0
rdc(:)       = 0.0
ci(:)        = 0.0
fo3_l(:)     = 0.0
flux_o3_l(:) = 0.0
fo3(:)       = 0.0
flux_o3(:)   = 0.0
o3mol(:)     = 0.0

isoprene(:)  = 0.0
terpene(:)   = 0.0
methanol(:)  = 0.0
acetone(:)   = 0.0


SELECT CASE ( can_rad_mod )
  CASE ( 4,5 )
    gc(:) = 0.0
END SELECT

!-----------------------------------------------------------------------
! Set the canopy CO2 concentration.
!-----------------------------------------------------------------------
IF ( l_co2_interactive ) THEN
!       Use full 3D CO2 field.
  DO m=1,veg_pts
    l = veg_index(m)
    j=(land_index(l)-1)/t_i_length + 1
    i = land_index(l) - (j-1)*t_i_length
    co2c(l) = co2_3d(i,j)
  END DO
ELSE
!       Use single CO2_MMR value.
  DO m=1,veg_pts
    l = veg_index(m)
    co2c(l) = co2
  END DO
END IF

!-----------------------------------------------------------------------
! Calculate the surface to level 1 humidity deficit and the surface
! density of the air
!-----------------------------------------------------------------------
! DEPENDS ON: qsat
CALL qsat(qs,tstar,pstar,land_pts)
DO m=1,veg_pts
  l = veg_index(m)
  j=(land_index(l)-1)/t_i_length + 1
  i = land_index(l) - (j-1)*t_i_length
  dq(l) = MAX(0.0,(qs(l) - q1(i,j)))
END DO

!-----------------------------------------------------------------------
! Calculate the PAR absorption factor
!-----------------------------------------------------------------------
IF ( can_rad_mod == 1 ) THEN
  DO m=1,veg_pts
    l = veg_index(m)
    fpar(l) = (1. - EXP(-kpar(ft)*lai(l))) / kpar(ft)
  END DO
END IF

!-----------------------------------------------------------------------
! Calculate the PAR absorbed by the top leaf and set leaf N value.
!-----------------------------------------------------------------------
DO m =1,veg_pts
  l = veg_index(m)
  j=(land_index(l)-1)/t_i_length + 1
  i = land_index(l) - (j-1)*t_i_length
  apar(l) = (1. - omega(ft)) * ipar(i,j)
  nl(l) = nl0(ft)
END DO

!-----------------------------------------------------------------------
! Calculate the LAI in each canopy layer.
!-----------------------------------------------------------------------
SELECT CASE ( can_rad_mod )
  CASE ( 2:5 )
    DO m =1,veg_pts
      l = veg_index(m)
      dlai(l) = lai(l) / FLOAT(ilayers)
    END DO
END SELECT

!-----------------------------------------------------------------------
! Convert O3 concentration from ppb to moles
!-----------------------------------------------------------------------
DO m =1,veg_pts
  l = veg_index(m)
  o3mol(l) = o3(l) * pstar(l) / (rmol * tstar(l))
END DO

!-----------------------------------------------------------------------
! Iterate to ensure that the canopy humidity deficit is consistent with
! the H2O flux. Ignore the (small) difference between the canopy and
! reference level CO2 concentration. Intially set the canopy humidity
! deficit using the previous value of GC.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! The new canopy light scheme with adjusted dark respiration and
! nitrogen profile required the iteration to be undertaken for each
! layer, unlike the other options.
! Therefore this option needs a differnet logical pathway.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! can_rad_mod=4 and 5 require iteration for each layer and hence follow
! a different path through the code.
!-----------------------------------------------------------------------

IF ( can_rad_mod == 4 ) THEN

!-----------------------------------------------------------------------
!       Varying N model+altered leaf respiration
!       Multiple canopy layers
!       N varies through canopy as exponential
!-----------------------------------------------------------------------

  DO m=1,veg_pts
    l = veg_index(m)
    wlitev(l) = wlite(l)
  END DO

  DO n=1,ilayers
!-----------------------------------------------------------------------
! Initialise GL and calculate the PAR absorbed in this layer.
!-----------------------------------------------------------------------
    DO m=1,veg_pts
      l = veg_index(m)
      gl(l) = 0.0
      nl_tmp(l)=nl(l)*EXP((n-1)/FLOAT(ilayers)*(-kn))
      faparv_layer(l,n) = faparv(l,n) * dlai(l)
    END DO

!-----------------------------------------------------------------------
! Iterate to ensure that the canopy humidity deficit is consistent with
! the H2O flux. Ignore the (small) difference between the canopy and
! reference level CO2 concentration. Intially set the canopy humidity
! deficit to using GL=0.
!-----------------------------------------------------------------------

    DO k=1,iter

!-----------------------------------------------------------------------
! Diagnose the canopy level humidity deficit and CO2 concentration
!-----------------------------------------------------------------------
      DO m=1,veg_pts
        l = veg_index(m)
        ra_rc(l) = ra(l) * gl(l)
        dqc(l) = dq(l) / (1.0 + ra_rc(l))
      END DO

!-----------------------------------------------------------------------
! Calculate the canopy resistance and photosynthesis
!-----------------------------------------------------------------------
      DO m =1,veg_pts
        l = veg_index(m)
        ca(l) = co2c(l) / epco2 * pstar(l)
        oa(l) = o2 / epo2 * pstar(l)
      END DO

!-----------------------------------------------------------------------
! Calculate the limiting factors for leaf photosynthesis
!-----------------------------------------------------------------------
! The NL*exp((N-1)/FLOAT(ILAYERS)*(-KN)) term is the non-uniform
! distribution for nitrogen in the canopy


! DEPENDS ON: leaf_limits
      CALL leaf_limits (land_pts,veg_pts,veg_index,ft             &
,                       nl_tmp                                    &
,                       dqc,apar,tstar,ca,oa,pstar,fsmc           &
,                       clos_pts,open_pts,clos_index,open_index   &
,                       ci,rd,wcarb,wexpt,wlite)

      DO m=1,open_pts
        l = veg_index(open_index(m))
        j=(land_index(l)-1)/t_i_length + 1
        i = land_index(l) - (j-1)*t_i_length
        wlitev(l)=wlite(l)/apar(l)*faparv(l,n)*ipar(i,j)
        acr(l) = apar(l) / conpar
        IF (acr(l)*1.e6 *faparv_layer(l,n) >  10.)                &
          rd(l)=( 0.5-0.05 *                                      &
                LOG(acr(l)*faparv_layer(l,n)*1.e6))*rd(l)
      END DO

! DEPENDS ON: leaf
      CALL leaf (land_pts,veg_pts,veg_index,ft                    &
,                clos_pts,open_pts,clos_index,open_index          &
,                o3mol,ra,flux_o3_l                               &
,                fsmc,tstar,ca,ci,rd,wcarb,wexpt,wlitev           &
,                gl,anetl,fo3_l)

    END DO                 ! K-ITER

    DO m=1,veg_pts
      l = veg_index(m)
      
      anetc(l) = anetc(l) + anetl(l) * dlai(l)
      gc(l) = gc(l) + gl(l) * dlai(l)
      rdc(l) = rdc(l) + rd(l) * dlai(l)
      
      flux_o3(l) = flux_o3(l) + flux_o3_l(l) * dlai(l)
      fo3(l) = fo3(l) + fo3_l(l) * dlai(l)
    END DO

  END DO                   ! N LAYERS


ELSE IF ( can_rad_mod == 5 ) THEN

!-----------------------------------------------------------------------
!       Sunlit and shaded leaves treated separately
!       Multiple canopy layers
!       N varies through canopy as exponential
!-----------------------------------------------------------------------

  DO n=1,ilayers

!-----------------------------------------------------------------------
! Initialise GL for this layer.
! We could initialise to gl(n-1) here, but simpler to use zero
! and seems to converge pretty quickly anyway.
! The NL*exp((N-1)/FLOAT(ILAYERS)*(-KN)) term is the non-uniform
! distribution for nitrogen in the canopy.
!-----------------------------------------------------------------------
    DO m=1,veg_pts
      l = veg_index(m)
      gl(l)=0.0
      nl_tmp(l)=nl(l)*EXP((n-1)/FLOAT(ilayers)*(-kn))
    END DO

!-----------------------------------------------------------------------
! Iterate to ensure that the canopy humidity deficit is consistent with
! the H2O flux. Ignore the (small) difference between the canopy and
! reference level CO2 concentration.
!-----------------------------------------------------------------------

    DO k=1,iter

!-----------------------------------------------------------------------
! Diagnose the canopy level humidity deficit and CO2 concentration
!-----------------------------------------------------------------------
      DO m=1,veg_pts
        l = veg_index(m)
        ra_rc(l) = ra(l) * gl(l)
        dqc(l) = dq(l) / (1.0 + ra_rc(l))
      END DO

!-----------------------------------------------------------------------
! Calculate the canopy resistance and photosynthesis
!-----------------------------------------------------------------------
      DO m =1,veg_pts
        l = veg_index(m)
        ca(l) = co2c(l) / epco2 * pstar(l)
        oa(l) = o2 / epo2 * pstar(l)
      END DO

!-----------------------------------------------------------------------
! Calculate the limiting factors for leaf photosynthesis
!-----------------------------------------------------------------------
! DEPENDS ON: leaf_limits
      CALL leaf_limits (land_pts,veg_pts,veg_index,ft             &
,                       nl_tmp                                    &
,                       dqc,apar,tstar,ca,oa,pstar,fsmc           &
,                       clos_pts,open_pts,clos_index,open_index   &
,                       ci,rd,wcarb,wexpt,wlite)
    
!-----------------------------------------------------------------------
! Set the sunlit and shaded respiration rates to the uninhibited, sunlit
! respiration rate from leaf_limits.
!-----------------------------------------------------------------------
      DO m=1,veg_pts        
        l = veg_index(m)    
        rd_sun(l) = rd(l)   
        rd_shd(l) = rd(l)   
      END DO

      DO m=1,open_pts
        l = veg_index(open_index(m))
        j=(land_index(l)-1)/t_i_length + 1
        i = land_index(l) - (j-1)*t_i_length
        wlitev_sun(l) = wlite(l)/apar(l)*fapar_sun(l,n)*ipar(i,j)
        wlitev_shd(l) = wlite(l)/apar(l)*fapar_shd(l,n)*ipar(i,j)

        acr(l) = ipar(i,j) / conpar
!-----------------------------------------------------------------------
! Introducing inhibition of leaf respiration in the light for sunlit and
! shaded leaves, from papers by Atkin et al. This is an improvement over
! the description used for can_rad_mod=4.    
!-----------------------------------------------------------------------
        IF ( fapar_sun(l,n) * acr(l) * fsun(l,n) * dlai(l) *           &
             1.0e6 >  10.0 ) rd_sun(l) = 0.7 * rd_sun(l)    
        IF ( fapar_shd(l,n) * acr(l) * (1.0-fsun(l,n)) * dlai(l) *     &
             1.0e6 >  10.0 ) rd_shd(l) = 0.7 * rd_shd(l)    
      END DO

!-----------------------------------------------------------------------
! Call leaf routine separately for sunlit and shaded leaves.
!-----------------------------------------------------------------------
! DEPENDS ON: leaf
      CALL leaf (land_pts,veg_pts,veg_index,ft                    &
,                clos_pts,open_pts,clos_index,open_index          &
,                o3mol,ra,flux_o3_l_sun                           &
,                fsmc,tstar,ca,ci,rd_sun,wcarb,wexpt,wlitev_sun   &
,                gl_sun,anetl_sun,fo3_l_sun)

! DEPENDS ON: leaf
      CALL leaf (land_pts,veg_pts,veg_index,ft                    &
,                clos_pts,open_pts,clos_index,open_index          &
,                o3mol,ra,flux_o3_l_shd                           &
,                fsmc,tstar,ca,ci,rd_shd,wcarb,wexpt,wlitev_shd   &
,                gl_shd,anetl_shd,fo3_l_shd)


!           Update layer conductance.
      DO m=1,veg_pts
        l = veg_index(m)
        gl(l) = fsun(l,n)*gl_sun(l) + (1.-fsun(l,n))*gl_shd(l)
        rd(l) = fsun(l,n)*rd_sun(l) + (1.-fsun(l,n))*rd_shd(l)
      END DO

    END DO                 ! K-ITER

    DO m=1,veg_pts
      l = veg_index(m)
      
      anetl(l) = fsun(l,n) * anetl_sun(l)                         &
               + ( 1.0-fsun(l,n) ) * anetl_shd(l)
      anetc(l) = anetc(l) + anetl(l) * dlai(l)
      
      gc(l) = gc(l) + gl(l) * dlai(l)
      rdc(l) = rdc(l)+ rd(l) * dlai(l)
      
      flux_o3_l(l) = fsun(l,n) * flux_o3_l_sun(l)                 &
                   + ( 1.0 - fsun(l,n) ) * flux_o3_l_shd(l)
      fo3_l(l) = fsun(l,n) * fo3_l_sun(l)                         &
               + ( 1.0 - fsun(l,n) ) * fo3_l_shd(l)
      
      flux_o3(l) = flux_o3(l) + flux_o3_l(l) * dlai(l)
      fo3(l) = fo3(l) + fo3_l(l) * dlai(l)
    END DO

  END DO                   ! N LAYERS

!-----------------------------------------------------------------------
ELSE                      ! CAN_RAD_MOD

!-----------------------------------------------------------------------
! can_rad_mod=1, 2 or 3
! Iterate to ensure that the canopy humidity deficit is consistent with
! the H2O flux. Ignore the (small) difference between the canopy and
! reference level CO2 concentration. Intially set the canopy humidity
! deficit using the previous value of GC.
!-----------------------------------------------------------------------

  DO k=1,iter

!-----------------------------------------------------------------------
! Diagnose the canopy level humidity deficit and CO2 concentration
!-----------------------------------------------------------------------
    DO m=1,veg_pts
      l = veg_index(m)
      ra_rc(l) = ra(l) * gc(l)
      dqc(l) = dq(l) / (1.0 + ra_rc(l))
    END DO

!-----------------------------------------------------------------------
! Calculate the canopy resistance and photosynthesis
!-----------------------------------------------------------------------
    DO m =1,veg_pts
      l = veg_index(m)
      ca(l) = co2c(l) / epco2 * pstar(l)
      oa(l) = o2 / epo2 * pstar(l)
    END DO

!-----------------------------------------------------------------------
! Calculate the limiting factors for leaf photosynthesis
!-----------------------------------------------------------------------
! DEPENDS ON: leaf_limits
    CALL leaf_limits (land_pts,veg_pts,veg_index,ft               &
,                     nl,dqc,apar,tstar,ca,oa,pstar,fsmc          &
,                     clos_pts,open_pts,clos_index,open_index     &
,                     ci,rd,wcarb,wexpt,wlite)

    DO m=1,veg_pts
      l = veg_index(m)
      apar_crit(l) = 0.0
    END DO

! Treat closed stomata as unlit when using can_rad_mod==3.
    DO j=1,clos_pts
      l = veg_index(clos_index(j))
      apar_crit(l) = 1.0E20
    ENDDO

!-----------------------------------------------------------------------
! Calculate leaf level quantities and scale up to canopy through user
! requested approach
!-----------------------------------------------------------------------
    IF ( can_rad_mod == 2 ) THEN

!-----------------------------------------------------------------------
!           Multiple canopy layers
!           N constant through canopy
!-----------------------------------------------------------------------

      DO m=1,veg_pts
        l = veg_index(m)
        anetc(l) = 0.0
        gc(l) = 0.0
        wlitev(l) = wlite(l)
      END DO

      DO n=1,ilayers
        DO m=1,open_pts
          l = veg_index(open_index(m))
          j=(land_index(l)-1)/t_i_length + 1
          i = land_index(l) - (j-1)*t_i_length

          wlitev(l)=wlite(l)/apar(l)*faparv(l,n)*ipar(i,j)
        END DO

! DEPENDS ON: leaf
        CALL leaf (land_pts,veg_pts,veg_index,ft                  &
,                  clos_pts,open_pts,clos_index,open_index        &
,                  o3mol,ra,flux_o3_l                             &
,                  fsmc,tstar,ca,ci,rd,wcarb,wexpt,wlitev         &
,                  gl,anetl,fo3_l)

        DO m=1,veg_pts
          l = veg_index(m)
          
          anetc(l) = anetc(l) + anetl(l) * dlai(l)
          gc(l) = gc(l) + gl(l) * dlai(l)
          
          flux_o3(l) = flux_o3(l) + flux_o3_l(l) * dlai(l)
          fo3(l) = fo3(l) + fo3_l(l) * dlai(l)
        END DO
      END DO

    ELSE IF ( can_rad_mod == 3 ) THEN

!-----------------------------------------------------------------------
!           sunlit and shaded leaves treated separately
!           N constant through canopy
!-----------------------------------------------------------------------
!           Initialise.
      DO m=1,veg_pts
        l = veg_index(m)
        anetc(l) = 0.0
        gc(l) = 0.0
        wlitev(l) = wlite(l)
        lai_lit(l) = 0.0
        lai_unlit(l) = 0.0
        apar_lit(l) = 0.0
        apar_unlit(l) = 0.0
      END DO

      DO j=1,open_pts
        l = veg_index(open_index(j))
! Not sure why 10^-10 was chosen here - will find out from Doug at CEH
! Probably something to do with physically realistic values
! Doug says: sort of. It's a recurring issue, don't know answer.
! We want to avoid having a "very small" number on the denominator.
! EPSILON is not that small for a 32-bit real. TINY is a bit too small.
        IF ( wlite(l) < 1.0e-10 .AND. wcarb(l) < 1.0e-10 ) THEN
          apar_crit(l) = 0.0
        ELSE
          apar_crit(l) = MIN(wcarb(l),wexpt(l)) * apar(l)/wlite(l)
        END IF
      END DO

!-----------------------------------------------------------------------
! Divide LAI into sunlit and shaded layers.  The sum,
!   lai_lit + lai_unlit = lai
! must always hold else plant-level GPP, NPP and respiration rates will
! be diagnosed incorrectly towards the end of this subroutine.
!-----------------------------------------------------------------------
      DO n=1,ilayers
        DO m=1,veg_pts
          l = veg_index(m)
          j=(land_index(l)-1)/t_i_length + 1
          i = land_index(l) - (j-1)*t_i_length

          aparv(l) = faparv(l,n) * ipar(i,j)

          IF (aparv(l) < apar_crit(l)) THEN
            apar_unlit(l) = apar_unlit(l)                         &
                        + aparv(l)*dlai(l)
            lai_unlit(l) = lai_unlit(l) + dlai(l)
          ELSE
            apar_lit(l) = apar_lit(l)                             &
                        + aparv(l)*dlai(l)
            lai_lit(l) = lai_lit(l) + dlai(l)
          END IF
        END DO
      END DO

! Calculate net photosynthesis and conductance for sunlit leaves.
      DO m=1,open_pts
        l = veg_index(open_index(m))
        IF (lai_lit(l) > 0.0)                                     &
          apar_lit(l) = apar_lit(l) / lai_lit(l)
        wlitev(l)=wlite(l)/apar(l)*apar_lit(l)
      END DO

! DEPENDS ON: leaf
      CALL leaf (land_pts,veg_pts,veg_index,ft                    &
,                clos_pts,open_pts,clos_index,open_index          &
,                o3mol,ra,flux_o3_l                               &
,                fsmc,tstar,ca,ci,rd,wcarb,wexpt,wlitev           &
,                gl,anetl,fo3_l)

      DO m=1,veg_pts
        l = veg_index(m)

        anetc(l) = anetc(l) + anetl(l) * lai_lit(l)
        gc(l) = gc(l) + gl(l) * lai_lit(l)

        flux_o3(l) = flux_o3(l) + flux_o3_l(l) * lai_lit(l)
        fo3(l) = fo3(l) + fo3_l(l) * lai_lit(l)
      END DO

! Calculate net photosynthesis and conductance for shaded leaves.
      DO m=1,open_pts
        l = veg_index(open_index(m))
        IF (lai_unlit(l) > 0.0)                                   &
          apar_unlit(l) = apar_unlit(l) / lai_unlit(l)
        wlitev(l)=wlite(l)/apar(l)*apar_unlit(l)
      END DO

! DEPENDS ON: leaf
      CALL leaf (land_pts,veg_pts,veg_index,ft                    &
,                clos_pts,open_pts,clos_index,open_index          &
,                o3mol,ra,flux_o3_l                               &
,                fsmc,tstar,ca,ci,rd,wcarb,wexpt,wlitev           &
,                gl,anetl,fo3_l)

      DO m=1,veg_pts
        l = veg_index(m)

        anetc(l) = anetc(l) + anetl(l) * lai_unlit(l)
        gc(l) = gc(l) + gl(l) * lai_unlit(l)

        flux_o3(l) = flux_o3(l) + flux_o3_l(l) * lai_unlit(l)
        fo3(l) = fo3(l) + fo3_l(l) * lai_unlit(l)
      END DO

    ELSE  !  can_rad_mod

!-----------------------------------------------------------------------
!           can_rad_mod = 1
!           "big leaf" model
!           N varies through canopy according to Beers Law
!-----------------------------------------------------------------------

! DEPENDS ON: leaf
      CALL leaf (land_pts,veg_pts,veg_index,ft                    &
,                clos_pts,open_pts,clos_index,open_index          &
,                o3mol,ra,flux_o3_l                               &
,                fsmc,tstar,ca,ci,rd,wcarb,wexpt,wlite            &
,                gl,anetl,fo3_l)

      DO m=1,veg_pts
        l = veg_index(m)

        anetc(l) = anetl(l) * fpar(l)
        gc(l) = fpar(l) * gl(l)
        rdc(l) = rd(l) * fpar(l)

        flux_o3(l)= flux_o3_l(l) * fpar(l)
        fo3(l)= fo3_l(l) * fpar(l)
      END DO

    END IF  !  can_rad_mod (= 2 or 3)

  END DO   ! End of iteration loop

END IF  ! can_rad_mod (= 4 or 5)

IF ( can_rad_mod == 3 ) THEN
!-----------------------------------------------------------------------
!       Impose glmin at closed points.
!-----------------------------------------------------------------------
  DO m=1,clos_pts
    l = veg_index(clos_index(m))
    gc(l) = glmin(ft)
  END DO
END IF

!-----------------------------------------------------------------------
!     Calculate plant level respiration, NPP and GPP
!-----------------------------------------------------------------------

SELECT CASE ( can_rad_mod )
  CASE ( 2,3 )
    DO m=1,veg_pts
      l = veg_index(m)
      rdc(l) = rd(l) * lai(l)
    END DO
END SELECT

DO m=1,veg_pts
  l = veg_index(m)

!-----------------------------------------------------------------------
! Assume that root biomass is equal to balanced growth leaf biomass
!-----------------------------------------------------------------------
  lai_bal(l) = (a_ws(ft)*eta_sl(ft)*ht(l)/a_wl(ft))             &
             **(1.0/(b_wl(ft)-1.0))
  root(l) = sigl(ft) * lai_bal(l)

!-----------------------------------------------------------------------
! Calculate the actual and balanced mean leaf nitrogen concentration
! assuming perfect light acclimation
!-----------------------------------------------------------------------
  nl(l) = nl0(ft)
  nl_bal(l) = nl0(ft)

!---------------------------------------------------------------
! Calculate the total nitrogen content of the leaf, root and stem
!-----------------------------------------------------------------------
  n_leaf(l) = nl(l) * sigl(ft) * lai(l)
  n_root(l) = nr_nl(ft) * nl_bal(l) * root(l)
  n_stem(l) = ns_nl(ft) * nl_bal(l) * eta_sl(ft) * ht(l) * lai(l)

!-----------------------------------------------------------------------
! Calculate the Gross Primary Productivity, the plant maintenance
! respiration rate, and the wood maintenance respiration rate
! in kg C/m2/sec
!-----------------------------------------------------------------------
  gpp(l) = cconu * (anetc(l) + rdc(l)*fsmc(l))
  resp_p_m(l) = cconu * rdc(l)                                  &
       * (n_leaf(l)*fsmc(l) + n_stem(l) + n_root(l)) / n_leaf(l)
  resp_w(l) = cconu * rdc(l) * n_stem(l) / n_leaf(l)

!-----------------------------------------------------------------------
! Calculate the total plant respiration and the Net Primary Productivity
!-----------------------------------------------------------------------
  resp_p_g(l) = r_grow(ft) * (gpp(l) - resp_p_m(l))
  resp_p(l) = resp_p_m(l) + resp_p_g(l)
  npp(l) = gpp(l) - resp_p(l)

END DO

!-----------------------------------------------------------------------
! Calculate BVOC emissions
!-----------------------------------------------------------------------
IF ( l_bvoc_emis )                                                        &
! DEPENDS ON: bvoc_emissions
  CALL bvoc_emissions(land_pts,veg_pts,ft,veg_index,ci,anetc,rdc,tstar,   &
                      isoprene,terpene,methanol,acetone)


IF (lhook) CALL dr_hook('SF_STOM',zhook_out,zhook_handle)
RETURN
END SUBROUTINE sf_stom

!***********************************************************************
! Calculates the canopy resistance, net photosynthesis and transpiration
! by scaling-up the leaf level response using the "Big-Leaf" approach
! of Sellers et al. (1994)

! Written by Peter Cox (May 1995)
!***********************************************************************
