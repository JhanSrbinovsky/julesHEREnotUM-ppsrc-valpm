! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

SUBROUTINE bvoc_emissions(land_pts,veg_pts,ft,veg_index,ci,al,rd,tstar,       &
                          isoprene,terpene,methanol,acetone)

  USE c_bvoc

  USE c_0_dg_c, ONLY : zerodegc

  USE surf_param, ONLY : o2,q10_leaf,fwe_c3,fwe_c4,beta1,beta2

  USE ccarbon, ONLY : epo2,epco2

  USE pftparm, ONLY : c3,neff,nl0,tupp,tlow,fd,f0,alpha,ief,tef,mef,aef,sigl

  USE parkind1, ONLY: jprb, jpim
  USE yomhook, ONLY: lhook, dr_hook

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Calculates emission of BVOCs (see Pacifico et al. (2011) Atm. Chem Phys.)
!
!   Also calculates emissions of (mono)terpenes, methanol and acetone.
!   Based on model developed by
!   Guenther et al. JGR 1995, V.100(D5) 8873-8892.
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Arguments
  INTEGER, INTENT(IN) ::                                                      &
    land_pts,                                                                 &
                   ! IN Total number of land points
    veg_pts,                                                                  &
                   ! IN Number of vegetated points
    veg_index(land_pts),                                                      &
                   ! IN Index of vegetated points
    ft
                   ! IN Plant functional type.

  REAL, INTENT(IN) ::                                                         &
    ci(land_pts),                                                             &
                   ! IN Internal CO2 pressure (Pa)
    al(land_pts),                                                             &
                   ! IN Net Leaf photosynthesis (mol CO2/m2/s)
    rd(land_pts),                                                             &
                   ! IN Dark respiration (mol CO2/m2/s)
    tstar(land_pts)
                   ! IN Leaf temperature (K) = Surface temperature (K)

  REAL, INTENT(OUT) ::                                                        &
    isoprene(land_pts),                                                       &
                   ! OUT Isoprene Emission Flux (kgC/m2/s)
    terpene(land_pts),                                                        &
                   ! OUT (Mono-)Terpene Emission Flux (kgC/m2/s)
    methanol(land_pts),                                                       &
                   ! OUT Methanol Emission Flux (kgC/m2/s)
    acetone(land_pts)
                   ! OUT Acetone Emission Flux (kgC/m2/s)


! Parameters
  REAL ::                                                                     &
    tdegc_st,                                                                 &
                    ! Temperature in standard conditions = 30 Celsius degrees
    tau_st,                                                                   &
    oa_st,                                                                    &
    ccp_st,                                                                   &
    vcmax_st,                                                                 &
    qtenf_st,                                                                 &
    denom_st,                                                                 &
    vcm_st,                                                                   &
    rd_st,                                                                    &
    co2_st,                                                                   &
    ca_st,                                                                    &
    ci_st,                                                                    &
    kc_st,                                                                    &
    ko_st,                                                                    &
    acr_st,                                                                   &
                    ! Absorbed PAR in standard conditions = 1000 micro mol photons/m2/s
    pstar_st,                                                                 &
                    ! Atmospheric pressure in standard conditions = 101325 Pa
    wcarb_st,                                                                 &
    wlite_st,                                                                 &
    wexpt_st,                                                                 &
    a1,                                                                       &
    a2,                                                                       &
    a3,                                                                       &
    b1,                                                                       &
    b2,                                                                       &
    b3,                                                                       &
    wp_st,                                                                    &
    wl_st,                                                                    &
    al_st


! Work variables
  REAL ::                                                                     &
    f_co2(land_pts),                                                          &
                   ! WORK CO2 scaling factor
    f_t_isop(land_pts),                                                       &
                   ! WORK Temperature scaling factor for isoprene
    f_t_terp(land_pts)
                   ! WORK Temperature scaling factor for isoprene

  INTEGER :: m,l  ! WORK Loop counters


  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle


!-----------------------------------------------------------------------------


  IF (lhook) CALL dr_hook('BVOC_EMISSIONS',zhook_in,zhook_handle)
  
!-----------------------------------------------------------------------------
! Initialise arrays
!-----------------------------------------------------------------------------
  f_co2(:)    = 0.0
  f_t_isop(:) = 0.0
  f_t_terp(:) = 0.0
  
  isoprene(:) = 0.0
  terpene(:)  = 0.0
  methanol(:) = 0.0
  acetone(:)  = 0.0


!-----------------------------------------------------------------------      
! Initialize parameter values for emission calculation
!
! Calculates AL_ST and RD_ST in standard conditions:
! T=30 Celsius degree, Light=ACR=1000 micro mol photons/m2/s, Pressure=101325 Pa,
! (Atm CO2 Conc)=CA=370ppmv, (Canopy level specific humidity deficit)=DQ=0
!-----------------------------------------------------------------------     
  acr_st   = 1.e-3
  pstar_st = 101325.0

  tdegc_st = t_ref - zerodegc
  tau_st   = 2600.0 * (0.57 ** (0.1 * (tdegc_st - 25.0)))
  oa_st    = o2 / epo2 * pstar_st
  ccp_st   = 0.5 * oa_st / tau_st * REAL(c3(ft))
  vcmax_st = neff(ft) * nl0(ft)
  qtenf_st = vcmax_st * (q10_leaf ** (0.1 * (tdegc_st - 25.0)))
  denom_st = (1 + EXP(0.3 * (tdegc_st - tupp(ft))))                           &
           * (1 + EXP(0.3 * (tlow(ft) - tdegc_st)))
  vcm_st   = qtenf_st / denom_st
  rd_st    = fd(ft) * vcm_st
  co2_st   = 370.0 * 44.0 / 28.97 * 1.e-6  ! convert CO2 from ppmv to mass mixing ratio
  ca_st    = co2_st / epco2 * pstar_st
  ci_st    = (ca_st - ccp_st) * f0(ft) + ccp_st

  IF ( c3(ft) == 1 ) THEN
    kc_st    = 30.0 * (2.1 ** (0.1 * (tdegc_st - 25.0)))
    ko_st    = 30000.0 * (1.2 ** (0.1 * (tdegc_st - 25.0)))
    wcarb_st = vcm_st * (ci_st - ccp_st)                                      &
             / (ci_st + kc_st * (1.0 + oa_st / ko_st))
    wlite_st = alpha(ft) * acr_st * (ci_st - ccp_st) / (ci_st + 2.0 * ccp_st)
    wexpt_st = fwe_c3 * vcm_st
  ELSE
    wcarb_st = vcm_st
    wlite_st = alpha(ft) * acr_st
    wexpt_st = fwe_c4 * vcm_st * ci_st / pstar_st
  ENDIF

  a1    = beta1
  a2    = -(wcarb_st + wlite_st)
  a3    = wcarb_st * wlite_st
  wp_st = -a2 / (2.0 * a1) - sqrt(a2 * a2 / (4.0 * a1 * a1) - a3 / a1)
  
  b1    = beta2
  b2    = -(wp_st + wexpt_st)
  b3    = wp_st * wexpt_st
  wl_st = -b2 / (2.0 * b1) - sqrt(b2 * b2 / (4.0 * b1 * b1) - b3 / b1)
  al_st = (wl_st - rd_st)


!-----------------------------------------------------------------------------
! Calculate BVOC emissions
!-----------------------------------------------------------------------------  
  DO m = 1,veg_pts
    l = veg_index(m)
    
!-----------------------------------------------------------------------------      
! CO2 factor           
!-----------------------------------------------------------------------------
    f_co2(l) = ci_st / ci(l)

!-----------------------------------------------------------------------------
! Temperature factors (eq A4b Arneth et al., 2007)
!-----------------------------------------------------------------------------
    f_t_isop(l) = MIN(f_tmax, EXP(atau * (tstar(l) - t_ref)))
    f_t_terp(l) =             EXP(btau * (tstar(l) - t_ref))

!-----------------------------------------------------------------------------
! Isoprene emission
!-----------------------------------------------------------------------------
    isoprene(l) = MAX(0.0, ief(ft) * sigl(ft) *                               &
                         ! convert kgc/m2 leaf to g_dw/m2 leaf
                           2.0e+03 *                                          &
                         ! convert micro gc to kgc
                           1.0e-09 *                                          &
                         ! convert hours to secs
                           1.0 / 3600.0 *                                     &
                           (al(l) + rd(l)) / (al_st + rd_st) *                &
                           f_t_isop(l) * f_co2(l))
                        
!-----------------------------------------------------------------------------
! (Mono-)Terpene emission (Niinemets et al.,1999;Arneth et al.,2007)
!-----------------------------------------------------------------------------
    terpene(l) = MAX(0.0, tef(ft) * sigl(ft) *                                &
                        ! convert kgc/m2 leaf to g_dw/m2 leaf
                          2.0e+03 *                                           &
                        ! convert micro gc to kgc
                          1.0e-09 *                                           &
                        ! convert hours to secs
                          1.0 / 3600.0 *                                      &
                          f_t_terp(l))

!-----------------------------------------------------------------------------
! Methanol emission
!-----------------------------------------------------------------------------
    methanol(l) = MAX(0.0, mef(ft) * sigl(ft) *                               &
                         ! convert kgc/m2 leaf to g_dw/m2 leaf
                           2.0e+03 *                                          &
                         ! convert micro gc to kgc
                           1.0e-09 *                                          &
                         ! convert hours to secs
                           1.0 / 3600.0 *                                     &
                           f_t_terp(l))

!-----------------------------------------------------------------------------
! Acetone emission
!-----------------------------------------------------------------------------
    acetone(l) = MAX(0.0, aef(ft) * sigl(ft) *                                &
                        ! convert kgc/m2 leaf to g_dw/m2 leaf
                          2.0e+03 *                                           &
                        ! convert micro gc to kgc
                          1.0e-09 *                                           &
                        ! convert hours to secs
                          1.0 / 3600.0 *                                      &
                          f_t_terp(l))
      
  END DO

  IF (lhook) CALL dr_hook('BVOC_EMISSIONS',zhook_out,zhook_handle)
      
  RETURN

END SUBROUTINE bvoc_emissions
