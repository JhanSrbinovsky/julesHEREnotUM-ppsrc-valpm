
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine ALLOCATE_JULES_ARRAYS
!
! Description: Routine that allocates memory to the JULES arrays
!     This assume that the values in the NSTYPES module have been set
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.
!
!   Code Owner: See Unified Model Code Owner's HTML page
!   This file belongs in section: Land

SUBROUTINE allocate_jules_arrays( land_field,ntiles,sm_levels,nice,nice_use)

USE nstypes

USE atm_fields_bounds_mod
USE theta_field_sizes, ONLY : t_i_length, t_j_length

USE switches, ONLY : l_phenol                                         &
                    ,l_triffid                                        &
                    ,l_flake_model                                    & 
                    ,l_spec_albedo                                    & 
                    ,l_albedo_obs 

USE ancil_info
USE prognostics
USE jules_mod, ONLY : snowdep_surf, albobs_scaling

USE c_z0h_z0m
USE c_elevate

USE soil_param, ONLY:                                                 &
    dzsoil

USE snow_param, ONLY:                                                 &
    cansnowtile,dzsnow,ds

USE surf_param, ONLY:                                                 &
    diff_frac

USE nvegparm
USE pftparm
USE trif

USE lake_mod

USE switches_urban, ONLY : l_urban2t

USE urban_param, ONLY :                                               &
   hgt, hwr, wrr, disp, ztm, albwl, albrd, emisw, emisr

USE fluxes, ONLY :                                                    &
   anthrop_heat, surf_ht_store, sw_sice, sw_sice_rts, alb_sice

USE ozone_vars, ONLY : flux_o3_ft, fo3_ft

USE bvoc_vars, ONLY : isoprene, isoprene_ft, terpene, terpene_ft,     &
                      methanol, methanol_ft, acetone, acetone_ft

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE PrintStatus_mod
IMPLICIT NONE

!-----------------------------------------------------------------------
! Input variables for dimensioning
!-----------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                &
  land_field,                                                         &
                       ! Number of land points
  ntiles,                                                             &
                       ! Number of surface tiles
  sm_levels,                                                          &
                       ! Number of soil layers
  nice,                                                               &
                       ! Number of sea ice categories
  nice_use  
                       ! Number of sea ice cats used in radiation and
                       !  explicit part of surface exchange

!-----------------------------------------------------------------------
! Local variables for error trapping
!-----------------------------------------------------------------------
INTEGER ::                                                            &
  check_stat = 0,                                                     &
                       ! Variable for trapping the error from each
                       ! individual call to allocate
  err_sum    = 0

                       ! Variable to track the sum of all errors
                       ! resulting from calls to allocate. Hence we
                       ! know that everything was successful if and
                       ! only if this is zero at the end

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


!
! START OF EXECUTABLE CODE
!

IF (lhook) CALL dr_hook('ALLOCATE_JULES_ARRAYS',zhook_in,zhook_handle)

!-----------------------------------------------------------------------
! Compute length of theta field in i and j directions.
! (The module keeps these values for future use, this is
!  earliest place in the code that they are needed.)
!-----------------------------------------------------------------------
t_i_length = tdims%i_end - tdims%i_start + 1
t_j_length = tdims%j_end - tdims%j_start + 1

!-----------------------------------------------------------------------
! Allocate space for arrays from c_z0h_z0m and c_elevate
!-----------------------------------------------------------------------
ALLOCATE( z0h_z0m(ntype), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( z0h_z0m_classic(ntype), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( surf_hgt(land_field,ntiles), stat = check_stat )
err_sum = err_sum + check_stat

!-----------------------------------------------------------------------
! Allocate space for arrays from soil_param
!-----------------------------------------------------------------------
ALLOCATE( dzsoil(sm_levels), stat = check_stat )
err_sum = err_sum + check_stat

!-----------------------------------------------------------------------
! Allocate space for arrays from surf_param
!-----------------------------------------------------------------------
ALLOCATE( diff_frac(t_i_length * t_j_length), stat = check_stat )
err_sum = err_sum + check_stat

!-----------------------------------------------------------------------
! Allocate arrays from nvegparm
!-----------------------------------------------------------------------
ALLOCATE( albsnc_nvg(nnvg), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( albsnf_nvgu(nnvg), STAT = check_stat ) 
err_sum = err_sum + check_stat 

ALLOCATE( albsnf_nvg(nnvg), STAT = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( albsnf_nvgl(nnvg), STAT = check_stat ) 
err_sum = err_sum + check_stat 

ALLOCATE( catch_nvg(nnvg), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( gs_nvg(nnvg), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( infil_nvg(nnvg), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( z0_nvg(nnvg), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( ch_nvg(nnvg), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( vf_nvg(nnvg), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( emis_nvg(nnvg), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( z0hm_nvg(nnvg), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( z0hm_classic_nvg(nnvg), stat = check_stat )
err_sum = err_sum + check_stat

!-----------------------------------------------------------------------
! Allocate space for arrays from pftparm
!-----------------------------------------------------------------------
ALLOCATE( c3(npft), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( orient(npft), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( a_wl(npft), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( a_ws(npft), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( albsnc_max(npft), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( albsnc_min(npft), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( albsnf_maxu(npft), STAT = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( albsnf_max(npft), STAT = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( albsnf_maxl(npft), STAT = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( alpha(npft), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( alniru(npft), STAT = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( alnir(npft), STAT = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( alnirl(npft), STAT = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( alparu(npft), STAT = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( alpar(npft), STAT = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( alparl(npft), STAT = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( b_wl(npft), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( catch0(npft), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( dcatch_dlai(npft), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( dgl_dm(npft), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( dgl_dt(npft), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( dqcrit(npft), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( dz0v_dh(npft), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( eta_sl(npft), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( fd(npft), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( fsmc_of(npft), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( f0(npft), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( glmin(npft), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( g_leaf_0(npft), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( infil_f(npft), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( kext(npft), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( kpar(npft), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( neff(npft), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( nl0(npft), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( nr_nl(npft), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( ns_nl(npft), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( omegau(npft), STAT = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( omega(npft), STAT = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( omegal(npft), STAT = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( omniru(npft), STAT = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( omnir(npft), STAT = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( omnirl(npft), STAT = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( r_grow(npft), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( rootd_ft(npft), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( sigl(npft), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( tleaf_of(npft), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( tlow(npft), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( tupp(npft), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( emis_pft(npft), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( z0hm_pft(npft), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( z0hm_classic_pft(npft), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( dust_veg_scj(npft), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( fl_o3_ct(npft), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( dfp_dcuo(npft), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( ief(npft), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( tef(npft), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( mef(npft), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( aef(npft), stat = check_stat )
err_sum = err_sum + check_stat

!-----------------------------------------------------------------------
! Allocate space for arrays for TRIFFID
! Only need to allocate space if at least phenology is enabled
!-----------------------------------------------------------------------
IF ( l_triffid .OR. l_phenol ) THEN
  ALLOCATE( crop(npft), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( g_area(npft), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( g_grow(npft), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( g_root(npft), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( g_wood(npft), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( lai_max(npft), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( lai_min(npft), stat = check_stat )
  err_sum = err_sum + check_stat
ELSE
  ALLOCATE( crop(1), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( g_area(1), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( g_grow(1), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( g_root(1), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( g_wood(1), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( lai_max(1), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( lai_min(1), stat = check_stat )
  err_sum = err_sum + check_stat

END IF

!-----------------------------------------------------------------------
! Allocate space for arrays from snow_param
!-----------------------------------------------------------------------
ALLOCATE( cansnowtile(ntiles), stat = check_stat )
err_sum = err_sum + check_stat

!-----------------------------------------------------------------------
! Allocate space for snow scheme variables
!-----------------------------------------------------------------------
ALLOCATE( nsnow(land_field,ntiles), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( snowdepth(land_field,ntiles), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( rho_snow_grnd(land_field,ntiles), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( snowdep_surf(land_field,ntiles), stat = check_stat )
err_sum = err_sum + check_stat

!-----------------------------------------------------------------------
! Snow scheme variables that can only be allocated if nsmax > 0
!-----------------------------------------------------------------------
IF ( nsmax > 0 ) THEN
  ALLOCATE( dzsnow(nsmax), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( ds(land_field,ntiles,nsmax), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( sice(land_field,ntiles,nsmax), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( sliq(land_field,ntiles,nsmax), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( tsnow(land_field,ntiles,nsmax), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( rgrainl(land_field,ntiles,nsmax), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( rho_snow(land_field,ntiles,nsmax), stat = check_stat )
  err_sum = err_sum + check_stat

ELSE

  ALLOCATE( dzsnow(1), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( ds(1,1,1), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( sice(1,1,1), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( sliq(1,1,1), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( tsnow(1,1,1), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( rgrainl(1,1,1), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( rho_snow(1,1,1), stat = check_stat )
  err_sum = err_sum + check_stat

END IF

!-----------------------------------------------------------------------
! Allocate space for sea/sea-ice arrays from ancil_info
!-----------------------------------------------------------------------
ALLOCATE( sea_frac(t_i_length * t_j_length), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( sice_frac(t_i_length * t_j_length), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( sice_frac_ncat(t_i_length * t_j_length,nice), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( ssi_index(t_i_length * t_j_length), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( fssi(t_i_length,t_j_length), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( sea_index(t_i_length * t_j_length), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( sice_index(t_i_length * t_j_length), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( sice_pts_ncat(nice), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( sice_index_ncat(t_i_length * t_j_length,nice), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( alb_sice(t_i_length * t_j_length, nice_use, 4), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( sw_sice_rts(t_i_length * t_j_length, nice_use), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( sw_sice(t_i_length * t_j_length, nice_use), stat = check_stat )
err_sum = err_sum + check_stat

!-----------------------------------------------------------------------
! Allocate space for MORUSES arrays from urban_param
!-----------------------------------------------------------------------

IF ( l_urban2t ) THEN
  IF ( printstatus > PrStatus_Normal ) THEN
    WRITE(6,*) 'Allocating URBAN-2T / MORUSES arrays'
  END IF

  ALLOCATE( hgt(land_field), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( hwr(land_field), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( wrr(land_field), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( disp(land_field), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( ztm(land_field), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( albwl(land_field), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( albrd(land_field), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( emisw(land_field), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( emisr(land_field), stat = check_stat )
  err_sum = err_sum + check_stat
ELSE

  ALLOCATE( hgt(1), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( hwr(1), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( wrr(1), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( disp(1), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( ztm(1), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( albwl(1), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( albrd(1), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( emisw(1), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( emisr(1), stat = check_stat )
  err_sum = err_sum + check_stat

END IF

!-----------------------------------------------------------------------
! Allocate space for arrays from fluxes
!-----------------------------------------------------------------------

ALLOCATE( anthrop_heat(land_field,ntiles), stat = check_stat )
err_sum = err_sum + check_stat
IF ( check_stat == 0 ) anthrop_heat(:,:) = 0.0

ALLOCATE( surf_ht_store(land_field,ntiles), stat = check_stat )
err_sum = err_sum + check_stat

!-----------------------------------------------------------------------
! Allocate space for FLake arrays
!-----------------------------------------------------------------------

IF ( l_flake_model ) THEN
  ALLOCATE( surf_ht_flux_lake(tdims%i_start:tdims%i_end  &
                             ,tdims%j_start:tdims%j_end), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( surf_ht_flux_lk(land_field), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( sw_down(land_field), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( coriolis_param(land_field), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( u_s_lake(land_field), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( lake_depth(land_field), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( lake_fetch(land_field), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( lake_albedo(land_field), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( lake_t_snow(land_field), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( lake_t_ice(land_field), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( lake_t_mean(land_field), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( lake_t_mxl(land_field), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( lake_shape_factor(land_field), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( lake_h_snow(land_field), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( lake_h_ice(land_field), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( lake_h_mxl(land_field), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( lake_t_sfc(land_field), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( ts1_lake(land_field), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( nusselt(land_field), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( g_dt(land_field), stat = check_stat )
  err_sum = err_sum + check_stat
ELSE
  ALLOCATE( surf_ht_flux_lake(1,1), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( surf_ht_flux_lk(1), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( sw_down(1), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( coriolis_param(1), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( u_s_lake(1), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( lake_depth(1), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( lake_fetch(1), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( lake_albedo(1), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( lake_t_snow(1), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( lake_t_ice(1), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( lake_t_mean(1), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( lake_t_mxl(1), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( lake_shape_factor(1), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( lake_h_snow(1), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( lake_h_ice(1), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( lake_h_mxl(1), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( lake_t_sfc(1), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( ts1_lake(1), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( nusselt(1), stat = check_stat )
  err_sum = err_sum + check_stat

  ALLOCATE( g_dt(1), stat = check_stat )
  err_sum = err_sum + check_stat

END IF
!-----------------------------------------------------------------------  
! Allocate space for the albedo scaling to obs  
!-----------------------------------------------------------------------  
IF ( l_albedo_obs .AND. l_spec_albedo) THEN  
  ALLOCATE( albobs_scaling(land_field,ntype,2), STAT = check_stat )  
  err_sum = err_sum + check_stat  
ELSE IF ( l_albedo_obs .AND. (.NOT.l_spec_albedo) ) THEN  
  ALLOCATE( albobs_scaling(land_field,ntype,1), STAT = check_stat )  
  err_sum = err_sum + check_stat  
END IF  

!-----------------------------------------------------------------------
! Allocate space for ozone diagnostic arrays
!-----------------------------------------------------------------------
ALLOCATE( flux_o3_ft(land_field,npft), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( fo3_ft(land_field,npft), stat = check_stat )
err_sum = err_sum + check_stat

!-----------------------------------------------------------------------
! Allocate space for BVOC diagnostic vars
!-----------------------------------------------------------------------
ALLOCATE( isoprene(land_field), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( isoprene_ft(land_field,npft), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( terpene(land_field), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( terpene_ft(land_field,npft), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( methanol(land_field), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( methanol_ft(land_field,npft), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( acetone(land_field), stat = check_stat )
err_sum = err_sum + check_stat

ALLOCATE( acetone_ft(land_field,npft), stat = check_stat )
err_sum = err_sum + check_stat

!-----------------------------------------------------------------------
! Write out an error if there was one
!-----------------------------------------------------------------------
IF ( err_sum > 0 )                                                    &
  WRITE(*,*) 'An error occured while allocating JULES arrays.'

IF (lhook) CALL dr_hook('ALLOCATE_JULES_ARRAYS',zhook_out,zhook_handle)
RETURN
END SUBROUTINE allocate_jules_arrays
