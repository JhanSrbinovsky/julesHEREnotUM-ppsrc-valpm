! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    SUBROUTINE HYDROL-------------------------------------------------

! Description:
!     Surface hydrology module which also updates the
!     sub-surface temperatures. Includes soil water phase
!     changes and the effect of soil water and ice on the
!     thermal and hydraulic characteristics of the soil.
!     This version is for use with MOSES II land surface
!     scheme.

! Documentation : UM Documentation Paper 25

! Subroutine Interface:
SUBROUTINE hydrol (                                               &
                   lice_pts,lice_index,soil_pts,soil_index,       &
                   nsnow,                                         &
                   npnts,nshyd,b,can_cpy,con_rain,                &
                   e_canopy,ext,hcap,hcon,ls_rain,                &
                   satcon,sathh,snowdepth,                        &
                   surf_ht_flux,timestep,                         &
                   v_sat,v_wilt,                                  &
                   can_wcnt,stf_hf_snow_melt,                     &
                   stf_sub_surf_roff,smcl,sthf,sthu,tsoil,        &
                   can_wcnt_gb,hf_snow_melt,smc,                  &
                   snow_melt,                                     &
                   sub_surf_roff,surf_roff,tot_tfall,             &
! add new inland basin variable
                   inlandout_atm,l_inland,                        &
! Additional variables for MOSES II
                   ntiles,tile_pts,tile_index,                    &
                   infil_tile,                                    &
                   melt_tile,tile_frac,                           &
! Additional variables required for large-scale hydrology:
                   l_top,l_pdm,fexp,gamtot,ti_mean,ti_sig,cs,     &
                   dun_roff,drain,fsat,fwetl,qbase,qbase_zw,      &
                   zw,sthzw,a_fsat,c_fsat,a_fwet,c_fwet,          &
                   fch4_wetl,dim_cs1,l_soil_sat_down,l_triffid,   &
 ltimer                                                           &
 )

USE c_topog, ONLY :                                               &
!      imported scalar parameters
 ti_max                                                           &
           !  Maximum topographic index considered
,zw_max  !  Maximum allowed water table depth (m)

USE soil_param, ONLY :                                            &
 dzsoil       !  Thicknesses of the soil layers (m)

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Subroutine arguments
!   Scalar arguments with intent(IN) :
INTEGER, INTENT(IN) ::                                            &
 lice_pts                                                         &
                     ! IN Number of land ice points.
,npnts                                                            &
                     ! IN Number of gridpoints.
,nshyd                                                            &
                     ! IN Number of soil moisture levels.
,soil_pts                                                         &
                     ! IN Number of soil points.
,ntiles                                                           &
                     ! IN Number of tiles
,dim_cs1             ! IN Number of soil carbon pools

REAL, INTENT(IN) ::                                               &
 timestep            ! IN Model timestep (s).

LOGICAL, INTENT(IN) :: ltimer ! Logical switch for TIMER diags

LOGICAL , INTENT(IN) ::                                           &
 stf_hf_snow_melt                                                 &
                     ! IN Stash flag for snowmelt heat flux.
!cxyz STF_HF_SNOW_MELT is not used in this version.
,stf_sub_surf_roff                                                &
                     ! IN Stash flag for sub-surface runoff.
,l_top                                                            &
             ! IN Flag for TOPMODEL-based hydrology.
,l_pdm                                                            &
             ! IN Flag for PDM hydrology.
,l_soil_sat_down                                                  &
             ! IN Switch controlling direction of movement of
             !    soil moisture in excess of saturation
,l_triffid   ! IN Switch to use TRIFFID.


!   Array arguments with intent(IN) :
INTEGER, INTENT(IN) ::                                            &
 lice_index(npnts)                                                &
                     ! IN Array of land ice points.
,soil_index(npnts)                                                &
                     ! IN Array of soil points.
,nsnow(npnts,ntiles) ! IN Number of snow layers

REAL, INTENT(IN) ::                                               &
 b(npnts,nshyd)                                                   &
                     ! IN Clapp-Hornberger exponent.
,can_cpy(npnts,ntiles)                                            &
                      !IN Canopy/surface capacity of
!                          !    land tiles (kg/m2).
,con_rain(npnts)                                                  &
                     ! IN Convective rain (kg/m2/s).
,e_canopy(npnts,ntiles)                                           &
!                          ! IN Canopy evaporation from
!                          !    land tiles (kg/m2/s).
,ext(npnts,nshyd)                                                 &
                     ! IN Extraction of water from each soil
!                          !    layer (kg/m2/s).
,hcap(npnts,nshyd)                                                &
                     ! IN Soil heat capacity (J/K/m3).
,hcon(npnts,0:nshyd)                                              &
                     ! IN Soil thermal conductivity (W/m/K).
,ls_rain(npnts)                                                   &
                     ! IN Large-scale rain (kg/m2/s).
,satcon(npnts,0:nshyd)                                            &
                     ! IN Saturated hydraulic conductivity
!                          !    (kg/m2/s).
,sathh(npnts,nshyd)                                               &
                     ! IN Saturated soil water pressure (m).
,snow_melt(npnts)                                                 &
                     ! IN Snowmelt (kg/m2/s).
,snowdepth(npnts,ntiles)                                          &
                     ! Snow depth (on ground) (m)
,surf_ht_flux(npnts)                                              &
                     ! IN Net downward surface heat flux (W/m2)
,v_sat(npnts,nshyd)                                               &
                     ! IN Volumetric soil moisture
!                          !    concentration at saturation
!                          !    (m3 H2O/m3 soil).
,v_wilt(npnts,nshyd)                                              &
                     ! IN Volumetric soil moisture
!                          !    concentration below which
!                          !    stomata close (m3 H2O/m3 soil).
,fexp(npnts)                                                      &
                     ! IN Decay factor in Sat. Conductivity
!                          !    in water table layer.
,gamtot(npnts)                                                    &
                     ! IN Integrated complete Gamma function.
,ti_mean(npnts)                                                   &
                     ! IN Mean topographic index.
,ti_sig(npnts)                                                    &
                     ! IN Standard dev. of topographic index.
,cs(npnts,dim_cs1)                                                &
                     ! IN Soil carbon (kg C/m2).
!                          !   For RothC (dim_cs1=4), the pools
!                          !    are DPM, RPM, biomass and humus.
,a_fsat(npnts)                                                    &
                     ! IN Fitting parameter for Fsat in LSH model
,c_fsat(npnts)                                                    &
                     ! IN Fitting parameter for Fsat in LSH model
,a_fwet(npnts)                                                    &
                     ! IN Fitting parameter for Fwet in LSH model
,c_fwet(npnts)
                     ! IN Fitting parameter for Fwet in LSH model


!   Array arguments with intent(INOUT) :

REAL, INTENT(INOUT) ::                                            &
 can_wcnt(npnts,ntiles)                                           &
!                          ! INOUT Canopy water content for
!                          !       land tiles (kg/m2).
,smcl(npnts,nshyd)                                                &
                     ! INOUT Soil moisture content of each
!                          !       layer (kg/m2).
,sthf(npnts,nshyd)                                                &
                     ! INOUT Frozen soil moisture content of
!                          !       each layer as a fraction of
!                          !       saturation.
,sthu(npnts,nshyd)                                                &
                     ! INOUT Unfrozen soil moisture content of
!                          !       each layer as a fraction of
!                          !       saturation.
,tsoil(npnts,nshyd)                                               &
                     ! INOUT Sub-surface temperatures (K).
,fsat(npnts)                                                      &
                      ! INOUT Surface saturation fraction.
,fwetl(npnts)                                                     &
                      ! INOUT Wetland fraction.
,zw(npnts)                                                        &
                      ! INOUT Water table depth (m).
,sthzw(npnts)         ! INOUT soil moist fract. in deep-zw layer.

! Arguments which are required for compatibility but not used.
! Shouldnt use OUT.
REAL              ::                                              &
 hf_snow_melt(npnts)
                      ! OUT Gridbox snowmelt heat flux (W/m2).


!   Array arguments with intent(OUT) :
REAL, INTENT(OUT) ::                                              &
 can_wcnt_gb(npnts)                                               &
                      ! OUT Gridbox canopy water content (kg/m2).
,smc(npnts)                                                       &
                      ! OUT Available soil moisture in a layer at
                      !     the surface (kg/m2)
,sub_surf_roff(npnts)                                             &
                      ! OUT Sub-surface runoff (kg/m2/s).
,surf_roff(npnts)                                                 &
                      ! OUT Surface runoff (kg/m2/s).
,tot_tfall(npnts)                                                 &
                      ! OUT Total throughfall (kg/m2/s).
,dun_roff(npnts)                                                  &
                      ! OUT Dunne part of sfc runoff (kg/m2/s).
,qbase(npnts)                                                     &
                      ! OUT Base flow (kg/m2/s).
,qbase_zw(npnts)                                                  &
                      ! OUT Base flow from ZW layer (kg/m2/s).
,drain(npnts)                                                     &
                      ! OUT Drainage out of nshyd'th level (kg/m2/s).
,fch4_wetl(npnts)     ! OUT Scaled wetland methane flux.
                      !     (10^-9 kg C/m2/s).


! Additional variables for MOSES II
INTEGER, INTENT(IN) ::                                            &
 tile_pts(ntiles)                                                 &
                     ! IN Number of tile points.
,tile_index(npnts,ntiles)
!                          ! IN Index of tile points.

REAL, INTENT(IN) ::                                               &
 infil_tile(npnts,ntiles)                                         &
!                          ! IN Maximum surface infiltration
,melt_tile(npnts,ntiles)                                          &
!                          ! IN Snowmelt on tiles (kg/m2/s).
,tile_frac(npnts,ntiles)                                          &
                     ! IN Tile fractions.

! Declare variable for inland basin outflow
,inlandout_atm(npnts)            ! IN TRIP INLAND BASIN
!                       OUTFLOW FOR LAND POINTS ONLY,kg/m2/s=mm
LOGICAL, INTENT(IN) ::                                            &
 l_inland                   ! IN True if re-routing inland
                            !   basin flow to soil moisture

! Local scalars:
INTEGER                                                           &
 i,j                                                              &
                      ! WORK Loop counters.
,n                    ! WORK Tile loop counter.

! Local arrays:

REAL                                                              &
 dsmc_dt(npnts)                                                   &
                      ! WORK Rate of change of soil moisture
!                           !      due to water falling onto the
!                           !      surface after surface runoff
!                           !      (kg/m2/s).
,w_flux(npnts,0:nshyd)                                            &
                      ! WORK Fluxes of water between layers
!                           !      (kg/m2/s).
,ksz(npnts,0:nshyd)                                               &
                      ! WORK Saturated hydraulic
!                           !      conductivity in layer (kg/m2/s).
,qbase_l(npnts,nshyd+1)                                           &
!                           ! WORK Base flow from each level (kg/m2/s).
,top_crit(npnts)                                                  &
                      ! WORK Critical TI when ZW <=0.0
,zdepth(0:nshyd)                                                  &
                      ! WORK Lower soil layer boundary depth (m).
,tsoil_d(npnts)                                                   &
                      ! WORK Soil temperature in the top metre
,wutot(npnts)         ! WORK Ratio of unfrozen to total soil
!                                             !    moisture at ZW.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! End of header--------------------------------------------------------

IF (lhook) CALL dr_hook('HYDROL',zhook_in,zhook_handle)

!----------------------------------------------------------------------
! Set up variables required for LSH scheme:
!----------------------------------------------------------------------
zdepth(:)=0.0

DO n=1,nshyd
   zdepth(n)=zdepth(n-1)+dzsoil(n)
END DO
!-----------------------------------------------------------------------
! Calculate throughfall and surface runoff, and update the canopy water
! content
!-----------------------------------------------------------------------
!! DEPENDS ON: surf_hyd
!CALL surf_hyd (npnts,ntiles,tile_pts,tile_index,                  &
!               can_cpy,e_canopy,tile_frac,infil_tile,con_rain,    &
!               ls_rain,melt_tile,snow_melt,timestep,              &
!               can_wcnt,can_wcnt_gb,dsmc_dt,                      &
!               l_top,l_pdm,nshyd,soil_pts,soil_index,             &
!               surf_roff,tot_tfall,                               &
!               dun_roff,fsat,v_sat,sthu,sthf)

!-----------------------------------------------------------------------
! Specify the reduction of hydraulic conductivity with depth:
! Initialiase base flow to zero:
!-----------------------------------------------------------------------

DO n=0,nshyd
!CDIR NODEP
  DO j=1,soil_pts
    i=soil_index(j)
    ksz(i,n)=satcon(i,n)
  END DO
END DO

DO n=1,nshyd
!CDIR NODEP
  DO j=1,soil_pts
    qbase_l(soil_index(j),n)=0.0
  END DO
END DO

DO i=1,npnts
  qbase(i)=0.0
  qbase_zw(i)=0.0
  wutot(i)=0.0
  drain(i)=0.0
END DO

!IF(l_top)THEN
!  IF (soil_pts /= 0) THEN
!! DEPENDS ON: calc_baseflow_jules
!    CALL calc_baseflow_jules(                                     &
!        soil_pts,soil_index,npnts,nshyd                           &
!       ,zdepth,ksz                                                &
!       ,b,fexp,ti_mean,zw,sthf,sthu                               &
!       ,wutot,top_crit,qbase,qbase_l                              &
! )
!  END IF
!END IF


IF(l_inland)THEN

 DO i=1,npnts

! Add inland basin outflow to change in soil moisture store

   dsmc_dt(i)=dsmc_dt(i)+inlandout_atm(i)

 END DO
END IF


!-----------------------------------------------------------------------
! Update the layer soil moisture contents and calculate the
! gravitational drainage.
!-----------------------------------------------------------------------
!IF (soil_pts /= 0) THEN
!! DEPENDS ON: soil_hyd
!  CALL soil_hyd (npnts,nshyd,soil_pts,soil_index,b,dzsoil,        &
!               ext,dsmc_dt,satcon,ksz,sathh,timestep,v_sat,       &
!               sub_surf_roff,smcl,sthu,surf_roff,w_flux,          &
!               stf_sub_surf_roff,                                 &
!               zw,sthzw,zdepth,qbase,qbase_l,                     &
!               dun_roff,drain,l_top,l_soil_sat_down)
!!-----------------------------------------------------------------------
!! Calculate surface saturation and wetland fractions:
!!-----------------------------------------------------------------------
!  IF(l_top)THEN
!
!    DO i=1,npnts
!      fsat(i)=0.0
!      fwetl(i)=0.0
!! Zero soil porosity over land ice:
!      IF(v_sat(i,nshyd) <= 0.0) zw(i)=zw_max
!    END DO
!
!    DO j=1,soil_pts
!      i=soil_index(j)
!      qbase_zw(i)=qbase_l(i,nshyd+1)
!!Now use fit for fsat and fwet:
!      fsat(i)=a_fsat(i)*EXP(-c_fsat(i)*top_crit(i))
!      fwetl(i)=a_fwet(i)*EXP(-c_fwet(i)*top_crit(i))
!      IF(top_crit(i) >= ti_max)THEN
!        fsat(i)=0.0
!        fwetl(i)=0.0
!      END IF
!    END DO
!  END IF  !  l_top
!
!ELSE
!
!!---------------------------------------------------------------------
!! If required by STASH flag and there are no soil points,
!! set sub-surface runoff to zero.
!!---------------------------------------------------------------------
!
!  IF(stf_sub_surf_roff) THEN
!    DO i=1,npnts
!      sub_surf_roff(i)=0.0
!    END DO
!  END IF
!
!END IF  !  soil_pts

!-----------------------------------------------------------------------
! Update the soil temperatures and the frozen moisture fractions
!-----------------------------------------------------------------------
!IF (soil_pts /= 0) THEN
!! DEPENDS ON: soil_htc
!  CALL soil_htc (npnts,nshyd,ntiles,soil_pts,soil_index,          &
!                 tile_pts,tile_index,nsnow,                       &
!                 b,dzsoil,tile_frac,hcap,hcon,                    &
!                 sathh,surf_ht_flux,timestep,v_sat,               &
!                 w_flux,smcl,snowdepth,sthu,sthf,tsoil)
!END IF
!
!!-----------------------------------------------------------------------
!! Update the sub-surface temperatures for land ice
!!-----------------------------------------------------------------------
!IF (lice_pts /= 0) THEN
!! DEPENDS ON: ice_htc
!  CALL ice_htc (npnts,nshyd,lice_pts,lice_index,dzsoil,           &
!                surf_ht_flux,timestep,                            &
!                tsoil)
!END IF

!-----------------------------------------------------------------------
! Diagnose the available soil moisture in a layer at the surface.
!-----------------------------------------------------------------------
! DEPENDS ON: soilmc
CALL soilmc ( npnts,nshyd,soil_pts,soil_index,                    &
              dzsoil,sthu,v_sat,v_wilt,smc )

!-----------------------------------------------------------------------
! Calculate mean soil temperature and scaled CH4 flux:
!-----------------------------------------------------------------------

DO i=1,npnts
  fch4_wetl(i)=0.0
END DO
!IF(l_top)THEN
!  IF (soil_pts /= 0) THEN
!! DEPENDS ON: soilt
!    CALL soilt(npnts,nshyd,soil_pts,soil_index                    &
!             ,dzsoil,tsoil,tsoil_d)
!! DEPENDS ON: ch4_wetl
!    CALL ch4_wetl(npnts,soil_pts,dim_cs1,soil_index,l_triffid     &
!      ,tsoil_d,cs,fwetl,fch4_wetl)
!  END IF
!END IF

!-----------------------------------------------------------------------
IF (lhook) CALL dr_hook('HYDROL',zhook_out,zhook_handle)
RETURN
END SUBROUTINE hydrol
