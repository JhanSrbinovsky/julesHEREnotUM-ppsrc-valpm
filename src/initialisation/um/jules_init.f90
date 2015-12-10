
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine JULES_INIT ----------------------------------
!
! Description: Initialisation of JULES.
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.
!
!   Code Owner: See Unified Model Code Owner's HTML page
!   This file belongs in section: Land


SUBROUTINE jules_init(  &
! ARGLNDM Constants for physics routines
        land_index, land_ice_index, soil_index,                         &
! ARGLNDM end
                        land_field,ntiles,sm_levels                   &
                       ,nice,nice_use,frac_land                       &
                       ,l_snow_albedo                                 &
                       ,rgrain,snow_tile,deep_soil_temp,snow_grnd     &
                       ,nsnow,rgrainl,rho_snow_grnd,sice,sliq         &
                       ,snowdepth,ds,tsnow )

USE nstypes, ONLY : ntype, lake

! JULES variables to be initialised
USE jules_mod, ONLY :  clapp_levs                                     &
                     , sathh_levs                                     &
                     ,  hcap_levs                                     &
                     ,  hcon_levs                                     &
                     ,satcon_levs                                     &
                     ,smvccl_levs                                     &
                     ,smvcwt_levs                                     &
                     ,smvcst_levs

! UM variables from which to initialise the JULES variables
USE atm_fields_mod, ONLY : clapp_horn                                 &
                          ,sat_soilw_suction                          &
                          ,sat_soil_cond                              &
                          ,therm_cap                                  &
                          ,therm_cond                                 &
                          ,vol_smc_crit                               &
                          ,vol_smc_wilt                               &
                          ,vol_smc_sat

! JULES module snow depth
USE prognostics, ONLY :                                               &
        snowdepth_jules => snowdepth

! JULES switches
      USE switches, ONLY : l_um_jules                                 &
                          ,l_aggregate                                &
                          ,l_flake_model

! JULES max no. of snow layers
  USE ancil_info, ONLY :  &
     nsmax

  USE dyn_coriolis_mod,  ONLY : f3_at_u
  USE theta_field_sizes, ONLY : t_i_length

  USE lake_mod, ONLY :    &
     coriolis_param       &
    ,nusselt              &
    ,nusselt_0

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Subroutine arguments

!   Scalar arguments with intent(in):
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

  REAL ::  frac_land( land_field,ntype)       &
!
          ,rgrain(    land_field,ntiles)      &
          ,snow_tile( land_field,ntiles)      &
          ,snow_grnd( land_field,ntiles)      &
!
          ,deep_soil_temp(land_field,sm_levels)  &
!
          ,snowdepth(    land_field,ntiles)          &
          ,nsnow(        land_field,ntiles)          &
          ,rho_snow_grnd(land_field,ntiles)          &
!
          ,tsnow(        land_field,ntiles,nsmax)    &
          ,rgrainl(      land_field,ntiles,nsmax)    &
          ,sice(         land_field,ntiles,nsmax)    &
          ,sliq(         land_field,ntiles,nsmax)    &
          ,ds(           land_field,ntiles,nsmax)    &
          ,rho_snow(     land_field,ntiles,nsmax)

! land_index
! TYPLNDM
! Formerly integral part of TYPCONA, the variables below have been
! separated from the rest of TYPCONA as they are required by some
! of the Ocean routines in the Ocean-Atmosphere configuration of
! the UM whilest TYPCONA is not.

      ! Primary Arrays
      INTEGER::land_points     ! No. of land points  (can be 0)
      INTEGER::land_ice_points ! Number of land ice points
      INTEGER::soil_points     ! Number of soil points

      ! Do not allow these arrays to have zero size
      INTEGER::land_index    (max(1,land_field)) ! set from land_sea_mask
      INTEGER::land_ice_index(max(1,land_field)) ! Array of land ice points.
      INTEGER::soil_index    (max(1,land_field)) ! Array of soil points.

      ! Gets some sizes transported around the model :
      COMMON /land_soil_dimensions/                                     &
     &  land_points , land_ice_points , soil_points

! TYPLNDM end

 LOGICAL :: l_snow_albedo

! WORK variables:
INTEGER :: i,j,k,l,n

  INTEGER ::                                &
    nsnow_integer(land_field,ntiles)


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!
! START OF EXECUTABLE CODE
!

! Dimension the JULES fields.
IF (lhook) CALL dr_hook('JULES_INIT',zhook_in,zhook_handle)
ALLOCATE( clapp_levs(land_field,  sm_levels))
ALLOCATE( sathh_levs(land_field,  sm_levels))
ALLOCATE(  hcap_levs(land_field,  sm_levels))
ALLOCATE(  hcon_levs(land_field,0:sm_levels))
ALLOCATE(satcon_levs(land_field,0:sm_levels))
ALLOCATE(smvccl_levs(land_field,  sm_levels))
ALLOCATE(smvcwt_levs(land_field,  sm_levels))
ALLOCATE(smvcst_levs(land_field,  sm_levels))

! Initialise JULES arrays.
! Soil properties on soil moisture levels.
DO i = 1, land_field
  DO j = 1, sm_levels
      clapp_levs(i,j) = clapp_horn(       i)
      sathh_levs(i,j) = sat_soilw_suction(i)
       hcap_levs(i,j) = therm_cap(        i)
     smvccl_levs(i,j) = vol_smc_crit(     i)
     smvcwt_levs(i,j) = vol_smc_wilt(     i)
     smvcst_levs(i,j) = vol_smc_sat(      i)
   END DO
   DO j = 0, sm_levels
       hcon_levs(i,j) = therm_cond(   i)
     satcon_levs(i,j) = sat_soil_cond(i)
  END DO
END DO

IF (l_um_jules) THEN

! snowdepth needed in AP1 for JULES radiation
!---------------------------------------------
  snowdepth_jules = snowdepth

END IF

! FLake model
!--------------
IF (     l_flake_model                       &
    .AND.(.NOT.l_aggregate)) THEN

! initialise the Nusselt number
    nusselt(:) = nusselt_0

    DO l=1,land_field

        j=(land_index(l)-1)/t_i_length + 1
        i = land_index(l) - (j-1)*t_i_length

! set the Coriolis parameter : ABSOLUTE VALUE
!
! To get the value at theta points,
! average the adjacent values at u points.
!
        coriolis_param(l) = ABS( (f3_at_u(i,j)+f3_at_u(i-1,j))/2.0 )

    END DO

END IF

IF (lhook) CALL dr_hook('JULES_INIT',zhook_out,zhook_handle)
RETURN

END SUBROUTINE jules_init
