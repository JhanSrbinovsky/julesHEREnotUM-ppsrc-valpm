! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description: Module containing the variables for the FLake lake scheme
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.
!
!   Code Owner: See Unified Model Code Owner's HTML page
!   This file belongs in section: Land
!
      MODULE LAKE_MOD
      IMPLICIT NONE

! parameters following values in FLake routine flake_parameters
      REAL, PARAMETER  :: h_snow_min_flk = 1.0E-5  & ! (m)
                         ,h_ice_min_flk  = 1.0E-9  & ! (m)
                         ,h_ice_max      = 3.0       ! (m)

! parameters following values in FLake routine flake_albedo_ref
      REAL, PARAMETER  :: albedo_whiteice_ref =  0.60 & ! White ice
                         ,albedo_blueice_ref  =  0.10 & ! Blue ice
                         ,c_albice_MR         = 95.6    ! Constant in the interpolation formula for 
                                                        ! the ice albedo (Mironov and Ritter 2004)

      REAL, DIMENSION(:,:), ALLOCATABLE :: SURF_HT_FLUX_LAKE
!                                        ! Net downward heat flux at surface over
!                                        ! lake fraction of gridbox, all points (W/m2)

      REAL, DIMENSION(:), ALLOCATABLE :: SURF_HT_FLUX_LK
!                                        ! Net downward heat flux at surface over
!                                        ! lake fraction of gridbox, land points (W/m2)
      REAL, DIMENSION(:), ALLOCATABLE :: U_S_LAKE
!                                        ! lake subsurface friction velocity (m/s)
      REAL, DIMENSION(:), ALLOCATABLE :: SW_DOWN
!                                        ! downwelling shortwave irradiance [W m^{-2}]
      REAL, DIMENSION(:), ALLOCATABLE :: CORIOLIS_PARAM
!                                        ! Coriolis parameter (s^-1)
      REAL, DIMENSION(:), ALLOCATABLE :: LAKE_DEPTH
!                                        ! lake depth (m)
      REAL, DIMENSION(:), ALLOCATABLE :: LAKE_FETCH
!                                        ! Typical wind fetch (m)
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: LAKE_ALBEDO
!                                        ! lake albedo
      REAL, DIMENSION(:), ALLOCATABLE :: LAKE_T_SNOW
!                                        ! temperature at the air-snow interface (K)
      REAL, DIMENSION(:), ALLOCATABLE :: LAKE_T_ICE
!                                        ! temperature at upper boundary of lake ice (K)
      REAL, DIMENSION(:), ALLOCATABLE :: LAKE_T_MEAN
!                                        ! lake mean temperature (K)
      REAL, DIMENSION(:), ALLOCATABLE :: LAKE_T_MXL
!                                        ! lake mixed-layer temperature (K)
      REAL, DIMENSION(:), ALLOCATABLE :: LAKE_SHAPE_FACTOR
!                                        ! thermocline shape factor (dimensionless?)
      REAL, DIMENSION(:), ALLOCATABLE :: LAKE_H_SNOW
!                                        ! snow thickness (m)
      REAL, DIMENSION(:), ALLOCATABLE :: LAKE_H_ICE
!                                        ! lake ice thickness (m)
      REAL, DIMENSION(:), ALLOCATABLE :: LAKE_H_MXL
!                                        ! lake mixed-layer thickness (m)
      REAL, DIMENSION(:), ALLOCATABLE :: LAKE_T_SFC
!                                        ! temperature (of water, ice or snow) at surface of lake (K)
      REAL, DIMENSION(:), ALLOCATABLE :: TS1_LAKE
!                                        ! "average" temperature of
!                                        ! lake-ice, lake and soil sandwich
!                                        ! (K).
      REAL, DIMENSION(:), ALLOCATABLE :: NUSSELT
!                                        ! Nusselt number
      REAL, DIMENSION(:), ALLOCATABLE :: G_DT
!                                        ! ground heat flux over delta T [W m-2 K-1]

! initial values
      REAL :: lake_depth_0  =    5.0  ! (m)
      REAL :: lake_fetch_0  =   25.0  ! (m)
      REAL :: lake_h_mxl_0  =    2.0  ! (m)
      REAL :: lake_shape_0  =    0.5  ! shape factor, dimensionless
      REAL :: nusselt_0     = 1000.0  ! Nusselt number, dimensionless
      REAL :: g_dt_0        = 1.0e-10 ! very small [W m-2 K-1]

! trap counters
      INTEGER :: trap_frozen   = 0
      INTEGER :: trap_unfrozen = 0

      END MODULE
