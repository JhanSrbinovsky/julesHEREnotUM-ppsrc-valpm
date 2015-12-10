! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: See Unified Model Code Owner's HTML page
! This file belongs in section: Land

MODULE urban_param

  IMPLICIT NONE

! Parameters for the MacDonald et al. (2008) formulation of displacement height
! and effective roughness length for momentum
  REAL, PARAMETER :: a       = 4.43
  REAL, PARAMETER :: cdz     = 1.2     ! Drag coefficient
  REAL, PARAMETER :: kappa2  = 0.16    ! Von Karman constant**2
  REAL, PARAMETER :: z0m_mat = 0.05    ! Material roughness length for momentum

! Note: z0m_mat has a value of 0.005 in original UM version, but 0.05 in
! original JULES version. z0m_mat = 0.05 was used in inter-comparison for VL92

  REAL, PARAMETER ::              &
     emiss       = 1.0,           & ! Emissivity sky
     omega_day   = 7.272e-5         ! Angular velocity of the Earth
                                    ! wrt sun( s-1 )

! At the moment set urban materials here. Could be re-written to be read from
! a look up table for different fabrics. If this was done these would need to
! be made into arrays(land_points) to store the values here and the code
! re-written to enable this.

  REAL, PARAMETER ::             &
     ! Road material = asphalt
     diffus_road  = 0.38e-6,     & ! Road: Thermal diffusivity (m2 s-1)
     cap_road     = 1.94e6,      & ! Road: Volumetric heat capacity (J K-1 m-3)

     ! Wall material = brick
     diffus_wall = 0.61e-6,      & ! Wall: Thermal diffusivity (m2 s-1)
     cap_wall    = 1.37e6,       & ! Wall: Volumetric heat capacity (J K-1 m-3)

     ! Roof material = clay
     diffus_roof = 0.47e-6,      & ! Roof: Thermal diffusivity (m2 s-1)
     cap_roof    = 1.77e6,       & ! Roof: Volumetric heat capacity (J K-1 m-3)




     dz_roof_p   = 0.02,         & ! Physical depth of roof as opposed to
                                   ! calculated damping depth

! These are required at the moment as these are not set through the STASH yet
     hgt_in    = 10.0,           & ! Building height
     hwr_in    =  1.0,           & ! Height-to-width ratio
     wrr_in    =  0.5,           & ! Width ratio
     disp_in   =  5.0,           & ! Displacement height
     ztm_in    =  1.0,           & ! Roughness length
     albwl_in  =  0.375,         & ! Albedo wall
     albrd_in  =  0.08,          & ! Albedo road
     emisw_in  =  0.875,         & ! Emissivity wall
     emisr_in  =  0.95             ! Emissivity road


  REAL, ALLOCATABLE, SAVE ::                                                  &
     hgt(:),        &  ! Building height
     hwr(:),        &  ! Height to width ratio
     wrr(:),        &  ! Width ratio
     disp(:),       &  ! Displacemnet height
     ztm(:),        &  ! Roughness length
     albwl(:),      &  ! Wall albedo
     albrd(:),      &  ! Road albedo
     emisw(:),      &  ! Wall emissivity
     emisr(:)          ! Road emissivity

  REAL ::                    &
     anthrop_heat_scale = 1.0 ! Scales anthropogenic heat source of roof &
                              ! canyon from being equally distributed (= 1.0)
                              ! to all being released in the canyon (= 0.0).
                              ! Takes a value between 0.0 - 1.0

  ! NVEGPARM URBAN-2T values. IF ( l_urban2T ) THEN set roof & canyon here
  REAL ::                    &
  ! The following are not used if appropriate MORUSES switch used
     ! ( .NOT. l_moruses_albedo; apart from albsnf_rf )
     ! Snow-covered albedos
     albsnc_c  = 0.40,       & ! Canyon
     albsnc_rf = 0.40,       & ! Roof
     ! Snow-free albedos
     albsnf_c  = 0.18,       & ! Canyon
     albsnf_rf = 0.18,       & ! Roof ( This is used regardless of switch )

     ! IF ( .NOT. l_moruses ) THEN ... ELSE ztm
     ! Roughness lengths for momentum
     ! NB These are not material roughness lengths
     z0_c      = 1.00,       & ! Canyon
     z0_rf     = 1.00,       & ! Roof
     ! Ratio of roughness length for heat to roughness length for momentum
     z0h_z0m_c  = 1.0e-7,    & ! Canyon
     z0h_z0m_rf = 1.0e-7,    & ! Roof

     ! ( .NOT. l_moruses_storage )
     ! "Canopy" heat capacity (J/K/m2)
     ch_c      = 0.28e6,     & ! Canyon
     ch_rf     = 0.04e6,     & ! Roof
     ! "Canopy" heat capacity (J/K/m2)
     vf_c      = 1.00,       & ! Canyon
     vf_rf     = 1.00,       & ! Roof

     ! ( .NOT. l_moruses_emissivity; apart from emis_rf )
     ! Emissivity
     emis_c      = 0.97,     & ! Canyon
     emis_rf     = 0.97,     & ! Roof ( This is used regardless of switch )

  ! The following are used regardless of MORUSES switch
     ! Canopy capacities (kg/m2)
     catch_c   = 0.5,        & ! Canyon
     catch_rf  = 0.5,        & ! Roof
     ! Surface conductance (m/s)
     gs_c      = 0.0,        & ! Canyon
     gs_rf     = 0.0,        & ! Roof
     ! Infiltration enhancement factors
     infil_c   = 0.1,        & ! Canyon
     infil_rf  = 0.1           ! Roof

!-----------------------------------------------------------------------
! Set up a namelist to allow URBAN-2T parameters to be set
!
! In standalone, all parameters except anthrop_heat_scale are populated
! from the non-veg parameters in init_nonveg
!-----------------------------------------------------------------------
  NAMELIST /urban2t_param/                                                   &

     albsnc_c,albsnc_rf,albsnf_c,albsnf_rf,catch_c,catch_rf,gs_c,gs_rf,      &
     infil_c,infil_rf,z0_c,z0_rf,z0h_z0m_c,z0h_z0m_rf,ch_c,ch_rf,vf_c,vf_rf, &
     emis_c,emis_rf,                                                         &

     anthrop_heat_scale

END MODULE urban_param
