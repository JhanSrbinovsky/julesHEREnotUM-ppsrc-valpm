MODULE switches_urban

! Description:
!   Module containing switches for the parametrisations of urban scheme MORUSES
!
! (c) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! -----------------------------------------------------------------------------

  IMPLICIT NONE

  LOGICAL ::                           &
! l_urban2T and l_moruses, when initially set, are mutually exclusive as they
! are used to set the other switches to their default values in init_urban.
!
     l_urban_empirical      = .TRUE.,  & ! Empirical relationships for urban
                                         ! geometry (WRR & HWR)
     l_moruses_macdonald    = .TRUE.,  & ! MacDonald formulation for
                                         ! displacement height and effective
                                         ! roughness length for momentum
     l_urban2T              = .FALSE., & ! Original URBAN-2T switch
     l_moruses              = .FALSE., & ! Indicates any MORUSES switch used

! Independent parametristaion switches
     l_moruses_albedo       = .TRUE.,  & ! SW canyon albedo
     l_moruses_emissivity   = .FALSE., & ! LW canyon emissivity
     l_moruses_rough        = .TRUE.,  & ! Heat transfer
     l_moruses_storage      = .TRUE.,  & ! Storage
     l_moruses_storage_thin = .TRUE.     ! Storage thin roof

!-----------------------------------------------------------------------
! Set up a namelist to allow switches to be set
! In the UM, l_urban2T is in the namelist, l_moruses is set by inspecting
!   other MORUSES switches (i.e. if l_moruses_albedo = T then l_moruses = T
!
! In standalone JULES, l_moruses is in the namelist, and l_urban2T is
!   set by looking at whether 2 urban tiles have been specified
!
! Either way, after initialisation, we have l_moruses = T => l_urban2T = T,
!   but not necessarily the other way round
!-----------------------------------------------------------------------
  NAMELIST /urban_switches/                                                   &

     l_urban2T,                                                               &



     l_moruses_albedo,l_moruses_emissivity,l_moruses_rough,l_moruses_storage, &
     l_moruses_storage_thin,l_moruses_macdonald,l_urban_empirical

END MODULE switches_urban
