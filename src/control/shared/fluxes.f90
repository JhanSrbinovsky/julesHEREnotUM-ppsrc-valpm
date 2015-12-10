! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module containing surface fluxes.
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.
!
! Code Owner: See Unified Model Code Owner's HTML page
! This file belongs in section: Land

MODULE fluxes

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! anthrop_heat is required by both the UM and standalone configurations
!-----------------------------------------------------------------------------
  REAL, ALLOCATABLE ::                                                        &
    anthrop_heat(:,:)      ! Additional heat source on tiles
!                          ! used for anthropgenic urban
!                          ! heat source (W/m2)

  REAL, ALLOCATABLE ::                                                        &
    surf_ht_store(:,:)     !   Diagnostic to store values
!                          !   of C*(dT/dt) during calculation
!                          !   of energy balance

  REAL, ALLOCATABLE, SAVE ::                                                  &
    sw_sice(:,:),                                                             &
                           !  Net SW on sea ice categories
    sw_sice_rts(:,:),                                                         &
                           !  Net SW on sea ice categories on the
!                          !     radiative timestep
    alb_sice(:,:,:)        !  Albedo of sea ice categories

!-----------------------------------------------------------------------------
! Everything else is required only in standalone configuration
!-----------------------------------------------------------------------------

END MODULE fluxes
