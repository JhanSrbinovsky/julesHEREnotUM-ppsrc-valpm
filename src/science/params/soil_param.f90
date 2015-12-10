! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module setting soil Parameter values.
! 
  
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.

MODULE soil_param

  USE max_dimensions, ONLY : sm_levels_max

  IMPLICIT NONE

!-----------------------------------------------------------------------
!     Parameter values for soil routines.
!-----------------------------------------------------------------------
! Thermal conductivities (Source: "The Frozen Earth" p.90).
  REAL, PARAMETER ::                                                &
   hcair = 0.025                                                    &
                   ! Thermal conductivity of air (W/m/K).
  ,hcice = 2.24                                                     &
                   ! Thermal conductivity of ice (W/m/K).
  ,hcwat = 0.56    ! Thermal conductivity of liquid water (W/m/K)

!-----------------------------------------------------------------------
! Parameters for soil moisture update.
!-----------------------------------------------------------------------
  REAL, PARAMETER ::                                                &
   gamma_w = 1.0          ! Forward timestep weighting.

!-----------------------------------------------------------------------
! Parameters for soil temperature update.
!-----------------------------------------------------------------------
  INTEGER, PARAMETER ::                                             &
   mmax=3   ! Maximum number of iterations on temperature.

  REAL, PARAMETER ::                                                &
   facur = 0.01                                                     &
                     ! Required flux conservation accuracy (W/m2).
  ,gamma_t = 1.0                                                    &
                     ! Forward timestep weighting.
  ,tacur = 0.00000   ! Required accuracy of temperature calculation
!                      (Celsius).

!-----------------------------------------------------------------------
! Other variables used in soil calculations
!-----------------------------------------------------------------------
  REAL ::                                                           &
   zsmc = 1.                                                        &
           !  depth of layer over which soil moisture diagnostic
           !  is averaged (m).
  ,zst = 1.                                                         &
           !  depth of layer over which soil temperature diagnostic
           !  is averaged (m).
  ,confrac = 0.3
           !  the fraction of the gridbox over which convective
           !  precipitation is assumed to fall

  REAL, ALLOCATABLE ::                                              &
   dzsoil(:)      !  Thicknesses of the soil layers (m)

!-----------------------------------------------------------------------
! Define fixed length counterparts to the arrays we want to set in IO
! using the namelist
!-----------------------------------------------------------------------
  REAL ::                                                           &
   dzsoil_io(sm_levels_max)

  NAMELIST /jules_soil_param/ zsmc,zst,confrac,dzsoil_io

END MODULE soil_param
