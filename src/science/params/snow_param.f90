! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module with Parameters for snow parameterisation.

MODULE snow_param

  USE max_dimensions, ONLY : snow_layers_max, npft_max

  IMPLICIT NONE

! Default values taken from example JULES input file
  REAL ::                                                           &
                  ! scalars
   rho_snow_const = 250.0                                           &
                  ! constant density of lying snow (kg per m**3)
                  ! This is used:
                  ! (a) as the snow density when nsmax=0
                  ! (b) with nsmax>0 and l_snowdep_surf=.TRUE. this
                  !     is the density of snow on the canopy
                  !     (not on the ground)
  ,rho_snow_fresh = 100.0                                           &
                  ! Density of fresh snow (kg per m**3).
                  ! Only used with nsmax>0.
  ,snow_hcon = 0.265                                                &
                  ! Thermal conductivity of lying snow
                  ! (Watts per m per K).
  ,snow_hcap = 0.63e6                                               &
                  ! Thermal capacity of lying snow (J/K/m3)
  ,snowliqcap = 0.05                                                &
                  ! Liquid water holding capacity of lying snow
                  ! as a fraction of snow mass.
  ,snowinterceptfact = 0.7                                          &
                  ! constant in relationship between mass of
                  ! intercepted snow and snowfall rate
  ,snowloadlai = 4.4                                                &
                  ! ratio of maximum canopy snow load to leaf
                  ! area index (kg m-2)
  ,snowunloadfact = 0.4
                  ! constant in relationship between canopy snow
                  ! unloading and canopy snow melt rate

  REAL, ALLOCATABLE ::                                              &
   dzsnow(:)                                                        &
                   !  Prescribed thickness of snow layers (m)
!                  !  This is the thickness of each snow layer
!                  !  when it is not the bottom layer.
!                  !  (Note that dzSnow(nsMax) is not used because
!                  !   that is always the bottom layer.)
  ,ds(:,:,:)       !  Snow layer thickness (m)

  LOGICAL, ALLOCATABLE ::                                           &
                    !  ARRAYS
   cansnowtile(:)   !  Switch for canopy snow model on each tile
!                      This can be TRUE only at PFT tiles.

!-----------------------------------------------------------------------
! Define fixed length counterparts for any arrays that we want to
! initialise in IO using a namelist
!-----------------------------------------------------------------------
  REAL ::                                                           &
    dzsnow_io(snow_layers_max)

  LOGICAL ::                                                        &
    cansnowpft(npft_max)

  DATA cansnowpft / npft_max * .FALSE. /


  NAMELIST /jules_snow_param/ rho_snow_const,rho_snow_fresh,        &
                              snow_hcon,snow_hcap,                  &
                              snowliqcap,snowinterceptfact,         &
                              snowloadlai,snowunloadfact,dzsnow_io, &
                              cansnowpft

END MODULE snow_param
