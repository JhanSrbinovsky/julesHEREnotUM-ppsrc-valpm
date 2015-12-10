! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! This module contains variables used for reading in triffid data
! and initialisations

  
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.


MODULE trif_io

  USE max_dimensions, ONLY:                                           &
    npft_max

  IMPLICIT NONE

!---------------------------------------------------------------------
! Set up variables to use in IO (a fixed size version of each array
! in trif that we want to initialise).
!---------------------------------------------------------------------
  INTEGER ::                                                          &
    crop_io(npft_max)

  REAL ::                                                             &
    g_area_io(npft_max),                                              &
    g_grow_io(npft_max),                                              &
    g_root_io(npft_max),                                              &
    g_wood_io(npft_max),                                              &
    lai_max_io(npft_max),                                             &
    lai_min_io(npft_max)

!---------------------------------------------------------------------
! Set up a namelist for reading and writing these arrays
!---------------------------------------------------------------------
  NAMELIST /jules_triffid/ crop_io,g_area_io,g_grow_io,g_root_io,     &
                           g_wood_io,lai_max_io,lai_min_io

END MODULE trif_io
