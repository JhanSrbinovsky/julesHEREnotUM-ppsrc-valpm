! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module for elevation values.

MODULE c_elevate

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.

! Declarations:

  USE max_dimensions, ONLY : ntiles_max

  IMPLICIT NONE








  REAL, DIMENSION(:,:), ALLOCATABLE ::                            &
   surf_hgt                   ! Height of tile above
!                               mean gridbox surface (m)
!-----------------------------------------------------------------------
! In JULES this is currently initialised as the same at every land point
! with one value for each tile
! So we create an array to use for IO which can hold a value for each
! tile
!-----------------------------------------------------------------------
  REAL ::                                                         &
   surf_hgt_io(ntiles_max)

! In example JULES control files, these are all 0, so initialise at that
  DATA surf_hgt_io / ntiles_max * 0.0 /

! Namelist used in the UM only
  NAMELIST /jules_elevate/ surf_hgt_io

END MODULE c_elevate
