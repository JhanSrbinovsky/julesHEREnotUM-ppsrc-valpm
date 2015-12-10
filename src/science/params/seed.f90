! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module sets the values of the variables FRAC_MIN and FRAC_SEED


! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.


MODULE seed

  IMPLICIT NONE

! These are set as parameters in UM include file - we just assign them
! default values.

! Minimum areal fraction for PFTs.
  REAL :: frac_min  = 1.0e-6

! "Seed" fraction for PFTs.
  REAL :: frac_seed = 0.01

  NAMELIST /jules_seed/ frac_min,frac_seed

END MODULE seed
