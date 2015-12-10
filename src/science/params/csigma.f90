! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module with UM setting of 
! Stefan-Boltzmann constant (W/m**2/K**4)
  
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.
    
MODULE csigma

IMPLICIT NONE

! CSIGMA start
      ! Stefan-Boltzmann constant (W/m**2/K**4).
      REAL, PARAMETER ::  SBCON=5.67E-8
! CSIGMA end

END MODULE csigma
