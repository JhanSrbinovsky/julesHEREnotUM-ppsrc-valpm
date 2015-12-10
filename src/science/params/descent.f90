! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module with UM setting of 
! 
  
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.
    
MODULE descent

IMPLICIT NONE

! DESCENT start

! Number of TRIFFID iterations for gradient descent to equilibrium.
      INTEGER,PARAMETER:: ITER_EQ = 10

! Minimum value for the denominator of the update equation. Ensures
! that gradient descent does not lead to an unstable solution.
      REAL,PARAMETER:: DENOM_MIN=1.0E-6

! Inverse timestep for gradient  descent to equilibrium (/360days).
      REAL,PARAMETER:: GAMMA_EQ = 1.0E-1

! DESCENT end

END MODULE descent
