! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module to set implicit weght coeffs

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.

MODULE c_gamma

  IMPLICIT NONE

!GAMMA holds the implicit weighting coeff. for up to 30 B.L. levels.
!It is only required for the the number of B.L. levels actually used,
!so it does not need to be set up to 30 when less BL levels are used.

  REAL gamma(30)       ! Max of 30 Boundary Layer levels assumed.

  DATA gamma / 2 * 2.0 , 1.5 , 27 * 1.0 /

END MODULE
