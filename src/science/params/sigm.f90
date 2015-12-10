! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module setting of Power in sigmoidal function.
! 
  
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.


MODULE sigm

  IMPLICIT NONE

! This is a parameter in the UM include file, whereas in JULES it can be
! changed. So we don't include the UM include file and just set it to a
! default value

  REAL :: pow = 20.0  ! Power in sigmoidal function.

  NAMELIST /jules_sigm/ pow

END MODULE sigm
