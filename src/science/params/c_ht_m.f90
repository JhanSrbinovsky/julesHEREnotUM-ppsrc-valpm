! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!

! Description:
!   Heights for 10m and 1.5m diagnotics.

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.

MODULE c_ht_m

IMPLICIT NONE

! C_HT_M constants for subroutine SF_EXCH

      ! height of 10m level for diagnostic calculations (m).
      REAL,PARAMETER:: Z10M  = 10.0

      ! height of 1.5m level for diagnostic calculations (m).
!CABLE? REAL,PARAMETER:: Z1P5M = 1.5
      REAL,PARAMETER:: Z1P5M = 2.0

! C_HT_M end

END MODULE
