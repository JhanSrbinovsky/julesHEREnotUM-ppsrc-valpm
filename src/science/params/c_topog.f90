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
 

MODULE c_topog

IMPLICIT NONE

! Use the include file if we are doing a UM run, else define variables



! C_TOPOG start
!
! Topographic index increment:
      REAL,PARAMETER :: DTI = 0.2
! Maximum topographic index considered:
      REAL,PARAMETER :: TI_MAX = 10.0
! Maximum allowed water table depth (m):
      REAL,PARAMETER :: ZW_MAX = 6.0
! Standard deviation of LOG(Ksat(0)):
      REAL,PARAMETER :: SIGMA_LOGK = 0.0
! Parameter to remove very high water tables
! from the calculated wetland fraction:
      REAL,PARAMETER :: TI_WETL = 1.5
!
! C_TOPOG end


END MODULE c_topog
