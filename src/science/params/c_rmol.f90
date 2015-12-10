! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module with UM setting of 
! molar universal gas constant (8.314 J K-1 MOL-1)
  
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.
  

MODULE c_rmol

IMPLICIT NONE

! C_RMOL start
!
! Description:
! This contains the molar universal gas constant (8.314 J K-1 MOL-1)
!
!
      REAL,PARAMETER:: RMOL = 8.314 ! molar gas const

! C_RMOL end

END MODULE c_rmol
