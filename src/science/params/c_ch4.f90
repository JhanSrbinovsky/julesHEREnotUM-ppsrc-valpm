! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module with UM setting of 
! params used in methane flux from wetlands.

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.

MODULE c_ch4

IMPLICIT NONE

! C_CH4 start
! 5.5 17/02/03    Required for large-scale hydrology L_TOP code.
!                 Used in calculation of methane flux from wetlands.
!
      REAL,PARAMETER :: T0_CH4 = 273.15               ! T0 value
      REAL,PARAMETER :: Q10_CH4 = 3.7                 ! Q10 value
! Scale factor when TRIFFID (and therefore ROTHC) is ON:
      REAL,PARAMETER :: CONST_CH4_trifid = 5.41E-10
! Scale factor with fixed vegetation:
      REAL,PARAMETER :: CONST_CH4_fixveg = 5.41E-12
! Scale factor with fixed vegetation:
      REAL,PARAMETER :: CONST_CH4        = 5.41E-12
!
! C_CH4 end

END MODULE c_ch4
