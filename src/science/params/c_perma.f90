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
              
MODULE c_perma

IMPLICIT NONE

! C_PERMA start

      ! Specific heat capacity of water vapour (J/kg/K)
      REAL,PARAMETER:: HCAPV=1850.0

      ! Specific heat capacity of water (J/kg/K)
      REAL,PARAMETER:: HCAPW=4180.0

      ! Specific heat capacity of ice (J/kg/K)
      REAL,PARAMETER:: HCAPI=2100.0

      ! Density of ice (kg/m3)
      REAL,PARAMETER:: RHO_ICE=917

      ! Rate of change of ice potential with temperature
      ! RHO_ICE*LF/ZERODEGC*1/(RHO_WATER*G) (m/K)
      REAL,PARAMETER:: DPSIDT=114.3

! C_PERMA end

END MODULE c_perma
