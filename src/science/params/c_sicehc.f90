! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module with UM setting of 
! reciprocal effective areal heat capacity of sea-ice 
  
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.
  
MODULE c_sicehc

IMPLICIT NONE

! C_SIECHC has constants for subroutine IMPL_CAL

      ! reciprocal effective areal heat capacity of sea-ice,
      !   ( 1 / (J per sq m per K)).
      REAL,PARAMETER:: AI  = 4.8E-6

! C_SIECHC end

END MODULE c_sicehc
