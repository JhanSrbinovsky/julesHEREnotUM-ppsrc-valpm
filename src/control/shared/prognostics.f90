! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module containing all of the prognostic variables 

MODULE prognostics

  IMPLICIT NONE

! Description
! Module containing all of the prognostic variables,
! i.e.those required to be kept from one timestep
! to another. Variables all appear in a model dump - NOT AT PRESENT!
! And some of these are not prognostics (eg smc)....

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.

! Declarations:


  INTEGER, ALLOCATABLE ::                                               &
    nsnow(:,:)      !  Number of snow layers on ground on tiles

  REAL, ALLOCATABLE ::                                                  &
    sice(:,:,:),                                                        &
                    ! Snow layer ice mass on tiles (Kg/m2)
    sliq(:,:,:),                                                        &
                    ! Snow layer liquid mass on tiles (Kg/m2)
    snowdepth(:,:),                                                     &
                    ! Snow depth on ground on tiles (m)
    tsnow(:,:,:),                                                       &
                    ! Snow layer temperature (K)
    rgrainl(:,:,:),                                                     &
                    ! Snow layer grain size on tiles (microns)
    rho_snow_grnd(:,:),                                                 &
                    !   Snowpack bulk density (kg/m3)
    rho_snow(:,:,:)
                    ! Snow layer densities (m)



!-----------------------------------------------------------------------
! If we are not in UM, define everything else that is needed
!-----------------------------------------------------------------------

END MODULE prognostics
