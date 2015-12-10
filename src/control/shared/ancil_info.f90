! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module containing size and dimension parameters, indexing variables
! ...and more.
!
! Most of these are not required in the UM implementation
!

MODULE ancil_info

IMPLICIT NONE

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.

! Declarations:

!-----------------------------------------------------------------------
! Define variables that are needed by both standalone and UM
!-----------------------------------------------------------------------
  INTEGER ::                                                            &
    nsmax = 0                                                           &
                      !  Maximum number of snow layers
   ,ssi_pts                                                             &
                      !  Number of sea or sea-ice points
   ,sea_pts                                                             &
                      !  Number of sea points
   ,sice_pts
                      !  Number of sea-ice points

  INTEGER, ALLOCATABLE ::                                               &
    ssi_index(:)                                                        &
                      !  index of sea and sea-ice points
   ,sea_index(:)                                                        &
                      !  index of sea points
   ,sice_index(:)                                                       &
                      !  index of sea-ice points
   ,sice_pts_ncat(:)                                                    &
                      !  Number of points for each sea-ice category
   ,sice_index_ncat(:,:)  !  index of points for each sea-ice category


  REAL, ALLOCATABLE ::                                                  &
    fssi(:,:)                                                           &
                      !  Fraction of gridbox covered by sea
!                     !  or sea-ice
   ,sea_frac(:)                                                         &
                      !  Fraction of gridbox covered by sea
!                     !  (converted to single vector array)
   ,sice_frac(:)                                                        &
                      !  Fraction of gridbox covered by sea-ice
!                     !  (converted to single vector array)
   ,sice_frac_ncat(:,:)     !  Fraction of gridbox covered by each
!                           !  sea-ice category
!                           !  (converted to single vector array)


!-----------------------------------------------------------------------
! If we are not in UM, define everything else that is needed
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Namelist for UM only - sm_levels and nsmax are read in the namelist
! jules_model_levels in standalone
!-----------------------------------------------------------------------
  NAMELIST /jules_ancil_info/ nsmax

END MODULE ancil_info
