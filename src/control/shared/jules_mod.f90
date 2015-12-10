! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************
!
!+ module of declarations and initialisation routines for JULES/UM

MODULE jules_mod

  ! Description:
  !   Module containing declarations of JULES-related variables
  !   and initialisation subroutines.
  !
  ! Method:
  !   Module contains declarations of variables needed for
  !   implementation of JULES routines in the UM, and also
  !   subroutines required to initialise these variables.
  !
  !   Since these variables will be in the UM whether JULES routines
  !   are implemented or not, the idea is to collect them here
  !   and then have different versions of this module for different
  !   physics versions. Specifically, the new variables will have
  !   different dimensions when JULES routines are not used.
  !
  ! Current Code Owner: Martin Best
  !
  ! Code Description:
  !   Language: FORTRAN 90
  !
  ! Declarations:

  IMPLICIT NONE

!========================
! extra JULES variables
!========================
  REAL , DIMENSION(:,:), ALLOCATABLE ::  CLAPP_LEVS
  REAL , DIMENSION(:,:), ALLOCATABLE ::  SATHH_LEVS
  REAL , DIMENSION(:,:), ALLOCATABLE ::   HCAP_LEVS
  REAL , DIMENSION(:,:), ALLOCATABLE ::   HCON_LEVS
  REAL , DIMENSION(:,:), ALLOCATABLE :: SATCON_LEVS
  REAL , DIMENSION(:,:), ALLOCATABLE :: SMVCCL_LEVS
  REAL , DIMENSION(:,:), ALLOCATABLE :: SMVCWT_LEVS
  REAL , DIMENSION(:,:), ALLOCATABLE :: SMVCST_LEVS

! snow variable
!--------------------
  REAL, ALLOCATABLE ::                                     &
    snowdep_surf(:,:)
                    ! Depth of snow (=snowdepth) for all   
                    ! surfaces except those using the 
                    ! snow canopy, for which it is the 
                    ! depth of snow in the canopy (m)
  REAL, ALLOCATABLE ::                                     & 
    albobs_scaling(:,:,:) 
                    ! scaling factor applied to match observed 
                    ! albedo, within limits, for each tile, for 
                    ! SW radn or VIS and NIR  

END MODULE jules_mod

