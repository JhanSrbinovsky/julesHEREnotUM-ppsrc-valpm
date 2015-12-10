! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!

! Module containing BVOC diagnostics

MODULE bvoc_vars

  IMPLICIT NONE

  REAL, ALLOCATABLE ::                                                        &
    isoprene(:),                                                              &
            ! Gridbox mean isoprene emission flux (kgC/m2/s)
    isoprene_ft(:,:),                                                         &
            ! Isoprene emission flux on PFTs (kgC/m2/s)
    terpene(:),                                                               &
            ! Gridbox mean (mono-)terpene emission flux (kgC/m2/s)
    terpene_ft(:,:),                                                          &
            ! (Mono-)Terpene emission flux on PFTs (kgC/m2/s)
    methanol(:),                                                              &
            ! Gridbox mean methanol emission flux (kgC/m2/s)
    methanol_ft(:,:),                                                         &
            ! Methanol emission flux on PFTs (kgC/m2/s)
    acetone(:),                                                               &
            ! Gridbox mean acetone emission flux (kgC/m2/s)
    acetone_ft(:,:)
            ! Acetone emission flux on PFTs (kgC/m2/s)

END MODULE bvoc_vars
