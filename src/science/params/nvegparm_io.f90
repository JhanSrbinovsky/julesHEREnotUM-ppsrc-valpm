! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module contains variables used for reading in nvegparm data
! and initialisations
  
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.


MODULE nvegparm_io

  USE max_dimensions, ONLY:                                           &
    nnvg_max

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Set up variables to use in IO (a fixed size version of each array
! in nvegparm that we want to initialise).
!-----------------------------------------------------------------------
  REAL ::                                                             &
    albsnc_nvg_io(nnvg_max),                                          &
    albsnf_nvgu_io(nnvg_max),                                          &
    albsnf_nvg_io(nnvg_max),                                          &
    albsnf_nvgl_io(nnvg_max),                                          &
    catch_nvg_io(nnvg_max),                                           &
    gs_nvg_io(nnvg_max),                                              &
    infil_nvg_io(nnvg_max),                                           &
    z0_nvg_io(nnvg_max),                                              &
    ch_nvg_io(nnvg_max),                                              &
    vf_nvg_io(nnvg_max),                                              &
    emis_nvg_io(nnvg_max),                                            &
    z0hm_nvg_io(nnvg_max),                                            &
    z0hm_classic_nvg_io(nnvg_max)





!-----------------------------------------------------------------------
! Set up a namelist for reading and writing these arrays
!-----------------------------------------------------------------------
  NAMELIST /jules_nvegparm/                                           &



                            albsnc_nvg_io,albsnf_nvgu_io,             &
                            albsnf_nvg_io, albsnf_nvgl_io,            &
                            catch_nvg_io,gs_nvg_io,infil_nvg_io,      &
                            z0_nvg_io,ch_nvg_io,vf_nvg_io,            &
                            emis_nvg_io,z0hm_nvg_io,                  &
                            z0hm_classic_nvg_io
                            
END MODULE nvegparm_io
