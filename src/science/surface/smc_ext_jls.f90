! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    SUBROUTINE SMC_EXT-----------------------------------------------

! Description:
!     Calculates the soil moisture availability factor and
!     the fraction of the transpiration which is extracted from each
!     soil layer.

! Documentation : UM Documentation Paper 25
!---------------------------------------------------------------------
! Subroutine Interface:
SUBROUTINE smc_ext (npnts,nshyd,tile_pts,tile_index               &
,                   f_root,sthu,v_crit,v_sat,v_wilt               &
,                   wt_ext,fsmc)


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

! Subroutine arguments:
!   Scalar arguments with intent(IN) :
INTEGER                                                           &
 npnts                                                            &
                      ! IN Number of gridpoints.
,nshyd                                                            &
                      ! IN Number of soil moisture layers.
,tile_pts                                                         &
                      ! IN Number of points containing the
!                           !    given surface type.
,tile_index(npnts)    ! IN Indices on the land grid of the
!                           !    points containing the given
!                           !    surface type.


!   Array arguments with intent(IN) :
REAL                                                              &
 f_root(nshyd)                                                    &
                      ! IN Fraction of roots in each soil
!                           !    layer.
,sthu(npnts,nshyd)                                                &
                      ! IN Unfrozen soil moisture content of
!                           !    each layer as a fraction of
!                           !    saturation.
,v_crit(npnts,nshyd)                                              &
                      ! IN Volumetric soil moisture
!                           !    concentration above which
!                           !    evapotranspiration is not sensitive
!                           !    to soil water (m3 H2O/m3 soil).
,v_sat(npnts,nshyd)                                               &
                      ! IN Volumetric soil moisture
!                           !    concentration at saturation
!                           !    (m3 H2O/m3 soil).
,v_wilt(npnts,nshyd)  ! IN Volumetric soil moisture
!                           !    concentration below which
!                           !    stomata close (m3 H2O/m3 soil).

!   Array arguments with intent(INOUT) :
REAL                                                              &
 wt_ext(npnts,nshyd)  ! OUT Cummulative fraction of transpiration
!                           !     extracted from each soil layer
!                           !     (kg/m2/s).


!   Array arguments with intent(OUT) :
REAL                                                              &
 fsmc(npnts)          ! OUT Soil moisture availability
!                           !     factor.

! Local scalars:
INTEGER                                                           &
 i,j,n                ! WORK Loop counters

! Local arrays:
REAL                                                              &
 fsmc_l(npnts,nshyd)  ! WORK Soil moisture availability
!                           !      factor for each soil layer.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!----------------------------------------------------------------------
! Initialisations
!----------------------------------------------------------------------
IF (lhook) CALL dr_hook('SMC_EXT',zhook_in,zhook_handle)
DO i=1,npnts
  fsmc(i)=0.0
END DO

!----------------------------------------------------------------------
! Calculate the soil moisture availability factor for each layer and
! weight with the root fraction to calculate the total availability
! factor.
!----------------------------------------------------------------------
DO n=1,nshyd
!CDIR NODEP
  DO j=1,tile_pts
    i=tile_index(j)

    IF ( abs(v_crit(i,n)-v_wilt(i,n)) > 0.0 ) THEN
      fsmc_l(i,n)=(sthu(i,n)*v_sat(i,n)-v_wilt(i,n))                &
                 /(v_crit(i,n)-v_wilt(i,n))
    ELSE
      fsmc_l(i,n)=0.0
    END IF

    fsmc_l(i,n)=MAX(fsmc_l(i,n),0.0)
    fsmc_l(i,n)=MIN(fsmc_l(i,n),1.0)

    fsmc(i)=fsmc(i)+f_root(n)*fsmc_l(i,n)

  END DO
END DO

!----------------------------------------------------------------------
! Calculate the fraction of the tranpiration which is extracted from
! each soil layer.
!----------------------------------------------------------------------
DO n=1,nshyd
  DO j=1,tile_pts
    i=tile_index(j)
    IF (fsmc(i)  >   0.0)                                         &
      wt_ext(i,n) = f_root(n)*fsmc_l(i,n)/fsmc(i)
  END DO
END DO

IF (lhook) CALL dr_hook('SMC_EXT',zhook_out,zhook_handle)
RETURN
END SUBROUTINE smc_ext
