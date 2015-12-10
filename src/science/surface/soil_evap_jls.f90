! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!**********************************************************************
! SUBROUTINE soil_evap

! Description:
! Subroutine to adjust canopy conductance and soil moisture extraction
! for soil evaporation beneath vegetation.

!**********************************************************************

SUBROUTINE soil_evap (npnts,nshyd,tile_pts,tile_index,            &
                      gsoil,lai,gs,wt_ext,fsoil)


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

INTEGER, INTENT(IN) ::                                            &
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

REAL, INTENT(IN) ::                                               &
 gsoil(npnts)                                                     &
                      ! IN Soil surface conductance (m/s).
,lai(npnts)           ! IN Leaf area index.

REAL, INTENT(INOUT) ::                                            &
 gs(npnts)                                                        &
                      ! INOUT Surface conductance (m/s).
,wt_ext(npnts,nshyd)  ! INOUT Fraction of evapotranspiration
!                           !       extracted from each soil layer.

REAL, INTENT(OUT) ::                                              &
 fsoil(npnts)         ! Fraction of ground below canopy
!                           ! contributing to evaporation.

INTEGER ::                                                        &
 j,k,l                ! Loop indices

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('SOIL_EVAP',zhook_in,zhook_handle)

!CDIR NODEP - order independent so can vectorise

! Initialisations
fsoil(:)=0.0

DO j=1,tile_pts
  l=tile_index(j)
  fsoil(l) = EXP(-0.5*lai(l))
END DO

DO k=2,nshyd
  DO j=1,tile_pts
    l=tile_index(j)
    wt_ext(l,k) = gs(l)*wt_ext(l,k)/(gs(l) + fsoil(l)*gsoil(l))
  END DO
END DO

!CDIR NODEP
DO j=1,tile_pts
  l=tile_index(j)
  wt_ext(l,1) = (gs(l)*wt_ext(l,1) + fsoil(l)*gsoil(l))           &
                 / (gs(l) + fsoil(l)*gsoil(l))
  gs(l) = gs(l) + fsoil(l)*gsoil(l)
END DO

IF (lhook) CALL dr_hook('SOIL_EVAP',zhook_out,zhook_handle)
RETURN
END SUBROUTINE soil_evap
