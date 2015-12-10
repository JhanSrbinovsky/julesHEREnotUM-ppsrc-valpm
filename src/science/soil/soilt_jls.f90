! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    SUBROUTINE SOILT--------------------------------------------------

! Description:
!     Diagnoses the soil temperature in a layer at the surface


! Subroutine Interface:
SUBROUTINE soilt(npnts,nshyd,soil_pts,soil_index                  &
  ,dz,tsoil,tsoil_d )

USE soil_param, ONLY :                                            &
 zst                  ! Depth of layer for soil moisture
!                           ! diagnostic (m).

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

INTEGER                                                           &
 npnts                                                            &
                      ! IN Number of gridpoints.
,nshyd                                                            &
                      ! IN Number of soil moisture levels.
,soil_pts                                                         &
                      ! IN Number of soil points.
,soil_index(npnts)    ! IN Array of soil points.

REAL                                                              &
 dz(nshyd)                                                        &
                      ! IN Thicknesses of the soil layers (m).
,tsoil(npnts,nshyd)   ! IN Soil temperature (K)

REAL                                                              &
 tsoil_d(npnts)       ! OUT Soil temperature in layer (K)


REAL                                                              &
 z1,z2                ! WORK Depth of the top and bottom of the
!                           !      soil layers (m).

INTEGER                                                           &
 i,j,n

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('SOILT',zhook_in,zhook_handle)

z2 = 0.

DO i=1,npnts
  tsoil_d(i) = 0.0
END DO

DO n=1,nshyd
  z1 = z2
  z2 = z2 + dz(n)
  IF ( z2 <  zst ) THEN
    DO j=1,soil_pts
      i=soil_index(j)
      tsoil_d(i) = tsoil_d(i) + dz(n) * tsoil(i,n)
    END DO
  ELSE IF ( z2 >= zst .AND. z1 <  zst ) THEN
    DO j=1,soil_pts
      i = soil_index(j)
      tsoil_d(i) = tsoil_d(i) + ( zst - z1 ) * tsoil(i,n)
    END DO
  END IF

END DO

IF (lhook) CALL dr_hook('SOILT',zhook_out,zhook_handle)
RETURN
END SUBROUTINE soilt
