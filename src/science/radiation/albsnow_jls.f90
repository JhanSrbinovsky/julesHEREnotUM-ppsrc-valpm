! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Routine to calculate spectral snow albedos for MOSES II, based on
! the Marshall (1989) parametrization of the Wiscombe and Warren (1980)
! model. Influence of contaminants in snow has not been included - see
! UM vn4.5 deck FTSA1A.

! *********************************************************************
SUBROUTINE albsnow (p_field,land_field,land_index,                &
                    ntiles,tile_index,tile_pts,l_aggregate,       &
                    cosz,rgrain,snowdepth,soot,alb_snow)

USE nstypes

USE rad_param, ONLY : r0, amax

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Subroutine arguments
!   Logical arguments with intent(in):
LOGICAL                                                           &
 l_aggregate                 ! Logical to set aggregate
                             ! surface scheme

!   Scalar arguments with intent(in):
INTEGER                                                           &
 p_field                                                          &
                             ! Total number of grid points.
,land_field                                                       &
                             ! Number of land points.
,ntiles                      ! Number of surface tiles.

!   Array arguments with intent(in):
INTEGER                                                           &
 land_index(land_field)                                           &
                             ! Index of land points.
,tile_pts(ntype)                                                  &
                             ! Number tile points.
,tile_index(land_field,ntype)! Index of tile points.

REAL                                                              &
 cosz(p_field)                                                    &
                             ! Zenith cosine.
,rgrain(land_field,ntiles)                                        &
                             ! Snow grain size (microns).
,snowdepth(land_field,ntiles)                                     &
                             ! Depth of snow on the ground (m).
,soot(p_field)               ! Snow soot content (kg/kg).

!   Array arguments with intent(out):
REAL                                                              &
 alb_snow(land_field,ntype,4)! Snow albedo.
!                              (*,*,1) - Direct beam visible
!                              (*,*,2) - Diffuse visible
!                              (*,*,3) - Direct beam near-IR
!                              (*,*,4) - Diffuse near-IR

! Local scalars:
REAL                                                              &
 reff                                                             &
                             ! Zenith effective grain size.
,sigma                       ! Scaled soot content

INTEGER                                                           &
 band,i,j,l,n                ! Loop counters.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


IF (lhook) CALL dr_hook('ALBSNOW',zhook_in,zhook_handle)
DO n=1,ntiles
  DO l=1,land_field
    i = land_index(l)
    IF (snowdepth(l,n)  >   0.) THEN
      reff = rgrain(l,n) * ( 1. + 0.77*(cosz(i)-0.65) )**2
      alb_snow(l,n,1) = amax(1) - 0.002*(SQRT(reff) - SQRT(r0))
      alb_snow(l,n,2) = amax(1) -                                 &
                        0.002*(SQRT(rgrain(l,n)) - SQRT(r0))
      alb_snow(l,n,3) = amax(2) - 0.09*ALOG(reff/r0)
      alb_snow(l,n,4) = amax(2) - 0.09*ALOG(rgrain(l,n)/r0)
    END IF
  END DO
END DO

! Adjust visible snow albedos for soot content
DO n=1,ntiles
  DO j=1,tile_pts(n)
    l = tile_index(j,n)
    IF ( snowdepth(l,n)  >   0. ) THEN
      i = land_index(l)
      sigma = soot(i) * rgrain(l,n) / 0.0017
      IF ( sigma  >   1. ) THEN
        alb_snow(l,n,1) = 0.07 +                                  &
                          0.5*(alb_snow(l,n,1) - 0.07)/(sigma**0.46)
        alb_snow(l,n,2) = 0.07 +                                  &
                          0.5*(alb_snow(l,n,2) - 0.07)/(sigma**0.46)
      ELSE
        alb_snow(l,n,1) = alb_snow(l,n,1) -                       &
                          0.5*(alb_snow(l,n,1) - 0.07)*(sigma**0.6)
        alb_snow(l,n,2) = alb_snow(l,n,2) -                       &
                          0.5*(alb_snow(l,n,2) - 0.07)*(sigma**0.6)
      END IF
    END IF
  END DO
END DO

! Adjust near-IR snow albedos for soot content
DO n=1,ntiles
  DO j=1,tile_pts(n)
    l = tile_index(j,n)
    IF ( snowdepth(l,n)  >   0. ) THEN
      i = land_index(l)
      sigma = soot(i) * rgrain(l,n) / 0.004
      IF ( sigma  >   1. ) THEN
        alb_snow(l,n,3) = 0.06 +                                  &
                          0.5*(alb_snow(l,n,3) - 0.06)/(sigma**0.6)
        alb_snow(l,n,4) = 0.06 +                                  &
                          0.5*(alb_snow(l,n,4) - 0.06)/(sigma**0.6)
      ELSE
        alb_snow(l,n,3) = alb_snow(l,n,3) -                       &
                          0.5*(alb_snow(l,n,3) - 0.06)*(sigma**0.7)
        alb_snow(l,n,4) = alb_snow(l,n,4) -                       &
                          0.5*(alb_snow(l,n,4) - 0.06)*(sigma**0.7)
      END IF
    END IF
  END DO
END DO

IF ( l_aggregate ) THEN
  DO band=1,4
    DO n=2,ntype
      DO j=1,tile_pts(n)
        l = tile_index(j,n)
        alb_snow(l,n,band) = alb_snow(l,1,band)
      END DO
    END DO
  END DO
END IF

IF (lhook) CALL dr_hook('ALBSNOW',zhook_out,zhook_handle)
RETURN
END SUBROUTINE albsnow
