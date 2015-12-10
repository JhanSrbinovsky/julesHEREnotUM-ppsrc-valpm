! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    SUBROUTINE HYD_CON_CH---------------------------------------------

! Description:
!     Calculates the hydraulic conductivity using Clapp and Hornberger
!     relationships.

! Documentation : UM Documentation Paper 25

! Subroutine Interface:
SUBROUTINE hyd_con_ch (npnts,soil_pts,soil_index,b,ks,thetak,k    &
,                      dk_dthk                                    &
! LOGICAL LTIMER
,ltimer                                                           &
)


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

! Subroutine arguments:
!   Scalar arguments with intent(IN) :
INTEGER                                                           &
 npnts                                                            &
                  ! IN points in grid
,soil_pts         ! IN Number of soil points.

!   Array arguments with intent(IN) :
INTEGER                                                           &
 soil_index(npnts)! IN Array of soil points.

REAL                                                              &
 b(npnts)                                                         &
                  ! IN Exponent in conductivity and soil water
!                       !    suction fits.
,ks(npnts)                                                        &
                  ! IN The saturated hydraulic conductivity (kg/m2
,thetak(npnts)    ! IN Fractional saturation.

LOGICAL ltimer    ! Logical switch for TIMER diags

!   Array arguments with intent(OUT) :
REAL                                                              &
 k(npnts)                                                         &
                  ! OUT The hydraulic conductivity (kg/m2/s).
,dk_dthk(npnts)   ! OUT The rate of change of K with THETAK
!                       !     (kg/m2/s).

! Local scalars:
INTEGER                                                           &
 i,j              ! WORK Loop counter.

REAL small_value

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('HYD_CON_CH',zhook_in,zhook_handle)
small_value=EPSILON(0.0)

!CDIR NODEP
DO j=1,soil_pts
  i=soil_index(j)

  dk_dthk(i)=0.0
  IF ( (thetak(i) >= small_value) .AND. (thetak(i) < 1.0) ) THEN
    k(i)=ks(i)*thetak(i)**(2*b(i)+3)
    dk_dthk(i)=(2*b(i)+3)*ks(i)*(thetak(i)**(2*b(i)+2))
  ELSE IF (thetak(i) < small_value) THEN
    k(i)=0.0
  ELSE
    k(i)=ks(i)
  END IF

END DO

IF (lhook) CALL dr_hook('HYD_CON_CH',zhook_out,zhook_handle)
RETURN
END SUBROUTINE hyd_con_ch
