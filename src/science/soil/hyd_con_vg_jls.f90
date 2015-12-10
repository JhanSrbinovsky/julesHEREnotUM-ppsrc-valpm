! *****************************COPYRIGHT*******************************
! (c) CROWN COPYRIGHT 2000, Met Office, All Rights Reserved.
! Please refer to file $UMDIR/vn$VN/copyright.txt for further details
! *****************************COPYRIGHT*******************************
!    SUBROUTINE HYD_CON_VG---------------------------------------------

! Description:
!     Calculates the hydraulic conductivity using Van Genuchten curves

! Documentation : UM Documentation Paper 25

! Subroutine Interface:
SUBROUTINE hyd_con_vg(npnts,soil_pts,soil_index,b,ks,thetak,k,    &
                    dk_dthk,ltimer )


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
REAL                                                              &
 bracket(npnts)                                                   &
                  ! WORK 1-S^(b+1)
,dbr_dthk(npnts)                                                  &
                  ! WORK The rate of change of BRACKET
!                       !       with THETAK.
,kred(npnts)                                                      &
                  ! WORK KSAT*S^L_WAG (kg/m2/s).
,sdum(npnts)      ! WORK Bounded THETAK value.

INTEGER                                                           &
 i,j              ! WORK Loop counter.

REAL                                                              &
 l_wag            ! Exponent in the van Mualem / Van Genuchten
!                       ! fit to the hydraulic conduactivity curve.
PARAMETER (l_wag=0.5)

REAL                                                              &
 theta_min                                                        &
                  ! Minimum value of THETAK for K calculation.
,theta_max        ! Maximum value of THETAK for K calculation.
PARAMETER (theta_min=0.05, theta_max=0.95)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('HYD_CON_VG',zhook_in,zhook_handle)

!CDIR NODEP
DO j=1,soil_pts
  i=soil_index(j)

  sdum(i)=MAX(thetak(i),theta_min)
  sdum(i)=MIN(sdum(i),theta_max)

  dk_dthk(i)=0.0
  bracket(i)=1-sdum(i)**(b(i)+1)
  kred(i)=ks(i)*sdum(i)**l_wag

  k(i)=kred(i)*(1-bracket(i)**(1.0/(b(i)+1)))**2

!----------------------------------------------------------------------
! To avoid blow-up of implicit increments approximate by piecewise
! linear functions
! (a) for THETA>THETA_MAX  (ensuring that K=KS at THETA=1)
! (b) for THETA<THETA_MIN  (ensuring that K=0 at THETA=THETA_MIN)
!----------------------------------------------------------------------
  IF (thetak(i).lt.theta_min) THEN
    dk_dthk(i)=k(i)/theta_min
    k(i)=k(i)+dk_dthk(i)*(MAX(thetak(i),0.0)-theta_min)
  ELSE IF (thetak(i).gt.theta_max) THEN
    dk_dthk(i)=(ks(i)-k(i))/(1.0-theta_max)
    k(i)=k(i)+dk_dthk(i)*(MIN(thetak(i),1.0)-theta_max)
  ELSE
    dbr_dthk(i)=-(b(i)+1)*sdum(i)**b(i)
    dk_dthk(i)=l_wag*k(i)/sdum(i)                                 &
            -2*kred(i)/(b(i)+1)                                   &
            *(1-bracket(i)**(1.0/(b(i)+1)))                       &
            *(bracket(i)**(-b(i)/(b(i)+1)))                       &
            *dbr_dthk(i)
  END IF

  IF ((thetak(i).gt.1.0).OR.(thetak(i).lt.0.0)) THEN
    dk_dthk(i)=0.0
    k(i)=MAX(k(i),0.0)
    k(i)=MIN(k(i),ks(i))
  END IF

END DO

IF (lhook) CALL dr_hook('HYD_CON_VG',zhook_out,zhook_handle)
RETURN
END SUBROUTINE hyd_con_vg
