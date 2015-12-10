! *****************************COPYRIGHT*******************************
! (c) CROWN COPYRIGHT 2000, Met Office, All Rights Reserved.
! Please refer to file $UMDIR/vn$VN/copyright.txt for further details
! *****************************COPYRIGHT*******************************

!  L  SUBROUTINE DEWPNT-------------------------------------------------

!    Purpose: Calculates the 1.5 metre dewpoint from 1.5 metre specific
!             humidity, 1.5 metre temperature and 1.5 metre pressure.

!    Suitable for single column usage.

!  Arguments:----------------------------------------------------------
SUBROUTINE dewpnt(q,p,t,p_field,td)

USE c_epslon
USE c_r_cp
USE c_lheat
USE c_0_dg_c

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

INTEGER p_field         ! IN Size of field arrays.

REAL, INTENT(IN) ::                                               &
 p(p_field)                                                       &
                        ! IN Pressure.
,q(p_field)                                                       &
                        ! IN Specific humidity.
,t(p_field)             ! IN Temperature.

REAL, INTENT(OUT) ::                                              &
 td(p_field)            ! OUT Dew point.


!----------------------------------------------------------------------
! Workspace - local variables
!----------------------------------------------------------------------
REAL ::                                                           &
 rv                                                               &
                        ! LOCAL Gas constant for water vapour.
,rl1                                                              &
                        ! LOCAL Latent heat of evaporation.
,rt                                                               &
                        ! LOCAL.
,p1(p_field)                                                      &
                        ! LOCAL Pressure.
!                                     j/Kg at 0 deg C.
,rl(p_field)                                                      &
                        ! LOCAL.
,q0(p_field)                                                      &
                        ! LOCAL local SH.
,es0                                                              &
                        ! LOCAL Saturated vapour pressure.
,v_pres(p_field)        ! LOCAL Vapour pressure.

INTEGER ::                                                        &
 i                      ! LOCAL loop variable.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

PARAMETER ( rv = r / repsilon )
PARAMETER ( rl1 = -2.73e3 )

!-----------------------------------------------------------------------
! EXTERNAL SUBROUTINES CALLED
!-----------------------------------------------------------------------
EXTERNAL  qsat_wat


!----------------------------------------------------------------------
!  Calculate P in HPa.
!----------------------------------------------------------------------
IF (lhook) CALL dr_hook('DEWPNT',zhook_in,zhook_handle)
DO i=1,p_field
  p1(i) = p(i) / 100.0
!----------------------------------------------------------------------
!  Calculate RL - The latent heat of evaporation.
!----------------------------------------------------------------------
  rl(i) = lc + rl1 * ( t(i) - tm )
!----------------------------------------------------------------------
!  Calculate Vapour pressure, and from that the dewpoint in Kelvins.
!----------------------------------------------------------------------
  v_pres(i) = q(i) * p1(i) / ( REPSILON + q(i))
END DO

! DEPENDS ON: qsat_wat
CALL qsat_wat(q0,t,p,p_field)

DO i=1,p_field
  IF (v_pres(i) .gt. 0.0) THEN
    es0 = (q0(i) * p1(i)) / (REPSILON + q0(i))
    rt = (1 / t(i)) - ( rv * ALOG(v_pres(i)/es0) )/rl(i)
    td(i) = 1.0 / rt
    IF (td(i) .gt. t(i)) td(i) = t(i)
  ELSE
    td(i)=0.0
  END IF
END DO
IF (lhook) CALL dr_hook('DEWPNT',zhook_out,zhook_handle)
RETURN
END SUBROUTINE dewpnt
