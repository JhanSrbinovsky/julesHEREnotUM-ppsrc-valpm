! *****************************COPYRIGHT*******************************

! (c) [University of Edinburgh] [2009]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the JULES collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC237]

! *****************************COPYRIGHT*******************************
!  SUBROUTINE TRIDAG---------------------------------------------------

! Description:
!     Solve tridiagonal matrix equation, numerical Recipes p43

! Subroutine Interface:
SUBROUTINE tridag(nvec,nmax,a,b,c,r,s)


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

! Scalar arguments with intent(in)
INTEGER, INTENT(IN) ::                                            &
 nvec                                                             &
                !  Vector length
,nmax           !  Maximum vector length

! Array arguments with intent(in)
REAL, INTENT(IN) ::                                               &
 a(nmax)                                                          &
                ! Below-diagonal matrix elements
,b(nmax)                                                          &
                ! Diagonal matrix elements
,c(nmax)                                                          &
                ! Above-diagonal matrix elements
,r(nmax)        ! Matrix equation rhs

! Array arguments with intent(out)
REAL, INTENT(OUT) ::                                              &
 s(nmax)        ! Solution vector

! Local scalars
INTEGER ::                                                        &
 n              ! loop counter

REAL ::                                                           &
 beta           ! work

! Local arrays
REAL ::                                                           &
 GAMMA(nvec)    ! work

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook('TRIDAG',zhook_in,zhook_handle)
beta = b(1)
s(1) = r(1) / beta

DO n=2,nvec
  GAMMA(n) = c(n-1) / beta
  beta = b(n) - a(n)*GAMMA(n)
  s(n) = ( r(n) - a(n)*s(n-1) ) / beta
END DO

DO n=nvec-1,1,-1
  s(n) = s(n) - GAMMA(n+1)*s(n+1)
END DO
IF (lhook) CALL dr_hook('TRIDAG',zhook_out,zhook_handle)
RETURN

END SUBROUTINE tridag
