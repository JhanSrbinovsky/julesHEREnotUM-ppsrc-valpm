! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    Subroutine VGRAV -------------------------------------------------

! Purpose: To calculate the gravitational sedimentation velocity of
!          tracer particles according to Stoke's law, including the
!          Cunningham correction factor.

!Documentation: Ref. Pruppacher & Klett
!                    Microphysics of clouds & ppn    1978,1980 edns.

!-----------------------------------------------------------------------

SUBROUTINE vgrav(                                                 &
                 nlevs,diam,rhop,p,t,vstokes,ccf,etaa )

USE atm_fields_bounds_mod

USE dust_param, ONLY: accf, bccf, cccf
USE c_sulchm
USE c_0_dg_c
USE c_g

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE


INTEGER nlevs              !IN number of levels

REAL diam                  !IN particle diameter
REAL rhop                  !IN particles density
REAL p(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nlevs) !IN pressure
REAL t(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nlevs) !IN temperature

REAL vstokes( tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nlevs) !OUT sedimentation velocity
REAL etaa(    tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nlevs) !OUT viscosity of air
REAL ccf(     tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nlevs) !OUT cunningham correction factor

! local variables

INTEGER ilev               !LOC loop counter for levels
INTEGER i                  !LOC loop counter
INTEGER j                  !LOC loop counter

REAL tc(       tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end) !LOC temperature in deg C
REAL lamdaa(   tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end) !LOC mean free path of particle
REAL alphaccf( tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end) !LOC

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! Calculate viscosity of air (Pruppacher & Klett p.323)
IF (lhook) CALL dr_hook('VGRAV',zhook_in,zhook_handle)

!$OMP  PARALLEL DEFAULT(SHARED) PRIVATE(tc, ilev, j, lamdaa, i,    &
!$OMP& alphaccf)

!$OMP DO SCHEDULE(STATIC)
DO ilev=1,nlevs
  DO j=tdims%j_start,tdims%j_end
    DO i = tdims%i_start,tdims%i_end
     tc(i,j)=t(i,j,ilev)-zerodegc
     IF (tc(i,j)  >=  0.) THEN
       etaa(i,j,ilev)=(1.718+0.0049*tc(i,j))*1.e-5
     ELSE
       etaa(i,j,ilev)=                                             &
            (1.718+0.0049*tc(i,j)-1.2e-5*tc(i,j)*tc(i,j))*1.e-5
     END IF

    END DO !ROW_LENGTH
  END DO !ROWS
END DO !NLEVS
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
DO ilev=1,nlevs
  DO j=tdims%j_start,tdims%j_end
    DO i=tdims%i_start,tdims%i_end

! Calculate mean free path of particle (Pruppacher & Klett p.323)
      lamdaa(i,j)=mfp_ref*pref_mfp*t(i,j,ilev)/                    &
      (p(i,j,ilev)*tref_mfp)
! Calculate Cunningham correction factor(Pruppacher & Klett p.361)
      alphaccf(i,j)=accf+bccf*EXP(cccf*diam*.5/lamdaa(i,j))
      ccf(i,j,ilev)=(1.+alphaccf(i,j)*lamdaa(i,j)/(.5*diam))
! Calculate sedimentation velocity (Pruppacher & Klett p.362)
      vstokes(i,j,ilev)=ccf(i,j,ilev)*(diam*diam*g*rhop)/          &
             (18.*etaa(i,j,ilev))

    END DO !ROW_LENGTH
  END DO !ROWS
END DO !NLEV
!$OMP END DO

!$OMP END PARALLEL

IF (lhook) CALL dr_hook('VGRAV',zhook_out,zhook_handle)
RETURN
END SUBROUTINE vgrav
