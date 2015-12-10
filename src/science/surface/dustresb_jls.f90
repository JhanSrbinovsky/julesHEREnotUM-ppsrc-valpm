! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! subroutine DUSTRESB

! Purpose:
!   To calculate the surface layer resistance for mineral dust

! Code Description:
!  Language: FORTRAN77 + common extensions
!  This code is written to UMDP3 v6 programming standards

! Documentation: "Modelling the atmospheric lifecycle..."
!                 Woodward, JGR106, D16, pp18155-18166
!---------------------------------------------------------------------

 SUBROUTINE dustresb(                                             &
  pstar,tstar,rhostar,aresist,vshr,cd_std_dust,                   &
  r_b_dust                                                        &
  )

USE atm_fields_bounds_mod

USE c_sulchm
USE dust_param, ONLY: ndiv, drep, rhop
USE c_r_cp
USE c_0_dg_c
USE c_pi
USE c_g

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

REAL                                                              &
     !IN
 pstar(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)       &
                              !IN surface pressure
,tstar(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)       &
                              !IN surface temperature
,rhostar(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)     &
                              !IN surface air density
,aresist(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)     &
                              !IN aerodynamic resistance
,vshr(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)        &
                              !IN surface to lowest lev windspeed
!                                   !   difference
,cd_std_dust(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                              !IN surface transfer coeffient for
!                                   !   momentum, excluding orographic
!                                   !   form drag

REAL                                                              &
     !OUT
 r_b_dust(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,ndiv)
                              !OUT surface layer resistance for
!                                     !    mineral dust

!     local variables

INTEGER                                                           &
 idiv                                                             &
      !loop counter, dust divisions
,i                                                                &
      !loop counter
,j                                                                &
      !loop counter
,lev1 !number of levels for vstokes calculation

REAL                                                              &
 nu(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)          &
                           !kinematic viscosity
,etaa(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)        &
                           !dynamic viscosity of air
,lamdaa(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)      &
                           !mean free path of air molecules
,vstokes1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)    &
                           !gravitational settling velocity, lev1
,nstokes(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)     &
                           !stokes number = VstokesVshrVshr/nu g
,nschmidt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)    &
                           !schmidt number = nu/diffusivit
,tc(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)          &
                           !temperature in deg C
,alphaccf(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)    &
                           !alpha in cunningham correction factor
,ccf(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)         &
                           !Cunningham correction factor
,fvsq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)        &
                           !friction velocity squared
,work(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)        &
                           !workspace
,stokes_exp                                                       &
                           !stokes term in R_B_DUST equation
,smallp                    !small +ve number, negligible compared to 1

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

 EXTERNAL vgrav

IF (lhook) CALL dr_hook('DUSTRESB',zhook_in,zhook_handle)

!... epsilon() is defined as almost negligible, so eps/100 is negligible
smallp = EPSILON(1.0) / 100.0

!...calc stokes number, schmidt number and finally resistance

lev1=1

DO idiv=1,ndiv
! DEPENDS ON: vgrav
  CALL vgrav(                                                     &
  lev1,drep(idiv),rhop,pstar,tstar,vstokes1,ccf,etaa              &
  )

!CDIR NOVECTOR
  DO j = tdims%j_start,tdims%j_end
    DO i= tdims%i_start,tdims%i_end
      nschmidt(i,j)=3.*pi*etaa(i,j)*etaa(i,j)*drep(idiv)/         &
       (rhostar(i,j)*boltzmann*tstar(i,j)*ccf(i,j))
      nstokes(i,j)=vstokes1(i,j)*cd_std_dust(i,j)*rhostar(i,j)*   &
       vshr(i,j)*vshr(i,j)/(etaa(i,j)*g)
      ! Avoid underflow in Stokes term by setting to zero if
      ! negligible compared to Schmidt term, i.e., if NSTOKES
      ! is too small.
      IF ( 3.0 / nstokes(i,j) <                                   &
           - LOG10( smallp *nschmidt(i,j)**(-2./3.) ) ) THEN
         stokes_exp = 10.**(-3./nstokes(i,j))
      ELSE
         stokes_exp = 0.0
      END IF
      r_b_dust(i,j,idiv)=1./( SQRT(cd_std_dust(i,j)) *            &
       (nschmidt(i,j)**(-2./3.)+stokes_exp) )
    END DO !ROW_LENGTH
  END DO !ROWS
END DO !NDIV

IF (lhook) CALL dr_hook('DUSTRESB',zhook_out,zhook_handle)
RETURN
END SUBROUTINE dustresb
