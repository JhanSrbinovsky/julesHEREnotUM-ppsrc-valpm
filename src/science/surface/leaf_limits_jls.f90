! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Purpose:
! Calculates the leaf resistance and net photosynthesis using:
!  (i) Collatz et al. (1992) C3 photosynthesis model
! (ii) Jacobs (1994) CI/CA closure.
!**********************************************************************
SUBROUTINE leaf_limits(land_field,veg_pts,veg_index,ft            &
,                      nl,dq,apar,tl,ca,oa,pstar,fsmc             &
,                      clos_pts,open_pts,clos_index,open_index    &
,                      ci,rd,wcarb,wexpt,wlite)

USE c_0_dg_c, ONLY : zerodegc
USE pftparm
USE surf_param

! Need to work out how to set this
!      REAL                                                              &
!     & FD(NPFT)                                                         &
!                              ! Dark respiration coefficient.
!     &,NEFF(NPFT)             ! Constant relating VCMAX and leaf N

!      DATA FD      /   0.015,  0.015,   0.015,     0.025, 0.015 /
!      DATA NEFF    /  0.8E-3 , 0.8E-3 , 0.8E-3 ,  0.4E-3, 0.8E-3  /

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

!  (mol/sec) / (watts) conversion for PAR:
REAL, PARAMETER      ::  conpar = 2.19e5

INTEGER                                                           &
 land_field                                                       &
                            ! IN Total number of land points.
,veg_pts                                                          &
                            ! IN Number of vegetated points.
,veg_index(land_field)                                            &
                            ! IN Index of vegetated points
                            !    on the land grid.
,ft                         ! IN Plant functional type.

INTEGER                                                           &
 clos_index(land_field)                                           &
                            ! OUT Index of land points
                            !     with closed stomata.
,clos_pts                                                         &
                            ! OUT Number of land points
                            !     with closed stomata.
,open_index(land_field)                                           &
                            ! OUT Index of land points
                            !     with open stomata.
,open_pts                   ! OUT Number of land points
!                                 !     with open stomata.

REAL                                                              &
 nl(land_field)                                                   &
                            ! IN Leaf nitrogen
!                                 !    concentration (kg N/kg C).
,dq(land_field)                                                   &
                            ! IN Canopy level specific humidity
!                                 !    deficit (kg H2O/kg air).
,apar(land_field)                                                 &
                            ! IN Absorbed PAR (W/m2)
,tl(land_field)                                                   &
                            ! IN Leaf temperature (K).
,ca(land_field)                                                   &
                            ! IN Canopy CO2 pressure (Pa).
,oa(land_field)                                                   &
                            ! IN Atmospheric O2 pressure (Pa).
,pstar(land_field)                                                &
                            ! IN Atmospheric pressure (Pa).
,fsmc(land_field)                                                 &
                            ! IN Soil water factor.
,ci(land_field)                                                   &
                            ! OUT Internal CO2 pressure (Pa).
,rd(land_field)                                                   &
                            ! OUT Dark respiration (mol CO2/m2/s).
,wcarb(land_field)                                                &
                            ! OUT Carboxylation, ...
,wlite(land_field)                                                &
                            !     ... Light, and ...
,wexpt(land_field)                                                &
                            !     ... export limited gross ...
!                                 !     ... photosynthetic rates ...
!                                 !     ... (mol CO2/m2/s).
,acr(land_field)                                                  &
                            ! WORK Absorbed PAR
!                                 !      (mol photons/m2/s).
,ccp(land_field)                                                  &
                            ! WORK Photorespiratory compensatory
!                                 !      point (mol/m3).
,denom(land_field)                                                &
                            ! WORK Denominator in equation for VCM
,kc(land_field)                                                   &
                            ! WORK Michaelis constant for CO2 (Pa)
,ko(land_field)                                                   &
                            ! WORK Michaelis constant for O2 (Pa).
,qtenf(land_field)                                                &
                            ! WORK Q10 function.
,tau(land_field)                                                  &
                            ! WORK CO2/O2 specificity ratio.
,tdegc(land_field)                                                &
                            ! WORK Leaf temperature (deg C).
,vcm(land_field)                                                  &
                            ! WORK Maximum rate of carboxylation
!                                 !      of Rubisco (mol CO2/m2/s).
,vcmax(land_field)          ! WORK Maximum rate of carboxylation
!                                 !      of Rubisco - without the
!                                 !      temperature factor
!                                 !      (mol CO2/m2/s).

INTEGER                                                           &
 j,l                        ! WORK Loop counters.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('LEAF_LIMITS',zhook_in,zhook_handle)

!----------------------------------------------------------------------
! Initialise counters
!----------------------------------------------------------------------
clos_pts = 0
open_pts = 0

DO j=1,veg_pts
  l = veg_index(j)

!----------------------------------------------------------------------
! Calculate the points with closed stomata
!----------------------------------------------------------------------

  IF (fsmc(l)==0.0 .OR. dq(l) >= dqcrit(ft)                       &
                     .OR. apar(l)==0.0) THEN
    clos_pts = clos_pts + 1
    clos_index(clos_pts) = j
  ELSE
    open_pts = open_pts + 1
    open_index(open_pts) = j
  END IF


END DO

!----------------------------------------------------------------------
! Calculate the photosynthetic parameters
!----------------------------------------------------------------------
!DIR$ IVDEP
DO j=1,open_pts
  l = veg_index(open_index(j))

  vcmax(l) = neff(ft) * nl(l)
  tdegc(l) = tl(l) - zerodegc

! TAU is the Rubisco specificity for CO2 relative to O2. The numbers
! in this equation are from Cox, HCTN 24, "Description ... Vegetation
! Model", equation 53.
  tau(l) = 2600.0 * (0.57 ** (0.1 * (tdegc(l) - 25.0)))
  ccp(l) = 0.5 * oa(l) / tau(l) * REAL(c3(ft))


  qtenf(l) = vcmax(l) * (q10_leaf ** (0.1 * (tdegc(l) - 25.0)))
  denom(l) = (1 + EXP (0.3 * (tdegc(l) - tupp(ft))))              &
           * (1 + EXP (0.3 * (tlow(ft) - tdegc(l))))
  vcm(l) = qtenf(l) / denom(l)  ! Cox, HCTN 24, equation 49.
  rd(l) = fd(ft) * vcm(l)

!----------------------------------------------------------------------
! Calculate the internal CO2 pressure (Jacobs, 1994).
!----------------------------------------------------------------------
  ci(l) = (ca(l) - ccp(l)) * f0(ft)                               &
        * (1 - dq(l) / dqcrit(ft)) + ccp(l)

!----------------------------------------------------------------------
! Convert absorbed PAR into mol PAR photons/m2/s
!----------------------------------------------------------------------
  acr(l) = apar(l) / conpar

END DO

!DIR$ IVDEP
DO j=1,clos_pts
  l = veg_index(clos_index(j))

  vcmax(l) = neff(ft) * nl(l)
  tdegc(l) = tl(l) - zerodegc

! TAU is the Rubisco specificity for CO2 relative to O2. The numbers
! in this equation are from Cox, HCTN 24, "Description ... Vegetation
! Model", equation 53.
  tau(l) = 2600.0 * (0.57 ** (0.1 * (tdegc(l) - 25.0)))
  ccp(l) = 0.5 * oa(l) / tau(l) * REAL( c3(ft) )

  qtenf(l) = vcmax(l) * (q10_leaf ** (0.1 * (tdegc(l) - 25.0)))
  denom(l) = (1 + EXP (0.3 * (tdegc(l) - tupp(ft))))             &
           * (1 + EXP (0.3 * (tlow(ft) - tdegc(l))))
  vcm(l) = qtenf(l) / denom(l)  ! Cox, HCTN 24, equation 49.

  rd(l) = fd(ft) * vcm(l)

!----------------------------------------------------------------------
! Calculate the internal CO2 pressure (Jacobs, 1994).
!----------------------------------------------------------------------
  ci(l) = (ca(l) - ccp(l)) * f0(ft)                               &
        * (1 - dq(l) / dqcrit(ft)) + ccp(l)

END DO


!----------------------------------------------------------------------
! Calculate the gross photosynthesis for RuBP-Carboxylase, Light and
! Export limited photosynthesis (Collatz et al., 1992).
!----------------------------------------------------------------------
IF (c3(ft) == 1) THEN

!DIR$ IVDEP
  DO j=1,open_pts
    l = veg_index(open_index(j))
! The numbers
! in these 2 equations are from Cox, HCTN 24, "Description ... Vegetation
! Model", equations 54 and 55.
    kc(l) = 30.0 * (2.1 ** (0.1 * (tdegc(l) - 25.0)))
    ko(l) = 30000.0 * (1.2 ** (0.1 * (tdegc(l) - 25.0)))

    wcarb(l) = vcm(l) * (ci(l) - ccp(l))                          &
             / (ci(l) + kc(l) * (1. + oa(l) / ko(l)))

    wlite(l) = alpha(ft) * acr(l) * (ci(l) - ccp(l))              &
             / (ci(l) + 2 * ccp(l))
    wlite(l) = MAX(wlite(l), TINY(1.0e0))

    wexpt(l) = fwe_c3 * vcm(l)
  END DO

ELSE

  DO j=1,open_pts

    l = veg_index(open_index(j))

    wcarb(l) = vcm(l)

    wlite(l) = alpha(ft) * acr(l)
    wlite(l) = MAX(wlite(l), TINY(1.0e0))

    wexpt(l) = fwe_c4 * vcm(l) * ci(l) / pstar(l)

  END DO

END IF

IF (lhook) CALL dr_hook('LEAF_LIMITS',zhook_out,zhook_handle)
RETURN
END SUBROUTINE leaf_limits
