! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Purpose:
! Subroutine to calculate gridbox mean values of surface conductance
! and carbon fluxes. Also returns net primary productivity, leaf
! turnover and wood respiration of each plant functional type for
! driving TRIFFID.

! Calculates the leaf resistance and net photosynthesis using:
!  (i) Collatz et al. (1992) C3 photosynthesis model, and the
!      Collatz et al. (1991) C4 photosynthesis model.
! (ii) Jacobs (1994) CI/CA closure.
! Written by Peter Cox (February 1996)
! Adapted for MOSES II tile model by Richard Essery (July 1997)
! Coded into UMVN 6.2  by Pete Falloon (July 2006)
!**********************************************************************
SUBROUTINE leaf (land_field,veg_pts,veg_index,ft                  &
,                clos_pts,open_pts,clos_index,open_index          &
,                o3mol,ra,flux_o3                                 &
,                fsmc,tl,ca,ci,rd,wcarb,wexpt,wlite               &
,                gl,al,fo3)

USE pftparm
USE c_rmol
USE surf_param, ONLY : beta1,beta2,ratio,ratio_o3

USE switches, ONLY : l_o3_damage

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

INTEGER                                                           &
 land_field                                                       &
                            ! IN Total number of land points.
,veg_pts                                                          &
                            ! IN Number of vegetated points.
,veg_index(land_field)                                            &
                            ! IN Index of vegetated points
                            !    on the land grid.
,ft                                                               &
                            ! IN Plant functional type.
,clos_index(land_field)                                           &
                            ! IN Index of land points
                            !    with closed stomata.
,clos_pts                                                         &
                            ! IN Number of land points
                            !    with closed stomata.
,open_index(land_field)                                           &
                            ! IN Index of land points
                            !    with open stomata.
,open_pts                   ! IN Number of land points
                            !    with open stomata.

REAL                                                              &
 fsmc(land_field)                                                 &
                            ! IN Soil water factor.
,tl(land_field)                                                   &
                            ! IN Leaf temperature (K).
,rd(land_field)                                                   &
                            ! IN Dark respiration (mol CO2/m2/s).
,ca(land_field)                                                   &
                            ! IN Canopy CO2 pressure (Pa).
,ci(land_field)                                                   &
                            ! IN Internal CO2 pressure (Pa).
,wcarb(land_field)                                                &
,wlite(land_field)                                                &
,wexpt(land_field)                                                &
                            ! IN Carboxylation,
                            !    Light, and
                            !    export limited gross
                            !    photosynthetic rates
                            !    (mol CO2/m2/s).
! Ozone-related variables
,o3mol(land_field)                                                &
                            ! IN Molar concentration of ozone
                            !    at reference level (nmol/m3).
,ra(land_field)                                                   &
                            ! IN Total aerodynamic+boundary layer resistance
                            !    between leaf surface and reference level (s/m).
,gl(land_field)                                                   &
                            ! OUT Leaf conductance for H2O (m/s).
,al(land_field)                                                   &
                            ! OUT Net Leaf photosynthesis
!                                 !     (mol CO2/m2/s).
,flux_o3(land_field)                                              &
                            ! OUT Flux of O3 to stomata (nmol O3/m2/s).
,fo3(land_field)                                                  &
                            ! OUT Ozone exposure factor.
,b1(land_field)                                                   &
,b2(land_field)                                                   &
,b3(land_field)                                                   &
                            ! WORK Coefficients of the quadratic.
,conv(land_field)                                                 &
                            ! WORK Factor for converting mol/m3
!                                 !      into Pa (J/mol).
,glco2(land_field)                                                &
                            ! WORK Leaf conductance for CO2 (m/s).
,wl(land_field)                                                   &
                            ! WORK Gross leaf phtosynthesis
!                                 !      (mol CO2/m2/s).
,wp(land_field)             ! WORK Smoothed minimum of
!                                 !      Carboxylation and Light
!                                 !      limited gross photosynthesis
!                                 !      (mol CO2/m2/s).

! Work variables for ozone calculations
REAL :: b,c


INTEGER                                                           &
 j,l                        ! WORK Loop counters.


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('LEAF',zhook_in,zhook_handle)

!----------------------------------------------------------------------
! Calculate the co-limited rate of gross photosynthesis
!----------------------------------------------------------------------
!DIR$ IVDEP
DO j=1,open_pts
  l = veg_index(open_index(j))

  b1(l) = beta1
  b2(l) = - (wcarb(l) + wlite(l))
  b3(l) = wcarb(l) * wlite(l)

  wp(l) = -b2(l)/(2*b1(l))                                        &
         - SQRT(b2(l)*b2(l)/(4*b1(l)*b1(l)) - b3(l)/b1(l))

END DO

!DIR$ IVDEP
DO j=1,open_pts
  l = veg_index(open_index(j))

  b1(l) = beta2
  b2(l) = - (wp(l) + wexpt(l))
  b3(l) = wp(l) * wexpt(l)

  wl(l) = -b2(l)/(2*b1(l))                                        &
         - SQRT(b2(l)*b2(l)/(4*b1(l)*b1(l)) - b3(l)/b1(l))

END DO

!----------------------------------------------------------------------
! Carry out calculations for points with open stomata
!----------------------------------------------------------------------
!DIR$ IVDEP
DO j=1,open_pts
  l = veg_index(open_index(j))

!----------------------------------------------------------------------
! Calculate the net rate of photosynthesis
!----------------------------------------------------------------------
  al(l) = (wl(l) - rd(l)) * fsmc(l)


!----------------------------------------------------------------------
! Calculate the factor for converting mol/m3 into Pa (J/m3).
!----------------------------------------------------------------------
  conv(l) = rmol * tl(l)

!----------------------------------------------------------------------
! Diagnose the leaf conductance
!----------------------------------------------------------------------
  glco2(l) = (al(l) * conv(l)) / (ca(l) - ci(l))
  gl(l) = ratio * glco2(l)

END DO

!----------------------------------------------------------------------
! Close stomata at points with negative or zero net photosynthesis
! or where the leaf resistance exceeds its maximum value.
!----------------------------------------------------------------------
!DIR$ IVDEP
DO j=1,open_pts
  l = veg_index(open_index(j))

  IF (gl(l) <= glmin(ft) .OR. al(l) <= 0.0) THEN
    gl(l) = glmin(ft)
    glco2(l) = gl(l) / ratio
    al(l) = -rd(l)*fsmc(l)
  END IF

END DO

IF ( l_o3_damage ) THEN
!-----------------------------------------------------------------------
! Modify the stomatal conductance and photosynthesis for ozone effects
! (Peter Cox, 12/11/04)
!-----------------------------------------------------------------------
  DO j=1,open_pts
    l = veg_index(open_index(j))

!-----------------------------------------------------------------------
! Flux of O3 without ozone effects (for use in analytical eqn)
!-----------------------------------------------------------------------
    flux_o3(l) = o3mol(l) / (ra(l) + (ratio_o3 / gl(l)))

!-----------------------------------------------------------------------
! Analytic solution for the ozone exposure factor
!-----------------------------------------------------------------------
! Use EPSILON to avoid overflow on division
    IF (ABS(ra(l)) < EPSILON(1.0)) THEN
      fo3(l) = (1 + dfp_dcuo(ft) * fl_o3_ct(ft))                  &
             / (1 + dfp_dcuo(ft) * flux_o3(l))
    ELSE
      b = ratio_o3 / (gl(l) * ra(l))                              &
        + dfp_dcuo(ft) * o3mol(l) / ra(l)                         &
        - (1 + dfp_dcuo(ft) * fl_o3_ct(ft))
      c = -ratio_o3 / (gl(l) * ra(l))                             &
        * (1 + dfp_dcuo(ft) * fl_o3_ct(ft))
      fo3(l) = -0.5 * b + 0.5 * SQRT(b * b - 4.0 * c)
    END IF

    fo3(l) = MIN(MAX(fo3(l),0.0),1.0)

!-----------------------------------------------------------------------
! Update the leaf conductance and photosynthesis
!-----------------------------------------------------------------------
    gl(l) = gl(l) * fo3(l)
    al(l) = al(l) * fo3(l)
  END DO

!----------------------------------------------------------------------
! Close stomata at points with negative or zero net photosynthesis
! or where the leaf resistance exceeds its maximum value.
!----------------------------------------------------------------------
!DIR$ IVDEP
  DO j=1,open_pts
    l = veg_index(open_index(j))

    IF (gl(l) <= glmin(ft) .OR. al(l) <= 0.0) THEN
      gl(l) = glmin(ft)
      glco2(l) = gl(l) / ratio
      al(l) = -rd(l) * fsmc(l) * fo3(l)
    END IF
  END DO

END IF ! o3 damage

!----------------------------------------------------------------------
! Define fluxes and conductances for points with closed stomata
!----------------------------------------------------------------------
!DIR$ IVDEP
DO j=1,clos_pts
  l = veg_index(clos_index(j))

  gl(l) = glmin(ft)
  glco2(l) = gl(l) / ratio
  al(l) = -rd(l) * fsmc(l)

! Add the contribution from ozone if enabled
  IF (l_o3_damage) al(l) = al(l) * fo3(l)
END DO

IF ( l_o3_damage ) THEN
!-----------------------------------------------------------------------
! Diagnose the ozone deposition flux on all points
!-----------------------------------------------------------------------
  DO j=1,veg_pts
    l = veg_index(j)
    flux_o3(l) = o3mol(l) / (ra(l) + (ratio_o3 / gl(l)))
    rd(l) = rd(l) * fo3(l)
  END DO
END IF

IF (lhook) CALL dr_hook('LEAF',zhook_out,zhook_handle)
RETURN
END SUBROUTINE leaf

