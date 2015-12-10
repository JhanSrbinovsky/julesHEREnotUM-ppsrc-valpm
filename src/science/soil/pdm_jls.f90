! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    SUBROUTINE PDM----------------------------------------------------

! Description:
!     Calculates the saturation excess runoff using PDM.
!     See Moore, 1985


! Subroutine Interface:
SUBROUTINE pdm(                                                   &
               npnts,soil_pts,soil_index,nshyd,                   &
               tot_tfall,snow_melt,infiltro,timestep,             &
               v_sat,satexro,sthu,sthf)

USE switches, ONLY :                                              &
!      imported scalar parameters
 b_pdm,dz_pdm

USE c_densty, ONLY :                                              &
!      imported scalar parameters
 rho_water  !  density of pure water (kg/m3)

USE soil_param, ONLY :                                            &
 dzsoil         !  Thicknesses of the soil layers (m)

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Subroutine arguments
!   Scalar arguments with intent(IN) :

INTEGER                                                           &
 npnts                                                            &
                     ! IN Number of gridpoints.
,nshyd                                                            &
                     ! IN Number of soil moisture levels.
,soil_pts            ! IN Number of soil points.

REAL                                                              &
 timestep            ! IN Model timestep (s).

!   Array arguments with intent(IN) :

INTEGER                                                           &
 soil_index(npnts)   ! IN Array of soil points.

REAL                                                              &
 tot_tfall(npnts)                                                 &
                      ! IN Cumulative canopy throughfall
                      !     (kg/m2/s).
,snow_melt(npnts)                                                 &
                      ! IN GBM snow melt (kg/m2/s).
,infiltro(npnts)                                                  &
                      ! IN Surface runoff by infiltration
                      !    excess (kg/m2/s).
,v_sat(npnts,nshyd)                                               &
                      ! IN Volumetric soil moisture
!                           !    concentration at saturation
!                           !    (m^3 water/m^3 soil).
,sthu(npnts,nshyd)                                                &
                      ! IN Unfrozen soil moisture content of
!                           !    each layer as fraction of saturation.
,sthf(npnts,nshyd)                                                &
                      ! IN Frozen soil moisture content of
!                           !    each layer as fraction of saturation.

!   Array arguments with intent(OUT) :
,satexro(npnts)       ! OUT Sat. excess runoff (kg/m2/s).

! Local scalars:
INTEGER                                                           &
 i,j,n                ! WORK Loop counters.

REAL                                                              &
 p_in                                                             &
                 ! WORK Water reaching soil surface (m^3/m^3).
,thcap_max                                                        &
                 ! WORK Used in equation for TH_STAR (m^3/m^3).
,th_star                                                          &
                 ! WORK Critical moisture storage (m^3/m^3).
,th_star_p                                                        &
                 ! WORK Crit. TH_STAR after rain input (m^3/m^3).
,sth_pdm                                                          &
                 ! WORK Soil moisture for PDM / sat. value
,dz_pdm_t                                                         &
                 ! WORK Soil depth temporary
,dz_pdm_b        ! WORK Soil depth temporary

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('PDM',zhook_in,zhook_handle)

DO j=1,soil_pts
  i=soil_index(j)

  p_in = (tot_tfall(i)+snow_melt(i)-infiltro(i))                  &
    *timestep/(dz_pdm*rho_water)

  IF (p_in  >   0.0) THEN

! DZ_PDM_T is the residual soil depth for the PDM below the top
! of the current layer.

    dz_pdm_t = dz_pdm
    sth_pdm = 0.0

    DO n=1,nshyd

      dz_pdm_b = dz_pdm_t - dzsoil(n)

! DZ_PDM_B is the residual soil depth for the PDM below the
! bottom of the current layer

      IF (dz_pdm_b  >=  0.0) THEN
        sth_pdm = sth_pdm + (sthu(i,n)+sthf(i,n))*dzsoil(n)
      ELSE
        sth_pdm = sth_pdm+(sthu(i,n)+sthf(i,n))*MAX(0.,dz_pdm_t)
      END IF

      dz_pdm_t = dz_pdm_b

    END DO ! Loop over soil layers

    sth_pdm = sth_pdm/dz_pdm

    thcap_max = (b_pdm+1.) * v_sat(i,1)
    IF(sth_pdm  >   1.0) sth_pdm = 1.0
    th_star = thcap_max * (1.-(1.-sth_pdm)**(1/(b_pdm+1.)))
    th_star_p = th_star + p_in

    IF(th_star_p  <   thcap_max) THEN
      satexro(i) = p_in - v_sat(i,1)*                             &
        (1. - (1. - th_star_p/thcap_max)**(b_pdm+1.) - sth_pdm)
    ELSE
      satexro(i) = p_in - v_sat(i,1)*(1. - sth_pdm)
    END IF

    IF (satexro(i)  <   0.) satexro(i) = 0.

    satexro(i) = satexro(i) / timestep * (dz_pdm*rho_water)

  END IF

END DO


IF (lhook) CALL dr_hook('PDM',zhook_out,zhook_handle)
RETURN
END SUBROUTINE pdm
