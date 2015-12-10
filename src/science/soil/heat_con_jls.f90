! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    SUBROUTINE HEAT_CON----------------------------------------------

! Description:
!    Calculates the soil thermal conductivity including the
!    effects of water and ice.

! Method == 1
!      Described in Cox et al (1999), Appendix B.
!      http://www.springerlink.com/content/9b459pyfhyjwk1ln/

! Method == 2
!      Simplified Johansen (1975).
!      See http://www-nwp/~frid/thermal_conductivity.pdf
!                                             (Dharssi, Jan 2008)

! Documentation : UM Documentation Paper 25

! Subroutine Interface:
SUBROUTINE heat_con(npnts,hcon,sthu,sthf,                         &
                    v_sat,hcons                                   &
)

USE switches, ONLY :                                              &
 soilhc_method
!        SOILHC_METHOD=1: Method of Cox et al (1999).
!        SOILHC_METHOD=2: Simplified Johansen (1975).

USE snow_param, ONLY :                                            &
!      imported scalars with intent(in)
   snow_hcon

USE soil_param, ONLY :                                            &
!      imported scalar parameters
   hcair,hcice,hcwat

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Subroutine arguments:
!   Scalar arguments with intent(IN) :
INTEGER                                                           &
 npnts              ! IN Number of gridpoints

REAL                                                              &
 hcon(npnts)                                                      &
                    ! IN Dry soil thermal conductivity (W/m/K).
,sthu(npnts)                                                      &
                    ! IN Fractional saturation of unfrozen water
!                         !    at layer boundaries.
,sthf(npnts)                                                      &
                    ! IN Fractional saturation of frozen water
!                         !    at layer boundaries.
,v_sat(npnts)       ! IN Volumetric soil moisture concentration
!                         !    at saturation (m3/m3 soil).

!   Array arguments with intent(OUT) :
REAL                                                              &
 hcons(npnts)       ! OUT The thermal conductivity between adjacen
!                         !     layers including effects of water and ic
!                         !     (W/m/K).

! Local scalars:
INTEGER                                                           &
 i                  ! WORK Loop counter.

! Local arrays:
REAL                                                              &
 hcsat(npnts)                                                     &
                    ! WORK The thermal conductivity of the
!                         !  saturated  soil at current ratio of ice to
!                         !      liquid water (W/m/K).
,sth(npnts)                                                       &
                    ! WORK Fractional saturation of water
!                         !     (liquid+ice) at layer boundaries.
,thice(npnts)                                                     &
                    ! WORK The concentration of ice at saturation
!                         !      for the current mass fraction of liquid
!                         !      water (m3 H2O/m3 soil).
,thwat(npnts)                                                     &
                    ! WORK The concentration of liquid water at
!                         !      saturation for the current mass
!                         !  fraction of liquid water  (m3 H2O/m3 soil).
,ke(npnts)          ! WORK Kersten number

!-----------------------------------------------------------------------
! Following local parameter values determined by linear regression,
! see Dharssi(2008).
!-----------------------------------------------------------------------
REAL                                                              &
 hcsat_max_allowed                                                &
,hcsat_min_allowed                                                &
,hcsat_gradient                                                   &
,hcsat_intercept
PARAMETER ( hcsat_max_allowed = 2.20, hcsat_min_allowed = 1.58    &
           ,hcsat_gradient    = 12.4, hcsat_intercept   = 0.25 )

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('HEAT_CON',zhook_in,zhook_handle)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)

!----------------------------------------------------------------------
! Initialise all points
!----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC) 
DO i=1,npnts
  IF (v_sat(i) >  0.0) THEN ! Soil points
    hcons(i)=hcon(i)
  ELSE ! Ice points
    hcons(i)=snow_hcon
  END IF
END DO
!$OMP END DO

IF (soilhc_method==1) THEN
!$OMP DO SCHEDULE(STATIC)
  DO i=1,npnts
!---------------------------------------------------------------
! Only do calculation for non land-ice pts
! V_SAT is set to zero for land-ice points
!---------------------------------------------------------------
    IF (v_sat(i) >  0.0) THEN

      IF (sthu(i) >  0.0) THEN
        thwat(i)=v_sat(i)*sthu(i)/(sthu(i)+sthf(i))
      ELSE
        thwat(i)=0.0
      END IF

      IF (sthf(i) >  0.0) THEN
        thice(i)=v_sat(i)*sthf(i)/(sthu(i)+sthf(i))
      ELSE
        thice(i)=0.0
      END IF

      sth(i)=sthu(i)+sthf(i)
      hcsat(i)=hcon(i)*(hcwat**thwat(i))*(hcice**thice(i))        &
                     /(hcair**v_sat(i))
      hcons(i)=(hcsat(i)-hcon(i))*sth(i)+hcon(i)
    END IF

  END DO
!$OMP END DO
END IF ! (SOILHC_METHOD==1)

IF (soilhc_method==2) THEN
!$OMP DO SCHEDULE(STATIC)
  DO i=1,npnts
    IF (v_sat(i) >  0.0) THEN

      IF (sthf(i) >  0.0) THEN
        thice(i)=v_sat(i)*sthf(i)/(sthu(i)+sthf(i))
      ELSE
        thice(i)=0.0
      END IF

      thwat(i)=v_sat(i)-thice(i)

      sth(i)=sthu(i)+sthf(i)

      hcsat(i)=hcsat_min_allowed                                  &
              +hcsat_gradient*(hcon(i)-hcsat_intercept)

      IF(hcsat(i) > hcsat_max_allowed) hcsat(i)=hcsat_max_allowed
      IF(hcsat(i) < hcsat_min_allowed) hcsat(i)=hcsat_min_allowed

      ! Adjust HCSAT for frozen soil water
      hcsat(i)=hcsat(i)*(hcwat**thwat(i))*(hcice**thice(i))       &
              /(hcwat**v_sat(i))

      IF(sth(i) <= 0.1) THEN
        ke(i)=0.0
      ELSE
        ke(i)=1.0 + LOG10(sth(i))
      END IF

      hcons(i)=(hcsat(i)-hcon(i))*ke(i)+hcon(i)
    END IF

  END DO
!$OMP END DO
END IF ! (SOILHC_METHOD==2)

!$OMP END PARALLEL 

IF (lhook) CALL dr_hook('HEAT_CON',zhook_out,zhook_handle)
RETURN
END SUBROUTINE heat_con
