! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Description:
! Calculates the soil respiration based on a simplified version of the
! model of Raich et al. (1991).

!***********************************************************************
SUBROUTINE microbe (land_pts,dim_cs1,dim_cs2,l_triffid,l_q10,cs,  &
                    sth_soil,v_sat,v_wilt,tsoil,resp_s,veg_frac)

USE surf_param, ONLY :                                            &
 kaps,kaps_roth,q10 => q10_soil

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

INTEGER                                                           &
 land_pts                                                         &
                            ! IN Number of land points to be
!                                 !    processed.
,dim_cs1, dim_cs2           ! IN soil carbon dimensions

LOGICAL                                                           &
        l_triffid                                                 &
                            ! TRUE if using TRIFFID
,       l_q10               ! TRUE if using Q10 for soil resp

REAL                                                              &
 cs(land_pts,dim_cs1)                                             &
                            ! IN Soil carbon (kg C/m2).
!                                 !    For RothC (dim_cs1=4), the pools
!                                 !    are DPM, RPM, biomass and humus.
,veg_frac(dim_cs2)                                                &
                            ! IN vegetated fraction.
,sth_soil(land_pts)                                               &
                            ! IN Top layer soil moisture as a
!                                 !    fraction of saturation (m3/m3).
,v_sat(land_pts)                                                  &
                            ! IN Volumetric soil moisture
!                                 !    concentration at saturation
!                                 !    (m3 H2O/m3 soil).
,v_wilt(land_pts)                                                 &
                            ! IN Volumetric soil moisture
!                                 !    concentration below which
!                                 !    stomata close (m3 H2O/m3 soil).
                            !    as a fraction of saturation.
,tsoil(land_pts)                                                  &
                            ! IN Soil temperature (K).
,resp_s(land_pts,dim_cs1)                                         &
                            ! OUT Soil respiration (kg C/m2/s).
,fsth(land_pts)                                                   &
                            ! WORK Factors describing the
,fprf(land_pts)                                                   &
                            !      influence of soil moisture,
,ftemp(land_pts)                                                  &
                            !      vegetation cover and
,sth_resp_min                                                     &
                            ! WORK soil moist giving min. resp
!                                 !      soil temperature respectively
!                                 !      on the soil respiration.
,sth_opt                                                          &
                            ! WORK Fractional soil moisture at
!                                 !      which respiration is maximum.
,sth_wilt                   ! WORK Wilting soil moisture as a
!                                 !      fraction of saturation.
INTEGER ::                                                        &
 l,n                        ! Loop counters

!-----------------------------------------------------------------------
! Local parameters
!-----------------------------------------------------------------------
REAL ::                                                           &
 min_factor = 1.7          ! FACTOR to scale WILT to get RESP_MIN
!                                ! at 25 deg C and optimum soil
!                                ! moisture (/s).

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('MICROBE',zhook_in,zhook_handle)
IF (.NOT. l_triffid) THEN
  l_q10 = .TRUE.
  min_factor = 1.0
END IF

! FSTH
DO l=1, land_pts
  fsth(l)=0.0
END DO
DO l=1,land_pts
  IF (v_sat(l)  >   0.0) THEN

    sth_wilt = v_wilt(l) / v_sat(l)
    sth_opt = 0.5 * (1 + sth_wilt)
    sth_resp_min = sth_wilt * min_factor

    fsth(l) = 0.2
    IF (sth_soil(l)  >   sth_resp_min .AND.                       &
            sth_soil(l)  <=  sth_opt) THEN
      fsth(l) = 0.2 + 0.8 * ((sth_soil(l) - sth_resp_min)         &
                          /  (sth_opt - sth_resp_min))
    ELSE IF (sth_soil(l)  >   sth_opt) THEN
      fsth(l) = 1 - 0.8 * (sth_soil(l) - sth_opt)
    END IF
  END IF
END DO

! FTEMP
IF (l_q10) THEN
! use original HadCM3LC Q10 formula
  DO l=1,land_pts
    ftemp(l) = q10 ** (0.1 * (tsoil(l) - 298.15))
  END DO
ELSE
! use RothC temperature formula (using TSOIL for now...)
  DO l=1,land_pts
    ftemp(l)=0.0
    IF (tsoil(l)  >   263.15)                                     &
        ftemp(l) = 47.9 / (1.0+EXP(106.0/(tsoil(l)-254.85)))
  END DO
END IF

! FPRF - plant retainment factor
!          =0.6 for fully vegetated
!          =1.0 for bare soil
!        only set for RothC runs, (i.e. L_TRIFFID=TRUE)

IF (l_triffid) THEN
  DO l=1,land_pts
    fprf(l) = 0.6 + 0.4*(1-veg_frac(l))
  END DO
END IF

! set 1-D or 4-D soil resp depending on whether using RothC or not
IF (l_triffid) THEN
  DO n=1,dim_cs1
    DO l=1,land_pts
      resp_s(l,n) = kaps_roth(n)*cs(l,n)*fsth(l)*ftemp(l)*fprf(l)
    END DO
  END DO
ELSE
  DO l=1,land_pts
    resp_s(l,1) = kaps * cs(l,1) * fsth(l) * ftemp(l)
  END DO
END IF

IF (lhook) CALL dr_hook('MICROBE',zhook_out,zhook_handle)
RETURN

END SUBROUTINE microbe
