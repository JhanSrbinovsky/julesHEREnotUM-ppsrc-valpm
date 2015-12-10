! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    SUBROUTINE CH4_WETL-----------------------------------------------

! Description:
!     Calculates methane emissions from wetland area.

! Subroutine Interface:
SUBROUTINE ch4_wetl(npnts,soil_pts,dim_cs1,soil_index,l_triffid   &
  ,tsoil_d,cs,f_wetl,fch4_wetl)

USE c_0_dg_c, ONLY :                                              &
!      imported scalar parameters
   tm

USE c_ch4, ONLY :                                                 &
!      imported scalar parameters
   const_ch4,q10_ch4,t0_ch4

USE surf_param, ONLY :                                            &
!      imported arrays with intent(in)
   kaps_roth

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

INTEGER, INTENT(IN) ::                                            &
 npnts                                                            &
                     ! IN Number of gridpoints.
,soil_pts                                                         &
                     ! IN Number of soil points.
,dim_cs1                                                          &
                     ! IN Number of soil carbon pools
,soil_index(npnts)   ! IN Array of soil points.

LOGICAL, INTENT(IN) ::                                            &
 l_triffid           ! TRUE if using TRIFFID

REAL, INTENT(IN) ::                                               &
 tsoil_d(npnts)                                                   &
                     ! IN Diagnosed soil temp to 1 metre (K).
,cs(npnts,dim_cs1)                                                &
                     ! IN Soil carbon (kg C/m2).
!                          !   For RothC (dim_cs1=4), the pools
!                          !    are DPM, RPM, biomass and humus.
,f_wetl(npnts)       ! IN Wetland fraction

REAL, INTENT(INOUT) ::                                            &
 fch4_wetl(npnts)    ! OUT Scaled methane flux (10^-9 kg C/m2/s)

REAL ::                                                           &
 q10t_ch4                                                         &
                     ! Q10 value at T
,const_tdep          ! T and Q10(0) dependent function

INTEGER ::                                                        &
 i,j,k

REAL ::                                                           &
 sumkaps             ! Sum of kaps_roth values.

!     Local aray variables.
REAL ::                                                           &
 cs_ch4(npnts)       ! Effective soil carbon (kg C/m2).

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------
!     Calculate an effective soil carbon for wetland methane emission.
!-----------------------------------------------------------------------
IF (lhook) CALL dr_hook('CH4_WETL',zhook_in,zhook_handle)
IF (l_triffid) THEN
  sumkaps = 0.0
  DO k=1,dim_cs1
    sumkaps = sumkaps + kaps_roth(k)
  END DO
!       Weight each pool by specific respiration rate.
  DO j=1,soil_pts
    i=soil_index(j)
    cs_ch4(i) = 0.0
    DO k=1,dim_cs1
      cs_ch4(i) = cs_ch4(i) + cs(i,k) * kaps_roth(k)
    END DO
  END DO
  cs_ch4(i) = cs_ch4(i) / sumkaps
ELSE
!       Use the single soil carbon pool.
  DO j=1,soil_pts
    i=soil_index(j)
    cs_ch4(i) = cs(i,1)
  END DO
END IF  !  l_triffid

!-----------------------------------------------------------------------
!     Calculate scaled wetland methane emission.
!-----------------------------------------------------------------------
const_tdep = t0_ch4 * LOG(q10_ch4)
DO j=1,soil_pts
  i=soil_index(j)
  IF( tsoil_d(i) > tm .AND. f_wetl(i) > 0.0 )THEN
    q10t_ch4 = EXP(const_tdep/tsoil_d(i))
    fch4_wetl(i) = 1.e9*const_ch4*cs_ch4(i)*                     &
                   f_wetl(i)*q10t_ch4**(0.1*(tsoil_d(i)-t0_ch4))
  END IF
END DO

IF (lhook) CALL dr_hook('CH4_WETL',zhook_out,zhook_handle)
RETURN
END SUBROUTINE ch4_wetl
