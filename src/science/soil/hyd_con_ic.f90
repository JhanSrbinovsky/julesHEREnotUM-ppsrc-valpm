! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    SUBROUTINE HYD_CON_IC---------------------------------------------

! Description:
!     Intermediate control routine to choose between the two
!     versions of the soil dynamics
!      i.e. Clapp and Hornburger or Van Genuchten

! Documentation : UM Documentation Paper 25

! Subroutine Interface:
SUBROUTINE hyd_con_ic( npnts,soil_pts,soil_index,b,ks,thetak      &
,                      k,dk_dthk,ltimer )

USE switches, ONLY :                                              &
 l_vg_soil

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

INTEGER, INTENT(IN) ::                                            &
                  !  IN scalars
 npnts                                                            &
                  ! IN points in grid
,soil_pts         ! IN Number of soil points.

INTEGER, INTENT(IN) ::                                            &
                  !  IN arrays
 soil_index(npnts)!  IN Array of soil points.

REAL, INTENT(IN) ::                                               &
                  !  IN arrays
 b(npnts)                                                         &
                  ! IN Exponent in conductivity and soil water
!                       !    suction fits.
,ks(npnts)                                                        &
                  ! IN The saturated hydraulic conductivity (kg/m2/s)
,thetak(npnts)    ! IN Fractional saturation.

LOGICAL ltimer        ! Logical switch for TIMER diags

REAL, INTENT(OUT) ::                                              &
                  !  OUT arrays
 k(npnts)                                                         &
                  ! OUT The hydraulic conductivity (kg/m2/s).
,dk_dthk(npnts)   ! OUT The rate of change of K with THETAK
!                       !     (kg/m2/s).

 INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
 INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
 REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook('HYD_CON_IC',zhook_in,zhook_handle)

IF ( l_vg_soil ) THEN
!       van Genuchten model
! DEPENDS ON: hyd_con_vg
  CALL hyd_con_vg( npnts,soil_pts,soil_index,b,ks,thetak,k        &
,                  dk_dthk,ltimer )
ELSE
!       Clapp and Hornberger model
! DEPENDS ON: hyd_con_ch
  CALL hyd_con_ch( npnts,soil_pts,soil_index,b,ks,thetak,k        &
,                  dk_dthk,ltimer )
END IF

IF (lhook) CALL dr_hook('HYD_CON_IC',zhook_out,zhook_handle)
RETURN
END SUBROUTINE hyd_con_ic
