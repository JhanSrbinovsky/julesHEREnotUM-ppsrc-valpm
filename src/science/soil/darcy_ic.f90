! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    SUBROUTINE DARCY_IC-----------------------------------------------

! Description:
!     Intermediate control routine to choose between the two
!     versions of the soil dynamics
!      i.e. Clapp and Hornburger or Van Genuchten

! Documentation : UM Documentation Paper 25

! Subroutine Interface:
SUBROUTINE darcy_ic( npnts,soil_pts,soil_index,b,ks,sathh         &
,                    sthu1,dz1,sthu2,dz2,wflux                    &
,                    dwflux_dsthu1,dwflux_dsthu2,ltimer )

USE switches, ONLY :                                              &
 l_vg_soil

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

INTEGER, INTENT(IN) ::                                            &
                      !  IN scalars
 npnts                                                            &
                      ! IN Number of gridpoints.
,soil_pts             ! IN Number of soil points.

INTEGER, INTENT(IN) ::                                            &
                      !  IN arrays
 soil_index(npnts)    ! IN Array of soil points.

REAL, INTENT(IN) ::                                               &
                      !  IN arrays
 b(npnts,2)                                                       &
                      ! IN Clapp-Hornberger exponent.
,dz1                                                              &
                      ! IN Thickness of the upper layer (m).
,dz2                                                              &
                      ! IN Thickness of the lower layer (m).
,ks(npnts)                                                        &
                      ! IN Saturated hydraulic conductivity
!                           !    (kg/m2/s).
,sathh(npnts,2)                                                   &
                      ! IN Saturated soil water pressure (m).
,sthu1(npnts)                                                     &
                      ! IN Unfrozen soil moisture content of upper
!                           !    layer as a fraction of saturation
,sthu2(npnts)         ! IN Unfrozen soil moisture content of lower
!                           !    layer as a fraction of saturation.

LOGICAL ltimer        ! Logical switch for TIMER diags

REAL, INTENT(OUT) ::                                              &
                      !   OUT arrays
 wflux(npnts)                                                     &
                      ! OUT The flux of water between layers
!                           !     (kg/m2/s).
,dwflux_dsthu1(npnts)                                             &
                      ! OUT The rate of change of the explicit
!                           !     flux with STHU1 (kg/m2/s).
,dwflux_dsthu2(npnts) ! OUT The rate of change of the explicit
!                           !     flux with STHU2 (kg/m2/s).

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook('DARCY_IC',zhook_in,zhook_handle)

IF ( l_vg_soil ) THEN
!       van Genuchten model
! DEPENDS ON: darcy_vg
  CALL darcy_vg ( npnts,soil_pts,soil_index,b,ks,sathh            &
,                 sthu1,dz1,sthu2,dz2,wflux                       &
,                 dwflux_dsthu1,dwflux_dsthu2,ltimer )
ELSE
!       Clapp and Hornberger model
! DEPENDS ON: darcy_ch
  CALL darcy_ch ( npnts,soil_pts,soil_index,b,ks,sathh            &
,                 sthu1,dz1,sthu2,dz2,wflux                       &
,                 dwflux_dsthu1,dwflux_dsthu2,ltimer )
END IF

IF (lhook) CALL dr_hook('DARCY_IC',zhook_out,zhook_handle)
RETURN
END SUBROUTINE darcy_ic
