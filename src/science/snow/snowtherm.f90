! *****************************COPYRIGHT*******************************

! (c) [University of Edinburgh] [2009]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the JULES collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC237]

! *****************************COPYRIGHT*******************************
!  SUBROUTINE SNOWTHERM------------------------------------------------

! Description:
!     Thermal properties of snow layers

! Subroutine Interface:
SUBROUTINE snowtherm ( land_pts,tile_pts,nsnow,tile_index,ds,     &
                       sice,sliq,csnow,ksnow )

USE ancil_info, ONLY :                                            &
!  imported scalars with intent(in)
 nsmax        !  Maximum number of snow layers

USE c_perma, ONLY :                                               &
!  imported scalar parameters
 hcapi                                                            &
              ! Specific heat capacity of ice (J/kg/K)
,hcapw        ! Specific heat capacity of water (J/kg/K)

USE c_rhowat, ONLY :                                              &
!  imported scalar parameters
 rho_water    ! density of water (kg/m3)


USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Scalar arguments with intent(in)
INTEGER,INTENT(IN) ::                                             &
 land_pts                                                         &
                     !  Number of land points
,tile_pts            !  Number of tile points

! Array arguments with intent(in)
INTEGER, INTENT(IN) ::                                            &
 nsnow(land_pts)                                                  &
                         !  Number of snow layers
,tile_index(land_pts)    !  Index of tile points

REAL, INTENT(IN) ::                                               &
 ds(land_pts,nsmax)                                               &
                         !  Snow layer thicknesses (m)
,sice(land_pts,nsmax)                                             &
                         !  Ice content of snow layers (kg/m2)
,sliq(land_pts,nsmax)    !  Liquid content of snow layers (kg/m2)

! Array arguments with intent(out)
REAL, INTENT(OUT) ::                                              &
 csnow(land_pts,nsmax)                                            &
                         !  Areal heat capacity of layers (J/K/m2)
,ksnow(land_pts,nsmax)   !  Thermal conductivity of layers (W/m/K)

! Local scalars
INTEGER ::                                                        &
 k                                                                &
                     !  Tile point index
,l                                                                &
                     !  Land point index
,n                   !  Snow layer index

REAL :: rho_snow     !  Snow layer density (kg/m3)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook('SNOWTHERM',zhook_in,zhook_handle)

DO k=1,tile_pts
  l = tile_index(k)
  DO n=1,nsnow(l)
    csnow(l,n) = sice(l,n)*hcapi + sliq(l,n)*hcapw
    rho_snow = (sice(l,n) + sliq(l,n)) / ds(l,n)
    ksnow(l,n) = 2.22*(rho_snow/rho_water)**1.88
  END DO
END DO
IF (lhook) CALL dr_hook('SNOWTHERM',zhook_out,zhook_handle)
RETURN

END SUBROUTINE snowtherm
