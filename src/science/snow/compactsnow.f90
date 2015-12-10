! *****************************COPYRIGHT*******************************

! (c) [University of Edinburgh] [2009]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the JULES collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC237]

! *****************************COPYRIGHT*******************************
!  SUBROUTINE COMPACTSNOW-----------------------------------------------

! Description:
!     Mechanical compaction of snow

! Subroutine Interface:
SUBROUTINE compactsnow ( land_pts,tile_pts,timestep,nsnow,        &
                         tile_index,sice,sliq,tsnow,rho_snow,ds )


USE ancil_info, ONLY :                                            &
!  imported scalars with intent(in)
 nsmax   !  Maximum possible number of snow layers

USE c_g, ONLY :                                                   &
!  imported scalar parameters
 g       !  mean acceleration due to gravity at
!              !  earth's surface (m s-2)

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Scalar arguments with intent(in)
INTEGER, INTENT(IN) ::                                            &
 land_pts                                                         &
                         ! Total number of land points
,tile_pts                ! Number of tile points

REAL, INTENT(IN) ::                                               &
 timestep                ! Timestep (s)

! Array arguments with intent(in)
INTEGER, INTENT(IN) ::                                            &
 nsnow(land_pts)                                                  &
                         ! Number of snow layers
,tile_index(land_pts)    ! Index of tile points

REAL, INTENT(IN) ::                                               &
 sice(land_pts,nsmax)                                             &
                         ! Ice content of snow layers (kg/m2)
,sliq(land_pts,nsmax)                                             &
                         ! Liquid content of snow layers (kg/m2)
,tsnow(land_pts,nsmax)   ! Snow layer temperatures (K)

! Array arguments with intent(inout)
REAL, INTENT(INOUT) ::                                            &
 rho_snow(land_pts,nsmax)  ! Snow layer densities (kg/m3)

! Array arguments with intent(out)
REAL, INTENT(OUT) ::                                              &
 ds(land_pts,nsmax)          ! Snow layer thicknesses (m)

! Local scalars
INTEGER ::                                                        &
 k                                                                &
                       ! Tile point index
,l                                                                &
                       ! Land point index
,n                     ! Snow layer index

REAL ::                                                           &
 mass                                                             &
                       ! Overlying mass of snow (kg/m2)
,rho                   ! Snow density (kg/m3)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!-------------------------------------------------------------------------------

IF (lhook) CALL dr_hook('COMPACTSNOW',zhook_in,zhook_handle)

DO k=1,tile_pts
  l = tile_index(k)

  mass = 0.
  ds(l,:) = 0.0

  DO n=1,nsnow(l)
    mass = mass + 0.5 * (sice(l,n) + sliq(l,n))
    rho = rho_snow(l,n)
    rho = rho + 0.5e-7 * rho * g * mass * timestep *              &
                EXP(14.643 - 4e3 / tsnow(l,n) - 0.02 * rho)
!     Note: mass (and hence rho) can be zero but nsnow>0 (likely 1!) if a very
!     shallow snowpack has been exhausted in this timestep.
    IF ( rho > EPSILON(rho) )                                     &
      ds(l,n) = (sice(l,n) + sliq(l,n)) / rho
    rho_snow(l,n) = rho
    mass = mass + 0.5 * (sice(l,n) + sliq(l,n))
  END DO
END DO
IF (lhook) CALL dr_hook('COMPACTSNOW',zhook_out,zhook_handle)
RETURN

END SUBROUTINE compactsnow
