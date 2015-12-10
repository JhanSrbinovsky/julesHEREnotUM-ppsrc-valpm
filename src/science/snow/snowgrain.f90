! *****************************COPYRIGHT*******************************

! (c) [University of Edinburgh] [2009]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the JULES collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC237]

! *****************************COPYRIGHT*******************************
!  SUBROUTINE SNOWGRAIN------------------------------------------------

! Description:
!     Calculate growth of snow grains.

! Subroutine Interface:
SUBROUTINE snowgrain ( land_pts,tile_pts,timestep,nsnow,          &
                       tile_index,sice,snowfall,snowmass,tsnow,  &
                       tstar_tile,rgrain,rgrainl,rgrain0 )

USE ancil_info, ONLY :                                            &
!  imported scalars with intent(in)
 nsmax        ! Maximum possible number of snow layers

USE c_0_dg_c, ONLY :                                              &
!  imported scalar parameters
 tm           ! Temperature at which fresh water freezes
!                   ! and ice melts (K)

USE c_pi, ONLY :                                                  &
!  imported scalar parameters
 pi

USE c_rmol, ONLY :                                                &
 rmol         ! universal gas constant

USE rad_param, ONLY :                                             &
!  imported scalars with intent(in)
 r0                                                               &
              !  Grain size for fresh snow (microns)
,rmax                                                             &
              !  Maximum snow grain size (microns)
!  imported arrays with intent(in)
,snow_ggr     !  snow grain area growth rates (microns**2 s-1).


USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Scalar arguments with intent(in)
INTEGER, INTENT(IN) ::                                            &
 land_pts                                                         &
                  ! Total number of land points
,tile_pts         ! Number of tile points

REAL, INTENT(IN) ::                                               &
 timestep         ! Timestep (s)

! Array arguments with intent(in)
INTEGER, INTENT(IN) ::                                            &
 nsnow(land_pts)                                                  &
                  ! Number of snow layers
,tile_index(land_pts)
                  ! Index of tile points

REAL, INTENT(IN) ::                                               &
 sice(land_pts,nsmax)                                             &
                  ! Ice content of snow layers (kg/m2)
,snowfall(land_pts)                                               &
                  ! snowfall rate (kg m-2 s-1)
,snowmass(land_pts)                                               &
                  ! Snow mass on ground on tile (kg m-2)
,tsnow(land_pts,nsmax)                                            &
                  ! Snow layer temperatures (K)
,tstar_tile(land_pts)
                  ! Tile surface temperature (K)

! Array arguments with intent(inout)
REAL, INTENT(INOUT) ::                                            &
 rgrain(land_pts)                                                 &
                  ! Snow grain size (microns)
,rgrainl(land_pts,nsmax)
                  ! Snow grain size in snow layers (microns)

! Array arguments with intent(out)
REAL, INTENT(OUT) ::                                              &
 rgrain0(land_pts)
                  ! Fresh snow grain size (microns)

! Local scalars.
INTEGER ::                                                        &
 i                                                                &
                  ! Land point index
,k                                                                &
                  ! Tile point index
,n                ! Snow layer index

REAL :: rate      ! Grain area growth rate (microns2/s)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook('SNOWGRAIN',zhook_in,zhook_handle)

DO k=1,tile_pts
  i = tile_index(k)
!-----------------------------------------------------------------------
! Set grain size for fresh snow.
!-----------------------------------------------------------------------
  rgrain0(i) = r0

  IF ( nsnow(i) == 0 ) THEN
!-----------------------------------------------------------------------
! No snow, or zero-layer model selected.
!-----------------------------------------------------------------------
    IF ( snowmass(i) > 0.0 ) THEN
      rate = snow_ggr(1)
      IF ( tstar_tile(i) < tm ) THEN
        IF ( rgrain(i) < 150. ) THEN
          rate = snow_ggr(2)
        ELSE
          rate = snow_ggr(3) *                                    &
                    EXP( -3.7e4 / (rmol * tstar_tile(i)) )
        END IF
      END IF
      rgrain(i) = SQRT( rgrain(i)**2 + (rate / pi) * timestep )   &
                    - (rgrain(i) - r0) * snowfall(i)/2.5
      rgrain(i) = MIN( rmax, rgrain(i) )
      rgrain(i) = MAX( r0, rgrain(i) )
    ELSE
! No snow. Set grain size to that for fresh snow (ready for next occurence).
      rgrain(i) = r0
    END IF  !  SNOWMASS

  ELSE

!-----------------------------------------------------------------------
!     NSNOW>0: one or more snow layers.
!-----------------------------------------------------------------------
    DO n=1,nsnow(i)
      IF ( sice(i,n) > 0 ) THEN
        rate = snow_ggr(1)
        IF ( tsnow(i,n) < tm ) THEN
          IF ( rgrainl(i,n) < 150.0 ) THEN
            rate = snow_ggr(2)
          ELSE
            rate = snow_ggr(3) *                                  &
                      EXP( -3.7e4 / (rmol * tsnow(i,n)) )
          END IF
        END IF
        rgrainl(i,n) = SQRT(rgrainl(i,n)**2 + (rate/pi)*timestep)
        rgrainl(i,n) = MIN( rmax, rgrainl(i,n) )
        rgrainl(i,n) = MAX( r0, rgrainl(i,n) )
      END IF
    END DO

!-----------------------------------------------------------------------
! For all currently empty layers, set grain size to that for fresh snow.
!-----------------------------------------------------------------------
    IF ( nsnow(i) < nsmax )                                       &
      rgrainl(i,nsnow(i)+1:) = rgrain0(i)

  END IF  !  NSNOW

END DO  !  k (tile points)
IF (lhook) CALL dr_hook('SNOWGRAIN',zhook_out,zhook_handle)
RETURN

END SUBROUTINE snowgrain
