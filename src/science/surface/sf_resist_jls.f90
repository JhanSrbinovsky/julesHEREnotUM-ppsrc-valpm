! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  L  SUBROUTINE SF_RESIST----------------------------------------------

!    Purpose: Calculate surface moisture flux resistance factors.

!    Arguments --------------------------------------------------------
SUBROUTINE sf_resist (                                            &
 land_pts,tile_pts,land_index,tile_index,                         &
 canopy,catch,ch,dq,epdt,flake,gc,snow_tile,vshr,                 &
 fraca,resfs,resft                                                &
 )

USE atm_fields_bounds_mod
USE theta_field_sizes, ONLY : t_i_length

USE switches, ONLY :                                              &
 frac_snow_subl_melt

USE rad_param, ONLY :                                             &
 maskd

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

INTEGER                                                           &
 land_pts                                                         &
                       ! IN Total number of land points.
,tile_pts                                                         &
                       ! IN Number of tile points.
,land_index(land_pts)                                             &
                       ! IN Index of land points.
,tile_index(land_pts)  ! IN Index of tile points.

REAL                                                              &
 canopy(land_pts)                                                 &
                     ! IN Surface water (kg per sq metre).  F642.
,catch(land_pts)                                                  &
                     ! IN Surface capacity (max. surface water)
!                          !    (kg per sq metre).  F6416.
,ch(land_pts)                                                     &
                     ! IN Transport coefficient for heat and
!                          !    moisture transport
,dq(land_pts)                                                     &
                     ! IN Sp humidity difference between surface
!                          !    and lowest atmospheric level (Q1 - Q*).
,epdt(land_pts)                                                   &
                     ! IN "Potential" Evaporation * Timestep.
!                          !    Dummy variable for first call to routine
,flake(land_pts)                                                  &
                     ! IN Lake fraction.
,gc(land_pts)                                                     &
                     ! IN Interactive canopy conductance
!                          !    to evaporation (m/s)
,snow_tile(land_pts)                                              &
                     ! IN Lying snow amount (kg per sq metre).
,vshr(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                     ! IN Magnitude of surface-to-lowest-level
!                          !     windshear

REAL                                                              &
 fraca(land_pts)                                                  &
                     ! OUT Fraction of surface moisture flux with
!                          !     only aerodynamic resistance.
,resfs(land_pts)                                                  &
                     ! OUT Combined soil, stomatal and aerodynamic
!                          !     resistance factor for fraction 1-FRACA.
,resft(land_pts)     ! OUT Total resistance factor
!                          !     FRACA+(1-FRACA)*RESFS.

! Workspace -----------------------------------------------------------
INTEGER                                                           &
 i,j                                                              &
             ! Horizontal field index.
,k                                                                &
             ! Tile field index.
,l           ! Land field index.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('SF_RESIST',zhook_in,zhook_handle)

!-----------------------------------------------------------------------
!     Evaporation over land surfaces without snow is limited by
!     soil moisture availability and stomatal resistance.
!     Set FRACA (= fA in the documentation) according to P243.68,
!     and RESFS (= fS) according to P243.75 and P243.61.
!-----------------------------------------------------------------------
!CDIR NODEP
DO k=1,tile_pts
  l = tile_index(k)
  j=(land_index(l)-1)/t_i_length + 1
  i = land_index(l) - (j-1)*t_i_length

!-----------------------------------------------------------------------
! Calculate the fraction of the flux with only aerodynamic resistance
! (canopy evaporation).
! Set to 1 for negative moisture flux or snow-covered land
! (no surface/stomatal resistance to condensation).
!-----------------------------------------------------------------------
  fraca(l) = 1.0
  IF (dq(l) <  0. .AND. snow_tile(l) <= 0.) fraca(l) = 0.0
  IF (dq(l) <  0. .AND. snow_tile(l) <= 0. .AND. catch(l) >  0.)  &
    fraca(l) = canopy(l) / ( epdt(l) + catch(l) )
  IF (snow_tile(l) > 0. .AND. frac_snow_subl_melt == 1) THEN
    fraca(l) = 1.0 - EXP(-maskd*snow_tile(l))
  END IF
  fraca(l) = MIN(fraca(l),1.0)

!-----------------------------------------------------------------------
! Calculate resistance factors for transpiration from vegetation tiles
! and bare soil evaporation from soil tiles.
!-----------------------------------------------------------------------
  resfs(l) = gc(l) / ( gc(l) + ch(l)*vshr(i,j) )
  resft(l) = flake(l) + (1. - flake(l)) *                         &
                        ( fraca(l) + (1. - fraca(l))*resfs(l) )

END DO

IF (lhook) CALL dr_hook('SF_RESIST',zhook_out,zhook_handle)
RETURN
END SUBROUTINE sf_resist
