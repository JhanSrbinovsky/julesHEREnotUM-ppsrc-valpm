
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Purpose:
! Calculates the anthropogenic contribution to surface heat fluxes for
! urban tiles by linear interpolation of monthly values. The value is
! then passed in anthrop_heat(n), that has a value of 0.0 except when
! n=6 (urban) and l_anthrop_heat=.true., and added to the surface
! heat fluxes in sf_expl, sf_impl and sf_impl_2.

! Original code from Martin Best and Peter Clark (December 2005).
! Updated for UM7.1 by Jorge Bornemann (May 2008)


SUBROUTINE generate_anthropogenic_heat(                           &
 val_year, val_day_number, val_hour, val_minute, val_second       &
, ntiles, land_points, frac, l_anthrop_heat_src                   &
 )

USE nstypes, ONLY : urban, ntype, urban_canyon, urban_roof

USE urban_param, ONLY : wrr, anthrop_heat_scale

USE switches_urban, ONLY : l_urban2t

USE switches, ONLY : l_aggregate, l_360

USE fluxes, ONLY : anthrop_heat

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook


USE PrintStatus_mod

IMPLICIT NONE

! IN time information for current timestep
INTEGER, INTENT(IN) ::                                            &
  val_year                                                        &
, val_day_number                                                  &
, val_hour                                                        &
, val_minute                                                      &
, val_second

! IN Number of tiles
INTEGER, INTENT(IN) :: &
   ntiles,             &
   land_points           ! No.of land points being processed, can be 0.

REAL, INTENT(IN)    :: &
   frac(land_points,ntype)            ! IN Fractions of surface types.

! IN Switch for  anthropogenic heat source
LOGICAL, INTENT(IN) :: l_anthrop_heat_src

! OUT

! Local Variables
INTEGER :: &
   tile_index(land_points,ntype),   & ! Index of tile points
   tile_pts(ntype)                    ! Number of tile points

REAL                                                                    &
   urban_month(12),                                                     &
                     ! Monthly values of anthropogenic
                     ! contribution to surface heat fluxes W/m2
   anthrop_heat_urban(land_points),                                     &
               ! For checking: Copy of anthrop_heat(:,urban)
   urban_agg,                                                           &
               ! For checking: anthrop_heat aggregated for urban tiles
   mm, dpm
INTEGER i,im,im1,n,k,l

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


! urban anthropogenic heat source is taken from the digest of energy
! statistics (1995 - 2003) monthly averaged to 9 years, converted
! to w/m2 and adjusted to fraction dissipated in urban areas

DATA urban_month /                                                &
  25.1447                                                         &
, 23.4529                                                         &
, 24.6175                                                         &
, 20.7802                                                         &
, 18.8677                                                         &
, 18.1567                                                         &
, 17.1881                                                         &
, 16.6732                                                         &
, 18.5490                                                         &
, 20.7435                                                         &
, 22.9870                                                         &
, 26.1500 /

IF (lhook) CALL dr_hook('GENERATE_ANTHROPOGENIC_HEAT',zhook_in,zhook_handle)

! Initialise arrays
anthrop_heat_urban(:) = 0.0

!! tile_pts isn't set until after ni_bl_ctl so need to call it here.
! DEPENDS ON: tilepts
CALL tilepts( land_points, frac, tile_pts, tile_index )

! Check whether anthropogenic heat is switched and make sure that
! we have an urban tile
IF ( .NOT. l_aggregate .AND.                                      &
   l_anthrop_heat_src .AND. ( urban > 0 .OR. urban_canyon > 0 ) ) THEN

! Check which urban tile is present to make sure that the following is
! excecuted regardless of urban scheme. "urban" may not be present if URBAN-2T
! i.e. "urban_canyon" is named instead
  IF ( urban_canyon > 0 .AND. urban < 0 ) urban = urban_canyon

! Make sure we have the correct average days per month depending on
! whether we are running with a 360 day year
  IF ( l_360 ) THEN
    dpm = 360.0 / 12.0
  ELSE
    dpm = 365.0 / 12.0
  END IF

  mm = val_day_number / dpm - 0.5

  im  = INT(mm)
  mm  = mm - im
  im  = im + 1
  im1 = im + 1

  IF (im == 0) THEN
    im  = 12
  ELSE IF (im1 == 13) THEN
    im1 = 1
  END IF

  anthrop_heat(:,urban)= mm * urban_month(im1)                    &
                      + (1.0 - mm) * urban_month(im)
  anthrop_heat_urban(:) = anthrop_heat(:,urban)

! URBAN-2T: Over-write anthrop_heat for roof and canyon doing roof at same
! time as canyon to avoid over-writting land-ice. Capability is now here to
! do spatially varying anthrop_heat.
  IF ( l_urban2t ) THEN
    n = urban_canyon
    DO k=1,tile_pts(n)
      l = tile_index(k,n)
! Here anthrop_heat(urban_roof) is a fraction (anthrop_heat_scale) of
! anthrop_heat(urban_canyon) where wrr is the canyon fraction. The total from
! the urban tile is conserved.
      anthrop_heat(l,n) = anthrop_heat(l,urban) /                 &
         ( anthrop_heat_scale * ( 1.0 - wrr(l) ) + wrr (l) )
      anthrop_heat(l,urban_roof) =                                &
         anthrop_heat_scale * anthrop_heat(l,n)
    END DO
  END IF

END IF

! Check that anthrop_heat is conserved
IF ( l_urban2t .AND. .NOT. l_aggregate ) THEN
  n = urban_canyon
  DO k = 1, tile_pts(n)
    l = tile_index(k,n)
    urban_agg = anthrop_heat(l,urban_roof) * ( 1.0 - wrr(l) ) +         &
       anthrop_heat(l,urban_canyon) * wrr(l)
    IF ( ABS( urban_agg - anthrop_heat_urban(l) ) > 1E-6 ) THEN

      IF ( printstatus > PrStatus_Normal ) THEN
        WRITE(6,*) 'WARNING: anthrop_heat not conserved on urban tiles for',&
           k, l
        WRITE(6,*) 'Aggregated', urban_agg, 'Urban', anthrop_heat_urban(l)
      END IF





    END IF
  END DO
END IF

IF (lhook) CALL dr_hook('GENERATE_ANTHROPOGENIC_HEAT',zhook_out,zhook_handle)
RETURN
END SUBROUTINE generate_anthropogenic_heat

