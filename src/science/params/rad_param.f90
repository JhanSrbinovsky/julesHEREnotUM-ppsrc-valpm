! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module sets parameter values for radiation routines.
  
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.


MODULE rad_param


  IMPLICIT NONE

!-----------------------------------------------------------------------
!  Parameters for snow grain size.
!-----------------------------------------------------------------------
! DEFAULTS TAKEN FROM JULES EXAMPLE CONTROL FILES
  REAL ::                                                           &
   r0 = 50.                                                         &
                  !  Grain size for fresh snow (microns)
  ,rmax = 2000.
                  !  Maximum snow grain size (microns).

!-----------------------------------------------------------------------
! Parameters for snow grain size growth (HCTN30.16).
! Values are for melting snow, cold fresh snow and cold aged snow
! respectively.
!-----------------------------------------------------------------------
  REAL ::                                                           &
   snow_ggr(3)    !  snow grain area growth
!                    rates (microns**2 s-1).
! DEFAULTS TAKEN FROM JULES EXAMPLE CONTROL FILES
  DATA snow_ggr / 0.6, 0.06, 0.23e6 /

!-----------------------------------------------------------------------
! Parameters for prognostic, spectral snow albedo.
!-----------------------------------------------------------------------
  REAL ::                                                           &
   amax(2)        !  Maximum albedo for fresh snow
!                    (values for VIS and NIR)
! DEFAULTS TAKEN FROM JULES EXAMPLE CONTROL FILES
  DATA amax / 0.98, 0.7 /

!-----------------------------------------------------------------------
!  Parameters for ( diagnostic) all-band snow albedo.
  REAL :: dtland = 2.
                  !  degrees Celsius below freezing point at which
!                    snow albedo equals cold deep snow albedo. This
!                    is 2 in HCTN30.4. Must not be zero!

  REAL :: kland_numerator = 0.3
                  !  KLAND is calculated as KLAND_NUMERATOR / DTLAND
!                    in jules_init once the namelist has been read

  REAL :: kland
                  !  used in snow-ageing effect on snow albedo.
!                    This is 0.3 in HCTN30.4, although
!                    note that the last term in that eqn should be
!                    divided by dtland.

  REAL :: maskd = 50.
                  !  used in the exponent of equation weighting
!                    snow and snow-free albedo to get tile
!                    albedo.
!          This is 0.2 in HCTN30.5, where it is used as maskd*snowMass,
!          assuming snow density=250kg/m3.
!          It is now used as maskd*snowDepth, so maskd=50 gives the same
!          relationship as HCTN30.5.

  REAL :: tcland
                  !  temperature below which snow albedo equals
!                    cold deep snow albedo (HCTN30.4)
!                    This is initialised as TM - DTLAND in
!                    jules_init once the namelist has been read
!-----------------------------------------------------------------------

  NAMELIST /jules_rad_param/ r0,rmax,snow_ggr,amax,maskd,dtland,    &
                             kland_numerator

END MODULE rad_param
