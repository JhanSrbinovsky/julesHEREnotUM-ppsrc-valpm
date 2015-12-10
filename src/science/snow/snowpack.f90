! *****************************COPYRIGHT*******************************

! (c) [University of Edinburgh] [2009]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the JULES collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC237]

! *****************************COPYRIGHT*******************************
!  SUBROUTINE SNOWPACK-------------------------------------------------

! Description:
!     Snow thermodynamics and hydrology

! Subroutine Interface:
SUBROUTINE snowpack ( land_pts,tile_pts,timestep,cansnowtile,     &
                      nsnow,tile_index,csnow,ei_tile,hcaps,hcons, &
                      ksnow,rho_snow_grnd,smcl1,snowfall,sthf1,   &
                      surf_htf_tile,tile_frac,v_sat1,ds,          &
                      melt_tile,sice,sliq,snomlt_sub_htf,         &
                      snowdepth,snowmass,tsnow,tsoil1,            &
                      snow_soil_htf,rho_snow,rho0,sice0,tsnow0 )

USE c_0_dg_c, ONLY :                                              &
!  imported scalar parameters
 tm              ! Temperature at which fresh water freezes
!                      ! and ice melts (K)

USE c_densty, ONLY :                                              &
!  imported scalar parameters
 rho_water       ! density of pure water (kg/m3)

USE c_lheat , ONLY :                                              &
!  imported scalar parameters
 lc                                                               &
                 ! latent heat of condensation of water
!                      ! at 0degC (J kg-1)
,lf              !  latent heat of fusion at 0degC (J kg-1)

USE c_perma, ONLY :                                               &
!  imported scalar parameters
 hcapi                                                            &
                 !  Specific heat capacity of ice (J/kg/K)
,hcapw           !  Specific heat capacity of water (J/kg/K)

USE soil_param, ONLY :                                            &
!  imported arrays with intent(in)
 dzsoil          !  Thicknesses of the soil layers (m)

USE ancil_info, ONLY :                                            &
!  imported scalars with intent(in)
 nsmax           ! Maximum possible number of snow layers

USE snow_param, ONLY :                                            &
 rho_snow_const                                                   &
                 !   constant density of lying snow (kg per m**3)
,rho_snow_fresh                                                   &
                 !   density of fresh snow (kg per m**3)
,snowliqcap      !   Liquid water holding capacity of lying snow
!                          as a fraction of snow mass.

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Scalar parameters
REAL, PARAMETER :: GAMMA = 0.5
                 ! Implicit timestep weighting

! Scalar arguments with intent(in)
INTEGER, INTENT(IN) ::                                            &
 land_pts                                                         &
                 ! Total number of land points
,tile_pts        ! Number of tile points

REAL, INTENT(IN) ::                                               &
 timestep        ! Timestep (s)

LOGICAL, INTENT(IN) ::                                            &
 cansnowtile     ! Switch for canopy snow model

! Array arguments with intent(in)
INTEGER, INTENT(IN) ::                                            &
 nsnow(land_pts)                                                  &
                 ! Number of snow layers
,tile_index(land_pts)
                 ! Index of tile points

REAL, INTENT(IN) ::                                               &
 csnow(land_pts,nsmax)                                            &
                 ! Areal heat capacity of layers (J/K/m2)
,ei_tile(land_pts)                                                &
                 ! Sublimation of snow (kg/m2/s)
,hcaps(land_pts)                                                  &
                 ! Heat capacity of soil surface layer (J/K/m3)
,hcons(land_pts)                                                  &
                 ! Thermal conductivity of top soil layer,
!                      ! including water and ice (W/m/K)
,ksnow(land_pts,nsmax)                                            &
                 ! Thermal conductivity of layers (W/m/K)
,rho_snow_grnd(land_pts)                                          &
                 ! Snowpack bulk density (kg/m3)
,smcl1(land_pts)                                                  &
                 ! Moisture content of surface soil layer (kg/m2)
,snowfall(land_pts)                                               &
                 ! Snow reaching the ground (kg/m2)
,sthf1(land_pts)                                                  &
                 ! Frozen soil moisture content of surface layer
!                      ! as a fraction of saturation.
,surf_htf_tile(land_pts)                                          &
                 ! Snow surface heat flux (W/m2)
,tile_frac(land_pts)                                              &
                 ! Tile fractions
,v_sat1(land_pts)
                 ! Surface soil layer volumetric
!                      ! moisture concentration at saturation

! Array arguments with intent(inout)
REAL, INTENT(INOUT) ::                                            &
 ds(land_pts,nsmax)                                               &
                 ! Snow layer depths (m)
,melt_tile(land_pts)                                              &
                 ! Surface snowmelt rate (kg/m2/s)
,sice(land_pts,nsmax)                                             &
                 ! Ice content of snow layers (kg/m2)
,sliq(land_pts,nsmax)                                             &
                 ! Liquid content of snow layers (kg/m2)
,snomlt_sub_htf(land_pts)                                         &
                 ! Sub-canopy snowmelt heat flux (W/m2)
,snowdepth(land_pts)                                              &
                 ! Snow depth (m)
,snowmass(land_pts)                                               &
                 ! Snow mass on the ground (kg/m2)
,tsnow(land_pts,nsmax)                                            &
                 ! Snow layer temperatures (K)
,tsoil1(land_pts)
                 ! Soil surface layer temperature(K)

! Array arguments with intent(out)
REAL, INTENT(OUT) ::                                              &
 snow_soil_htf(land_pts)                                          &
                 ! Heat flux into the uppermost subsurface layer
!                      !        (W/m2)
!                      ! i.e. snow to ground, or into snow/soil
!                      !      composite layer
,rho_snow(land_pts,nsmax)                                         &
                 ! Density of snow layers (kg/m3)
,rho0(land_pts)                                                   &
                 ! Density of fresh snow (kg/m3)
!                         ! Where NSNOW=0, rho0 is the density
!                         ! of the snowpack.
,sice0(land_pts)                                                  &
                 ! Ice content of fresh snow (kg/m2)
!                      ! Where NSNOW=0, SICE0 is the mass of
!                      ! the snowpack.
,tsnow0(land_pts)
                 ! Temperature of fresh snow (K)

! Local scalars
INTEGER ::                                                        &
 i                                                                &
                 ! land point index
,k                                                                &
                 ! Tile point index
,n               ! Snow layer index

REAL ::                                                           &
 asoil                                                            &
                 ! 1 / (dz*hcap) for surface soil layer
,can_melt                                                         &
                 ! Melt of snow on the canopy (kg/m2/s)
,coldsnow                                                         &
                 ! layer cold content (J/m2)
,dsice                                                            &
                 ! Change in layer ice content (kg/m2)
,g_snow_surf                                                      &
                 ! Heat flux at the snow surface (W/m2)
,sliqmax                                                          &
                 ! Maximum liquid content for layer (kg/m2)
,submelt                                                          &
                 ! Melt of snow beneath canopy (kg/m2/s)
,smclf                                                            &
                 ! Frozen soil moisture content of
!                      ! surface layer (kg/m2)
,win             ! Water entering layer (kg/m2)

! Local arrays
REAL ::                                                           &
 asnow(nsmax)                                                     &
                 ! Effective thermal conductivity (W/m2/k)
,a(nsmax)                                                         &
                 ! Below-diagonal matrix elements
,b(nsmax)                                                         &
                 ! Diagonal matrix elements
,c(nsmax)                                                         &
                 ! Above-diagonal matrix elements
,dt(nsmax)                                                        &
                 ! Temperature increments (k)
,r(nsmax)        ! Matrix equation rhs

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------
! External routines called:
!-----------------------------------------------------------------------
EXTERNAL tridag


!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook('SNOWPACK',zhook_in,zhook_handle)
DO k=1,tile_pts
  i = tile_index(k)

  g_snow_surf = surf_htf_tile(i)
!-----------------------------------------------------------------------
! Add melt to snow surface heat flux, unless using the snow canopy model
!-----------------------------------------------------------------------
  IF ( .NOT. cansnowtile )                                        &
    g_snow_surf = g_snow_surf + lf*melt_tile(i)

  IF ( nsnow(i) == 0 ) THEN

! Add snowfall (including canopy unloading) to ground snowpack.
    snowmass(i) = snowmass(i) + snowfall(i)

    IF ( .NOT. cansnowtile ) THEN
!-----------------------------------------------------------------------
! Remove sublimation and melt from snowpack.
!-----------------------------------------------------------------------
      snowmass(i) = snowmass(i) -                                 &
                        ( ei_tile(i) + melt_tile(i) ) * timestep

    ELSE IF ( tsoil1(i) > tm ) THEN
!-----------------------------------------------------------------------
! For canopy model, calculate melt of snow on ground underneath canopy.
!-----------------------------------------------------------------------
      smclf = rho_water * dzsoil(1) * v_sat1(i) * sthf1(i)
      asoil = 1./ ( dzsoil(1) * hcaps(i) +                        &
                      hcapw * (smcl1(i) - smclf) + hcapi*smclf )
      submelt = MIN( snowmass(i) / timestep,                      &
                     (tsoil1(i) - tm) / (lf * asoil * timestep) )
      snowmass(i) = snowmass(i) - submelt * timestep
      tsoil1(i) = tsoil1(i) -                                     &
                   tile_frac(i) * asoil * timestep * lf * submelt
      melt_tile(i) = melt_tile(i) + submelt
      snomlt_sub_htf(i) = snomlt_sub_htf(i) +                     &
                             tile_frac(i) * lf * submelt

    END IF

! Set flux into uppermost snow/soil layer (after melting).
    snow_soil_htf(i) = surf_htf_tile(i)

! Diagnose snow depth.
    snowdepth(i) = snowmass(i) / rho_snow_grnd(i)

! Set values for the surface layer. These are only needed for nsmax>0.
    IF ( nsmax > 0 ) THEN
      rho0(i) = rho_snow_grnd(i)
      sice0(i) = snowmass(i)
      tsnow0(i) = MIN( tsoil1(i), tm )
    END IF

! Note with nsmax>0 the density of a shallow pack (nsnow=0) does not
! evolve with time (until it is exhausted). We could consider updating
! the density using a mass-weighted mean of the pack density and that
! of fresh snow, so that a growing pack would reach nsnow=1 more quickly.
! This would only affect a pack that grows from a non-zero state, and
! is not an issue if the pack grows from zero, because in that case the
! the density was previously set to the fresh snow value.

  ELSE

!-----------------------------------------------------------------------
! There is at least one snow layer. Calculate heat conduction between
!   layers and temperature increments.
!-----------------------------------------------------------------------

! Save rate of melting of snow on the canopy.
    IF ( cansnowtile ) THEN
      can_melt = melt_tile(i)
    ELSE
      can_melt = 0.0
    END IF

    IF ( nsnow(i) == 1 ) THEN

! Single layer of snow.
      asnow(1) = 2.0 /                                            &
                 ( snowdepth(i)/ksnow(i,1) + dzsoil(1)/hcons(i) )
      snow_soil_htf(i) = asnow(1) * ( tsnow(i,1) - tsoil1(i) )
      dt(1) = ( g_snow_surf - snow_soil_htf(i) ) * timestep /     &
                     ( csnow(i,1) + GAMMA * asnow(1) * timestep )
      snow_soil_htf(i) = asnow(1) *                                 &
                         ( tsnow(i,1) + GAMMA*dt(1) - tsoil1(i) )
      tsnow(i,1) = tsnow(i,1) + dt(1)

    ELSE

! Multiple snow layers.
      DO n=1,nsnow(i)-1
        asnow(n) = 2.0 /                                          &
                   ( ds(i,n)/ksnow(i,n) + ds(i,n+1)/ksnow(i,n+1) )
      END DO
      n = nsnow(i)
      asnow(n) = 2.0 / ( ds(i,n)/ksnow(i,n) + dzsoil(1)/hcons(i) )
      a(1) = 0.
      b(1) = csnow(i,1) + GAMMA*asnow(1)*timestep
      c(1) = -GAMMA * asnow(1) * timestep
      r(1) = ( g_snow_surf - asnow(1)*(tsnow(i,1)-tsnow(i,2)) )   &
                                                      * timestep
      DO n=2,nsnow(i)-1
        a(n) = -GAMMA * asnow(n-1) * timestep
        b(n) = csnow(i,n) + GAMMA * ( asnow(n-1) + asnow(n) )     &
                                                       * timestep
        c(n) = -GAMMA * asnow(n) * timestep
        r(n) = asnow(n-1)*(tsnow(i,n-1) - tsnow(i,n) ) * timestep &
               +  asnow(n)*(tsnow(i,n+1) - tsnow(i,n)) * timestep
      END DO
      n = nsnow(i)
      a(n) = -GAMMA * asnow(n-1) * timestep
      b(n) = csnow(i,n) + GAMMA * (asnow(n-1)+asnow(n)) * timestep
      c(n) = 0.
      r(n) = asnow(n-1)*( tsnow(i,n-1) - tsnow(i,n) ) * timestep  &
              + asnow(n) * ( tsoil1(i) - tsnow(i,n) ) * timestep

!-----------------------------------------------------------------------
! Call the tridiagonal solver.
!-----------------------------------------------------------------------
! DEPENDS ON: tridag
      CALL tridag( nsnow(i),nsmax,a,b,c,r,dt )

      n = nsnow(i)
      snow_soil_htf(i) = asnow(n) *                               &
                         ( tsnow(i,n) + GAMMA*dt(n) - tsoil1(i) )
      DO n=1,nsnow(i)
        tsnow(i,n) = tsnow(i,n) + dt(n)
      END DO
    END IF  !  NSNOW

!-----------------------------------------------------------------------
! Melt snow in layers with temperature exceeding melting point
!-----------------------------------------------------------------------
    DO n=1,nsnow(i)
      coldsnow = csnow(i,n)*(tm - tsnow(i,n))
      IF ( coldsnow < 0 ) THEN
        tsnow(i,n) = tm
        dsice = -coldsnow / lf
        IF ( dsice > sice(i,n) ) dsice = sice(i,n)
        ds(i,n) = ( 1.0 - dsice/sice(i,n) ) * ds(i,n)
        sice(i,n) = sice(i,n) - dsice
        sliq(i,n) = sliq(i,n) + dsice
      END IF
    END DO
! Melt still > 0? - no snow left

!-----------------------------------------------------------------------
! Remove snow by sublimation unless snow is beneath canopy
!-----------------------------------------------------------------------
    IF ( .NOT. cansnowtile ) THEN
      dsice = MAX( ei_tile(i), 0. ) * timestep
      IF ( dsice > 0.0 ) THEN
        DO n=1,nsnow(i)
          IF ( dsice > sice(i,n) ) THEN
! Layer sublimates completely
            dsice = dsice - sice(i,n)
            sice(i,n) = 0.
            ds(i,n) = 0.
          ELSE
! Layer sublimates partially
            ds(i,n) = (1.0 - dsice/sice(i,n))*ds(i,n)
            sice(i,n) = sice(i,n) - dsice
            EXIT    !   sublimation exhausted
          END IF
        END DO
      END IF  !  DSICE>0
    END IF  !  CANSNOWTILE

!-----------------------------------------------------------------------
! Move liquid water in excess of holding capacity downwards or refreeze.
!-----------------------------------------------------------------------
    win = 0.0
    DO n=1,nsnow(i)
      sliq(i,n) = sliq(i,n) + win
      win = 0.0
      sliqmax = snowliqcap * rho_water * ds(i,n)
      IF (sliq(i,n) > sliqmax) THEN
! Liquid capacity exceeded
        win = sliq(i,n) - sliqmax
        sliq(i,n) = sliqmax
      END IF
      coldsnow = csnow(i,n)*(tm - tsnow(i,n))
      IF (coldsnow > 0) THEN
! Liquid can freeze
        dsice = MIN(sliq(i,n), coldsnow / lf)
        sliq(i,n) = sliq(i,n) - dsice
        sice(i,n) = sice(i,n) + dsice
        tsnow(i,n) = tsnow(i,n) + lf*dsice/csnow(i,n)
      END IF
    END DO

!-----------------------------------------------------------------------
! The remaining liquid water flux is melt.
! Include any separate canopy melt in this diagnostic.
!-----------------------------------------------------------------------
    melt_tile(i) = ( win / timestep ) + can_melt

!-----------------------------------------------------------------------
! Diagnose layer densities
!-----------------------------------------------------------------------
    DO n=1,nsnow(i)
      rho_snow(i,n) = 0.
      IF ( ds(i,n) > EPSILON(ds) )                                &
        rho_snow(i,n) = (sice(i,n) + sliq(i,n)) / ds(i,n)
    END DO

!-----------------------------------------------------------------------
! Add snowfall and frost as layer 0.
!-----------------------------------------------------------------------
    sice0(i) = snowfall(i)
    IF ( .NOT. cansnowtile )                                      &
      sice0(i) = snowfall(i) - MIN(ei_tile(i), 0.) * timestep
    tsnow0(i) = tsnow(i,1)
    rho0(i) = rho_snow_fresh

!-----------------------------------------------------------------------
! Diagnose total snow depth and mass
!-----------------------------------------------------------------------
    snowdepth(i) = sice0(i) / rho0(i)
    snowmass(i) = sice0(i)
    DO n=1,nsnow(i)
      snowdepth(i) = snowdepth(i) + ds(i,n)
      snowmass(i) = snowmass(i) + sice(i,n) + sliq(i,n)
    END DO

  END IF    !  NSNOW

END DO  !  k (points)
IF (lhook) CALL dr_hook('SNOWPACK',zhook_out,zhook_handle)
RETURN

END SUBROUTINE snowpack
