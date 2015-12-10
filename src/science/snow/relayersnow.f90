! *****************************COPYRIGHT*******************************

! (c) [University of Edinburgh] [2009]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the JULES collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC237]

! *****************************COPYRIGHT*******************************
!  SUBROUTINE RELAYERSNOW-----------------------------------------------

! Description:
!     Redivide snowpack after changes in depth, conserving mass and energy

! Subroutine Interface:
SUBROUTINE relayersnow( land_pts,tile_pts,tile_index,             &
                        rgrain0,rho0,sice0,snowfall,snowmass,     &
                        tsnow0,nsnow,ds,rgrain,rgrainl,sice,      &
                        rho_snow_grnd,sliq,tsnow,rho_snow,        &
                        snowdepth )

USE ancil_info, ONLY :                                            &
!  imported scalars with intent(in)
 nsmax           !  Maximum possible number of snow layers

USE c_0_dg_c, ONLY :                                              &
!  imported scalar parameters
 tm   !   temperature at which fresh water freezes and ice melts (K)

USE c_perma, ONLY :                                               &
!  imported scalar parameters
 hcapi                                                            &
            !  Specific heat capacity of ice (J/kg/K)
,hcapw      !  Specific heat capacity of water (J/kg/K)

USE snow_param, ONLY :                                            &
!  imported scalars with intent(in)
 rho_snow_fresh  !  density of fresh snow (kg per m**3)

USE rad_param, ONLY :                                             &
!  imported scalars with intent(in)
 r0         !  Grain size for fresh snow (microns)

USE switches, ONLY :                                              &
!  imported scalars with intent(in)
 l_snow_albedo                                                    &
            ! Switch for prognostic snow albedo
,l_rho_snow_corr
            ! Switch for snowpack density correction

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Scalar arguments with intent(in)
INTEGER, INTENT(IN) ::                                            &
 land_pts                                                         &
                    ! Total number of land points
,tile_pts           ! Number of tile points

! Array arguments with intent(in)
INTEGER, INTENT(IN) ::                                            &
 tile_index(land_pts)     ! Index of tile points

REAL, INTENT(IN) ::                                               &
 rgrain0(land_pts)                                                &
                    ! Fresh snow grain size (microns)
,rho0(land_pts)                                                   &
                    ! Density of fresh snow (kg/m3)
!                         ! Where NSNOW=0, rho0 is the density
!                         ! of the snowpack.
,sice0(land_pts)                                                  &
                    ! Ice content of fresh snow (Kg/m2)
!                         ! Where NSNOW=0, SICE0 is the mass
!                         ! of the snowpack.
,snowfall(land_pts)                                               &
                    ! Snowfall reaching the ground (kg/m2)
,snowmass(land_pts)                                               &
                    ! Snow mass on the ground (kg/m2)
,tsnow0(land_pts)   ! Temperature of fresh snow (K)

! Array arguments with intent(inout)
INTEGER, INTENT(INOUT) ::                                         &
 nsnow(land_pts)       ! Number of snow layers

REAL, INTENT(INOUT) ::                                            &
 ds(land_pts,nsmax)                                               &
                       ! Snow layer thicknesses (m)
,rgrain(land_pts)                                                 &
                       ! Snow surface grain size (microns)
,rgrainl(land_pts,nsmax)                                          &
                       ! Snow grain size (microns)
,rho_snow_grnd(land_pts)                                          &
                       ! Snowpack bulk density (kg/m3)
,sice(land_pts,nsmax)                                             &
                       ! Ice content of snow layers (kg/m2)
,sliq(land_pts,nsmax)                                             &
                       ! Liquid content of snow layers (kg/m2)
,tsnow(land_pts,nsmax) ! Snow layer temperatures (K)

! Array arguments with intent(out)
REAL, INTENT(OUT) ::                                              &
 rho_snow(land_pts,nsmax)                                         &
                       ! Snow layer densities (kg/m3)
,snowdepth(land_pts)   ! Snow depth (m)

! Local scalar parameters
REAL, PARAMETER :: thin_snow_limit = 1.0e-6
                  ! Maximum snow thickness (m) that is neglected
                  ! during relayering. All contributions
                  ! (mass, energy etc) from that snow are
                  ! neglected.

! Local scalars
INTEGER ::                                                        &
 i                                                                &
                       ! Land point index
,iznew                                                            &
                       ! layer index
,izz                                                              &
                       ! layer index
,k                                                                &
                       ! Tile point index
,n                                                                &
                       ! Snow layer index
,new                                                              &
                       ! layer index
,old                   ! layer index

REAL ::                                                           &
 csnow                                                            &
                       ! Areal heat capacity of layer (J/K/m2)
,oldremains                                                       &
                       ! remaining depth in an old layer (m)
,wt                    ! weight given to a layer value

! Local arrays
INTEGER ::                                                        &
 nold(land_pts)        ! Number of layers before adjustment

REAL ::                                                           &
 d0(land_pts,0:nsmax)                                             &
                        ! Layer thicknesses before adjustment (m)
!                             ! D0(:,0) represents new snow if NSNOW>0,
!                             ! otherwise it is all snow.
,e(0:nsmax)                                                       &
                     ! Internal energy before adjustment (J/m2)
,newremains(nsmax)                                                &
                     ! available (unfilled) depth in new layer (m)
,r(0:nsmax)                                                       &
                     ! Grain size before adjustment (kg/m2)
,s(0:nsmax)                                                       &
                     ! Ice content before adjustment (kg/m2)
,w(0:nsmax)                                                       &
                     ! Liquid content before adjustment (kg/m2)
,u(nsmax)            ! Layer energy contents

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------
! External routines used
!-----------------------------------------------------------------------
EXTERNAL layersnow

IF (lhook) CALL dr_hook('RELAYERSNOW',zhook_in,zhook_handle)

!-----------------------------------------------------------------------
! Initialise grain variables if they aren't used.
!-----------------------------------------------------------------------
IF (.NOT.l_snow_albedo) THEN
  rgrain(:)  = r0
  rgrainl(:,:) = r0
END IF

!-----------------------------------------------------------------------
! Initialise snowdepth with value that is retained where tile frac=0.
!-----------------------------------------------------------------------
snowDepth(:) = 0.0

!-----------------------------------------------------------------------
! Store previous layer thicknesses.
!-----------------------------------------------------------------------
d0(:,0) = 0.0
nold(:) = nsnow(:)
DO k=1,tile_pts
  i = tile_index(k)
  IF ( sice0(i) > 0.0 ) d0(i,0) = sice0(i) / rho0(i)
  d0(i,1:nsnow(i)) = ds(i,1:nsnow(i))
END DO

!-----------------------------------------------------------------------
! Calculate snowdepth
!-----------------------------------------------------------------------
DO k=1,tile_pts
  i = tile_index(k)
  snowdepth(i) = d0(i,0)
  DO n=1,nsnow(i)
    snowdepth(i) = snowdepth(i) + d0(i,n)
  END DO
END DO

!-----------------------------------------------------------------------
! Divide snowpack into new layers
!-----------------------------------------------------------------------
! DEPENDS ON: layersnow
CALL layersnow( land_pts,tile_pts,tile_index,snowdepth,nsnow,ds )

DO k=1,tile_pts
  i = tile_index(k)

  IF ( nsnow(i) > 0 ) THEN
!-----------------------------------------------------------------------
! Store previous snow layer energy and mass contents
!-----------------------------------------------------------------------
    csnow = sice0(i) * hcapi
    e(0) = csnow * ( tsnow0(i) - tm )
    r(0) = rgrain0(i)
! If the previous timestep had no layers (but could still have snow),
! estimate R(0) using average of previous grain size and that of
! fresh snow. In fact fresh snow would be on top and so more important!
    IF ( nold(i) == 0 )                                           &
      r(0) = ( 1.0-snowfall(i)/snowmass(i)) * rgrain(i) +         &
                      snowfall(i)/snowmass(i) * rgrain0(i)
    s(0) = sice0(i)
    w(0) = 0.0
    DO n=1,nold(i)
      csnow = sice(i,n)*hcapi + sliq(i,n)*hcapw
      e(n) = csnow * ( tsnow(i,n) - tm )
      r(n) = rgrainl(i,n)
      s(n) = sice(i,n)
      w(n) = sliq(i,n)
    END DO

!-----------------------------------------------------------------------
! Initialise accumulations for new layer values.
!-----------------------------------------------------------------------
    u(:) = 0.
    sice(i,:) = 0.
    sliq(i,:) = 0.
    rgrainl(i,:) = 0.

!-----------------------------------------------------------------------
! Set the state of the new layers.
!-----------------------------------------------------------------------
!    Initialise with all new layers empty.
    newremains(1:nsnow(i)) = ds(i,1:nsnow(i))
!     Start by filling top new layer.
    iznew = 1

!     Loop over the old layers.
    DO old=0,nold(i)

!       All of this old layer remains to be reassigned to new layer(s).
      oldremains = d0(i,old)

!       Point to first new layer with remaining space.
      izz = iznew

!       Loop over new layers with remaining space.
      DO new=izz,nsnow(i)

        IF ( oldremains > newremains(new) ) THEN
!-----------------------------------------------------------------------
! The remaining depth in the new layer will be exhausted by some or
! all of the remaining depth from the old layer.
!-----------------------------------------------------------------------

!           Decrement old layer by the remaining space in new layer.
          oldremains = oldremains - newremains(new)

!           Add properties from old layer to accumulation for new layer.
!           Note that wt is <= 1 since here we have oldRemains>newRemains,
!           and oldRemains <= d0.
          IF ( d0(i,old) > thin_snow_limit ) THEN
            wt =  newremains(new) / d0(i,old)
            u(new) = u(new) + e(old) * wt
            sice(i,new) = sice(i,new) + s(old) * wt
            sliq(i,new) = sliq(i,new) + w(old) * wt
            rgrainl(i,new) = rgrainl(i,new) +                     &
                             r(old) * newremains(new)
          END IF

!           Update the pointer to the next new layer with space.
          izz = new + 1

        ELSE

!-----------------------------------------------------------------------
! The old layer will be exhausted by this increment.
!-----------------------------------------------------------------------
!           Decrement available space in the new layer.
          newremains(new) = newremains(new) - oldremains
!           Add properties from old layer to accumulation for new layer.
          IF ( d0(i,old) > thin_snow_limit ) THEN
            wt = oldremains /  d0(i,old)
            u(new) = u(new) + e(old) * wt
            sice(i,new) = sice(i,new) + s(old) * wt
            sliq(i,new) = sliq(i,new) + w(old) * wt
            rgrainl(i,new) = rgrainl(i,new) + r(old) * oldremains
          END IF
!           Proceed to the next old layer by exiting from the new layer loop.
          EXIT
        END IF
      END DO  !  new layers
!       Update pointer to the next new layer with space.
      iznew = izz
    END DO  !  old layers

!-----------------------------------------------------------------------
! Diagnose layer temperatures and densities.
!-----------------------------------------------------------------------
    DO n=1,nsnow(i)
      csnow = sice(i,n)*hcapi + sliq(i,n)*hcapw
      tsnow(i,n) = tm + u(n) / csnow
      rho_snow(i,n) = ( sice(i,n) + sliq(i,n) ) / ds(i,n)
      rgrainl(i,n) = rgrainl(i,n) / ds(i,n)
    END DO

!-----------------------------------------------------------------------
! Snow surface grain size for radiative calculations
!-----------------------------------------------------------------------
    rgrain(i) = rgrainl(i,1)

!-----------------------------------------------------------------------
! Diagnose bulk density of pack.
!-----------------------------------------------------------------------
    rho_snow_grnd(i) = snowmass(i) / snowdepth(i)

  ELSE

!-----------------------------------------------------------------------
! Set bulk density of pack to a constant value if there is (effectively)
! no snow. This density is then used the next time a shallow pack forms.
! We could also recalculate snowdepth, to make it consistent
! with the revised density, but we're not bothering as depths are tiny!
! Note that because the bulk density is not calculated for nsnow=0 and
! snowmass>0, it remains constant (until the pack is exhausted or it
! grows to nsnow=1).
!-----------------------------------------------------------------------
    IF ( snowmass(i) < 1.0e-9 ) THEN
      rho_snow_grnd(i) = rho_snow_fresh
    ELSE IF ( l_rho_snow_corr ) THEN
!-----------------------------------------------------------------------
! Diagnose bulk density of pack. On most occasions this density is
! unchanged from input value, but the calculation is required for
! timesteps when nsnow changes to zero.
!-----------------------------------------------------------------------
      rho_snow_grnd(i) = snowmass(i) / snowdepth(i)
    END IF

  END IF   !  nsnow

!-----------------------------------------------------------------------
! Set values for unused snow layers.
! Note: not needed for algorithm, but clearer to follow.
!-----------------------------------------------------------------------
  IF (  nsnow(i) < nsmax ) THEN
    n = nsnow(i) + 1
    rgrainl(i,n:) = r0
    rho_snow(i,n:) = 0.0
    sice(i,n:) = 0.
    sliq(i,n:) = 0.
    tsnow(i,n:) = tm
  END IF

END DO  !  K (points)
IF (lhook) CALL dr_hook('RELAYERSNOW',zhook_out,zhook_handle)
RETURN

END SUBROUTINE relayersnow
