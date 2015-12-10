! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  SUBROUTINE SURF_HYD-----------------------------------------------

!  PURPOSE : TO CARRY OUT CANOPY AND SURFACE HYDROLOGY CALCULATIONS

!            CANOPY WATER CONTENT IS DEPRECIATED BY EVAPORATION

!            SNOWMELT IS RUNOFF THE SURFACE WITHOUT INTERACTING
!            WITH THE CANOPY

!            THE CANOPY INTERCEPTION AND SURFACE RUNOFF OF
!            LARGE-SCALE RAIN IS CALCUALTED

!            THE CANOPY INTERCEPTION AND SURFACE RUNOFF OF
!            CONVECTIVE RAIN IS CALCUALTED


!  SUITABLE FOR SINGLE COLUMN MODEL USE

!  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER NO 25

!--------------------------------------------------------------------

!    ARGUMENTS---------------------------------------------------------

SUBROUTINE surf_hyd (npnts,ntiles,tile_pts,tile_index,            &
                     can_cpy,e_canopy,tile_frac,infil_tile,       &
                     con_rain,ls_rain,melt_tile,snow_melt,        &
                     timestep,can_wcnt,can_wcnt_gb,dsmc_dt,       &
                     l_top,l_pdm,nshyd,soil_pts,soil_index,       &
                     surf_roff,tot_tfall,dun_roff,fsat,           &
                     v_sat,sthu,sthf)

USE soil_param, ONLY :                                            &
 confrac   !  the fraction of the gridbox over which convective
           !  precipitation is assumed to fall

USE switches, ONLY :                                              &
 l_point_data         ! Switch for using point rainfall data

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

INTEGER                                                           &
 npnts                                                            &
                      ! IN Total number of land points.
,ntiles                                                           &
                      ! IN Number of tiles.
,nshyd                                                            &
                     ! IN Number of soil moisture levels.
,soil_pts                                                         &
                     ! IN Number of soil points.
,tile_pts(ntiles)                                                 &
                      ! IN Number of tile points.
,tile_index(npnts,ntiles)                                         &
!                           ! IN Index of tile points.
,soil_index(npnts)   ! IN Array of soil points.

REAL                                                              &
 can_cpy(npnts,ntiles)                                            &
                      ! IN Canopy capacity for land tiles (kg/m2).
,e_canopy(npnts,ntiles)                                           &
                       !IN Canopy evaporation (kg/m2/s).
,tile_frac(npnts,ntiles)                                          &
                      ! IN Tile fractions.
,infil_tile(npnts,ntiles)                                              &
                      ! IN Tile infiltration rate (kg/m2/s).
,con_rain(npnts)                                                  &
                      ! IN Convective rain (kg/m2/s).
,ls_rain(npnts)                                                   &
                      ! IN Large-scale rain (kg/m2/s).
,melt_tile(npnts,ntiles)                                          &
!                           ! IN Snow melt on tiles (kg/m2/s).
,snow_melt(npnts)                                                 &
                      ! IN GBM snow melt (kg/m2/s).
,fsat(npnts)                                                      &
                     ! IN Surface saturation fraction.
,sthf(npnts,nshyd)                                                &
                     ! INOUT Frozen soil moisture content of
!                          !       each layer as a fraction of
!                          !       saturation.
,sthu(npnts,nshyd)                                                &
                     ! INOUT Unfrozen soil moisture content of
!                          !       each layer as a fraction of
!                          !       saturation.
,v_sat(npnts,nshyd)                                               &
                     ! IN Volumetric soil moisture
!                          !    concentration at saturation
!                          !    (m3 H2O/m3 soil).

,timestep             ! IN Timestep (s).


LOGICAL                                                           &
 l_top                                                            &
                      ! IN TOPMODEL-based hydrology logical.
,l_pdm                ! IN PDM logical.
REAL                                                              &
 can_wcnt(npnts,ntiles)!INOUT Tile canopy water contents (kg/m2).

REAL                                                              &
 can_wcnt_gb(npnts)                                               &
                      ! OUT Gridbox canopy water content (kg/m2).
,dsmc_dt(npnts)                                                   &
                      ! OUT Rate of change of soil moisture
                      !     content (kg/m2/s).
,surf_roff(npnts)                                                 &
                      ! OUT Cumulative surface runoff (kg/m2/s).
,tot_tfall(npnts)                                                 &
                      ! OUT Cumulative canopy throughfall
                      !     (kg/m2/s).
,dun_roff(npnts)      ! OUT Cumulative Dunne sfc runoff (kg/m2/s).

!  Workspace -----------------------------------------------------------
REAL                                                              &
 r                                                                &
                      ! WORK Total downward water flux (i.e. rain
                      ! + condensation + snowmelt) kg/m**2/s
,tfall                                                            &
                      ! WORK Cumulative canopy throughfall
                      !      on a tile (kg/m2/s).
,s_roff               ! WORK Cumulative surface runoff
                      !      on a tile (kg/m2/s).
REAL                                                              &
 can_cond(npnts)      ! Canopy condensation (kg/m2/s).

INTEGER                                                           &
 i                                                                &
                      ! Loop counter (land points).
,j                                                                &
                      ! Loop counter (tile points).
,n                    ! Loop counter (tiles).

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!---------------------------------------------------------------------
!  External routines called :-

EXTERNAL frunoff,sieve,pdm


!-----------------------------------------------------------------------
! Zero cumulative stores
IF (lhook) CALL dr_hook('SURF_HYD',zhook_in,zhook_handle)
DO i=1,npnts
  can_wcnt_gb(i) = 0.0
  tot_tfall(i) = 0.0
  surf_roff(i) = 0.0
  dsmc_dt(i)   = 0.0
  dun_roff(i) = 0.0
END DO

! Reduce canopy water content by evaporation
DO n=1,ntiles
  DO j=1,tile_pts(n)
  i = tile_index(j,n)
    IF (e_canopy(i,n)  >   0.0)                                   &
      can_wcnt(i,n) =                                             &
                 MAX( can_wcnt(i,n) - e_canopy(i,n)*timestep, 0. )
  END DO
END DO

IF(l_point_data) THEN
!-----------------------------------------------------------------------
!       Using point precipitation data.
!-----------------------------------------------------------------------
  DO n=1,ntiles
    DO j=1,tile_pts(n)
      i = tile_index(j,n)
!           Calculate total downward flux of liquid water.
      r = ls_rain(i) + con_rain(i)
!           Add any condensation
      IF ( e_canopy(i,n) < 0.0 ) r = r - e_canopy(i,n)
!           Calculate throughfall.
      IF ( r <= can_cpy(i,n) / timestep .AND.                     &
           can_cpy(i,n) > 0 ) THEN
        tfall = r * can_wcnt(i,n) / can_cpy(i,n)
      ELSE
        tfall = r - ( can_cpy(i,n) - can_wcnt(i,n) ) / timestep
      END IF
!           Update canopy water content.
      can_wcnt(i,n) = can_wcnt(i,n) + ( r - tfall ) * timestep
!           Add melt to throughfall before calculation of runoff.
      tfall = tfall + melt_tile(i,n)
!           Calculate surface runoff.
      IF ( tfall > infil_tile(i,n) ) THEN
        s_roff = tfall - infil_tile(i,n)
      ELSE
        s_roff = 0.0
      END IF

!           Add to gridbox accumulations.
!           Don't include melt in throughfall.
      tot_tfall(i) = tot_tfall(i) +                               &
                 ( tfall - melt_tile(i,n) ) * tile_frac(i,n)
      surf_roff(i) = surf_roff(i) + s_roff * tile_frac(i,n)
      can_wcnt_gb(i) = can_wcnt_gb(i) +                           &
                           can_wcnt(i,n) * tile_frac(i,n)
    END DO
  END DO

ELSE
!-----------------------------------------------------------------------
!       Using area-average precipitation data. Assume spatial distribution.
!-----------------------------------------------------------------------

  DO n=1,ntiles

! Surface runoff of snowmelt, assumed to cover 100% of tile
!   DEPENDS ON: frunoff
    CALL frunoff (npnts,tile_pts(n),tile_index(:,n),1.,           &
                can_cpy(:,n),can_cpy(:,n),infil_tile(:,n),        &
                melt_tile(:,n),tile_frac(:,n),timestep,           &
                surf_roff)

! Define canopy condensation when evaporation is negative
    DO j=1,tile_pts(n)
    i = tile_index(j,n)
      IF ( e_canopy(i,n)  <   0. ) THEN
       can_cond(i) = - e_canopy(i,n)
      ELSE
       can_cond(i) = 0.
      END IF
    END DO

! Canopy interception, throughfall and surface runoff for condensation,
! assumed to cover 100% of gridbox
! DEPENDS ON: sieve
    CALL sieve (npnts,tile_pts(n),tile_index(:,n),1.,             &
                can_cpy(:,n),can_cond,tile_frac(:,n),timestep,    &
                can_wcnt(:,n),tot_tfall)
! DEPENDS ON: frunoff
    CALL frunoff (npnts,tile_pts(n),tile_index(:,n),1.,           &
                  can_cpy(:,n),can_wcnt(:,n),infil_tile(:,n),     &
                  can_cond,tile_frac(:,n),timestep,               &
                  surf_roff)

! Canopy interception, throughfall and surface runoff for large-scale
! rain, assumed to cover 100% of gridbox
! DEPENDS ON: sieve
    CALL sieve (npnts,tile_pts(n),tile_index(:,n),1.,             &
                can_cpy(:,n),ls_rain,tile_frac(:,n),timestep,     &
                can_wcnt(:,n),tot_tfall)
! DEPENDS ON: frunoff
    CALL frunoff (npnts,tile_pts(n),tile_index(:,n),1.,           &
                  can_cpy(:,n),can_wcnt(:,n),infil_tile(:,n),     &
                  ls_rain,tile_frac(:,n),timestep,                &
                  surf_roff)

! Canopy interception, throughfall and surface runoff for convective
! rain, assumed to cover fraction CONFRAC of gridbox
! DEPENDS ON: sieve
    CALL sieve (npnts,tile_pts(n),tile_index(:,n),confrac,        &
                can_cpy(:,n),con_rain,tile_frac(:,n),timestep,    &
                can_wcnt(:,n),tot_tfall)
! DEPENDS ON: frunoff
    CALL frunoff (npnts,tile_pts(n),tile_index(:,n),confrac,      &
                  can_cpy(:,n),can_wcnt(:,n),infil_tile(:,n),     &
                  con_rain,tile_frac(:,n),timestep,               &
                  surf_roff)

    DO j=1,tile_pts(n)
      i = tile_index(j,n)
      can_wcnt_gb(i) = can_wcnt_gb(i)                             &
                       + tile_frac(i,n)*can_wcnt(i,n)
    END DO

  END DO
END IF   !   l_point_data

!-----------------------------------------------------------------------
! Calculate Saturation excess runoff through PDM:
!-----------------------------------------------------------------------

IF(soil_pts >  0.AND.l_pdm)                                       &
! DEPENDS ON: pdm
  CALL pdm(                                                       &
           npnts,soil_pts,soil_index,nshyd,                       &
          tot_tfall,snow_melt,surf_roff,timestep,                 &
          v_sat,dun_roff,sthu,sthf)

IF(l_top)THEN
  DO i=1,npnts
    dun_roff(i)=fsat(i)*(tot_tfall(i) + snow_melt(i)              &
                  - surf_roff(i))
  END DO
END IF
DO i=1,npnts
  IF(l_top.OR.l_pdm)surf_roff(i) = surf_roff(i) + dun_roff(i)
  dsmc_dt(i) = tot_tfall(i) + snow_melt(i) - surf_roff(i)
END DO

IF (lhook) CALL dr_hook('SURF_HYD',zhook_out,zhook_handle)
RETURN
END SUBROUTINE surf_hyd
