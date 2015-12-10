! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    SUBROUTINE SOIL_HTC------------------------------------------------

! Description:
!     Updates deep soil temperatures, frozen and unfrozen
!     frozen soil water content.

! Documentation : UM Documentation Paper 25

! Subroutine Interface:
SUBROUTINE soil_htc (                                             &
 npnts,nshyd,ntiles,soil_pts,soil_index,tile_pts,tile_index       &
,nsnow,bexp,dz,tile_frac,hcap,hcon,sathh,surf_ht_flux,timestep    &
,v_sat,w_flux,smcl,snowdepth,sthu,sthf,tsoil                      &
,ltimer                                                           &
)

USE c_0_dg_c, ONLY :                                              &
!      imported scalar parameters
   zerodegc

USE c_perma, ONLY :                                               &
!      imported scalar parameters
   hcapi,hcapw,dpsidt

USE c_lheat, ONLY :                                               &
!      imported scalar parameters
   lf

USE c_densty, ONLY :                                              &
!      imported scalar parameters
   rho_water

USE snow_param, ONLY :                                            &
!      imported scalars with intent(in)
   snow_hcon

USE soil_param, ONLY :                                            &
!      imported scalar parameters
   facur,GAMMA=>gamma_t,mmax,tacur

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Subroutine arguments
!   Scalar arguments with intent(IN) :
INTEGER                                                           &
 npnts                                                            &
                      ! IN Number of gridpoints.
,nshyd                                                            &
                      ! IN Number of soil moisture levels.
,ntiles                                                           &
                      ! IN Number of tiles.
,soil_pts             ! IN Number of soil points.

REAL                                                              &
 timestep             ! IN Model timestep (s).


!   Array arguments with intent(IN) :
INTEGER                                                           &
 soil_index(npnts)                                                &
                      ! IN Array of soil points.
,tile_pts(ntiles)                                                 &
                      ! IN Number of tile points.
,tile_index(npnts,ntiles)                                         &
!                           ! IN Index of tile points.
,nsnow(npnts,ntiles)  ! IN Number of snow layerss

REAL                                                              &
 bexp(npnts,nshyd)                                                &
                      ! IN Clapp-Hornberger exponent.
,dz(nshyd)                                                        &
                      ! IN Thicknesses of the soil layers (m).
,tile_frac(npnts,ntiles)                                          &
                      ! IN Tile fractions.
,hcap(npnts,nshyd)                                                &
                      ! IN Soil heat capacity (J/K/m3).
,hcon(npnts,0:nshyd)                                              &
                      ! IN Soil thermal conductivity (W/m/K).
,sathh(npnts,nshyd)                                               &
                      ! IN Saturated soil water pressure (m).
,smcl(npnts,nshyd)                                                &
                      ! IN Soil moisture content of each
!                           !    layer (kg/m2).
,snowdepth(npnts,ntiles)                                          &
!                           ! IN Snow depth (on ground) (m)
,surf_ht_flux(npnts)                                              &
                      ! IN Net downward surface heat flux (W/m2).
,v_sat(npnts,nshyd)                                               &
                      ! IN Volumetric soil moisture
!                           !    concentration at saturation
!                           !    (m3 H2O/m3 soil).
,w_flux(npnts,0:nshyd)! IN The fluxes of water between layers
!                           !    (kg/m2/s).

LOGICAL ltimer        ! Logical switch for TIMER diags

!   Array arguments with intent(OUT) :

!   Array arguments with intent(INOUT) :
REAL                                                              &
 sthf(npnts,nshyd)                                                &
                      ! INOUT Frozen soil moisture content of
!                           !       each layer as a fraction of
!                           !       saturation.
,sthu(npnts,nshyd)                                                &
                      ! INOUT Unfrozen soil moisture content of
!                           !       each layer as a fraction of
!                           !       saturation.
,tsoil(npnts,nshyd)   ! INOUT Sub-surface temperatures (K).

! Local scalars:
INTEGER                                                           &
 i,j,m,n                                                          &
                      ! WORK Loop counters.
,iter_pts             ! WORK Number of soil points which require
!                           !      iteration.
REAL                                                              &
 si_tile                                                          &
                      ! WORK Tile snow insulation factor
,small_value                                                      &
                      ! WORK
,tiny_0               ! WORK

! Local arrays:
INTEGER                                                           &
 iter_index(npnts)    ! WORK Array of soil points which require
!                           !      iteration.

REAL                                                              &
 ceacur(npnts)                                                    &
                      ! WORK Flux conservation accuracy of the
!                           !      calculation (W/m2)
,dhsl0(npnts,nshyd)                                               &
                      ! WORK Total heat increment to the layer
!                           !      (J/m2/timestep)
,dhsl(npnts,nshyd)                                                &
                      ! WORK The heat available to update the
!                           !      layer temperature (J/m2/timestep)
,dsmclf(npnts,nshyd)                                              &
                      ! WORK The increment to the layer frozen
!                           !      soil moisture (kg/m2/timestep).
,dthu(npnts,nshyd)                                                &
                      ! WORK Rate of change of volumetric unfrozen
!                           !      soil moisture concentration with
!                           !      temperature (m3 liquid H2O/m3 soil/K)
,dtsl(npnts,nshyd)                                                &
                      ! WORK The increment to the layer temperatur
!                           !      (K/timestep).
,dtslmax(npnts,nshyd)                                             &
                      ! WORK Maximum value of DTSL (K/timestep).
,dtslmin(npnts,nshyd)                                             &
                      ! WORK Minimum value of DTSL (K/timestep).
,hcapt(npnts,nshyd)                                               &
                     ! WORK The total volumetric heat capacity
!                           !      (soil+water) of the layer (J/m3/K).
,hc(npnts,nshyd)                                                  &
                      ! WORK The thermal conductivity of each
!                           !      layer (W/m/K).
,hcons(npnts)                                                     &
                      ! WORK The thermal conductivity between
!                           !      adjacent soil layers (W/m/K).
,h_flux(npnts,0:nshyd)                                            &
                       !WORK The fluxes of heat between layers
!                           !      (W/m2).
,hadv(npnts,nshyd)                                                &
                      ! WORK Heat flux due to moisture advection
!                           !      (W/m2).
,sifact(npnts)                                                    &
                      ! WORK Snow insulation factor.
,smclf(npnts,nshyd)                                               &
                      ! WORK Frozen moisture content of each
!                           !      soil layer (kg/m2).
,smclf0(npnts,nshyd)                                              &
                      ! WORK Previous value of SMCLF (kg/m2).
,smclsat(npnts,nshyd)                                             &
                      ! WORK The saturation moisture content of
!                           !      each layer (kg/m2).
,smclu(npnts,nshyd)                                               &
                      ! WORK Unfrozen moisture content of each
!                           !      soil layer (kg/m2).
,smclu0(npnts,nshyd)                                              &
                      ! WORK Previous value of SMCLU (kg/m2).
,smcfu(npnts)                                                     &
                      ! WORK Fractional saturation (unfrozen water
!                           !      at layer boundaries.
,smcff(npnts)                                                     &
                      ! WORK Fractional saturation (frozen water)
!                           !      at layer boundaries.
,tmax(npnts,nshyd)                                                &
                      ! WORK Temperature above which all water is
!                           !      unfrozen (Celsius)
,tsl(npnts,0:nshyd)                                               &
                      ! WORK Soil layer temperatures (Celsius)
!                           !      TSL(0) temperature of incoming water.
!                           !      TSL(1:NSHYD) sub-surface soil
!                           !      temperatures .
,tsl0(npnts,0:nshyd)                                              &
                      ! WORK Previous value of TSL (Celsius).
,v_satk(npnts)                                                    &
                      ! WORK Saturated volumetric soil moisture
!                           !    at the layer boundary (m3/m3)
,work1(npnts,nshyd)   ! WORK for Tmax

LOGICAL                                                           &
 iter(npnts)          ! WORK .T. on points requiring iterations.

!-----------------------------------------------------------------------
! Variables required for the implicit calculation.
!-----------------------------------------------------------------------
REAL                                                              &
 dhflux_dtsl1(npnts,0:nshyd),dhflux_dtsl2(npnts,0:nshyd)          &
,dhadv_dtsl0(npnts,nshyd),dhadv_dtsl1(npnts,nshyd)                &
,dhadv_dtsl2(npnts,nshyd)                                         &
                          ! WORK Rate of change of the explicit
!                           ! fluxes with the layer temperatures
!                           ! (W/m2/K).
,a(npnts,nshyd),b(npnts,nshyd),c(npnts,nshyd),d(npnts,nshyd)      &
!                           ! WORK Matrix elements.
,gamcon               ! WORK Forward timestep weighting constant.

INTEGER :: jj
INTEGER, PARAMETER :: j_block = 4 ! this specifies the block length  
                                  ! which enables parallelisation over
                                  ! the j dimension. Must be of sufficient
                                  ! granularity to avoid inner loop overhead
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('SOIL_HTC',zhook_in,zhook_handle)

!-----------------------------------------------------------------------
! Initialisations
!-----------------------------------------------------------------------
tiny_0=TINY(0.0)
small_value=EPSILON(0.0)

DO j=1,soil_pts
  i=soil_index(j)
  tsl(i,0)=tsoil(i,1)-zerodegc
END DO

DO n=1,nshyd
  DO j=1,soil_pts
    i=soil_index(j)
!-----------------------------------------------------------------------
! Define soil layer temperatures TSL (in celsius).
!-----------------------------------------------------------------------
    tsl(i,n)=tsoil(i,n)-zerodegc
  END DO
END DO

DO n=1,nshyd
! CDIR$ IVDEP here would force vectorization but changes results!
  DO j=1,soil_pts
    i=soil_index(j)
!-----------------------------------------------------------------------
! Diagnose the frozen and unfrozen water.
!-----------------------------------------------------------------------

    smclsat(i,n)=rho_water*dz(n)*v_sat(i,n)
    smclf(i,n)=smclsat(i,n)*sthf(i,n)
    smclu(i,n)=smcl(i,n)-smclf(i,n)
  END DO                                  !  J=1,SOIL_PTS
END DO                                    ! N=1,NSHYD

!-----------------------------------------------------------------------
! Initialise the array of points for which calculations are required.
!-----------------------------------------------------------------------
DO i=1,npnts
  iter(i)=.FALSE.
END DO

DO j=1,soil_pts
  i=soil_index(j)
  iter(i)=.TRUE.
END DO

!----------------------------------------------------------------------
! Calculate the heat conductivity between adjacent layers
!----------------------------------------------------------------------

v_satk(:)=0.0

DO n=1,nshyd-1
  DO j=1,soil_pts
    i=soil_index(j)
    smcfu(i)=(dz(n+1)*sthu(i,n)+dz(n)*sthu(i,n+1))                &
              /(dz(n+1)+dz(n))
    smcff(i)=(dz(n+1)*sthf(i,n)+dz(n)*sthf(i,n+1))                &
              /(dz(n+1)+dz(n))
!!!
!!! THIS LINE IS REPLACED WITH THE ONE FOLLOWING
!!! IN ORDER TO OBTAIN BIT-COMPARABILITY
!!!          V_SATK(I)=(DZ(N+1)*V_SAT(I,N)+DZ(N)*V_SAT(I,N+1))             &
!!!     &              /(DZ(N+1)+DZ(N))
    v_satk(i)=v_sat(i,n)
  END DO
! DEPENDS ON: heat_con
  CALL heat_con (npnts,hcon(:,n)                                  &
,                smcfu(:),smcff(:)                                &
,                v_satk(:),hcons(:),ltimer)
  DO j=1,soil_pts
    i=soil_index(j)
    hc(i,n)=hcons(i)
  END DO
END DO

!--------------------------------------------------------------------
! Calculate the snow insulation factor
!--------------------------------------------------------------------
DO i=1,npnts
  sifact(i) = 0.
END DO

!$OMP PARALLEL  DEFAULT(SHARED) PRIVATE(i, j, si_tile, n)
DO n=1,ntiles
!$OMP DO SCHEDULE(STATIC)
  DO j=1,tile_pts(n)
    i = tile_index(j,n)
    si_tile = 1.
    IF ( v_sat(i,1) > EPSILON(v_sat(i,1)) .AND.                   &
         nsnow(i,n) == 0 ) THEN
      IF ( snowdepth(i,n) <= 0.5*dz(1) ) THEN
        si_tile = 1. / ( 1. + 2.*snowdepth(i,n)/(dz(1) + dz(2)) )
      ELSE
        si_tile =(dz(1) + dz(2)) /                                &
                 ( hc(i,1)*(2.*snowdepth(i,n) - dz(1))/snow_hcon  &
                   + 2.*dz(1) + dz(2) )
      END IF
    END IF
    sifact(i) = sifact(i) + tile_frac(i,n)*si_tile
  END DO
!$OMP END DO
! implicit barrier to ensure lockstep on n
END DO
!$OMP END PARALLEL

!--------------------------------------------------------------------
! Calculate heat fluxes across layer boundaries
!--------------------------------------------------------------------

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(jj, j, i, n)

! Blocking is needed as the nshyd dimension has a spatial 
! dependence (n +/-1). Blocking enables parallelising over the 
! j (soil points) dimension.
!$OMP DO SCHEDULE(STATIC)
DO jj=1, soil_pts, j_block
  DO n=1,nshyd-1    
    DO j=jj,MIN(jj+j_block-1,soil_pts)
      i=soil_index(j)
      h_flux(i,n)=-hc(i,n)*2.0*(tsl(i,n+1)-tsl(i,n))                &
           /(dz(n+1)+dz(n))
      dhflux_dtsl1(i,n)=hc(i,n)*2.0/(dz(n+1)+dz(n))
      dhflux_dtsl2(i,n)=-hc(i,n)*2.0/(dz(n+1)+dz(n))
    END DO
  END DO
END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
!DIR$ IVDEP
!CDIR NODEP
DO j=1,soil_pts
  i=soil_index(j)
  h_flux(i,nshyd)=0.0
  h_flux(i,0) = surf_ht_flux(i)
  h_flux(i,1) = sifact(i)*h_flux(i,1)
  dhflux_dtsl1(i,nshyd)=0.0
  dhflux_dtsl2(i,nshyd)=0.0
  dhflux_dtsl1(i,0)=0.0
  dhflux_dtsl2(i,0)=0.0
  dhflux_dtsl1(i,1)=sifact(i)*dhflux_dtsl1(i,1)
  dhflux_dtsl2(i,1)=sifact(i)*dhflux_dtsl2(i,1)
END DO
!$OMP END DO

!--------------------------------------------------------------------
! Calculate the advection of heat by moisture fluxes
!--------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
DO jj=1, soil_pts, j_block
  DO n=2,nshyd-1
    DO j=jj,MIN(jj+j_block-1,soil_pts)
      i=soil_index(j)
      hadv(i,n)=hcapw*dz(n)*                                        &
           (w_flux(i,n-1)*(tsl(i,n-1)-tsl(i,n))/(dz(n)+dz(n-1))     &
           +w_flux(i,n)*(tsl(i,n)-tsl(i,n+1))/(dz(n)+dz(n+1)))
      dhadv_dtsl0(i,n)=hcapw*dz(n)*w_flux(i,n-1)/(dz(n)+dz(n-1))
      dhadv_dtsl1(i,n)=hcapw*dz(n)*                                 &
           (-w_flux(i,n-1)/(dz(n)+dz(n-1))                          &
           +w_flux(i,n)/(dz(n)+dz(n+1)))
      dhadv_dtsl2(i,n)=-hcapw*dz(n)*w_flux(i,n)/(dz(n)+dz(n+1))
    END DO
  END DO
END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
!DIR$ IVDEP
!CDIR NODEP
DO j=1,soil_pts
  i=soil_index(j)
  hadv(i,1)=hcapw*dz(1)*                                          &
   (w_flux(i,0)*(tsl(i,0)-tsl(i,1))/dz(1)                         &
   +w_flux(i,1)*(tsl(i,1)-tsl(i,2))/(dz(1)+dz(2)))
  dhadv_dtsl0(i,1)=0.0
  dhadv_dtsl1(i,1)=hcapw*dz(1)*                                   &
    (-w_flux(i,0)/dz(1)+w_flux(i,1)/(dz(1)+dz(2)))
  dhadv_dtsl2(i,1)=-hcapw*dz(1)*w_flux(i,1)/(dz(1)+dz(2))
  hadv(i,nshyd)=hcapw*dz(nshyd)*                                  &
               w_flux(i,nshyd-1)*(tsl(i,nshyd-1)-tsl(i,nshyd))    &
               /(dz(nshyd)+dz(nshyd-1))
  dhadv_dtsl0(i,nshyd)=hcapw*dz(nshyd)*w_flux(i,nshyd-1)          &
                       /(dz(nshyd)+dz(nshyd-1))
  dhadv_dtsl1(i,nshyd)=-hcapw*dz(nshyd)*w_flux(i,nshyd-1)         &
                       /(dz(nshyd)+dz(nshyd-1))
  dhadv_dtsl2(i,nshyd)=0.0
END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
DO jj=1, soil_pts, j_block
  DO n=1,nshyd
!DIR$ IVDEP
!CDIR NODEP
    DO j=jj,MIN(jj+j_block-1,soil_pts)
      i=soil_index(j)
!-----------------------------------------------------------------------
! Calculate TMAX, the temperature above which all soil water is
! unfrozen.
!-----------------------------------------------------------------------
      IF ( (v_sat(i,n) > tiny_0) .AND. (smcl(i,n) > small_value)) THEN
        work1(i,n)=(smcl(i,n)/smclsat(i,n))**(bexp(i,n))
        IF ( work1(i,n) > small_value ) THEN
          tmax(i,n)=-sathh(i,n)/(dpsidt*work1(i,n))
          tmax(i,n)=MAX(tmax(i,n),-zerodegc)
        ELSE
          tmax(i,n)=-zerodegc
        END IF
      ELSE
        tmax(i,n)=-zerodegc
      END IF
      
      dhsl0(i,n)=timestep*(h_flux(i,n-1)-h_flux(i,n)+hadv(i,n))
      
      dhsl(i,n)=dhsl0(i,n)
      
    END DO  !  j (points)
  END DO  !  n (layers)
END DO
!$OMP END DO

!$OMP END PARALLEL 
!-----------------------------------------------------------------------
! Iteration loop
!-----------------------------------------------------------------------
DO m=1,mmax

!-----------------------------------------------------------------------
! Define the array of points which fail to meet the flux criterion.
!-----------------------------------------------------------------------

  iter_pts=0

  DO j=1,soil_pts
    i=soil_index(j)

    IF (iter(i)) THEN

      iter_pts=iter_pts+1
      iter_index(iter_pts)=i

    END IF
    iter(i)=.FALSE.

  END DO




  IF ( iter_pts == 0 ) EXIT

!-----------------------------------------------------------------------
! Update calculations at these points.
!-----------------------------------------------------------------------

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(jj, j, i, gamcon, n)

!$OMP DO SCHEDULE(STATIC) 
  DO jj = 1, iter_pts
    DO n=1,nshyd
      
!CDIR NODEP
      DO j=jj,MIN(jj+j_block-1,iter_pts)
        i=iter_index(j)
        
        tsl0(i,n)=tsl(i,n)
        smclf0(i,n)=smclf(i,n)
        smclu0(i,n)=smclu(i,n)
        
        dtslmax(i,n)=1.0e4-tsl(i,n)
        dtslmin(i,n)=-zerodegc-tsl(i,n)
        
        IF (tsl(i,n) >  tmax(i,n)) THEN         ! All water unfrozen
          dthu(i,n)=0.0
          dtslmin(i,n)=tmax(i,n)-tsl(i,n)
        ELSE IF (tsl(i,n) == tmax(i,n).AND.                       &
             ! Remains unfrozen
             dhsl(i,n) >= 0.0) THEN
          dthu(i,n)=0.0
        ELSE                                    ! Phase changes
          dthu(i,n)=dpsidt*smclsat(i,n)                           &
               /(bexp(i,n)*sathh(i,n)*rho_water*dz(n))            &
               *(-dpsidt*tsl(i,n)/sathh(i,n))**(-1.0/bexp(i,n)-1.0)
          dtslmax(i,n)=tmax(i,n)-tsl(i,n)
        END IF
        
        hcapt(i,n)=hcap(i,n)+(hcapw-hcapi)*smclu(i,n)/dz(n)       &
             +hcapi*smcl(i,n)/dz(n)                               &
             +rho_water*dthu(i,n)*((hcapw-hcapi)*tsl(i,n)+lf)
        
        
!-----------------------------------------------------------------------
! Calculate the matrix elements required for the implicit update.
!-----------------------------------------------------------------------
        gamcon=GAMMA*timestep/(hcapt(i,n)*dz(n))
        a(i,n)=-gamcon*(dhflux_dtsl1(i,n-1)+dhadv_dtsl0(i,n))
        b(i,n)=1.0-gamcon*(dhflux_dtsl2(i,n-1)-dhflux_dtsl1(i,n)  &
             +dhadv_dtsl1(i,n))
        c(i,n)=gamcon*(dhflux_dtsl2(i,n)+dhadv_dtsl2(i,n))
        d(i,n)=1.0/(hcapt(i,n)*dz(n))*dhsl(i,n)
        
      END DO
    END DO
  END DO
!$OMP END DO

!$OMP MASTER
!-----------------------------------------------------------------------
! Solve the triadiagonal matrix equation.
!-----------------------------------------------------------------------
! DEPENDS ON: gauss
  CALL gauss(nshyd,npnts,iter_pts,iter_index,a,b,c,d              &
,            dtslmin,dtslmax,dtsl)
!$OMP END MASTER 
!$OMP BARRIER

!-----------------------------------------------------------------------
! Diagnose the implicit DHSL
!-----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC) 
  DO jj = 1, iter_pts, j_block
    DO n=2,nshyd-1
      DO j=jj,MIN(jj+j_block-1,iter_pts)
        i=iter_index(j)
        dhsl(i,n)=dhsl(i,n)-dz(n)*hcapt(i,n)*(a(i,n)*dtsl(i,n-1)  &
             +(b(i,n)-1)*dtsl(i,n)+c(i,n)*dtsl(i,n+1))
      END DO
    END DO
  END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
  DO j=1,iter_pts
    i=iter_index(j)
    dhsl(i,1)=dhsl(i,1)-dz(1)*hcapt(i,1)*(                        &
                    +(b(i,1)-1)*dtsl(i,1)+c(i,1)*dtsl(i,2))
    dhsl(i,nshyd)=dhsl(i,nshyd)-dz(nshyd)*hcapt(i,nshyd)*         &
    (a(i,nshyd)*dtsl(i,nshyd-1)+(b(i,nshyd)-1)*dtsl(i,nshyd))
  END DO
!$OMP END DO

!-----------------------------------------------------------------------
! Update the layer temperatures
!-----------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC) 
  DO jj = 1, iter_pts, j_block
    DO n=1,nshyd
!CDIR NODEP
      DO j=jj,MIN(jj+j_block-1,iter_pts)
        i=iter_index(j)
        
        tsl(i,n)=tsl(i,n)+dtsl(i,n)
        
!-----------------------------------------------------------------------
! If the temperature increment is small and frozen water exists
! assume that the excess energy goes into phase change
!-----------------------------------------------------------------------
        IF (ABS(dtsl(i,n)) <  tacur.AND.                          &
             tsl(i,n) <= tmax(i,n)) THEN
          dsmclf(i,n)=-dhsl(i,n)/((hcapw-hcapi)*tsl(i,n)+lf)
          dsmclf(i,n)=MAX(dsmclf(i,n),-smclf(i,n))
          dsmclf(i,n)=MIN(dsmclf(i,n),smclu(i,n))
          smclu(i,n)=smclu(i,n)-dsmclf(i,n)
          smclf(i,n)=smclf(i,n)+dsmclf(i,n)
        END IF
        
!-----------------------------------------------------------------------
! Diagnose unfrozen and frozen water contents
!-----------------------------------------------------------------------
        IF (tsl(i,n) >= tmax(i,n)) THEN
          smclu(i,n)=smcl(i,n)
          smclf(i,n)=0.0
        ELSE
          smclu(i,n)=smclsat(i,n)                                   &
               *(-dpsidt*tsl(i,n)/sathh(i,n))**(-1.0/bexp(i,n))
          smclf(i,n)=smcl(i,n)-smclu(i,n)
        END IF

!-----------------------------------------------------------------------
! Calculate the error in heat conservation
!-----------------------------------------------------------------------
        dsmclf(i,n)=smclf(i,n)-smclf0(i,n)
        dhsl(i,n)=dhsl(i,n)-(hcap(i,n)*dz(n)+hcapw*smclu0(i,n)      &
             +hcapi*smclf0(i,n))*dtsl(i,n)                          &
             -dsmclf(i,n)*((hcapi-hcapw)*tsl0(i,n)-lf)
        
!-----------------------------------------------------------------------
! Calculate the error in flux conservation
!-----------------------------------------------------------------------
        ceacur(i)=ABS(dhsl(i,n))/timestep

        IF (ceacur(i)  >   facur) THEN
          iter(i)=.TRUE.
        END IF
        
      END DO
    END DO
  END DO
!$OMP END DO 

!$OMP END PARALLEL
!-----------------------------------------------------------------------
! End of iteration loop
!-----------------------------------------------------------------------
END DO

!-----------------------------------------------------------------------
! Diagnose soil temperatures (K) and fractional values of unfrozen and
! frozen water.
!-----------------------------------------------------------------------

!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(jj, n,  &
!$OMP& j, i) SHARED(tsoil, soil_index, sthu, sthf, soil_pts,      &
!$OMP& nshyd, tsl, smclu, smclsat, smclf)
DO jj = 1, soil_pts, j_block
  DO n=1,nshyd
    DO j=jj,MIN(jj+j_block-1,soil_pts)
      i=soil_index(j)
      tsoil(i,n)=tsl(i,n)+zerodegc
      sthu(i,n)=smclu(i,n)/smclsat(i,n)
      sthf(i,n)=smclf(i,n)/smclsat(i,n)
    END DO
  END DO
END DO
!$OMP END PARALLEL DO

IF (lhook) CALL dr_hook('SOIL_HTC',zhook_out,zhook_handle)
RETURN
END SUBROUTINE soil_htc
