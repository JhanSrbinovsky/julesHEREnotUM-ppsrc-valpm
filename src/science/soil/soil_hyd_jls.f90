! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    SUBROUTINE SOIL_HYD-----------------------------------------------

! Description:
!     Increments the layer soil moisture contents and calculates
!     calculates gravitational runoff.

! Documentation : UM Documentation Paper 25


! Subroutine Interface:
SUBROUTINE soil_hyd (npnts,nshyd,soil_pts,soil_index              &
,                    bexp,dz,ext,fw,ks,ksz,sathh,timestep,v_sat   &
,                    slow_runoff,smcl,sthu,surf_roff,w_flux       &
,                    stf_slow_runoff                              &
,                    zw,sthzw,zdepth,qbase,qbase_l                &
,                    dun_roff,drain,l_top,l_soil_sat_down         &
,                    ltimer                                       &
)

USE c_topog, ONLY :                                               &
!      imported scalar parameters
   zw_max

USE c_densty, ONLY :                                              &
!      imported scalar parameters
   rho_water  !  density of pure water (kg/m3)

USE soil_param, ONLY :                                            &
!      imported scalar parameters
   GAMMA=>gamma_w


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
,soil_pts             ! IN Number of soil points.

REAL                                                              &
 timestep             ! IN Model timestep (s).

LOGICAL                                                           &
 stf_slow_runoff                                                  &
                      ! IN Stash flag for sub-surface runoff.
,l_top                                                            &
                      ! IN Flag for TOPMODEL-based hydrology.
,l_soil_sat_down      ! IN Direction of super-sat soil moisture



!   Array arguments with intent(IN) :
INTEGER                                                           &
 soil_index(npnts)    ! IN Array of soil points.

REAL                                                              &
 bexp(npnts,nshyd)                                                &
                      ! IN Clapp-Hornberger exponent.
,dz(nshyd)                                                        &
                      ! IN Thicknesses of the soil layers (m).
,ext(npnts,nshyd)                                                 &
                      ! IN Extraction of water from each soil
!                           !    layer (kg/m2/s).
,fw(npnts)                                                        &
                      ! IN Throughfall from canopy plus snowmelt
!                           !    minus surface runoff (kg/m2/s).
,sathh(npnts,nshyd)                                               &
                      ! IN Saturated soil water pressure (m).
,v_sat(npnts,nshyd)                                               &
                      ! IN Volumetric soil moisture
!                           !    concentration at saturation
!                           !    (m3 H2O/m3 soil).
,ks(npnts)                                                        &
                      ! IN Saturated hydraulic
!                           !    conductivity at surface (kg/m2/s).
,ksz(npnts,nshyd)                                                 &
                      ! IN Saturated hydraulic conductivity
!                           !    in each soil layer (kg/m2/s).
,qbase(npnts)                                                     &
                      ! IN Base flow (kg/m2/s).
,qbase_l(npnts,nshyd+1)                                           &
!                           ! IN Base flow from each level (kg/m2/s).
!                           !    (kg/m2/s).
,zdepth(0:nshyd)                                                  &
                      ! IN Soil layer depth at lower boundary (m).
,zw(npnts)            ! IN Mean water table depth (m).


LOGICAL ltimer        ! Logical switch for TIMER diags

!   Array arguments with intent(OUT) :
REAL                                                              &
 slow_runoff(npnts)                                               &
                      ! OUT Drainage from the base of the
!                           !     soil profile (kg/m2/s).
,smclsat(npnts,nshyd)                                             &
                      ! OUT The saturation moisture content of
!                           !     each layer (kg/m2).
,w_flux(npnts,0:nshyd)                                            &
                      ! OUT The fluxes of water between layers
!                           !     (kg/m2/s).
,dun_roff(npnts)                                                  &
                      ! OUT Cumm. Dunne surface runoff (kg/m2/s).
,drain(npnts)         ! OUT Drainage out of nshyd'th level (kg/m2/s).

!   Array arguments with intent(INOUT) :
REAL                                                              &
 smcl(npnts,nshyd)                                                &
                      ! INOUT Total soil moisture contents
!                           !       of each layer (kg/m2).
,sthu(npnts,nshyd)                                                &
                      ! INOUT Unfrozen soil moisture content of ea
!                           !       layer as a fraction of saturation.
,surf_roff(npnts)                                                 &
                      ! INOUT Surface runoff (kg/m2/s).
,sthzw(npnts)         ! INOUT soil moist fraction in deep layer.


! Local scalars:
INTEGER                                                           &
 i,j,n                ! WORK Loop counters.

REAL                                                              &
 gamcon                                                           &
                      ! WORK Constant (s/mm).
,dw                                                               &
                      ! WORK
,dwzw                 ! WORK

! Local arrays:
REAL                                                              &
 a(npnts,nshyd)                                                   &
                      ! WORK Matrix elements corresponding
!                           !      to the coefficients of DSTHU(n-1).
,b(npnts,nshyd)                                                   &
                      ! WORK Matrix elements corresponding
!                           !      to the coefficients of DSTHU(n).
,c(npnts,nshyd)                                                   &
                      ! WORK Matrix elements corresponding
!                           !      to the coefficients of DSTHU(n+1).
,d(npnts,nshyd)                                                   &
                      ! WORK Matrix elements corresponding
!                           !      to the RHS of the equation.
,dsmcl(npnts,nshyd)                                               &
                      ! WORK Soil moisture increment
!                           !      (kg/m2/timestep).
,dsthu(npnts,nshyd)                                               &
                      ! WORK Increment to STHU (/timestep).
,dsthumin(npnts,nshyd)                                            &
                      ! WORK Minimum value of DSTHU.
,dsthumax(npnts,nshyd)                                            &
                      ! WORK Maximum value of DSTHU.
,dwflux_dsthu1(npnts,nshyd)                                       &
                            ! WORK The rate of change of the expli
!                           !      flux with STHU1 (kg/m2/s).
,dwflux_dsthu2(npnts,nshyd)                                       &
                            ! WORK The rate of change of the expli
!                           !      flux with STHU2 (kg/m2/s).
,smclu(npnts,nshyd)                                               &
                      ! WORK Unfrozen soil moisture contents
!                           !      of each layer (kg/m2).
,smclzw(npnts)                                                    &
                      ! WORK moisture content in deep layer.
,smclsatzw(npnts)     ! WORK moisture content in deep layer
!                           !      at saturation.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('SOIL_HYD',zhook_in,zhook_handle)

!-----------------------------------------------------------------------
! Calculate the unfrozen soil moisture contents and the saturation
! total soil moisture for each layer.
!-----------------------------------------------------------------------
DO n=1,nshyd
  DO j=1,soil_pts
    i=soil_index(j)
    smclsat(i,n)=rho_water*dz(n)*v_sat(i,n)
    smclu(i,n)=sthu(i,n)*smclsat(i,n)
    dsthumin(i,n)=-sthu(i,n)
    dsthumax(i,n)=1.0-smcl(i,n)/smclsat(i,n)
    dwflux_dsthu1(i,n)=0.0
    dwflux_dsthu2(i,n)=0.0
  END DO
END DO

!-----------------------------------------------------------------------
! Top boundary condition and moisture in deep layer:
!-----------------------------------------------------------------------
DO j=1,soil_pts
  i=soil_index(j)
  w_flux(i,0)=fw(i)
  smclsatzw(i)=rho_water*v_sat(i,nshyd)                           &
       *(zw_max-zdepth(nshyd))
  smclzw(i)=sthzw(i)*smclsatzw(i)

END DO

!-----------------------------------------------------------------------
! Calculate the Darcian fluxes and their dependencies on the soil
! moisture contents.
!-----------------------------------------------------------------------

! If L_VG_SOIL is T then Van Genuchten formulation is used otherwise
! Clapp Hornberger using Cosby parameters is used

! DEPENDS ON: hyd_con_ic
CALL hyd_con_ic (npnts,soil_pts,soil_index,bexp(:,nshyd)          &
,                ksz(:,nshyd),sthu(:,nshyd)                       &
,                w_flux(:,nshyd),dwflux_dsthu1(:,nshyd),ltimer)

DO n=2,nshyd
! DEPENDS ON: darcy_ic
  CALL darcy_ic (npnts,soil_pts,soil_index,bexp(:,n-1:n)          &
,                ksz(:,n-1),sathh(:,n-1:n)                        &
,                sthu(:,n-1),dz(n-1),sthu(:,n),dz(n)              &
,                w_flux(:,n-1)                                    &
,                dwflux_dsthu1(:,n-1),dwflux_dsthu2(:,n-1),ltimer)
END DO

!-----------------------------------------------------------------------
! Limit the explicit fluxes to prevent supersaturation in deep layer.
! Note that this w_flux can be overwritten by l_soil_sat_down loop.
!-----------------------------------------------------------------------
IF(l_top)THEN
   DO j=1,soil_pts
      i=soil_index(j)
      dwzw=(smclzw(i)-smclsatzw(i))/timestep+                     &
           (w_flux(i,nshyd)-qbase_l(i,nshyd+1))
      IF(dwzw >  0.0 ) THEN
         w_flux(i,nshyd)=w_flux(i,nshyd)-dwzw
      END IF
   END DO
END IF

!-----------------------------------------------------------------------
! Calculate the explicit increments.
! This depends on the direction in which moisture in excess of
! saturation is pushed (down if L_SOIL_SAT_DOWN, else up).
!-----------------------------------------------------------------------
IF (l_soil_sat_down) THEN
!-----------------------------------------------------------------------
! Moisture in excess of saturation is pushed down.
!-----------------------------------------------------------------------
  DO n=1,nshyd
    DO j=1,soil_pts
      i=soil_index(j)
      dsmcl(i,n)=(w_flux(i,n-1)-w_flux(i,n)-ext(i,n))*timestep
      IF(l_top) dsmcl(i,n)=dsmcl(i,n)-qbase_l(i,n)*timestep

!-----------------------------------------------------------------------
! Limit the explicit fluxes to prevent supersaturation.
!-----------------------------------------------------------------------
      IF (dsmcl(i,n) >  (smclsat(i,n)-smcl(i,n))) THEN
        dsmcl(i,n)=smclsat(i,n)-smcl(i,n)
        w_flux(i,n)=w_flux(i,n-1)-dsmcl(i,n)/timestep-ext(i,n)
        IF(l_top) w_flux(i,n)=w_flux(i,n)-qbase_l(i,n)
      END IF
!-----------------------------------------------------------------------
! Limit the explicit fluxes to prevent negative soil moisture.
! The flux out of a layer is reduced if necessary.
! This MAY no longer be required!!!!
!-----------------------------------------------------------------------
      IF ( l_top ) THEN
        IF ( smcl(i,n)+dsmcl(i,n) <  0.0 ) THEN
          dsmcl(i,n) = -smcl(i,n)
          w_flux(i,n)=w_flux(i,n-1)-ext(i,n)-qbase_l(i,n)         &
                     -dsmcl(i,n)/timestep
        END IF
      END IF

    END DO
  END DO

ELSE  !   NOT l_soil_sat_down
!-----------------------------------------------------------------------
! Moisture in excess of saturation is pushed up.
!-----------------------------------------------------------------------
  DO n=nshyd,1,-1
    DO j=1,soil_pts
      i=soil_index(j)
      dsmcl(i,n)=(w_flux(i,n-1)-w_flux(i,n)-ext(i,n))*timestep
      IF(l_top) dsmcl(i,n)=dsmcl(i,n)-qbase_l(i,n)*timestep

!-----------------------------------------------------------------------
! Limit the explicit fluxes to prevent supersaturation.
!-----------------------------------------------------------------------
      IF (dsmcl(i,n) >  (smclsat(i,n)-smcl(i,n))) THEN
        dsmcl(i,n)=smclsat(i,n)-smcl(i,n)
        w_flux(i,n-1)=dsmcl(i,n)/timestep+w_flux(i,n)+ext(i,n)
        IF(l_top) w_flux(i,n-1)=w_flux(i,n-1)+qbase_l(i,n)
      END IF
!-----------------------------------------------------------------------
! Limit the explicit fluxes to prevent negative soil moisture.
! The flux into a layer is increased if necessary.
! Note that we don't do this for N=1 because we can't increase the supply
! at the soil surface.
! This MAY no longer be required!!!!
!-----------------------------------------------------------------------
      IF ( l_top ) THEN
        IF ( smcl(i,n)+dsmcl(i,n) < 0.0 .AND. n>1 ) THEN
          dsmcl(i,n) = -smcl(i,n)
          w_flux(i,n-1)=w_flux(i,n)+ext(i,n)+qbase_l(i,n)         &
                       +dsmcl(i,n)/timestep
        END IF
      END IF

    END DO  !  j (points)
  END DO  !  n (layers)

END IF  !  l_soil_sat_down

!-----------------------------------------------------------------------
! Calculate the matrix elements required for the implicit update.
!-----------------------------------------------------------------------
DO j=1,soil_pts
  i=soil_index(j)
  gamcon=GAMMA*timestep/smclsat(i,1)
  a(i,1)=0.0
  b(i,1)=1.0+gamcon*dwflux_dsthu1(i,1)
  c(i,1)=gamcon*dwflux_dsthu2(i,1)
  d(i,1)=dsmcl(i,1)/smclsat(i,1)
END DO

DO n=2,nshyd
  DO j=1,soil_pts
    i=soil_index(j)
    gamcon=GAMMA*timestep/smclsat(i,n)
    a(i,n)=-gamcon*dwflux_dsthu1(i,n-1)
    b(i,n)=1.0-gamcon*(dwflux_dsthu2(i,n-1)-dwflux_dsthu1(i,n))
    c(i,n)=gamcon*dwflux_dsthu2(i,n)
    d(i,n)=dsmcl(i,n)/smclsat(i,n)
  END DO
END DO

!-----------------------------------------------------------------------
! Solve the triadiagonal matrix equation.
!-----------------------------------------------------------------------
! DEPENDS ON: gauss
CALL gauss(nshyd,npnts,soil_pts,soil_index,a,b,c,d                &
,          dsthumin,dsthumax,dsthu)

!-----------------------------------------------------------------------
! Diagnose the implicit fluxes.
!-----------------------------------------------------------------------
DO n=1,nshyd
  DO j=1,soil_pts
    i=soil_index(j)
    dsmcl(i,n)=dsthu(i,n)*smclsat(i,n)
    w_flux(i,n)=w_flux(i,n-1)-ext(i,n)-dsmcl(i,n)/timestep
    IF(l_top)w_flux(i,n)=w_flux(i,n)-qbase_l(i,n)
  END DO
END DO

IF(l_top)THEN
   DO j=1,soil_pts
      i=soil_index(j)
!-----------------------------------------------------------------------
! Limit implicit fluxes to prevent negative moisture in bottom layer.
!-----------------------------------------------------------------------
      IF(smcl(i,nshyd)+dsmcl(i,nshyd) <  0.0)THEN
         dsmcl(i,nshyd) = -smcl(i,nshyd)
         w_flux(i,nshyd) = w_flux(i,nshyd-1)                      &
                         - ext(i,nshyd)-qbase_l(i,nshyd)          &
                         - dsmcl(i,nshyd)/timestep
      END IF

!-----------------------------------------------------------------------
! Limit the implicit fluxes to prevent supersaturation in the deep layer.
! Adjust drainage flux out of layer above.
!-----------------------------------------------------------------------
      dwzw=(smclzw(i)-smclsatzw(i))/timestep+                     &
           (w_flux(i,nshyd)-qbase_l(i,nshyd+1))
      IF(dwzw >  0.0 ) THEN
         w_flux(i,nshyd) = w_flux(i,nshyd)-dwzw
         dsmcl(i,nshyd) = dsmcl(i,nshyd) + dwzw*timestep
      END IF
   END DO

!-----------------------------------------------------------------------
! Limit the implicit fluxes to prevent supersaturation in soil layers.
! Note that the form here effectively assumes l_soil_sat_down=FALSE.
!-----------------------------------------------------------------------
   DO n=nshyd,1,-1
      DO j=1,soil_pts
         i=soil_index(j)
         dw=(smcl(i,n)+dsmcl(i,n)-smclsat(i,n))/timestep
         IF (dw >= 0.0) THEN
            dsmcl(i,n)=smclsat(i,n)-smcl(i,n)
            w_flux(i,n-1)=w_flux(i,n-1)-dw
            IF(n /= 1)dsmcl(i,n-1)=dsmcl(i,n-1)+dw*timestep
         END IF
      END DO
   END DO
END IF  !  l_top

!-----------------------------------------------------------------------
! Update the prognostic variables.
!-----------------------------------------------------------------------
DO n=1,nshyd
  DO j=1,soil_pts
    i=soil_index(j)
    smclu(i,n)=smclu(i,n)+dsmcl(i,n)
    smcl(i,n)=smcl(i,n)+dsmcl(i,n)
    sthu(i,n)=smclu(i,n)/smclsat(i,n)
  END DO
END DO

IF(l_top)THEN

!-----------------------------------------------------------------------
! Diagnose the new water table depth.
! Assume local equilibrium psi profile
!-----------------------------------------------------------------------

  DO j=1,soil_pts
    i=soil_index(j)
    smclzw(i)=(smclzw(i)-(qbase_l(i,nshyd+1)-w_flux(i,nshyd))*    &
         timestep)
! Update prognostic deep layer soil moisture fraction:
    sthzw(i)=smclzw(i)/smclsatzw(i)
  END DO

! DEPENDS ON: calc_zw
  CALL calc_zw(npnts,nshyd,soil_pts,soil_index                    &
   ,bexp,sathh,smcl,smclzw,smclsat,smclsatzw,v_sat,zw)

!-----------------------------------------------------------------------
! Dont allow negative base flows:
!-----------------------------------------------------------------------
  DO j=1,soil_pts
    i=soil_index(j)
    qbase(i)=0.0
    DO n=1,nshyd+1
       qbase_l(i,n)=MAX(qbase_l(i,n),0.0)
       qbase(i)=qbase(i)+qbase_l(i,n)
    END DO
  END DO

END IF
!-----------------------------------------------------------------------
! Output slow runoff (drainage) diagnostic.
!-----------------------------------------------------------------------
IF (stf_slow_runoff) THEN
  DO i=1,npnts
    slow_runoff(i)=0.0
  END DO
  DO j=1,soil_pts
    i=soil_index(j)
! Ensure correct field is output as deep runoff: dependant on L_TOP
    IF(l_top)THEN
      slow_runoff(i)=qbase(i)
    ELSE
      slow_runoff(i)=w_flux(i,nshyd)
    END IF
    drain(i) = w_flux(i,nshyd)
  END DO
END IF

!-----------------------------------------------------------------------
! Update surface runoff diagnostic.
!-----------------------------------------------------------------------
DO j=1,soil_pts
  i=soil_index(j)
  surf_roff(i)=surf_roff(i)+(fw(i)-w_flux(i,0))
END DO

IF (lhook) CALL dr_hook('SOIL_HYD',zhook_out,zhook_handle)
RETURN
END SUBROUTINE soil_hyd
