! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    SUBROUTINE ICE_HTC------------------------------------------------

! Description:
!     Updates deep soil temperatures for ice. No external subroutines
!     are called.

! Documentation : UM Documentation Paper 25

! Subroutine Interface:
SUBROUTINE ice_htc (                                              &
 npnts,nshyd,lice_pts,lice_index,dz                               &
,surf_ht_flux,timestep                                            &
,tsoil                                                            &
,ltimer                                                           &
)

USE snow_param, ONLY :                                            &
!      imported scalars with intent(in)
 snow_hcap,snow_hcon
USE switches, ONLY :                                              &
 l_land_ice_imp
USE soil_param, ONLY :                                            &
!      imported scalar parameters
   GAMMA=>gamma_t

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Subroutine arguments
!   Scalar arguments with intent(IN) :
INTEGER                                                           &
 lice_pts                                                         &
                      ! IN Number of land ice points.
,npnts                                                            &
                      ! IN Number of gridpoints.
,nshyd                ! IN Number of soil moisture levels.

REAL                                                              &
 timestep             ! IN Model timestep (s).


!   Array arguments with intent(IN) :
INTEGER                                                           &
 lice_index(npnts)    ! IN Array of ice points.

REAL                                                              &
 dz(nshyd)                                                        &
                      ! IN Thicknesses of the soil layers (m).
,surf_ht_flux(npnts)  ! IN Net downward surface heat flux (W/m2).

LOGICAL ltimer        ! Logical switch for TIMER diags

!   Array arguments with intent(INOUT) :
REAL                                                              &
 tsoil(npnts,nshyd)   ! INOUT Sub-surface temperatures (K).

! Local scalars:
INTEGER                                                           &
 i,j,n                ! WORK Loop counters.

! Local arrays:
REAL                                                              &
 h_flux(npnts,0:nshyd)! WORK The fluxes of heat between layers
!                           !      (W/m2).
!-----------------------------------------------------------------------
! Variables required for the implicit calculation.
!-----------------------------------------------------------------------
REAL                                                              &
 dhflux_dtsl1(npnts,0:nshyd)                                      &
                          ! WORK Rate of change of the explicit
!                         ! downward flux at the base of the layer
!                         ! with the temperature of the layer
!                         ! (W/m2/K).
,dhflux_dtsl2(npnts,0:nshyd)          &
                          ! WORK Rate of change of the explicit
!                         ! downward flux at the base of the layer
!                         ! with the temperature of the lower layer
!                         ! (W/m2/K).
,a(npnts,nshyd),b(npnts,nshyd),c(npnts,nshyd),d(npnts,nshyd)      &
!                         ! WORK Matrix elements.
,gamcon                                                           &
!                         ! WORK Forward timestep weighting constant.
,dtsoil(npnts,nshyd)                                                &
!                           ! WORK The increment to the ice 
!                           ! temperature (K/timestep).
,dtsoilmax(npnts,nshyd)                                             &
!                         ! WORK Maximum allowed increment 
!                         ! to soil temperature (K/timestep).
,dtsoilmin(npnts,nshyd)
!                         ! WORK Minimum allowed increment 
!                         ! to soil temperature (K/timestep).



INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('ICE_HTC',zhook_in,zhook_handle)


IF (l_land_ice_imp) THEN
!
! Implicit scheme modelled on soil_htc
!
!--------------------------------------------------------------------
! Calculate heat fluxes across layer boundaries and the derivatives
! of the fluxes with respect to the temperatures.
!--------------------------------------------------------------------
  DO n=1,nshyd-1
    DO j=1,lice_pts
      i=lice_index(j)
      h_flux(i,n)=-snow_hcon*2.0*(tsoil(i,n+1)-tsoil(i,n))              &
                              /(dz(n+1)+dz(n))
      dhflux_dtsl1(i,n)=snow_hcon*2.0/(dz(n+1)+dz(n))
      dhflux_dtsl2(i,n)=-snow_hcon*2.0/(dz(n+1)+dz(n))
    END DO
  END DO
!
  DO j=1,lice_pts
    i=lice_index(j)
    h_flux(i,0) = surf_ht_flux(i)
    dhflux_dtsl1(i,0)=0.0
    dhflux_dtsl2(i,0)=0.0
    h_flux(i,nshyd)=0.0
    dhflux_dtsl1(i,nshyd)=0.0
    dhflux_dtsl2(i,nshyd)=0.0
  END DO
!
!-----------------------------------------------------------------------
! Calculate the matrix elements required for the implicit update.
!-----------------------------------------------------------------------
  DO n=1, nshyd
    DO j=1, lice_pts
      i=lice_index(j)
      gamcon=GAMMA*timestep/(snow_hcap*dz(n))
      a(i,n)=-gamcon*dhflux_dtsl1(i,n-1)
      b(i,n)=1.0-gamcon*(dhflux_dtsl2(i,n-1)-dhflux_dtsl1(i,n))
      c(i,n)=gamcon*dhflux_dtsl2(i,n)
      d(i,n)=(-h_flux(i,n)+h_flux(i,n-1))*timestep/(snow_hcap*dz(n))
    END DO
  END DO

!-----------------------------------------------------------------------
! Solve the triadiagonal matrix equation.
!-----------------------------------------------------------------------
!
! Allow wide limits as we should not need them.
  dtsoilmin(:,:)=-1.0e4
  dtsoilmax(:,:)=1.0e4
!
! DEPENDS ON: gauss
  CALL gauss(nshyd,npnts,lice_pts,lice_index,a,b,c,d              &
,            dtsoilmin,dtsoilmax,dtsoil)

!-----------------------------------------------------------------------
! Update the layer temperatures
!-----------------------------------------------------------------------
  DO n=1,nshyd
!CDIR NODEP
    DO j=1,lice_pts
      i=lice_index(j)
      tsoil(i,n)=tsoil(i,n)+dtsoil(i,n)
    END DO
  END DO
!
ELSE
!
! Original explicit scheme
!
!--------------------------------------------------------------------
! Calculate heat fluxes across layer boundaries
!--------------------------------------------------------------------
  DO n=1,nshyd-1
    DO j=1,lice_pts
      i=lice_index(j)
      h_flux(i,n)=-snow_hcon*2.0*(tsoil(i,n+1)-tsoil(i,n))          &
                             /(dz(n+1)+dz(n))
    END DO
  END DO

!DIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
  DO j=1,lice_pts
    i=lice_index(j)
    h_flux(i,nshyd)=0.0
    h_flux(i,0)=surf_ht_flux(i)
  END DO

!--------------------------------------------------------------------
! Update the sub-surface temperatures
!--------------------------------------------------------------------
  DO n=1,nshyd
!CDIR NOVECTOR
!   CDIR$ IVDEP here would force vectorization but changes results!
    DO j=1,lice_pts
      i=lice_index(j)
  
      tsoil(i,n)=tsoil(i,n)                                         &
       +1.0/(snow_hcap*dz(n))*(h_flux(i,n-1)                        &
       -h_flux(i,n))*timestep

    END DO
  END DO

END IF

IF (lhook) CALL dr_hook('ICE_HTC',zhook_out,zhook_handle)
RETURN
END SUBROUTINE ice_htc
