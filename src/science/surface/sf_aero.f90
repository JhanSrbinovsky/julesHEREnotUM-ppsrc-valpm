! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!   SUBROUTINE SF_AERO------------------------------------------------

!  Purpose: Calculate coefficients required for CLASSIC aerosol code

!  Suitable for Single Column use.

!  Documentation: UM Documentation Paper No 24, section P243.
!                 See especially sub-section (ix).
!---------------------------------------------------------------------

! Arguments :-
SUBROUTINE sf_aero (                                              &
 land_pts,ntiles,land_index,tile_index,tile_pts,                  &
 l_aero_classic,l_dust,l_dust_diag,                               &
 flandg,tile_frac,pstar,rhostar,tstar,vshr_land,vshr_ssi,         &
 cd_ssi,ch_ssi,cd_std,ch_tile,rhokh_gb,                           &
 rho_aresist,aresist,resist_b,rho_aresist_tile,aresist_tile,      &
 resist_b_tile,r_b_dust,cd_std_dust,rhokh_mix                     &
 )


USE atm_fields_bounds_mod
USE theta_field_sizes, ONLY : t_i_length

USE dust_param, ONLY: ndiv

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

INTEGER, INTENT(IN) ::                                            &
 land_pts                                                         &
                       ! IN No of land points being processed.
,ntiles                                                           &
                       ! IN Number of land tiles per land point.
,land_index(land_pts)                                             &
                       ! IN Index of land points.
,tile_index(land_pts,ntiles)                                      &
                       ! IN Index of tile points.
,tile_pts(ntiles)      ! IN Number of tile points.


LOGICAL, INTENT(IN) ::                                            &
 l_aero_classic                                                   &
                       ! IN switch for using CLASSIC aerosol scheme
,l_dust                                                           &
                       ! IN switch for prognostic mineral dust
,l_dust_diag           ! IN switch for diagnostic mineral dust
                       !    lifting

REAL, INTENT(IN) ::                                               &
 flandg(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)      &
                       ! IN Land fraction on all tiles.
,tile_frac(land_pts,ntiles)                                       &
                       ! IN Tile fractions.
,pstar(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)       &
                       ! IN Surface pressure (Pascals).
,rhostar(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)     &
                       ! IN Surface air density
,tstar(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)       &
                       ! IN Gridbox Mean Surface Temperature (K)
,vshr_land(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)   &
                       ! IN Magnitude of land sfc-to-lowest-level
!                            !    wind shear
,vshr_ssi(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)    &
                       ! IN Mag. of mean sea sfc-to-lowest-level
!                            !    wind shear
,cd_ssi(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)      &
                       ! IN Bulk transfer coefficient for
!                            !      momentum over sea mean.
,ch_ssi(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)      &
                       ! IN Bulk transfer coefficient for heat
!                            !    and/or moisture over sea mean.
,cd_std(land_pts,ntiles)                                          &
                       ! IN Local drag coefficient for calc
!                            !    of interpolation coefficient
,ch_tile(land_pts,ntiles)                                         &
                       ! IN Transfer coefficient for heat and
!                            !    moisture
,rhokh_gb(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                       ! IN Grid-box surface exchange
!                            !     coefficients


!  Output variables.

REAL, INTENT(OUT) ::                                              &
 rho_aresist(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end) &
                       ! OUT RHOSTAR*CD_STD*VSHR 
                       !     for CLASSIC aerosol scheme
,aresist(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)     &
                       ! OUT 1/(CD_STD*VSHR) 
                       !     for CLASSIC aerosol scheme
,resist_b(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)    &
                       ! OUT (1/CH-1/CD_STD)/VSHR 
                       !     for CLASSIC aerosol scheme
,rho_aresist_tile(land_pts,ntiles)                                &
                       ! OUT RHOSTAR*CD_STD*VSHR on land tiles
                       !     for CLASSIC aerosol scheme
,aresist_tile(land_pts,ntiles)                                    &
                       ! OUT 1/(CD_STD*VSHR) on land tiles
                       !     for CLASSIC aerosol scheme
,resist_b_tile(land_pts,ntiles)                                   &
                       ! OUT (1/CH-1/CD_STD)/VSHR on land tiles
                       !     for CLASSIC aerosol scheme
,r_b_dust(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,ndiv)&
                       ! OUT surf layer res for dust
,cd_std_dust(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end) &
                       ! OUT Bulk transfer coef. for momentum,
                       !     excluding orographic effects
,rhokh_mix(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                       ! OUT Exchange coeffs for moisture.


! Local workspace
REAL ::                                                           &
 rho_aresist_land(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)&
                       ! Land mean of rho_aresist_tile
,vshr(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                       ! friction velocity passed to DUSTRESB

! Local scalars.
INTEGER ::                                                        &
 i,j                                                              &
             ! Loop counter (horizontal field index).
,k                                                                &
             ! Loop counter (tile field index).
,l                                                                &
             ! Loop counter (land point field index).
,n                                                                &
             ! Loop counter (tile index).
,idiv        ! Loop counter (dust division).

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('SF_AERO',zhook_in,zhook_handle)

!-----------------------------------------------------------------------
!!  Calculate surface layer resistance for CLASSIC aerosol
!-----------------------------------------------------------------------

DO j=tdims%j_start,tdims%j_end
  DO i=tdims%i_start,tdims%i_end
    rho_aresist(i,j) = 0.
    rho_aresist_land(i,j)=0.0
    aresist(i,j) = 0.
    resist_b(i,j) = 0.
  END DO
END DO

IF (l_aero_classic) THEN 
  DO j=tdims%j_start,tdims%j_end
    DO i=tdims%i_start,tdims%i_end
      IF ( flandg(i,j) < 1.0 ) THEN
        rho_aresist(i,j) = rhostar(i,j)*cd_ssi(i,j)*vshr_ssi(i,j)
        aresist(i,j) =  1. / (cd_ssi(i,j) * vshr_ssi(i,j))
        resist_b(i,j)= (cd_ssi(i,j)/ch_ssi(i,j) - 1.0) *          &
                        aresist(i,j)
      END IF
    END DO
  END DO
END IF !(l_aero_classic)

! Land tiles
DO n=1,ntiles
  DO l=1,land_pts
    rho_aresist_tile(l,n) = 0.
    aresist_tile(l,n) = 0.
    resist_b_tile(l,n) = 0.
  END DO
END DO

IF (l_aero_classic) THEN 
  DO n=1,ntiles
    DO k=1,tile_pts(n)
      l = tile_index(k,n)
      j=(land_index(l)-1)/t_i_length + 1
      i = land_index(l) - (j-1)*t_i_length
      rho_aresist_tile(l,n) = rhostar(i,j) * cd_std(l,n)          &
                  * vshr_land(i,j)
      aresist_tile(l,n) = 1. / ( cd_std(l,n) * vshr_land(i,j) )
      resist_b_tile(l,n) = ( cd_std(l,n)/ch_tile(l,n) - 1.0 ) *   &
                                                 aresist_tile(l,n)
      IF (resist_b_tile(l,n) < 0.) resist_b_tile(l,n) = 0.
      rho_aresist_land(i,j) = rho_aresist_land(i,j) +             &
                       tile_frac(l,n)*rho_aresist_tile(l,n)
    END DO
  END DO
  
  DO l=1,land_pts
    j=(land_index(l)-1)/t_i_length + 1
    i = land_index(l) - (j-1)*t_i_length
    rho_aresist(i,j) = flandg(i,j)*rho_aresist_land(i,j) +        &
         (1.0-flandg(i,j))*rho_aresist(i,j)
    aresist(i,j) = rhostar(i,j) / rho_aresist(i,j)
  END DO
END IF !(l_aero_classic)

!-----------------------------------------------------------------------
!!  Calculate surface layer resistance for mineral dust
!-----------------------------------------------------------------------

IF (l_dust.OR.l_dust_diag) THEN

  DO j = tdims%j_start,tdims%j_end
    DO i = tdims%i_start,tdims%i_end
      cd_std_dust(i,j)=(1.-flandg(i,j))*cd_ssi(i,j)
    END DO
  END DO
  
  DO n=1,ntiles
    DO k=1,tile_pts(n)
      l = tile_index(k,n)
      j=(land_index(l)-1)/t_i_length + 1
      i = land_index(l) - (j-1)*t_i_length
      cd_std_dust(i,j) = cd_std_dust(i,j) +                      &
           flandg(i,j)*tile_frac(l,n)*cd_std(l,n)
    END DO
  END DO
  
  DO j = tdims%j_start,tdims%j_end
    DO i = tdims%i_start,tdims%i_end
      DO idiv = 1,ndiv
        r_b_dust(i,j,idiv) = 0.
      END DO !IDIV
      vshr(i,j)=(1.0 - flandg(i,j)) * vshr_ssi(i,j) +            &
           flandg(i,j) * vshr_land(i,j)
    END DO !I
  END DO !J
  
! DEPENDS ON: dustresb
  CALL dustresb (                                                &
       pstar,tstar,rhostar,aresist,vshr,cd_std_dust,                &
       r_b_dust                                                     &
       )
  
END IF !(L_DUST.OR.L_DUST_DIAG)


!-----------------------------------------------------------------------
! 10 Set RHOKH, the coefficients required for tracer mixing.
!    Required 5B and after due to change in contents of RHOKH in rest
!    of routine.
!-----------------------------------------------------------------------

DO j=tdims%j_start,tdims%j_end
  DO i=tdims%i_start,tdims%i_end
    rhokh_mix(i,j) = rhokh_gb(i,j)
  END DO
END DO

IF (lhook) CALL dr_hook('SF_AERO',zhook_out,zhook_handle)
RETURN
END SUBROUTINE sf_aero
