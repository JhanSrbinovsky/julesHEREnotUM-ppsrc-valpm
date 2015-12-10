! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  SUBROUTINE SICE_HTF-----------------------------------------------
!
!  Purpose: Updates sea-ice surface layer temperature.
!
! Arguments:---------------------------------------------------------
SUBROUTINE sice_htf (                                             &
 row_length,rows,flandg,simlt,nice,nice_use,                      &
 di_ncat,ice_fraction,ice_fract_ncat,surf_ht_flux_sice_ncat,      &
 tstar_sice_cat,timestep,                                         &
 ti,ti_gb,sice_mlt_htf,sea_ice_htf,l_sice_heatflux,               &
 ltimer)

USE atm_fields_bounds_mod

USE c_mdi
USE c_sicehc
USE c_kappai
USE c_0_dg_c

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

LOGICAL ltimer

INTEGER                                                           &
 row_length                                                       &
                      ! IN Number of X points?
,rows                                                             &
                      ! IN Number of Y points?
,nice                                                             &
                      ! IN Number of ice categories
,nice_use             ! IN Number of ice categories used fully in
                      !    surface exchange

LOGICAL                                                           &
 simlt                      ! IN STASH flag for sea-ice melting
!                                 !    ht flux.

REAL                                                              &
 flandg(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)      &
                              ! IN Land fraction.
,di_ncat(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice)&
                               ! IN "Equivalent thickness" of
!                                   !    sea-ice (m).
,ice_fraction(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)&
                              ! IN Fraction of gridbox covered by
!                                   !    sea-ice.
,ice_fract_ncat(tdims%i_start:tdims%i_end,                        &
                tdims%j_start:tdims%j_end,nice)                   &
                                     ! IN Fraction of ice in
!                                 ! gridbox covered by each ice catagory
,surf_ht_flux_sice_ncat(tdims%i_start:tdims%i_end,                &
                        tdims%j_start:tdims%j_end,nice)           &
!                                   ! IN Net downward heat flux at
!                                   !    sea-ice surface W/m2
,tstar_sice_cat(tdims%i_start:tdims%i_end,                        &
                tdims%j_start:tdims%j_end,nice_use)               &
                              ! INOUT Sea-ice surface temperature (K).
,timestep                     ! IN Timestep (s).

REAL                                                              &
 ti(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice)     &
                          ! INOUT Sea-ice surface layer
!                                   !       temperature(K) Set to TSTAR
!                                   !       for unfrozen sea, missing
!                                   !       data for land.
,sice_mlt_htf(tdims%i_start:tdims%i_end,                          &
              tdims%j_start:tdims%j_end,nice)                     &
                              ! INOUT Heat flux due to melt
!                                   !       of sea-ice (W/m2).
,ti_gb(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)       &
                              ! OUT GBM ice surface temperature (K)
,sea_ice_htf(tdims%i_start:tdims%i_end,                           &
             tdims%j_start:tdims%j_end,nice)
                              ! OUT Heat flux through sea-ice
!                                   !   (W per sq m, positive downwards)

LOGICAL                                                           &
 l_sice_heatflux          ! IN T: semi-implicit update of TI

REAL                                                              &
 tsave                    ! Temporary temperature

REAL                                                              &
 ti_gb_local(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                          ! local equivalent to TI_GB

INTEGER i,j,n             ! Loop Counter; horizontal field index.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('SICE_HTF',zhook_in,zhook_handle)

DO j=tdims%j_start,tdims%j_end
  DO i=tdims%i_start,tdims%i_end
    IF (flandg(i,j) == 1.0) THEN
      DO n=1,nice
        sea_ice_htf(i,j,n) = 0.0
        ti(i,j,n) = rmdi
      END DO
      ti_gb_local(i,j) = rmdi
    ELSE IF (ice_fraction(i,j) <= 0.0) THEN
      DO n=1,nice
        sea_ice_htf(i,j,n) = 0.0
        ti(i,j,n) = tfs
      END DO
!      IF (nice .NE. nice_use) tstar_sice_cat(i,j,1) = tfs  ! set in sf_imp2
      ti_gb_local(i,j) = tfs
    ELSE
      ti_gb_local(i,j) = 0.0
      DO n=1,nice
        IF (ice_fract_ncat(i,j,n) > 0.0) THEN
          IF (l_sice_heatflux) THEN
            ! Semi-implicit update of TI
            tsave=ti(i,j,n)
            ti(i,j,n)=( ti(i,j,n) + ai * timestep * (             &
             surf_ht_flux_sice_ncat(i,j,n) +                      &
             ((tfs-tsave*0.5)*kappai/di_ncat(i,j,n)) )  )  /      &
             ( 1.0+kappai*ai*timestep*0.5/di_ncat(i,j,n) )
            sea_ice_htf(i,j,n) = kappai *                         &
                      ((ti(i,j,n)+tsave)*0.5 - tfs)/di_ncat(i,j,n)
          ELSE
            ! Original explicit update of TI
            ! (unstable for small ice thicknesses)
            sea_ice_htf(i,j,n) = kappai *                         &
                                 (ti(i,j,n) - tfs)/di_ncat(i,j,n)
            ti(i,j,n) = ti(i,j,n) + ai * timestep *               &
               (surf_ht_flux_sice_ncat(i,j,n) - sea_ice_htf(i,j,n))
          END IF
          IF ( ti(i,j,n) > tm ) THEN
            IF (simlt) sice_mlt_htf(i,j,n) = sice_mlt_htf(i,j,n)+ &
                        (ti(i,j,n) - tm)/(ai * timestep)
            ti(i,j,n) = tm
          END IF
          ti_gb_local(i,j) = ti_gb_local(i,j)                     &
        + (ice_fract_ncat(i,j,n) / ice_fraction(i,j)) * ti(i,j,n)
        ELSE
          sea_ice_htf(i,j,n) = 0.0
          ti(i,j,n) = tfs
        END IF
      END DO
    END IF
  END DO
END DO

!-----------------------------------------------------------------------
! When NICE=1 the pointers for TI_GB and TI are identical and TI=TI_GB
! Thus we need to ensure that neither is reassigned until the end of the
! summing of the gridbox mean. This does not occur when NICE  /=  1 as
! then the pointers are not equivalent.
!-----------------------------------------------------------------------
DO j=tdims%j_start,tdims%j_end
  DO i=tdims%i_start,tdims%i_end
! check for missed settings, possible in last bit of calc above?
    IF(ti_gb_local(i,j)==0.0)THEN
      PRINT*,'SICE_HTF: TIDYING UP MISSED VALUES OF TI_GB_LOCAL'
      ti_gb_local(i,j)=tfs
    END IF
    ti_gb(i,j) = ti_gb_local(i,j)
  END DO
END DO

IF (lhook) CALL dr_hook('SICE_HTF',zhook_out,zhook_handle)
RETURN
END SUBROUTINE sice_htf
