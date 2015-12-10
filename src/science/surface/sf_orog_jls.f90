! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  SUBROUTINE SF_OROG------------------------------------------------
!
!  Purpose: Calculate roughness lengths, blending height and wind
!           profile factor
!
!--------------------------------------------------------------------
SUBROUTINE sf_orog(                                               &
 land_pts,tile_pts,land_index,tile_index,                         &
 fd_stab_dep,orog_drag_param,                                     &
 ho2r2_orog,rib,sil_orog,z0m,z1,                                  &
 wind_profile_factor,z0m_eff                                      &
 )

USE atm_fields_bounds_mod
USE theta_field_sizes, ONLY : t_i_length

USE c_surf, ONLY : ri_crit
USE c_vkman
USE c_mdi
USE surf_param, ONLY : h_blend_max,h_blend_min

USE bl_option_mod, ONLY : on

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

INTEGER                                                           &
 land_pts                                                         &
                      ! IN Number of land points.
,tile_pts                                                         &
                      ! IN Number of tile points.
,land_index(land_pts)                                             &
                      ! IN Index of land points.
,tile_index(land_pts)                                             &
                      ! IN Index of tile points
,fd_stab_dep          ! IN Switch to implement stability
!                           !    dependence of orographic form drag

REAL                                                              &
 ho2r2_orog(land_pts)                                             &
                      !IN Peak to trough height of unresolved
!                           !   orography divided by 2SQRT(2) (m).
,rib(land_pts)                                                    &
                      ! IN Bulk Richardson number for lowest layer
,sil_orog(land_pts)                                               &
                      ! IN Silhouette area of unresolved orography
!                           !    per unit horizontal area
,orog_drag_param                                                  &
!                           ! IN Drag coefficient for orographic
!                           !    form drag
,z0m(land_pts)                                                    &
                      ! IN Roughness length for momentum (m).
,z1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
!                     ! IN Height of lowest atmospheric level (m).

!  Output variables.

REAL                                                              &
 wind_profile_factor(land_pts)                                    &
!                           ! OUT For transforming effective surface
!                           !     transfer coefficients to those
!                           !     excluding form drag.
,z0m_eff(land_pts)    ! OUT Effective roughness length for
!                                 momentum (m)

!  Work Variables

INTEGER                                                           &
 i,j                                                              &
              ! Horizontal field index
,k                                                                &
              ! Tile field index
,l            ! Land field index

REAL                                                              &
 h_blend_orog                                                     &
              ! Blending height
,rib_fn                                                           &
              ! Interpolation coefficient for 0 < RIB < RI_CRIT
,zeta1                                                            &
              ! Work space
,zeta2                                                            &
              ! More work space
,zeta3        ! Even more work space

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('SF_OROG',zhook_in,zhook_handle)

! Set blending height, effective roughness length and
! wind profile factor at land points.
!CDIR NODEP
DO k=1,tile_pts
  l = tile_index(k)
  j=(land_index(l)-1)/t_i_length + 1
  i = land_index(l) - (j-1)*t_i_length

  zeta1 = 0.5 * orog_drag_param * sil_orog(l)
  zeta2 = LOG ( 1.0 + z1(i,j)/z0m(l) )
  zeta3 = 1.0 / SQRT ( zeta1/(vkman*vkman) + 1.0/(zeta2*zeta2) )
  zeta2 = 1.0 / EXP(zeta3)
  h_blend_orog = MAX ( z1(i,j) / (1.0 - zeta2) ,                  &
                     ho2r2_orog(l) * SQRT(2.0) )
  h_blend_orog = MIN ( h_blend_max, h_blend_orog )

! Apply simple stability correction to form drag if RIB is stable

  IF (sil_orog(l)  ==  rmdi) THEN
    zeta1 = 0.0
  ELSE
    IF (fd_stab_dep == on) THEN
      rib_fn =  ( 1.0 - rib(l) / ri_crit )
      IF (rib_fn >  1.0) rib_fn = 1.0
      IF (rib_fn <  0.0) rib_fn = 0.0
      zeta1 = 0.5 * orog_drag_param * sil_orog(l) * rib_fn
    ELSE
      zeta1 = 0.5 * orog_drag_param * sil_orog(l)
    END IF
  END IF

  z0m_eff(l) = h_blend_orog / EXP ( vkman / SQRT ( zeta1 +        &
               (vkman / LOG ( h_blend_orog / z0m(l) ) ) *         &
               (vkman / LOG ( h_blend_orog / z0m(l) ) ) ) )
  IF ( rib(l) >  ri_crit .AND.                                    &
       fd_stab_dep == on )  z0m_eff(l)=z0m(l)

  wind_profile_factor(l) = LOG( h_blend_orog / z0m_eff(l) ) /     &
                           LOG( h_blend_orog / z0m(l) )

END DO

IF (lhook) CALL dr_hook('SF_OROG',zhook_out,zhook_handle)
RETURN
END SUBROUTINE sf_orog
