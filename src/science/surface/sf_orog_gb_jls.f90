! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  SUBROUTINE SF_OROG_GB --------------------------------------------
!
!  Purpose: Calculate effective roughness length and blending height
!
!--------------------------------------------------------------------
SUBROUTINE sf_orog_gb(                                            &
 land_pts,land_index,                                             &
 fd_stab_dep,orog_drag_param,                                     &
 ho2r2_orog,rib,sil_orog,z0m,z1,                                  &
 h_blend_orog,z0m_eff,sz0heff,z0h,z0h_eff                         &
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
                      ! IN Number of land points to be processed.
,land_index(land_pts)                                             &
                      ! IN Index of land points.
,fd_stab_dep          ! IN Switch to implement stability
!                           !    dependence of orographic form drag

LOGICAL                                                           &
 sz0heff
                      ! IN .TRUE. for Z0H_EFF diagnostic

REAL                                                              &
 ho2r2_orog(land_pts)                                             &
                      !IN Peak to trough height of unresolved
!                           !   orography divided by 2SQRT(2) (m).
,rib(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)         &
                      ! IN GBM Bulk Richardson number for lowest
!                           !    layer
,sil_orog(land_pts)                                               &
                      ! IN Silhouette area of unresolved orography
!                           !    per unit horizontal area
,orog_drag_param                                                  &
!                           ! IN Drag coefficient for orographic
!                           !    form drag
,z0m(land_pts)                                                    &
                      ! IN GBM Roughness length for momentum (m).
,z0h(land_pts)                                                    &
                      ! IN GBM Roughness length for heat (m)
,z1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                      ! IN Height of lowest atmospheric level (m).

!  Output variables.

REAL                                                              &
 h_blend_orog(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)&
!                           ! OUT Blending height
,z0m_eff(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)     &
,z0h_eff(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
!                           ! OUT Effective roughness lengths for
!                                 momentum, and heat and moisture (m)

!  Work Variables

INTEGER                                                           &
 i,j                                                              &
              ! Horizontal field index
,l            ! Land field index

REAL                                                              &
 rib_fn                                                           &
              ! Interpolation coefficient for 0 < RIB < RI_CRIT
,zeta1                                                            &
              ! Work space
,zeta2                                                            &
,zeta3                                                            &
,zeta4
              ! More work space

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('SF_OROG_GB',zhook_in,zhook_handle)

! Set blending height, effective roughness length and
! wind profile factor at land points.

DO l = 1,land_pts
  j=(land_index(l)-1)/t_i_length + 1
  i = land_index(l) - (j-1)*t_i_length

  h_blend_orog(i,j) = h_blend_min
  z0m_eff(i,j) = z0m(l)

  zeta1 = 0.5 * orog_drag_param * sil_orog(l)
  zeta2 = LOG ( 1.0 + z1(i,j)/z0m(l) )
  zeta3 = 1.0 / SQRT ( zeta1/(vkman*vkman) + 1.0/(zeta2*zeta2) )
  zeta2 = 1.0 / EXP(zeta3)
  h_blend_orog(i,j) = MAX ( z1(i,j) / (1.0 - zeta2) ,             &
                     ho2r2_orog(l) * SQRT(2.0) )
  h_blend_orog(i,j) = MIN ( h_blend_max, h_blend_orog(i,j) )

! Apply simple stability correction to form drag if RIB is stable

  IF (sil_orog(l)  ==  rmdi) THEN
    zeta1 = 0.0
  ELSE
    IF (fd_stab_dep == on) THEN
      rib_fn =  ( 1.0 - rib(i,j) / ri_crit )
      IF (rib_fn >  1.0) rib_fn = 1.0
      IF (rib_fn <  0.0) rib_fn = 0.0
      zeta1 = 0.5 * orog_drag_param * sil_orog(l) * rib_fn
    ELSE
      zeta1 = 0.5 * orog_drag_param * sil_orog(l)
    END IF
  END IF

  z0m_eff(i,j) = h_blend_orog(i,j) / EXP ( vkman / SQRT ( zeta1   &
           + (vkman / LOG ( h_blend_orog(i,j) / z0m(l) ) ) *      &
             (vkman / LOG ( h_blend_orog(i,j) / z0m(l) ) ) ) )
  IF ( rib(i,j) >  ri_crit .AND.                                  &
       fd_stab_dep == on ) z0m_eff(i,j) = z0m(l)

END DO

!     ! Calculate effective Z0H, if required

IF (sz0heff) THEN

  DO l = 1,land_pts
    j=(land_index(l)-1)/t_i_length + 1
    i = land_index(l) - (j-1)*t_i_length

    IF (sil_orog(l)  ==  rmdi) THEN
      zeta1 = 0.0
      zeta4 = 0.0
    ELSE
      IF (fd_stab_dep == on) THEN
        rib_fn =  ( 1.0 - rib(i,j) / ri_crit )
        IF (rib_fn >  1.0) rib_fn = 1.0
        IF (rib_fn <  0.0) rib_fn = 0.0
      ELSE
        rib_fn = 1.0
      END IF
      zeta1 = 0.5 * orog_drag_param * sil_orog(l) * rib_fn
      zeta4 = SQRT(1.0 + zeta1*LOG( h_blend_orog(i,j) / z0m(l) )  &
                              *LOG( h_blend_orog(i,j) / z0m(l) )  &
                              /(vkman*vkman) )
!           ! If the Hewer and Wood (1998) modification of scalar flux
!           ! parametrization were included:
!             ZETA4 = ZETA4 * ( 1.0 - 2.2*RIB_FN*SIL_OROG(L) )
    END IF

!         ! Rearranging Eq 148 of UM6.3 documentation gives:
    z0h_eff(i,j) = h_blend_orog(i,j) /                            &
           EXP ( zeta4 * LOG( h_blend_orog(i,j) / z0h(l) ) )
    IF ( rib(i,j) >  ri_crit .AND.                                &
         fd_stab_dep == on ) z0h_eff(i,j) = z0h(l)

  END DO

END IF

IF (lhook) CALL dr_hook('SF_OROG_GB',zhook_out,zhook_handle)
RETURN
END SUBROUTINE sf_orog_gb
