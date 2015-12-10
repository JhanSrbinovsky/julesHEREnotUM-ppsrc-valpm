! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module for picking up BL, and other, options

! Description:
!   Permissible settings for BL options.
!   This module replaces blopt8a.h for settings required by land surface.
!   It also includes variables from bl_options_mod, bl_diags_mod 
!   and swapable_field_mod

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.

MODULE blopt8a

    IMPLICIT NONE

! Start blopt8a

! Description:
!   Permissible settings for BL options.


      INTEGER, PARAMETER :: off = 0  ! Switch disabled
      INTEGER, PARAMETER :: on  = 1  ! Switch enabled

!     Options for non-gradient stress following
      INTEGER, PARAMETER :: BrownGrant97 = 1
      INTEGER, PARAMETER :: BrownGrant97_limited = 2
!       Brown and Grant (1997), version 2 including a limit on its size

!     Options for flux gradient formulation
      INTEGER, PARAMETER :: Locketal2000   = 0
!       Flux gradients as in Lock et al. (2000)
      INTEGER, PARAMETER :: HoltBov1993 = 1
!       Flux gradients as in Lock et al (2000) but using
!       coefficients from Holtslag and Boville (1993)
      INTEGER, PARAMETER :: LockWhelan2006 = 2
!       Flux gradients as in Lock and Whelan (2006)

!     Options for entrainment enhancement in Sc over Cu
      INTEGER, PARAMETER :: Buoyrev_feedback = 1

!     Options for form drag
      INTEGER, PARAMETER :: No_drag         = 0
      INTEGER, PARAMETER :: Effective_z0    = 1
      INTEGER, PARAMETER :: Explicit_stress = 2

!     Options for marine boundary layers
      INTEGER, PARAMETER :: Fixed_Z0T = 0
!       Stanard flixed value of thermal roughness length over sea
      INTEGER, PARAMETER :: SurfDivZ0T = 1
!       Thermal roughness length over sea defined from surface
!       divergence theory
      INTEGER, PARAMETER :: DynDiag_ZL = 1
      INTEGER, PARAMETER :: DynDiag_ZL_corrn = 2
!       The ratio of the height of the inversion to the surface
!       Obukhov length is used as a dynamic criterion in the
!       diagnosis of BL types: version 2 includes changes to
!       cope with BL_LEVELS >> 3km
      INTEGER, PARAMETER :: DynDiag_ZL_CuOnly = 3
!       As 2 but only applied to points diagnosed with Cumulus 
!       and strictly for sea points (fland<0.01, cf 0.5)
      INTEGER, PARAMETER :: DynDiag_Ribased = 4
!       As 3 but also overrides Cumulus diagnosis if 
!          ZH(Ri) > ZLCL+zhloc_depth_fac*(ZHPAR-ZLCL)
!       Note that here Ri accounts for gradient adjustment by the 
!       non-local scheme.

!     Options for surface exchange
      INTEGER, PARAMETER :: Use_Correct_Ustar = 2
!       Option under the COR_MO_ITER switch for the dust scheme
!       to use the correct ustar
      INTEGER, PARAMETER :: Limit_ObukhovL = 3
!       Option under the COR_MO_ITER switch for imposing
!       lower limit on L in very stable conditions.
      INTEGER, PARAMETER :: Limit_expl_ustar = 2
!       Option under the COR_UST switch to limit the magnitude of the
!       explicitly calculated ustar
      INTEGER, PARAMETER :: IP_SrfExWithCnv = 1
!       Option to include deep convective gustiness in the surface
!       transfer

!     Options for convective boundary layers
      INTEGER, PARAMETER ::  UM_std     = 0
      INTEGER, PARAMETER ::  neut_cbl   = 1
      INTEGER, PARAMETER ::  LEM_conven = 2
      INTEGER, PARAMETER ::  LEM_std    = 3

!     Options for stable boundary layers
      INTEGER, PARAMETER ::  Long_tails           = 0
      INTEGER, PARAMETER ::  Sharpest             = 1
      INTEGER, PARAMETER ::  Sharp_sea_long_land  = 2
      INTEGER, PARAMETER ::  Mes_tails            = 3
      INTEGER, PARAMETER ::  Louis_tails          = 4
      INTEGER, PARAMETER ::  Depth_based          = 5
      INTEGER, PARAMETER ::  Sharp_sea_mes_land   = 6
      INTEGER, PARAMETER ::  LEM_stability        = 7
      INTEGER, PARAMETER ::  Sharp_sea_Louis_land = 8

!     Options for Prandtl number (in local Ri scheme)
      INTEGER, PARAMETER ::  Constant_SBL = 0
      INTEGER, PARAMETER ::  LockMailhot2004 = 1

! End blopt8a

!---------------------------------------------------------------
! Duplicating later BL options from bl_options_mod
! for use in the standalone JULES code.
!
! NOTE THAT DEFAULT VALUES SHOULD BE SUPPLIED HERE
!

  INTEGER :: ISrfExCnvGust = 0
!                      ! Switch to include the effect of convective
!                      ! downdraughts on surface exchange
!                      ! OFF (=0) => not used: only boundary-layer
!                      !   gustiness is considered (oiginal version)
!                      ! IP_SrfExWithCnv (=1) the impact of gustiness
!                      !   due to boundary layer eddies is reduced
!                      !   relative to the above, but eddies driven
!                      !   by convective downdraughts are included
  INTEGER :: Fric_heating = 0
!                      ! Switch to apply heating source from turbulence 
!                      ! dissipation
!                      ! OFF (=0) => not used
!                      ! ON  (=1) => used
!

  REAL :: Max_Stress_Grad = 0.05
!                      ! Maximum implied stress gradient across the 
!                      ! boundary layer, used to limit the explicit 
!                      ! stress applied in non-local scalings (m/s2)
!---------------------------------------------------------------
! From bl_diags_mod: only need dTfric which appears in imp_solver
      TYPE Strnewbldiag
        Logical :: L_dTfric
        Real, Pointer :: dTfric(:, :, :)
!                    188      Heating increment from turbulence dissipation
      END TYPE Strnewbldiag
!----------------------------------------------------
! From swapable_field_mod
  TYPE swapable_field_pointer_type
    INTEGER :: field_type
    INTEGER :: levels
    INTEGER :: rows
    LOGICAL :: vector
    REAL, POINTER :: field(:, :, :)
    REAL, POINTER :: field_2d(:, :)
  END TYPE swapable_field_pointer_type


END MODULE blopt8a
