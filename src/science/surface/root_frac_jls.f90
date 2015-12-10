! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    SUBROUTINE ROOT_FRAC---------------------------------------------

! Subroutine Interface:
SUBROUTINE root_frac (nshyd,dz,rootd,f_root)


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

! Description:
!     Calculates the fraction of the total plant roots within each
!     soil layer.

! Documentation : UM Documentation Paper 25


! Subroutine arguments:
!   Scalar arguments with intent(IN) :
INTEGER                                                           &
 nshyd                ! IN Number of soil moisture layers.

REAL                                                              &
 dz(nshyd)                                                        &
                      ! IN Soil layer thicknesses (m).
,rootd                ! IN Rootdepth (m).

!   Array arguments with intent(OUT) :
REAL                                                              &
 f_root(nshyd)        ! OUT Fraction of roots in each soil
                      !     layer.
! Local scalars:
INTEGER                                                           &
 n                    ! WORK Loop counters

REAL                                                              &
 ftot                                                             &
                      ! WORK Normalisation factor.
,ztot                                                             &
                      ! WORK Total depth of soil (m).
,z1,z2                ! WORK Depth of the top and bottom of the
                      !      soil layers (m).

! Local parameters:
REAL                                                              &
 p                    ! WORK Power describing depth dependence
                      !           of the root density profile.
PARAMETER (p=1.0)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('ROOT_FRAC',zhook_in,zhook_handle)
z2=0.0
ztot=0.0

DO n=1,nshyd
  z1=z2
  z2=z2+dz(n)
  ztot=ztot+dz(n)
  f_root(n)=EXP(-p*z1/rootd)-EXP(-p*z2/rootd)
END DO

ftot=1.0-EXP(-p*ztot/rootd)
DO n=1,nshyd
  f_root(n)=f_root(n)/ftot
END DO

IF (lhook) CALL dr_hook('ROOT_FRAC',zhook_out,zhook_handle)
RETURN
END SUBROUTINE root_frac
