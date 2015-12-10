MODULE c_bvoc

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Parameter values for BVOC emission routine
! See Pacifico et al., 2011, Atmos. Chem. Phys.
!-----------------------------------------------------------------------
  REAL, PARAMETER ::                                                          &
    F_TMAX = 2.3,                                                             &
               ! Empirical estimate that makes isoprene emission level off
               ! at temperature close to 40 C (Arneth personal communication)
    ATAU = 0.1,                                                               &
               ! Scaling parameter (isoprene) (eq A4b Arneth et al., 2007)
    BTAU = 0.09,                                                              &
               ! Scaling parameter (mono-terpene, methanol, acetone)
               !  (eq 1 Schurges et al., 2009)
    T_REF = 303.15
               ! Temperature reference value (K) (eq A4b Arneth et al., 2007)

END MODULE c_bvoc
