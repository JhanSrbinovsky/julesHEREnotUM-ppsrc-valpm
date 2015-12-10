SUBROUTINE urbanemis(hwr, emisr, emiss, emisw, sigmalw)

! Description:
!   Calculates the emissivity of the urban canyon
!
! (c) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! -----------------------------------------------------------------------------


  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim


  USE PrintStatus_mod

  IMPLICIT NONE

! Subroutine arguments
  REAL, INTENT(IN) ::                                                         &
     hwr,     & ! Canyon height-to-width ratio
     emisr,   & ! Emissivity road
     emiss,   & ! Emissivity sky
     emisw      ! Emissivity wall

  REAL, INTENT(OUT) ::                                                        &
     sigmalw    ! Effective emissivity

! Local declarations:

  INTEGER, PARAMETER :: arraydim = 3 ! Dimension of arrays RF and RFINV
  INTEGER ::  i, j                   ! Looping variables

  REAL ::                                                                     &
     sigmat,      &  ! Effective emissivity
     psiroad,     &  ! Factors for net_L
     psiwall,     &
     psisky,      &
     aa,bb,cc,dd, &   ! Elements of view factors array
     emisrratio,  &
     emiswratio

! Work variables - arrays
  REAL ::                                                                     &
     sf(arraydim,arraydim),     & ! Array of view factors
     delta(arraydim,arraydim),  & ! Delta function
     gammau(arraydim,arraydim)    ! Exchange cofficients

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle


!-----------------------------------

  IF (lhook) CALL dr_hook('URBANEMIS',zhook_in,zhook_handle)
  aa = (1.0 + hwr**2.0)**(0.5) - hwr
  bb = (1.0 + (1.0/hwr)**2.0)**(0.5) - (1.0/hwr)
  cc = (1.0 - aa)
  dd = (1.0 - bb) / 2.0

  !first the shape factors
  sf(1,1) = 0.0
  sf(1,2) = cc
  sf(1,3) = aa
  sf(2,1) = dd
  sf(2,2) = bb
  sf(2,3) = dd
  sf(3,1) = aa
  sf(3,2) = cc
  sf(3,3) = 0.0

  DO i = 1,3
    DO j = 1,3
      IF (i == 1) THEN
        sf(i,j) = sf(i,j) * (1.0 - emisr)
      END IF
      IF (i == 2) THEN
        sf(i,j) = sf(i,j) * (1.0 - emisw)
      END IF
      IF (i == 3) THEN
        sf(i,j) = sf(i,j) * (1.0 - emiss)
      END IF
    END DO
  END DO

!------------------------------------------------------------------------------
  DO i = 1,3
    DO j = 1,3
      IF (i == j) THEN
        delta(i,j) = 1.0
      ELSE
        delta(i,j) = 0.0
      END IF
    END DO
  END DO
!------------------------------------------------------------------------------

  gammau(:,:) = delta(:,:) - sf(:,:)

!------------------------------------------------------------------------------
! DEPENDS ON: matinv
  CALL matinv(gammau,arraydim)

! gamma is then the view factors
!------------------------------------------------------------------------------

  emisrratio = emisr/ (1.0 - emisr)
  emiswratio = emisw/ (1.0 - emisw)

  sigmalw = (emisrratio * gammau(1,3)) +                                      &
     (2.0 * hwr * emiswratio * gammau(2,3))
  sigmat = emisrratio * ((gammau(1,1) * emisr) +                              &
     (gammau(1,2) * emisw) - 1.0)
  sigmat = sigmat + (2.0 * hwr * emiswratio *                                 &
     ((gammau(2,1) * emisr) + (gammau(2,2) * emisw) - 1.0))

  psiroad = gammau(3,1)
  psiwall = gammau(3,2)
  psisky =  gammau(3,3)


  IF ( printstatus > PrStatus_Normal ) THEN
    IF ( ABS( sigmat + sigmalw ) > 1e-06 ) THEN
      WRITE(6,*) 'Problem balancing emissivity', ABS( sigmat + sigmalw )
      WRITE(6,'(a10,f8.3)') 'sigmat  = ',sigmat
      WRITE(6,'(a10,f8.3)') 'sigmalw = ',sigmalw
    END IF
  END IF








  IF (lhook) CALL dr_hook('URBANEMIS',zhook_out,zhook_handle)
  RETURN

END SUBROUTINE urbanemis
