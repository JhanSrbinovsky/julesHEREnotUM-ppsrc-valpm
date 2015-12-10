

SUBROUTINE init_urban (land_pts, frac, hgt_p, hwr_p, wrr_p, disp_p, ztm_p,    &
                       albwl_p, albrd_p, emisw_p, emisr_p)
! Description:
!   Routine to initialize URBAN parameters.
!   This is intended for testing only as it initialises the urban arrays using
!     set values
!
! (c) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! -----------------------------------------------------------------------------
!!!  6.1   10/01/07   First written. Peter Clark and Aurore Porson
!******************************************************************************
!
! Code Owner: See Unified Model Code Owner's HTML page
! This file belongs in section: Land

  USE urban_param, ONLY :                                                     &
     hgt, hwr, wrr, disp, ztm, albwl, albrd, emisw, emisr,                    &
     hgt_in, hwr_in, wrr_in, disp_in, ztm_in, albwl_in, albrd_in,             &
     emisw_in, emisr_in, a, cdz, kappa2, z0m_mat

  USE switches_urban, ONLY :                                                  &
     l_urban2T, l_moruses, l_moruses_albedo, l_moruses_emissivity,            &
     l_moruses_rough, l_moruses_storage, l_moruses_storage_thin,              &
     l_moruses_macdonald, l_urban_empirical

  USE nstypes, ONLY : urban, ice, urban_canyon, urban_roof, ntype

  USE parkind1, ONLY: jprb, jpim
  USE yomhook, ONLY: lhook, dr_hook
  USE PrintStatus_mod  
  IMPLICIT NONE

! Arguments:
  INTEGER, INTENT(IN) ::                                                      &
     land_pts                      ! Number of land points to be processed.
                                   ! include the nth surface type.

  REAL, INTENT(INOUT) ::                                                      &
     frac(land_pts,ntype),       & ! Fractional cover of each surface type.
     hgt_p(land_pts),            & ! From d1: Building height
     hwr_p(land_pts),            & ! From d1: Height to width
     wrr_p(land_pts),            & ! From d1: Width ratio
     albwl_p(land_pts),          & ! From d1: Wall albedo
     albrd_p(land_pts),          & ! From d1: Road albedo
     emisw_p(land_pts),          & ! From d1: Wall emmissivity
     emisr_p(land_pts),          & ! From d1: Road emmissivity
     disp_p(land_pts),           & ! From d1: Displacemnet height
     ztm_p(land_pts)               ! From d1: Roughness length

! Local declarations:
  REAL ::                                                                     &
     sc_hwr(land_pts),       & ! working variable
     d_h(land_pts),          & ! working variable
     lambdaf, lambdap          ! Frontal and planar area index

  INTEGER ::                                                                  &
     l                         ! WORK Loop counters

  REAL :: urban_fraction        ! Temporary store

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle


!----------------------------------------------------------------------
! Set parameters for urban morphology
!----------------------------------------------------------------------

  IF (lhook) CALL dr_hook('INIT_URBAN',zhook_in,zhook_handle)

! Initialising arrays
  hgt(:)   = 0.0
  hwr(:)   = 0.0
  wrr(:)   = 0.0
  albwl(:) = 0.0
  albrd(:) = 0.0
  emisw(:) = 0.0
  emisr(:) = 0.0
  ztm(:)   = 0.0
  disp(:)  = 0.0

! Check logic

! Initialise MORUSES switch
  l_moruses = .FALSE.

  IF ( l_moruses_albedo .OR. l_moruses_emissivity                       &
     .OR. l_moruses_rough .OR. l_moruses_storage ) THEN
    ! Turn on l_moruses if any of the independent switches used
    l_moruses = .TRUE. ! This does not mean that all moruses switches are true
                       ! This is used to set z0 to ztm in sparm.
    l_urban2T = .TRUE. ! MORUSES must be used with URBAN-2T
  END IF

  IF ( printstatus > PrStatus_Normal ) THEN
    WRITE(6,*) 'init_urban'
    WRITE(6,'(/,a)') 'Urban switches used'
    WRITE(6,*) 'l_urban2T             ', l_urban2T
    WRITE(6,*) 'l_moruses             ', l_moruses
    WRITE(6,*) 'l_moruses_albedo      ', l_moruses_albedo
    WRITE(6,*) 'l_moruses_emissivity  ', l_moruses_emissivity
    WRITE(6,*) 'l_moruses_rough       ', l_moruses_rough
    WRITE(6,*) 'l_moruses_storage     ', l_moruses_storage
    WRITE(6,*) 'l_moruses_storage_thin', l_moruses_storage_thin
    WRITE(6,*) 'l_moruses_macdonald   ', l_moruses_macdonald
    WRITE(6,*) 'l_urban_empirical     ', l_urban_empirical
    WRITE(6,*)

    ! Issue warnings on logic
    IF ( l_moruses_storage_thin .AND. .NOT. l_moruses_storage ) THEN
      WRITE(6,'(/,a)') 'WARNING: MORUSES storage parametrisation not used.'
      WRITE(6,*) 'l_moruses_storage      = .FALSE. when'
      WRITE(6,*) 'l_moruses_storage_thin = .TRUE.'
    END IF
    IF ( l_moruses_macdonald .AND. .NOT. l_moruses ) THEN
      WRITE(6,'(/,a)') 'WARNING: MORUSES is not switched on.'
      WRITE(6,*) 'l_moruses              = .FALSE. when'
      WRITE(6,*) 'l_moruses_macdonald    = .TRUE.'
    END IF
  END IF

! Ensure D1 array initialised then set if needed. D1 array not currently set
! through STASH
  wrr_p(:)    = 0.0
  hgt_p(:)    = 0.0
  hwr_p(:)    = 0.0
  disp_p(:)   = 0.0
  ztm_p(:)    = 0.0
  albwl_p(:)  = 0.0
  albrd_p(:)  = 0.0
  emisw_p(:)  = 0.0
  emisr_p(:)  = 0.0

  IF ( l_urban2t ) THEN
    wrr_p(:)      = wrr_in
    IF ( l_moruses ) THEN
      hgt_p(:)    = hgt_in
      hwr_p(:)    = hwr_in
      disp_p(:)   = disp_in
      ztm_p(:)    = ztm_in
      albwl_p(:)  = albwl_in
      albrd_p(:)  = albrd_in
      emisw_p(:)  = emisw_in
      emisr_p(:)  = emisr_in
    END IF
  END IF

  IF ( l_urban2T ) THEN

    IF ( printstatus > PrStatus_Normal ) THEN
      WRITE(6,*) "Setting URBAN-2T parameters"
    END IF

! Empirical relationships derived from correlating CEH urban fraction and
! LUCID urban geometry data for London. Obtained from collaboration with the
! University of Reading. See:
!     Bohnenstengel, S.I., Evans, S., Clark, P., Belcher, S.E. (2010);
!     Simulations of the London urban heat island, Q.J.R.Meteorol. Soc., to
!     be submitted.
! for more information

    IF ( l_urban_empirical ) THEN
      IF ( printstatus > PrStatus_Normal ) THEN
        WRITE(6,*) 'Using empirical relationships for urban geometry: wrr'
        IF ( l_moruses ) WRITE(6,*) 'Using empirical relationships for urban geometry: hwr'
      END IF
      DO l = 1, land_pts
        IF (frac(l,urban) > 0.0 .AND. frac(l,ice) == 0.0 ) THEN
          lambdap = 22.878*frac(l,urban)**6 - 59.473*frac(l,urban)**5      &
             + 57.749*frac(l,urban)**4 - 25.108*frac(l,urban)**3           &
             + 4.3337*frac(l,urban)**2 + 0.1926*frac(l,urban)              &
             + 0.036
          lambdaf = 16.412*frac(l,urban)**6 - 41.855*frac(l,urban)**5      &
             + 40.387*frac(l,urban)**4 - 17.759*frac(l,urban)**3           &
             + 3.2399*frac(l,urban)**2 + 0.0626*frac(l,urban)              &
             + 0.0271
          hwr(l) = 4.0 * atan(1.0)/2.0 * lambdaf / ( 1.0 - lambdap )
          wrr(l) = 1.0 - lambdap
        END IF
      END DO
      wrr_p(:) = wrr(:)
    ELSE
      IF ( printstatus > PrStatus_Normal ) THEN
        WRITE(6,'(/,a)') 'WARNING: Ancillary file used instead of empirical relationships for urban geometry (wrr)'
      END IF
! Fill allocatable arrays so that these can be passed in a module.
      wrr(:)= wrr_p(:)
    END IF

    IF ( l_moruses ) THEN
      IF ( printstatus > PrStatus_Normal ) THEN
        WRITE(6,*) "Setting MORUSES parameters"
      END IF

! First set MORUSES parameters with no other alternative parameterisations
      albwl(:) = albwl_p(:)
      albrd(:) = albrd_p(:)
      emisw(:) = emisw_p(:)
      emisr(:) = emisr_p(:)

      IF ( l_urban_empirical ) THEN
        IF ( printstatus > PrStatus_Normal ) THEN
          WRITE(6,*) 'Using empirical relationships for urban geometry: hgt'
        END IF
        DO l = 1, land_pts
          IF (frac(l,urban) > 0.0 .AND. frac(l,ice) == 0.0 ) THEN
            hgt(l) =                                                          &
                 167.409  * frac(l,urban)**5 - 337.853  * frac(l,urban)**4    &
               + 247.813  * frac(l,urban)**3 -  76.3678 * frac(l,urban)**2 &
               +  11.4832 * frac(l,urban)    +   4.48226
          END IF
        END DO
        hwr_p(:) = hwr(:)
        hgt_p(:) = hgt(:)
      ELSE
        IF ( printstatus > PrStatus_Normal ) THEN
          WRITE(6,*) 'WARNING: Ancillary file used instead of empirical relationship for urban geometry (hwr & hgt)'
        END IF
        hwr(:)   = hwr_p(:)
        hgt(:)   = hgt_p(:)
      END IF

      IF ( l_moruses_macdonald ) THEN
        !       Macdonald Formulation
        IF ( printstatus > PrStatus_Normal ) THEN
          WRITE(6,*) 'Using MacDonald formulation'
        END IF
        sc_hwr(:) = 0.5 * ( hwr(:) / (2.0 * ATAN(1.0)) )
        d_h(:)    = 1.0 - wrr(:) * ( a**(wrr(:) - 1.0) )
        disp(:)   = d_h(:) * hgt(:)
        DO l = 1, land_pts
          IF ( wrr(l) > 0.0 ) THEN
            ztm(l)    = (cdz * (1.0 - d_h(l)) *                               &
               sc_hwr(l) * wrr(l) / kappa2)**(-0.5)
            ztm(l)    = (1.0 - d_h(l))*EXP(-ztm(l))
            ztm(l)    = ztm(l) * hgt(l)
            ztm(l)    = MAX(ztm(l),z0m_mat)
          END IF
        END DO
        ztm_p(:)  = ztm(:)
        disp_p(:) = disp(:)
      ELSE
        IF ( printstatus > PrStatus_Normal ) THEN
          WRITE(6,'(/,a)') 'WARNING: Ancillary file used instead of MacDonald formulation'
        END IF
        ztm(:)  = ztm_p(:)
        disp(:) = disp_p(:)
      END IF


    END IF ! l_moruses

!------------------------------------------------------------------------------
! Expand urban tile to two tiles based on WRR
!------------------------------------------------------------------------------

    IF ( printstatus > PrStatus_Normal ) THEN
      WRITE(6,*) "Splitting urban tile to canyon/ roof & checking for land ice"
    END IF
    DO l = 1, land_pts
      IF (frac(l,urban) > 0.0 .AND. frac(l,ice) == 0.0 ) THEN
        IF ( printstatus > PrStatus_Normal ) THEN
          WRITE(6,'(a7,i9,9(2x,f5.3))') 'BEFORE:', l,frac(l,:)
        END IF
        urban_fraction  = frac(l,urban)
        frac(l,urban_canyon)   = urban_fraction * wrr(l)
        frac(l,urban_roof)     = urban_fraction - frac(l,urban_canyon)
        IF ( printstatus > PrStatus_Normal ) THEN
          WRITE(6,'(a7,i9,9(2x,f5.3))') 'AFTER :', l,frac(l,:)
        END IF
      ELSE IF ( frac(l,urban) > 0.0 .AND. frac(l,ice) > 0.0 ) THEN
        IF ( printstatus > PrStatus_Normal ) THEN
          WRITE(6,*) "Warning ice and urban co-exist for land point", l
        END IF
      END IF
    END DO

  ELSE    ! .NOT. l_urban2T

    IF ( printstatus > PrStatus_Normal ) THEN
      WRITE(6,*) 'URBAN-2T OR MORUSES not used'
      WRITE(6,*) 'All associated parameters initialised to zero'
    END IF

  END IF   ! l_urban2T

  IF (lhook) CALL dr_hook('INIT_URBAN',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE init_urban

