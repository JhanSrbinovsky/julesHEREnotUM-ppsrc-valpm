SUBROUTINE urbanz0( n, z1_uv, z1_tq, hgt, hwr, disp, z0m, z0m_snow, ztm, zth )

! Description:
!   Calculates the effective roughness lengths for the urban tiles.
!   Note that the meaning of the input z0m for urban areas has changed
!   z0m in the ancillary is now the material roughness length
!
! (c) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.

  USE urban_param, ONLY : kappa2


  USE nstypes, ONLY : urban_canyon, urban_roof


  USE um_types
  USE parkind1, ONLY: jprb, jpim
  USE yomhook, ONLY: lhook, dr_hook


  USE PrintStatus_mod

  IMPLICIT NONE

! This routine considers that wind speeds tend to 0 when z is equal to d+z0
! this suggests that the interpolation log functions are
! written as log((z+z0)/z0) or log((z+ztm)/ztm)
! for respectively transfer from surface to internal
! boundary-layer and from internal to 1st atm. layer

  ! Subroutine arguments with intent(in)
  INTEGER ::           &
     n                   !index of tile

  REAL, INTENT(IN) ::  &
     z1_uv,            & ! Height of the level 1 winds
     z1_tq,            & ! Height of the level 1 t-q
     hgt,              & ! Height of buildings at each land point
     hwr,              & ! Canyon aspect ratio at each land point
     disp,             & ! Displacement height
     z0m,              & ! Material r.l. for momentum
                         ! This z0m is unaffected by snow and used
                         ! only on walls
     z0m_snow,         & ! This z0m can be affected by snow
     ztm                 ! Bulk r.l. for momentum

  REAL, INTENT(OUT) :: &
     zth                 ! Bulk r.l. for heat

  !Work variables - scalars
  INTEGER ::           & ! Looping variables
     callnum             ! Call number
  DATA  callnum /0/      ! Initialise to zero

  REAL ::             &
!     ures1_can,       & ! Resistances to transport in canyon
!     ures2_can,       &
!     ures3_can,       &
     ures1_roof,      & ! Resistances to transport in roof
     ures3_roof,      &
     urest_can,       & ! Bulk resistance to transport into bl
     urest_roof,      &
     u_us,            & ! Scaling for wind speed
                        ! (normalised reciprcal of friction velocity)
!    fzct,            &
     fz1_roof,        &
!    fz1_road,        &
!    fz2_wall,        &
     zref,            & ! Thickness of IBL formed along facets (in get_us also)
     roof_log           ! Precision for log infinity
!    canyon_log         ! Precision for log infinity

  REAL ( kind = real64 ) ::  & ! Do need the double precision (Aurore Porson)
     zth_h                     ! Normalised effective r.l. for total heat

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

! May need z0h_z0m = 0.1 as well

!------------------------------------------------------------------------------

  IF (lhook) CALL dr_hook('URBANZ0',zhook_in,zhook_handle)
  IF ( z0m - (1.0e-6) <= ztm ) THEN ! i.e. change if increases

    zref = 0.1 * hgt

    u_us = LOG( z1_uv / ztm + 1.0 ) / SQRT( kappa2 )

    IF ( n == urban_roof ) THEN

      roof_log = MAX( ( 1.1 * hgt - disp ), ztm )
      fz1_roof = LOG( roof_log / ztm ) / LOG( z1_uv / ztm + 1.0)

      ures1_roof = LOG( zref / z0m + 1.0 ) *                                  &
         LOG( ( zref + z0m ) / ( 0.1 * z0m ) )
      ures1_roof = ures1_roof / ( kappa2 * fz1_roof )
      ures3_roof = ( 1.0 - fz1_roof ) * ( u_us**2.0 )

      urest_roof = ures1_roof + ures3_roof

      zth_h = kappa2 * urest_roof / LOG( z1_uv / ztm + 1.0 )
      zth_h = ( z1_tq + ztm ) / ( EXP( zth_h ) * hgt )

      zth = zth_h * hgt
      zth = MAX( zth, 1.0e-30 )
      IF ( hwr < (1.0/3.0) ) THEN
        zth = MIN( zth, 0.1 * z0m )
      END IF

    ELSE IF ( n == urban_canyon ) THEN

      ! DEPENDS ON: get_us
      CALL get_us( 1.0, 0.0, hgt, hwr, z0m, ztm, z1_uv, disp,    &
         u_us, urest_can)

      zth_h = kappa2 * urest_can / LOG( z1_uv / ztm + 1.0)
      zth_h = ( z1_tq + ztm ) / (hgt * EXP(zth_h) )

      zth = zth_h * hgt
      zth = MAX( zth, 1.0e-30 )
      IF ( hwr < (1.0/3.0) ) THEN
        zth = MIN( zth, 0.1*z0m )
      END IF

    ELSE

      IF ( printstatus > PrStatus_Normal ) THEN
        WRITE(6,*) 'WARNING: Call to urbanz0 occuring not just for urban tiles'
      END IF



    END IF

!------------------------------------------------------------------------------
! Print warning to screen on first calling
!------------------------------------------------------------------------------

    IF ( callnum == 0 ) THEN

      IF ( printstatus > PrStatus_Normal ) THEN
        WRITE(6,*) 'MORUSES urbanz0: Altered roughness lengths for heat for urban tiles'
      END IF



      callnum = 1
    END IF

  ELSE


    IF ( printstatus > PrStatus_Normal ) THEN
      WRITE(6,*) "ERROR - wrong calling of resistance network"
      WRITE(6,*) 'z0m', z0m, 'ztm', ztm
    END IF




    IF (lhook) CALL dr_hook('URBANZ0',zhook_out,zhook_handle)
    RETURN
  END IF

  IF (lhook) CALL dr_hook('URBANZ0',zhook_out,zhook_handle)
  RETURN

END SUBROUTINE urbanz0
