! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    SUBROUTINE CALC_BASEFLOW_jules------------------------------------

! Description:
!     Calculates subsurface runoff (aka baseflow).

! Subroutine Interface:
SUBROUTINE calc_baseflow_jules(                                   &
 soil_pts,soil_index,npnts,nshyd                                  &
,zdepth,ksz                                                       &
,b,fexp,ti_mean,zw,sthf,sthu                                      &
,wutot,top_crit,qbase,qbase_l                                     &
 )

USE c_topog, ONLY :                                               &
! imported scalar parameters
    ti_max,zw_max

USE soil_param, ONLY :                                            &
! imported arrays with intent(in)
    dzsoil

USE switches, ONLY :                                              &
! imported scalars with intent(in)
    l_vg_soil, l_baseflow_corr

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

!-----------------------------------------------------------------------
! Scalar arguments with intent(in)
!-----------------------------------------------------------------------
INTEGER, INTENT(in) ::                                            &
 npnts                                                            &
                     ! Number of gridpoints.
,nshyd                                                            &
                     ! Number of soil moisture levels.
,soil_pts            ! Number of soil points.

!-----------------------------------------------------------------------
!   Array arguments with intent(IN) :
!-----------------------------------------------------------------------
INTEGER, INTENT(in) ::                                            &
 soil_index(npnts)   ! Array of soil points.

REAL, INTENT(in) ::                                               &
 b(npnts,nshyd)                                                   &
                     ! Clapp-Hornberger exponent.
,fexp(npnts)                                                      &
                     ! Decay factor in Sat. Conductivity
!                    ! in water table layer.
,ti_mean(npnts)                                                   &
                     ! Mean topographic index.
,zw(npnts)                                                        &
                     ! Water table depth (m).
,ksz(npnts,0:nshyd)  
                     ! Saturated hydraulic conductivity
!                    ! for each layer (kg/m2/s).

!-----------------------------------------------------------------------
! Array arguments with intent(inout)
!-----------------------------------------------------------------------
REAL, INTENT(inout) ::                                            &
 sthf(npnts,nshyd)                                                &
                     ! Frozen soil moisture content of each layer
!                    ! as a fraction of saturation.
,sthu(npnts,nshyd)   ! Unfrozen soil moisture content of each
!                    ! layer as a fraction of saturation.


!-----------------------------------------------------------------------
! Array arguments with intent(out)
!-----------------------------------------------------------------------
REAL, INTENT(out) ::                                              &
 qbase(npnts)                                                     &
                      ! Total base flow (kg/m2/s).
,qbase_l(npnts,nshyd+1)                                           &
!                     ! Base flow from each layer (kg/m2/s).
,top_crit(npnts)                                                  &
                      ! Critical topographic index required
!                     !     to calc surface saturation frac.
,wutot(npnts)         ! UNFROZEN to TOTAL fraction at ZW.

!-----------------------------------------------------------------------
! Local arrays:
!-----------------------------------------------------------------------
REAL ::                                                           &
 bracket(npnts,nshyd+1)                                           &
!                     ! work: 1-(1-S^(b+1))^(1/(b+1))
,ksfz(npnts,nshyd+1)                                              &
                      ! Function of sat. hydraulic conductivity,
!                     ! frozen soil and mean topographic index.
,qbase_max(npnts)                                                 &
                      ! Max possible base flow (kg/m2/s).
,qbase_max_l(npnts,nshyd+1)                                       &
!                     ! Max possible base flow
!                     ! from each level (kg/m2/s).
,zdepth(0:nshyd)                                                  &
                      ! Lower soil layer boundary depth (m).
,qbase_min(npnts)     ! Residual base flow at zw_max (kg/m2/s).

!-----------------------------------------------------------------------
! Local scalars:
!-----------------------------------------------------------------------
INTEGER ::                                                        &
 i,j                                                              &
                      ! Loop counters.
,n                                                                &
                      ! Tile loop counter.
,errorstatus          ! Error status for ereport.

! variables to limit printout
INTEGER fail_count, max_print
LOGICAL l_printout

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------
! Initialise TOP_CRIT to maximum value.
! Initialise baseflow components to zero.
!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook('CALC_BASEFLOW_JULES',zhook_in,zhook_handle)

DO j=1,soil_pts
  i=soil_index(j)
  top_crit(i) = ti_max
  qbase(i)    = 0.0
  qbase_max(i)= 0.0
  DO n=1,nshyd+1
    qbase_l(i,n)     = 0.0
    qbase_max_l(i,n) = 0.0
  END DO
END DO

DO j=1,soil_pts
  i=soil_index(j)

!-----------------------------------------------------------------------
! Initialise print options.
!-----------------------------------------------------------------------
  fail_count =  0
  max_print  = 10
  l_printout = .TRUE.

!-----------------------------------------------------------------------
! Calculate a layer-dependent variable which is dependent on
! effective saturated conductivity:
!-----------------------------------------------------------------------
  IF ( l_vg_soil .AND. l_baseflow_corr ) THEN

    DO n = 1 , nshyd
      IF ( sthf(i,n) <  1.0 ) THEN
        ! bracket = ( 1 - S^(b+1) ) ^ ( 1/(b+1) ) 
        bracket(i,n) = ( 1.0 - ( 1.0 - sthf(i,n) )**( b(i,n) + 1.0 ) ) &
                       ** ( 1.0/(b(i,n)+1.0) )
        ksfz(i,n) = 0.5 * ( ksz(i,n-1) + ksz(i,n) )                    &
            * ( (1.0 - sthf(i,n))**0.5 )                               &
            * ( 1.0 - bracket(i,n) ) * ( 1.0 - bracket(i,n) )          &
            * EXP( -ti_mean(i) )
      ELSE
        ksfz(i,n) = 0.0
      END IF
    END DO

    IF ( sthf(i,nshyd) <  1.0 ) THEN
      ! Using the params from the lowest layer, so bracket has already been
      ! calculated in the last step of the previous loop
      ksfz(i,nshyd+1) = ( ksz(i,nshyd) / fexp(i) )                     &
            * ( (1.0 - sthf(i,nshyd))**0.5 )                             &
            * ( 1.0 - bracket(i,nshyd) ) * ( 1.0 - bracket(i,nshyd) )  &
            * EXP( -ti_mean(i) )
    ELSE
      ksfz(i,nshyd+1) = 0.0
    END IF

  ELSE

!   l_vg_soil = .FALSE.
    DO n=1,nshyd
      IF ( sthf(i,n) <  1.0 ) THEN
        ksfz(i,n) = 0.5 * (ksz(i,n-1)+ksz(i,n))                   &
            * (1.0-sthf(i,n)) ** (2.0*b(i,n)+3.0) * EXP(-ti_mean(i))
      ELSE
        ksfz(i,n) = 0.0
      END IF
    END DO
    IF ( sthf(i,nshyd) <  1.0 ) THEN
      ksfz(i,nshyd+1) = ksz(i,nshyd) / fexp(i)                    &
                     * (1.0-sthf(i,nshyd)) ** (2.0*b(i,nshyd)+3.0)&
                     * EXP(-ti_mean(i))
    ELSE
      ksfz(i,nshyd+1) = 0.0
    END IF
  END IF  !  l_vg_soil

  IF ( EXP( -fexp(i)*(zw_max-zdepth(nshyd) ) ) >  0.05 ) THEN
    fail_count = fail_count + 1
    IF (l_printout) THEN
    WRITE(6,*)'CB_J: maximum water table depth is too low!!!'
    WRITE(6,*)'at soil point',i,fexp(i),zw_max,zdepth(nshyd)
    WRITE(6,*)EXP(-fexp(i)*(zw_max-zdepth(nshyd)))
    WRITE(6,*)'try zw_max>',-LOG(0.05)/fexp(i)+zdepth(nshyd)
    errorstatus=1000
!            CALL EREPORT('CALC_BASEFLOW', ERRORSTATUS,                  &
!     &        'ERROR ZW_MAX is TOO SMALL')
    END IF
  END IF
  IF (fail_count >= max_print) THEN
    l_printout = .FALSE.
  END IF

!-----------------------------------------------------------------------
! Calculate base flow between maximum allowed water table depth and
! "infinity":
!-----------------------------------------------------------------------
  qbase_min(i) = ksfz(i,nshyd+1)                                  &
              * EXP( -fexp(i) * (zw_max-zdepth(nshyd)) )

!-----------------------------------------------------------------------
! Calculate maximum possible and actual base flow for each layer:
!-----------------------------------------------------------------------

  DO n=1,nshyd

    qbase_max_l(i,n) = ksfz(i,n) * ( zdepth(n) - zdepth(n-1) )

    IF ( zw(i) <= zdepth(n-1) )                                   &
      qbase_l(i,n) = qbase_max_l(i,n)

    IF ( zw(i) <  zdepth(n) .AND. zw(i) >  zdepth(n-1) ) THEN
      qbase_l(i,n) = ksfz(i,n) * ( zdepth(n) - zw(i) )
      IF ( sthu(i,n)+sthf(i,n) >  0.0 )                           &
        wutot(i) = sthu(i,n) / ( sthu(i,n) + sthf(i,n) )
    END IF
    IF ( n == 1 .AND. zw(i) <  zdepth(n) ) THEN
      qbase_l(i,n) = ksfz(i,n) * ( zdepth(n) - zw(i) )
      IF ( sthu(i,n)+sthf(i,n) >  0.0 )                           &
        wutot(i) = sthu(i,n) / ( sthu(i,n) + sthf(i,n) )
    END IF

  END DO

  qbase_max_l(i,nshyd+1) = ksfz(i,nshyd+1) - qbase_min(i)

  IF ( zw(i) <= zdepth(nshyd) ) THEN
    qbase_l(i,nshyd+1) = qbase_max_l(i,nshyd+1)
  ELSE
    qbase_l(i,nshyd+1) = ksfz(i,nshyd+1)                          &
                      * EXP( -fexp(i) * (zw(i)-zdepth(nshyd)) )   &
                      - qbase_min(i)
    IF ( sthu(i,nshyd)+sthf(i,nshyd) >  0.0 )                     &
      wutot(i) = sthu(i,nshyd) / ( sthu(i,nshyd) + sthf(i,nshyd) )
  END IF

!-----------------------------------------------------------------------
! Calculate total possible and actual base flow:
!-----------------------------------------------------------------------
  DO n=1,nshyd+1
    qbase_l(i,n) = MAX( 0.0, qbase_l(i,n) )
    qbase(i)     = qbase(i) + qbase_l(i,n)
    qbase_max(i) = qbase_max(i) + qbase_max_l(i,n)
  END DO

!-----------------------------------------------------------------------
! Calculate critical topographic index.
!-----------------------------------------------------------------------
  IF(qbase(i) >  qbase_max(i)) qbase(i) = qbase_max(i)

! Check that QBASE_MAX(I)/QBASE(I) will not underflow.
  IF ( qbase_max(i) >  EPSILON(qbase_max(i)) .AND.                &
     qbase(i) >  qbase_max(i)*(EPSILON(qbase(i))) )               &
    top_crit(i) = LOG( qbase_max(i) / qbase(i) )

END DO

 IF (fail_count > max_print) THEN
   WRITE(6,*)'CB_J: ZW_MAX point-by-point warnings terminated.'
   WRITE(6,*)'CB_J: Total pts with ZW_MAX too small = ',fail_count
 END IF

IF (lhook) CALL dr_hook('CALC_BASEFLOW_JULES',zhook_out,zhook_handle)
RETURN
END SUBROUTINE calc_baseflow_jules
