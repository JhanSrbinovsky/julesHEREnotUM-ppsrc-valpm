! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    SUBROUTINE CALC_FSAT----------------------------------------------

! Description:
!     Calculates the surface saturation fraction.

! Subroutine Interface:
SUBROUTINE calc_fsat(l_gamtot,soil_pts,soil_index,                &
   npnts,ti_mean,ti_sig,wutot,top_crit,gamtot,fsat,fwetl)

USE c_topog

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

INTEGER                                                           &
 npnts                                                            &
                !IN No. of land points.
,soil_pts                                                         &
,soil_index(soil_pts)

REAL                                                              &
 ti_mean(npnts)                                                   &
                !IN Gridbox mean topographic index.
,ti_sig(npnts)                                                    &
                !IN Standard deviation in topographic index.
,gamtot(npnts)                                                    &
                !IN Integrated complete Gamma function.
!                     !   Unless L_GAMTOT=TRUE, then INOUT
,top_crit(npnts)                                                  &
                !IN Critical topographic index required
!                     !   to calculate the surface saturation fraction.
,wutot(npnts)   !IN UNFROZEN to TOTAL fraction at ZW.

LOGICAL                                                           &
 l_gamtot        !IN TRUE if need to calculate GAMTOT

REAL                                                              &
 fsat(npnts)                                                      &
                !INOUT Surface saturation fraction.
,fwetl(npnts)   !INOUT Wetland fraction.

REAL                                                              &
 ti_sig_use(npnts)                                                &
                 !WORK Standard deviation in topographic index
!                !     with a minimum value set.
,ti_max_use(npnts)                                                &
!                !WORK TI_MAX dependent on ti_sig
,dti_use(npnts)                                                   &
                 !WORK increment which dependent on ti_sig
,alf                                                              &
                 !WORK Parameter in incomplete Gamma function.
,alf_ksat                                                         &
                 !WORK Parameter in incomplete Gamma function
!                      !     including horizontal ksat variability.
,cum                                                              &
                 !WORK Integrated incomplete Gamma fn.
,ti_sc                                                            &
                 !WORK Incremented Topographic index.
!                      !     Scaled by ti_mean
,dti_sc                                                           &
                 !WORK Scale increment in TI
,ticr_sc                                                          &
                 !WORK Critical topographic index. (Any value
!                      !     above this gives surface saturation).
!                      !     Scaled by TI_MEAN
,ticr_sc_w                                                        &
                 !WORK As above, but for wetland fraction calc.
,fbox_s                                                           &
                 !WORK Fraction of box for integration.
,fbox_w          !WORK Fraction of box for integration.

INTEGER                                                           &
 i,j                                                              &
                 ! Loop counter for points.
,nti                                                              &
                 ! Loop counter for topographic index.
,mti                                                              &
                 ! Max loop count for topographic index.
,mti_w           ! Max loop count for topographic index.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!----------------------------------------------------------------------
! Set up maximum topographic index to integrate to and increment
! factor. These are topography dependent to maximise accuracy and
! minimise runtime:
!----------------------------------------------------------------------

IF (lhook) CALL dr_hook('CALC_FSAT',zhook_in,zhook_handle)

DO j=1,soil_pts
  i=soil_index(j)
! Set a minimum possible value for TI_SIG for safety.
  ti_sig_use(i)=MAX(ti_sig(i),0.5)
  ti_max_use(i)=ti_mean(i)+5.*ti_sig_use(i) 
  dti_use(i)=dti/(ti_sig_use(i)*ti_sig_use(i))
END DO

!----------------------------------------------------------------------
! Calculate the total integral under the Gamma function. Carried
! out via the reconfiguration:
!----------------------------------------------------------------------
IF(l_gamtot)THEN

  DO j=1,soil_pts
    i=soil_index(j)

    mti=NINT(ti_max_use(i)/dti_use(i))
    dti_sc=dti_use(i)/ti_mean(i)

    alf=(ti_mean(i)/ti_sig_use(i))**2
    alf_ksat=1./(1./alf+(sigma_logk/ti_mean(i))**2)

    DO nti=1,mti
      ti_sc=(nti-0.5)*dti_sc
      gamtot(i)=gamtot(i)                                         &
        +ti_sc**(alf_ksat-1.0)*EXP(-alf_ksat*ti_sc)
    END DO
  END DO
ELSE

!----------------------------------------------------------------------
! Calculate the integral under the incomplete Gamma function for the
! saturation surface fraction:
!----------------------------------------------------------------------
  DO j=1,soil_pts
    i=soil_index(j)

    IF(top_crit(i) <  ti_max_use(i))THEN

      ticr_sc=top_crit(i)/ti_mean(i)+1.0
      ticr_sc_w=ticr_sc+ti_wetl/ti_mean(i)

      IF(ticr_sc*ti_mean(i) <  ti_max_use(i))THEN

        alf=(ti_mean(i)/ti_sig_use(i))**2
        alf_ksat=1./(1./alf+(sigma_logk/ti_mean(i))**2)

        dti_sc=dti_use(i)/ti_mean(i)
        mti=INT(ticr_sc/dti_sc)

        cum=0.0
        DO nti=1,mti
          ti_sc=(nti-0.5)*dti_sc
          cum=cum+ti_sc**(alf_ksat-1.0)                           &
            *EXP(-alf_ksat*ti_sc)
        END DO

! Include the fractional increment:
        fbox_s=(ticr_sc-mti*dti_sc)/dti_sc
        ti_sc=(mti+0.5*fbox_s)*dti_sc
        cum=cum+ti_sc**(alf_ksat-1.0)                             &
          *EXP(-alf_ksat*ti_sc)*fbox_s

        fsat(i)=1.0-cum/gamtot(i)


!----------------------------------------------------------------------
! Calculate the integral under the incomplete Gamma function for the
! wetland fraction:
!----------------------------------------------------------------------
        IF(ticr_sc_w*ti_mean(i) <  ti_max_use(i))THEN

          mti_w=INT(ticr_sc_w/dti_sc)

! Include the fractional increment:
          IF(mti_w == mti)THEN
            fbox_w=(ticr_sc_w-mti*dti_sc)/dti_sc
          ELSE
            fbox_w=1.0
          END IF

          ti_sc=(mti+0.5*fbox_w)*dti_sc
          cum=cum+ti_sc**(alf_ksat-1.0)                           &
            *EXP(-alf_ksat*ti_sc)*(fbox_w-fbox_s)

          DO nti=mti+2,mti_w
            ti_sc=(nti-0.5)*dti_sc
              cum=cum+ti_sc**(alf_ksat-1.0)                       &
            *EXP(-alf_ksat*ti_sc)
          END DO

! Include the fractional increment:
          fbox_w=(ticr_sc_w-mti_w*dti_sc)/dti_sc
          ti_sc=(mti_w+0.5*fbox_w)*dti_sc
          cum=cum+ti_sc**(alf_ksat-1.0)                           &
          *EXP(-alf_ksat*ti_sc)*fbox_w

          fwetl(i)=-1.0+cum/gamtot(i)+fsat(i)
        END IF

        IF(fwetl(i) <  0.0)fwetl(i)=0.0
        IF(fsat(i) <  0.0)fsat(i)=0.0
        IF(fwetl(i) >  fsat(i))fwetl(i)=fsat(i)

!----------------------------------------------------------------------
! Assume that in partially frozen water, flow is inhibited, therefore
! the critical flow for wetland area is no longer valid:
!----------------------------------------------------------------------
        fwetl(i)=fwetl(i)*wutot(i)+fsat(i)*(1.0-wutot(i))

      END IF
    END IF

  END DO

END IF

IF (lhook) CALL dr_hook('CALC_FSAT',zhook_out,zhook_handle)
RETURN
END SUBROUTINE calc_fsat

!-----------------------------------------------------------------------
