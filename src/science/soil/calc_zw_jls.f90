
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!   Subroutine CALC_ZW-------------------------------------------------

!   Purpose: To calculate the mean water table depth from the soil
!            moisture deficit as described in Koster et al., 2000.,
!            using the Newton-Raphson method.

! Documentation: UNIFIED MODEL DOCUMENTATION PAPER NO 25

SUBROUTINE calc_zw(npnts,nshyd,soil_pts,soil_index                &
  ,bexp,sathh,smcl,smclzw,smclsat,smclsatzw,v_sat,zw)

USE c_densty
USE c_topog

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Subroutine arguments
!   Scalar arguments with intent(IN) :
INTEGER                                                           &
 npnts                                                            &
                      ! IN Number of gridpoints.
,nshyd                                                            &
                      ! IN Number of soil moisture levels.
,soil_pts             ! IN Number of soil points.

!   Array arguments with intent(IN) :
INTEGER                                                           &
 soil_index(npnts)    ! IN Array of soil points.

REAL                                                              &
 bexp(npnts)                                                      &
                      ! IN Clapp-Hornberger exponent.
,sathh(npnts)                                                     &
                      ! IN Saturated soil water pressure (m).
,smcl(npnts,nshyd)                                                &
                      ! IN Total soil moisture contents
!                           !      of each layer (kg/m2).
,smclsat(npnts,nshyd)                                             &
                      ! IN Soil moisture contents of each
!                           !      layer at saturation (kg/m2).
,smclzw(npnts)                                                    &
                      ! IN moisture content in deep layer.
,smclsatzw(npnts)                                                 &
                      ! IN moisture content in deep layer
!                           !      at saturation.
,v_sat(npnts)         ! IN Volumetric soil moisture.


!   Array arguments with intent(INOUT) :
REAL                                                              &
 zw(npnts)            ! INOUT Water table depth (m).

REAL                                                              &
 zw_old                                                           &
            ! WORK zw from last timestep
,fn                                                               &
            ! WORK
,dfn                                                              &
            ! WORK
,zwest                                                            &
            ! WORK
,smd                                                              &
            ! WORK soil moisture deficit
,psisat                                                           &
            ! WORK
,diffinter  ! WORK

INTEGER                                                           &
 i,j,n                                                            &
            ! Loop counters
,it                                                               &
            ! Loop counters
,niter                                                            &
            ! Number of iterations
,maxinter
PARAMETER(niter=3)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('CALC_ZW',zhook_in,zhook_handle)
DO j=1,soil_pts
i=soil_index(j)
  zw_old=zw(i)

! Calculate soil moisture deficit:
  smd=0.0
  DO n=1,nshyd
     smd=smd+(smclsat(i,n)-smcl(i,n))/rho_water
  END DO
  smd=smd+(smclsatzw(i)-smclzw(i))/rho_water
  psisat=-sathh(i)
  zwest=zw_old

  maxinter=0
  DO it=1,niter
     IF(zwest <  0)zwest=0.0
     IF(zwest >  zw_max)zwest=zw_max
     zw_old=zwest
     zwest=zw_old

! Newton-Raphson. zw(next)=zw-f(zw)/f'(zw).

     fn = zwest*v_sat(i) - smd -                                  &
        v_sat(i) * bexp(i)/(bexp(i)-1.0) * psisat *               &
        ( 1.0 - ( 1.0-zwest/psisat ) **                           &
                ( 1.0 - 1.0/bexp(i) ) )       !   f(zw)
     dfn =  v_sat(i) +v_sat(i) *                                  &
        ( 1.0-zwest/psisat ) ** (-1.0/bexp(i))

     IF(ABS(dfn) >  1.e-10)zwest=zwest-fn/dfn
     IF(zwest <  0)zwest=0.0
     IF(zwest >  zw_max)zwest=zw_max
     diffinter=ABS(zw_old-zwest)
  END DO

  IF(diffinter >  0.01)maxinter=maxinter+1
  zw(i)=zwest
  IF(zw(i) <  0)zw(i)=0.0
  IF(zw(i) >  zw_max)zw(i)=zw_max
END DO
IF (lhook) CALL dr_hook('CALC_ZW',zhook_out,zhook_handle)
RETURN

END SUBROUTINE calc_zw

