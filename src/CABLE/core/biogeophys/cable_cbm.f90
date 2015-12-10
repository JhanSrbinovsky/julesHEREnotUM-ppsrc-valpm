!==============================================================================
! This source code is part of the 
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CABLE Academic User Licence Agreement 
! (the "Licence").
! You may not use this file except in compliance with the Licence.
! A copy of the Licence and registration form can be obtained from 
! http://www.cawcr.gov.au/projects/access/cable
! You need to register and read the Licence agreement before use.
! Please contact cable_help@nf.nci.org.au for any questions on 
! registration and the Licence.
!
! Unless required by applicable law or agreed to in writing, 
! software distributed under the Licence is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the Licence for the specific language governing permissions and 
! limitations under the Licence.
! ==============================================================================
!
! Purpose: Calls CABLE routines including define_air, surface_albedo, 
!          define_canopy, soilsnow, carbon
!          Note that cbm is called once per timestep in the offline case but
!          twice per timestep in the ACCESS case. Not all parts of cbm 
!          are executed in each of the ACCESS calls.
!          
! Called from: cable_driver for offline version
!              cable_explicit_driver, cable_implicit_driver for ACCESS
!
! Contact: Yingping.Wang@csiro.au
!
! History: Calling sequence changes for ACCESS compared to v1.4b
!
!
! ==============================================================================

!#define NO_CASA_YET 1

MODULE cable_cbm_module
   
   USE cable_canopy_module
   USE cable_albedo_module
  
   IMPLICIT NONE
  
   PRIVATE
   PUBLIC cbm 

CONTAINS

   SUBROUTINE cbm( dels, air, bgc, canopy, met,                                &
                   bal, rad, rough, soil,                                      &
                   ssnow, sum_flux, veg )
    
   USE cable_common_module
   USE cable_carbon_module
   USE cable_soil_snow_module
   USE cable_def_types_mod
   USE cable_roughness_module
   USE cable_radiation_module
   USE cable_air_module

   USE casadimension,     only : icycle ! used in casa_cnp

   USE cable_data_module, ONLY : icbm_type, point2constants 
   
   !ptrs to local constants 
   TYPE( icbm_type ) :: C
   ! CABLE model variables
   TYPE (air_type),       INTENT(INOUT) :: air
   TYPE (bgc_pool_type),  INTENT(INOUT) :: bgc
   TYPE (canopy_type),    INTENT(INOUT) :: canopy
   TYPE (met_type),       INTENT(INOUT) :: met
   TYPE (balances_type),  INTENT(INOUT) :: bal
   TYPE (radiation_type), INTENT(INOUT) :: rad
   TYPE (roughness_type), INTENT(INOUT) :: rough
   TYPE (soil_snow_type), INTENT(INOUT) :: ssnow
   TYPE (sum_flux_type),  INTENT(INOUT) :: sum_flux
    
   TYPE (soil_parameter_type), INTENT(INOUT)   :: soil 
   TYPE (veg_parameter_type),  INTENT(INOUT)    :: veg  

   REAL, INTENT(IN)               :: dels ! time setp size (s)
    
   INTEGER :: k,kk,j  






   ! assign local ptrs to constants defined in cable_data_module
   CALL point2constants(C)    

   IF( cable_runtime%um ) THEN
      
      cable_runtime%um_radiation = .FALSE.
      
      IF( cable_runtime%um_explicit ) THEN
         CALL ruff_resist(veg, rough, ssnow, canopy)
         met%tk = met%tk + C%grav/C%capp*(rough%zref_tq + 0.9*rough%z0m)
      ENDIF
      
      CALL define_air (met, air)
   
   ELSE
      call ruff_resist(veg, rough, ssnow, canopy)
   ENDIF


   CALL init_radiation(met,rad,veg, canopy) ! need to be called at every dt

   print *, "jhan:cbm: cable_runtime%um ", cable_runtime%um
      print *, "jhan:cbm: cable_runtime%explicit ", cable_runtime%um_explicit
   IF( cable_runtime%um ) THEN
      
      IF( cable_runtime%um_explicit ) THEN
         CALL surface_albedo(ssnow, veg, met, rad, soil, canopy)
      ENDIF
   
   ELSE
      CALL surface_albedo(ssnow, veg, met, rad, soil, canopy)
   ENDIf
    
   ! Calculate canopy variables:
   CALL define_canopy(bal,rad,rough,air,met,dels,ssnow,soil,veg, canopy)

   ssnow%otss_0 = ssnow%otss
   ssnow%otss = ssnow%tss
   ! RML moved out of following IF after discussion with Eva
   ssnow%owetfac = ssnow%wetfac

   IF( cable_runtime%um ) THEN
      
     IF( cable_runtime%um_implicit ) THEN
         CALL soil_snow(dels, soil, ssnow, canopy, met, bal,veg)
      ENDIF

   ELSE
      call soil_snow(dels, soil, ssnow, canopy, met, bal,veg)
   ENDIF

   ssnow%deltss = ssnow%tss-ssnow%otss
   ! correction required for energy balance in online simulations
   IF( cable_runtime%um ) THEN
   
      canopy%fhs = canopy%fhs + ( ssnow%tss-ssnow%otss ) * ssnow%dfh_dtg
      
      canopy%fhs_cor = canopy%fhs_cor + ( ssnow%tss-ssnow%otss ) * ssnow%dfh_dtg
      
      canopy%fh = canopy%fhv + canopy%fhs

   canopy%fes = canopy%fes + ( ssnow%tss-ssnow%otss ) *                        &
                ( ssnow%dfe_ddq * ssnow%ddq_dtg )
                !( ssnow%cls * ssnow%dfe_ddq * ssnow%ddq_dtg )
   
   canopy%fes_cor = canopy%fes_cor + ( ssnow%tss-ssnow%otss ) *                &
                    ( ssnow%cls * ssnow%dfe_ddq * ssnow%ddq_dtg )

   ENDIF

   ! need to adjust fe after soilsnow
   canopy%fev  = canopy%fevc + canopy%fevw
  
   ! Calculate total latent heat flux:
   canopy%fe = canopy%fev + canopy%fes

   ! Calculate net radiation absorbed by soil + veg
   canopy%rnet = canopy%fns + canopy%fnv

   ! Calculate radiative/skin temperature:
   rad%trad = ( ( 1.-rad%transd ) * canopy%tv**4 +                             &
              rad%transd * ssnow%tss**4 )**0.25

   ! rml 17/1/11 move all plant resp and soil resp calculations here            
   ! from canopy. in UM only call on implicit step.
   ! put old and new soil resp calculations into soilcarb subroutine
   ! make new plantcarb subroutine
   IF (.not.cable_runtime%um_explicit .AND. icycle == 0) THEN

      !calculate canopy%frp
      CALL plantcarb(veg,bgc,met,canopy)
     
      !calculate canopy%frs
      CALL soilcarb(soil, ssnow, veg, bgc, met, canopy)

      CALL carbon_pl(dels, soil, ssnow, veg, canopy, bgc)

      canopy%fnpp = -1.0* canopy%fpn - canopy%frp
      canopy%fnee = canopy%fpn + canopy%frs + canopy%frp

   ENDIF

!    do k=1,mp,20
!    print 51,k,(soil%albsoil(j,1),j=k,k+19)
!51  format('albsoil1',1x,i5,20(1x,e11.4))
!    enddo
!    do k=1,mp,20
!    print 52,k,(soil%albsoil(j,2),j=k,k+19)
!52  format('albsoil2',1x,i5,20(1x,e11.4))
!    enddo
!    do k=1,mp,20
!    print 53,k,(ssnow%albsoilsn(j,1),j=k,k+19)
!53  format('albsoilsn1',1x,i5,20(1x,e11.4))
!    enddo
!    do k=1,mp,20
!    print 54,k,(ssnow%albsoilsn(j,2),j=k,k+19)
!54  format('albsoilsn2',1x,i5,20(1x,e11.4))
!    enddo
!    do k=1,mp,20
!    print 55,k,(ssnow%snowd(j),j=k,k+19)
!55  format('snowd',1x,i5,20(1x,e11.4))
!    enddo
!    do k=1,mp,20
!    print 56,k,(ssnow%wb(j,1),j=k,k+19)
!56  format('wb1',1x,i5,20(1x,e11.4))
!    enddo
!    do k=1,mp,20
!    print 57,k,(ssnow%wb(j,2),j=k,k+19)
!57  format('wb2',1x,i5,20(1x,e11.4))
!    enddo
!    do k=1,mp,20
!    print 58,k,(ssnow%wb(j,3),j=k,k+19)
!58  format('wb3',1x,i5,20(1x,e11.4))
!    enddo
!    do k=1,mp,20
!    print 59,k,(ssnow%wb(j,4),j=k,k+19)
!59  format('wb4',1x,i5,20(1x,e11.4))
!    enddo
!    do k=1,mp,20
!    print 60,k,(ssnow%wb(j,5),j=k,k+19)
!60  format('wb5',1x,i5,20(1x,e11.4))
!    enddo
!    do k=1,mp,20
!    print 61,k,(ssnow%wb(j,6),j=k,k+19)
!61  format('wb6',1x,i5,20(1x,e11.4))
!    enddo
!   do k=1,mp,20
!    print 62,k,(ssnow%wbice(j,1),j=k,k+19)
!62  format('wbice1',1x,i5,20(1x,e11.4))
!    enddo
!    do k=1,mp,20
!    print 63,k,(ssnow%wbice(j,2),j=k,k+19)
!63  format('wbice2',1x,i5,20(1x,e11.4))
!    enddo
!    do k=1,mp,20
!    print 64,k,(ssnow%wbice(j,3),j=k,k+19)
!64  format('wbice3',1x,i5,20(1x,e11.4))
!    enddo
!    do k=1,mp,20
!    print 65,k,(ssnow%wbice(j,4),j=k,k+19)
!65  format('wbice4',1x,i5,20(1x,e11.4))
!    enddo
!    do k=1,mp,20
!    print 66,k,(ssnow%wbice(j,5),j=k,k+19)
!66  format('wbice5',1x,i5,20(1x,e11.4))
!    enddo
!    do k=1,mp,20
!    print 67,k,(ssnow%wbice(j,6),j=k,k+19)
!67  format('wbice6',1x,i5,20(1x,e11.4))
!    enddo
!   do k=1,mp,20
!    print 68,k,(ssnow%tgg(j,1),j=k,k+19)
!68  format('tgg1',1x,i5,20(1x,e11.4))
!    enddo
!    do k=1,mp,20
!    print 69,k,(ssnow%tgg(j,2),j=k,k+19)
!69  format('tgg2',1x,i5,20(1x,e11.4))
!    enddo
!    do k=1,mp,20
!    print 70,k,(ssnow%tgg(j,3),j=k,k+19)
!70  format('tgg3',1x,i5,20(1x,e11.4))
!    enddo
!    do k=1,mp,20
!    print 71,k,(ssnow%tgg(j,4),j=k,k+19)
!71  format('tgg4',1x,i5,20(1x,e11.4))
!    enddo
!    do k=1,mp,20
!    print 72,k,(ssnow%tgg(j,5),j=k,k+19)
!72  format('tgg5',1x,i5,20(1x,e11.4))
!    enddo
!    do k=1,mp,20
!    print 73,k,(ssnow%tgg(j,6),j=k,k+19)
!73  format('tgg6',1x,i5,20(1x,e11.4))
!    enddo
!   do k=1,mp,20
!    print 74,k,(rough%z0m(j),j=k,k+19)
!74  format('z0m',1x,i5,20(1x,e11.4))
!    enddo
!    do k=1,mp,20
!    print 75,k,(rough%z0soil(j),j=k,k+19)
!75  format('z0soil',1x,i5,20(1x,e11.4))
!    enddo
!    do k=1,mp,20
!    print 76,k,(veg%hc(j),j=k,k+19)
!76  format('hc',1x,i5,20(1x,e11.4))
!    enddo
!    do k=1,mp,20
!    print 77,k,(veg%vlai(j),j=k,k+19)
!77  format('vlai',1x,i5,20(1x,e11.4))
!    enddo
!    do k=1,mp,20
!    print 78,k,(canopy%vlaiw(j),j=k,k+19)
!78  format('vlaiw',1x,i5,20(1x,e11.4))
!    enddo
!    do k=1,mp,20
!    print 79,k,(soil%ssat(j),j=k,k+19)
!79  format('ssat',1x,i5,20(1x,e11.4))
!    enddo

  
END SUBROUTINE cbm

END MODULE cable_cbm_module


