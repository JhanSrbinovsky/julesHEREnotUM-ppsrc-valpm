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
! Purpose: Updates CABLE variables (as altered by first pass through boundary 
!          layer and convection scheme), calls cbm, passes CABLE variables back 
!          to UM. 'Implicit' is the second call to cbm in each UM timestep.
!
! Called from: UM/JULES sf_impl
!
! Contact: Jhan.Srbinovsky@csiro.au
!
! History: Developed for CABLE v1.8
!
! ==============================================================================

subroutine cable_implicit_driver( LS_RAIN, CON_RAIN, LS_SNOW, CONV_SNOW,       &
                                  DTL_1,DQW_1, TSOIL, TSOIL_TILE, SMCL,        &
                                  SMCL_TILE, timestep, SMVCST,STHF, STHF_TILE, &
                                  STHU,STHU_TILE, snow_tile, SNOW_RHO1L,       &
                                  ISNOW_FLG3L, SNOW_DEPTH3L, SNOW_MASS3L,      &
                                  SNOW_RHO3L, SNOW_TMP3L, SNOW_COND,           &
                                  FTL_1, FTL_TILE, FQW_1, FQW_TILE,    &
                                  TSTAR_TILE, &
                                  SURF_HT_FLUX_LAND, ECAN_TILE, ESOIL_TILE,    &
                                  EI_TILE, RADNET_TILE, TOT_ALB, SNAGE_TILE,   &
                                  CANOPY_TILE, GS, T1P5M_TILE, Q1P5M_TILE,     &
                                  ! r935
                                  !CANOPY_TILE, GS,GS_TILE,T1P5M_TILE, Q1P5M_TILE,&
                                  CANOPY_GB, FLAND, MELT_TILE, DIM_CS1,        &
                                  DIM_CS2, NPP, NPP_FT, GPP, GPP_FT, RESP_S,   &
                                  RESP_S_TOT, RESP_S_TILE, RESP_P, RESP_P_FT,  &
                     ! r825 added casa vars after G_LEAF, but we need
                    ! need. vars here satisfy _hydrol CALL that is now from _impl
                    !idoy added r1164+
                                  G_LEAF, & 
                                  LYING_SNOW, SURF_ROFF, SUB_SURF_ROFF,  &
                                  TOT_TFALL )
                                  !G_LEAF, TRANSP_TILE, CPOOL_TILE, NPOOL_TILE, &
                                  !PPOOL_TILE, GLAI, PHENPHASE, NPP_FT_ACC,     &
                                  !RESP_W_FT_ACC, idoy )

   USE cable_def_types_mod, ONLY : mp
   USE cable_data_module,   ONLY : PHYS
   USE cable_data_mod,      ONLY : cable
   USE cable_um_tech_mod,   ONLY : um1, conv_rain_prevstep, conv_snow_prevstep,&
                                  air, bgc, canopy, met, bal, rad, rough, soil,&
                                  ssnow, sum_flux, veg, basic_diag
   USE cable_common_module, ONLY : cable_runtime, cable_user, l_casacnp,       &
                                   l_vcmaxFeedbk, knode_gl, ktau_gl, kend_gl,  &
                                   knode_gl
   USE cable_um_init_subrs_mod, ONLY : um2cable_rr
   USE cable_cbm_module,    ONLY : cbm
   USE casavariable
   USE phenvariable
   USE casa_types_mod
   !USE casa_cable
   USE casa_um_inout_mod

   IMPLICIT NONE
   
   character(len=*), parameter :: subr_name = "cable_implicit_driver"
        
   REAL, DIMENSION(um1%ROW_LENGTH,um1%ROWS) ::                                 &
      LS_RAIN,  & ! IN Large scale rain
      LS_SNOW,  & ! IN Large scale snow
      CON_RAIN, & ! IN Convective rain
      CONV_SNOW,& ! IN Convective snow
      DTL_1,    & ! IN Level 1 increment to T field 
      DQW_1       ! IN Level 1 increment to q field 

   REAL :: timestep

   INTEGER ::                                                                  &
      DIM_CS1, DIM_CS2 

   REAL, DIMENSION(um1%land_pts) ::                                            &
      GS,      &  ! OUT "Stomatal" conductance to
      SMVCST,  &  ! IN Volumetric saturation point
      FLAND       ! IN Land fraction on land tiles
   
   REAL, DIMENSION(um1%ROW_LENGTH,um1%ROWS) ::                                 &
      !--- Net downward heat flux at surface over land.
      !--- fraction of gridbox (W/m2).
      SURF_HT_FLUX_LAND,                                                       &   
      !--- Moisture flux between layers. (kg/m^2/sec).
      !--- FQW(,1) is total water flux from surface, 'E'.
      FQW_1,                                                                   &  
      !--- FTL(,K) =net turbulent sensible heat flux into layer K
      !--- from below; so FTL(,1) = surface sensible heat, H.(W/m2)
      FTL_1         

   REAL, DIMENSION(um1%land_pts,um1%ntiles) ::                              &
      !___Surface FTL, FQL for land tiles
      FTL_TILE, FQW_TILE,                                                   &
      
      !___(tiled) latent heat flux, melting, stomatatal conductance
     LE_TILE, MELT_TILE, GS_TILE,                                           &
     
     !___ INOUT Surface net radiation on tiles (W/m2)
     RADNET_TILE, &
     TOT_ALB,     & ! total albedo
     EI_TILE,     & ! OUT EI for land tiles.
     ECAN_TILE,   & ! OUT ECAN for snow-free land tiles
     ESOIL_TILE     ! evapotranspiration from soil moisture store (kg/m2/s) 

   !___ soil prognostics: moisture, frozen, unfrozen content, soil temp.
   !___ runoff ??
   REAL, dimension(um1%land_pts,um1%sm_levels) ::                           &
      SMCL,       & ! 
      STHF,       & !
      STHU,       & !
      TSOIL         !

   !___(tiled) soil prognostics: as above 
   REAL, dimension(um1%land_pts,um1%ntiles,um1%sm_levels) ::                &
      SMCL_TILE, & !
      STHU_TILE, & !
      TSOIL_TILE,& !
      STHF_TILE    !

   !___flag for 3 layer snow pack
   INTEGER :: ISNOW_FLG3L(um1%LAND_PTS,um1%NTILES)
   
   !___(tiled, 3 layer) Snow depth (m), mass, density, temp., conductivity
   REAL, dimension(um1%land_pts,um1%ntiles,3) :: &
      SNOW_DEPTH3L,  & ! 
      SNOW_MASS3L,   & !
      SNOW_RHO3L,    & !
      SNOW_TMP3L,    & !
      SNOW_COND        !

   REAL, DIMENSION(um1%land_pts,um1%ntiles) ::                              &
!     FRS_TILE,   & ! Local
!     NEE_TILE,   & ! Local
!     NPP_TILE,   & ! Local
!     GPP_TILE,   & ! Local
!     GLEAF_TILE, & ! Local, kdcorbin, 10/10
!     FRP_TILE,   &
      NPP_FT,     &
!     NPP_FT_old, &
      GPP_FT      
!     GPP_FT_old       

   REAL, DIMENSION(um1%land_pts) ::                                         &
      SNOW_GRD,    & !
      CANOPY_GB,   & !
      RESP_P,      & !
      NPP,         & !
      GPP            !
      
   REAL, DIMENSION( um1%land_pts,um1%ntiles ) ::                               &
      SNOW_TILE,     &
      SNOW_RHO1L,    &  ! Mean snow density
      SNAGE_TILE,    &
      CANOPY_TILE,   &
      T1P5M_TILE,    &
      Q1P5M_TILE,    &
      TSTAR_TILE,    &
      RESP_S_TILE,   & 
      RESP_P_FT,     &
!     RESP_P_FT_old, &
      G_LEAF,        &
      TRANSP_TILE

   REAL ::                                                                     &
      RESP_S(um1%LAND_PTS,DIM_CS1),     &
!     RESP_S_old(um1%LAND_PTS,DIM_CS1), &
      RESP_S_TOT(DIM_CS2)    

   REAL, DIMENSION(um1%LAND_PTS,um1%NTILES,10) ::                              &
      CPOOL_TILE, &
      NPOOL_TILE     
   REAL, DIMENSION(um1%LAND_PTS,um1%NTILES,12) ::                              &
      PPOOL_TILE
   REAL, DIMENSION(um1%LAND_PTS,um1%NTILES) ::                                 &
      GLAI, &
   !INTEGER, DIMENSION(um1%LAND_PTS,um1%NTILES) ::                              &
      PHENPHASE

   ! Lestevens 23apr13
   REAL, DIMENSION(um1%LAND_PTS,um1%NTILES) ::                                 &
      NPP_FT_ACC, &
      RESP_W_FT_ACC

   INTEGER ::     &
      ktauday,    &  ! day counter for CASA-CNP
      idoy           ! day of year (1:365) counter for CASA-CNP
   INTEGER, SAVE :: &
      kstart = 1

   REAL, DIMENSION(mp) ::                                                      & 
      vfrac,&
      dtlc, & 
      dqwc
   
   REAL, DIMENSION(um1%LAND_PTS) ::                               &
      LYING_SNOW,    & ! OUT Gridbox snowmass (kg/m2)        
      SUB_SURF_ROFF, & !
      SURF_ROFF,     & !
      TOT_TFALL        !

   REAL, POINTER :: TFRZ

   !___ 1st call in RUN (!=ktau_gl -see below) 
   LOGICAL, SAVE :: first_cable_call = .TRUE.

   type ProgBank
      
      real, dimension(:,:,:), allocatable ::                &
         TSOIL, SMCL, STHF 
         
      real, dimension(:,:,:), allocatable ::                &
         SNOW_DEPTH,                                                           &
         SNOW_MASS, SNOW_TMP, SNOW_RHO
      
      real, dimension(:,:), allocatable ::                &
         SNOW_RHO1L, SNOW_AGE 
      
      integer, dimension(:,:), allocatable ::                &
         SNOW_FLG3L

   End type ProgBank

   type (ProgBank), dimension(cable% um% numcycles) :: PB
   !type (ProgBank), dimension(2) :: PB
   
   integer :: cpb, ipb, i, j
   integer :: flpt,ftile,flev
   integer :: lpt,tile,lev
      !decs to write text files mapping mp points to lat/lon, i/j 
      character(len=*), parameter :: hcomp="mype"
      character(len=*), parameter :: hcompa ="i"
      character(len=*), parameter :: hcompb ="j"
      character(len=*), parameter :: hcompc ="l"
      character(len=*), parameter :: hcompd ="n"
      character(len=*), parameter :: hcomp1 ="mp "
      character(len=*), parameter :: hcomp2 ="lat "
      character(len=*), parameter :: hcomp3 ="lon"
      character(len=*), parameter :: hcomp4 ="frac"
      character(len=*), parameter :: hcomp5 ="iveg"
      character(len=*), parameter :: footer1 =""
      character(len=*), parameter :: footer2 = &
                  "---------------------------------------------------------------------------------"
      character(len=*), parameter :: hfmt1 = &
                  '(A8, 4X, A8, 4X, A8, 4X, A8, 4X, A8, 4X, A8, 4X, A8,   4X, A8,   4X, A8,    4X, A8)'
      character(len=*), parameter :: dfmt1 = &
                  '(I8, 4X, I8, 4X, I8, 4X, I8, 4X, I8, 4X, I8, 4X, F8.3, 4X, F8.3, 4X, ES8.2, 4X, I8)'
      character(len=30) :: chnode
      character(len=12) :: filename
      character(len=9), parameter :: basename="im_tsoil"
 

  
   flpt  = um1%land_pts
   ftile = um1%ntiles
   flev  = um1%sm_levels
   flpt  = min(1, flpt)
   ftile = min(1, ftile)
   flev  = min(1, flev)
   cpb = cable% um% numcycles
   if (.NOT. allocated(PB(1) %tsoil) ) then
      do ipb = 1, cpb 
         allocate( PB(ipb)%TSOIL(um1%land_pts,um1%ntiles,um1%sm_levels) ) 
         allocate( PB(ipb)%SMCL(um1%land_pts,um1%ntiles,um1%sm_levels) ) 
         allocate( PB(ipb)%STHF(um1%land_pts,um1%ntiles,um1%sm_levels) ) 
         allocate( PB(ipb)%snow_depth(um1%land_pts,um1%ntiles,3) )
         allocate( PB(ipb)%snow_mass(um1%land_pts,um1%ntiles,3) )
         allocate( PB(ipb)%snow_tmp(um1%land_pts,um1%ntiles,3) )
         allocate( PB(ipb)%snow_rho(um1%land_pts,um1%ntiles,3) )
         allocate( PB(ipb)%snow_rho1l(um1%land_pts,um1%ntiles) )
         allocate( PB(ipb)%snow_age(um1%land_pts,um1%ntiles) )
         allocate( PB(ipb)%snow_flg3l(um1%land_pts,um1%ntiles) )
      enddo
   endif   

   ipb = cable% um% cycleno
   print *,"cable% um% cycleno ",  cable% um% cycleno

   PB(ipb)%tsoil     = 0.  
   PB(ipb)%smcl      = 0. 
   PB(ipb)%sthf      = 0. 
   PB(ipb)%snow_depth= 0. 
   PB(ipb)%snow_mass = 0. 
   PB(ipb)%snow_tmp  = 0. 
   PB(ipb)%snow_rho  = 0. 
   PB(ipb)%snow_rho1l= 0. 
   PB(ipb)%snow_age  = 0. 
   PB(ipb)%snow_flg3l= 0 

!print *, "jh:shape(tsoil_soil ", shape(tsoil_tile)
!print *, "jh:shape(PB(1)%tsoil ", shape(PB(1)%tsoil)
!print *, "jh:shape(PB(2)%tsoil ", shape(PB(2)%tsoil)
!print *, "jh:(tsoil_soil ", (tsoil_tile)
            ! ----------------------------------------------------------------------------------
            write(chnode,10) knode_gl
   10       format(I3.3)   
            filename=trim(trim(basename)//trim(chnode))
            
            open(unit=12517,file=filename,status="unknown", &
                  action="write", form="formatted",position='append' )
               
               write (12517, *) "shape(tsoil_soil) "
               write (12517, *) shape(tsoil_tile)
               write (12517, *) " "
               
               write (12517, *) "shape(PB(1)%tsoil) "
               write (12517, *) shape(PB(1)%tsoil)
               write (12517, *) " "

               write (12517, *) "shape(PB(2)%tsoil) "
               write (12517, *) shape(PB(2)%tsoil)
               write (12517, *) " "

               do i=1,flpt
               do j=1,ftile
                  write (12517, *) "tsoil_soil(",i,",",j,"1) ", tsoil_tile(i,j,1)
                  write (12517, *) "tsoil_soil(",i,",",j,"2) ", tsoil_tile(i,j,2)
                  write (12517, *) "tsoil_soil(",i,",",j,"3) ", tsoil_tile(i,j,3)
                  write (12517, *) "tsoil_soil(",i,",",j,"4) ", tsoil_tile(i,j,4)
                  write (12517, *) "tsoil_soil(",i,",",j,"5) ", tsoil_tile(i,j,5)
                  write (12517, *) "tsoil_soil(",i,",",j,"6) ", tsoil_tile(i,j,6)
                  write (12517, *) "ipb ", ipb, PB(ipb)%tsoil(i,j,:)
               enddo   
               enddo   
               write (12517, *) " "
            
               PB(ipb)%tsoil     = tsoil_tile 
            close(12517)



!if(um1%land_pts >0 ) then
   PB(ipb)%tsoil     = tsoil_tile 
   PB(ipb)%smcl      = smcl_tile
   PB(ipb)%sthf      = sthf_tile 
   PB(ipb)%snow_depth= snow_depth3l
   PB(ipb)%snow_mass = snow_mass3l
   PB(ipb)%snow_tmp  = snow_tmp3l
   PB(ipb)%snow_rho  = snow_rho3l
   PB(ipb)%snow_rho1l= snow_rho1l
   PB(ipb)%snow_age  = snage_tile
   PB(ipb)%snow_flg3l= isnow_flg3l
!endif

if (ipb == cable% um% numcycles) then
   tsoil_tile     = PB(1)%tsoil  
   smcl_tile      = PB(1)%smcl
   sthf_tile      = PB(1)%sthf
   snow_depth3l   = PB(1)%snow_depth
   snow_mass3l    = PB(1)%snow_mass
   snow_tmp3l     = PB(1)%snow_tmp
   snow_rho3l     = PB(1)%snow_rho
   snow_rho1l     = PB(1)%snow_rho1l
   snage_tile     = PB(1)%snow_age
   isnow_flg3l    = PB(1)%snow_flg3l

   open(unit=12511,file='c_impl_progs_tropics',status="unknown", &
                  action="write", form="formatted",position='append' )
      WRITE(12511,*) , "" 

      write(12511,*) "tropics"!, !(where - insert inices in WRITE) 
      lpt = 1
      tile = 17
      lev = 1
       
      write(12511,*) "tsoil ", PB(1)%tsoil(lpt,tile,lev), PB(2)%tsoil(lpt,tile,lev)
      write(12511,*) "smcl  ", PB(1)%smcl(lpt,tile,lev), PB(2)%smcl(lpt,tile,lev)        
      write(12511,*) "sthf  ", PB(1)%sthf(lpt,tile,lev), PB(2)%sthf(lpt,tile,lev)
         
      write(12511,*) "depth ", PB(1)%snow_depth(lpt,tile,lev), PB(2)%snow_depth(lpt,tile,lev)     
      write(12511,*) "mass  ", PB(1)%snow_mass(lpt,tile,lev) , PB(2)%snow_mass(lpt,tile,lev)    
      write(12511,*) "temp  ", PB(1)%snow_tmp(lpt,tile,lev)  , PB(2)%snow_tmp(lpt,tile,lev)    
      write(12511,*) "rho3l ", PB(1)%snow_rho(lpt,tile,lev)  , PB(2)%snow_rho(lpt,tile,lev)    
         
      write(12511,*) "rho1l ", PB(1)%snow_rho1l(lpt,tile), PB(2)%snow_rho1l(lpt,tile)  
      write(12511,*) "snage ", PB(1)%snow_age(lpt,tile), PB(2)%snow_age(lpt,tile)    
         
      write(12511,*) "flg3L ", PB(1)%snow_flg3l(lpt,tile), PB(2)%snow_flg3l(lpt,tile)  

   close(12511)

   open(unit=12511,file='c_impl_progs_sahara',status="unknown", &
                  action="write", form="formatted",position='append' )
      WRITE(12511,*) , "" 

      write(12511,*) "sahara"!, !(where - insert inices in WRITE) 
      lpt = 1
      tile = 17
      lev = 1
       
      write(12511,*) "tsoil ", PB(1)%tsoil(lpt,tile,lev), PB(2)%tsoil(lpt,tile,lev)
      write(12511,*) "smcl  ", PB(1)%smcl(lpt,tile,lev), PB(2)%smcl(lpt,tile,lev)        
      write(12511,*) "sthf  ", PB(1)%sthf(lpt,tile,lev), PB(2)%sthf(lpt,tile,lev)
         
      write(12511,*) "depth ", PB(1)%snow_depth(lpt,tile,lev), PB(2)%snow_depth(lpt,tile,lev)     
      write(12511,*) "mass  ", PB(1)%snow_mass(lpt,tile,lev) , PB(2)%snow_mass(lpt,tile,lev)    
      write(12511,*) "temp  ", PB(1)%snow_tmp(lpt,tile,lev)  , PB(2)%snow_tmp(lpt,tile,lev)    
      write(12511,*) "rho3l ", PB(1)%snow_rho(lpt,tile,lev)  , PB(2)%snow_rho(lpt,tile,lev)    
         
      write(12511,*) "rho1l ", PB(1)%snow_rho1l(lpt,tile), PB(2)%snow_rho1l(lpt,tile)  
      write(12511,*) "snage ", PB(1)%snow_age(lpt,tile), PB(2)%snow_age(lpt,tile)    
         
      write(12511,*) "flg3L ", PB(1)%snow_flg3l(lpt,tile), PB(2)%snow_flg3l(lpt,tile)  

   close(12511)
endif

!if (cable% um% cycleno == 1) then
!if (cable% um% cycleno == cable% um% numcycles) then

      !DO J=1,um1%sm_levels
      !   ssnow%tgg(:,J) = PACK(TSOIL_TILE(:,:,J),um1%l_tile_pts)
      !ENDDO 
      !IF( first_cable_call) ssnow%tggav = 0.
      !DO k = 1, um1%sm_levels
      !   ssnow%tggav = ssnow%tggav  + soil%zse(k)*ssnow%tgg(:,k) / 4.596
      !END DO
      !ijctr=0 
      !do i=1, um1%land_pts
      !   do j=1,um1%sm_levels
      !      if( tsoil(i,j) < 230. )then 
      !         ijctr=ijctr+1
      !         !print *, "_impli:tsoil ", tsoil(i,j)
      !      endif   
      !      do k=1,um1%ntiles
      !         if(k==17) then 
      !            print *, "_impli:tsoil17 ", tsoil_tile(i,k,j)
      !            tsoil_tile(i,k,j) = tsoil(i,j)
      !         endif   
      !      enddo
      !   enddo
      !enddo
      !print *, "# times tsoil < 230 ", ijctr
      !
      !DO J=1,um1%sm_levels
      !   ssnow%tgg(:,J) = PACK(TSOIL_TILE(:,:,J),um1%l_tile_pts)
      !ENDDO



      IF(cable_user%run_diag_level == "BASIC")                                    &     
         CALL basic_diag(subr_name, "Called.") 

      TFRZ => PHYS%TFRZ
   
      ! FLAGS def. specific call to CABLE from UM
      cable_runtime%um_explicit = .FALSE.
      !cable_runtime%um_implicit = .TRUE.
   
      dtlc = 0. ; dqwc = 0.

      !--- All these subrs do is pack a CABLE var with a UM var.
      !-------------------------------------------------------------------
      !--- UM met forcing vars needed by CABLE which have UM dimensions
      !---(rowlength,rows)[_rr], which is no good to CABLE. These have to be 
      !--- re-packed in a single vector of active tiles. Hence we use 
      !--- conditional "mask" l_tile_pts(land_pts,ntiles) which is .true.
      !--- if the land point is/has an active tile
      !--- generic format:
      !--- um2cable_rr( UM var, default value for snow tile, CABLE var, mask )
      !--- where mask tells um2cable_rr whether or not to use default value 
      !--- for snow tile 
      !-------------------------------------------------------------------
      CALL um2cable_rr( (LS_RAIN+CON_RAIN)*um1%TIMESTEP, met%precip)
      CALL um2cable_rr( (LS_SNOW+CONV_SNOW)*um1%TIMESTEP, met%precip_sn)
      CALL um2cable_rr( dtl_1, dtlc)
      CALL um2cable_rr( dqw_1, dqwc)
      
      !--- conv_rain(snow)_prevstep are added to precip. in explicit call
      CALL um2cable_rr( (CON_RAIN)*um1%TIMESTEP, conv_rain_prevstep)
      CALL um2cable_rr( (CONV_snow)*um1%TIMESTEP, conv_snow_prevstep)
      
      met%precip   =  met%precip + met%precip_sn
      met%tk = met%tk + dtlc
      met%qv = met%qv + dqwc
      met%tvair = met%tk
      met%tvrad = met%tk
 
      canopy%cansto = canopy%oldcansto

      CALL cbm(real(TIMESTEP), air, bgc, canopy, met, bal,  &
           rad, rough, soil, ssnow, sum_flux, veg)

      ! Lestevens - temporary ?
      ktauday = int(24.0*3600.0/TIMESTEP)
      IF(idoy==0) idoy =365
  
! Lestevens Sept2012 - Call CASA-CNP
!      if (l_casacnp) then
!      CALL bgcdriver(ktau_gl,kstart,kend_gl,TIMESTEP,met,ssnow,canopy,veg,soil, &
!                     casabiome,casapool,casaflux,casamet,casabal,phen,          &
!                     .FALSE., .FALSE., ktauday, idoy, .FALSE., .FALSE. )
!                     spinConv, spinup, ktauday, idoy, cable_user%casa_dump_read,&
!                     cable_user%casa_dump_write )
!      endif

!      CALL sumcflux(ktau_gl,kstart,kend_gl,TIMESTEP,bgc,canopy,soil,ssnow,      &
!                    sum_flux,veg,met,casaflux,l_vcmaxFeedbk)

      !cable_runtime%um_implicit = .FALSE.
!end if ! cycleno

!if (cable% um% cycleno == cable% um% numcycles) then
      CALL implicit_unpack( TSOIL, TSOIL_TILE, SMCL, SMCL_TILE,                &
                            SMVCST, STHF, STHF_TILE, STHU, STHU_TILE,          &
                            snow_tile, SNOW_RHO1L ,ISNOW_FLG3L, SNOW_DEPTH3L,  &
                            SNOW_MASS3L, SNOW_RHO3L, SNOW_TMP3L, SNOW_COND,    &
                            FTL_1, FTL_TILE, FQW_1,  FQW_TILE, TSTAR_TILE,     &
                            SURF_HT_FLUX_LAND, ECAN_TILE,        &
                            ESOIL_TILE, EI_TILE, RADNET_TILE, TOT_ALB,         &
                            SNAGE_TILE, CANOPY_TILE, GS, T1P5M_TILE,           &
                            !r935 SNAGE_TILE, CANOPY_TILE, GS,GS_TILE, T1P5M_TILE,        &
                            Q1P5M_TILE, CANOPY_GB, FLAND, MELT_TILE, DIM_CS1,  &
                            DIM_CS2, NPP, NPP_FT, GPP, GPP_FT, RESP_S,         &
                            RESP_S_TOT, RESP_S_TILE, RESP_P, RESP_P_FT, G_LEAF,& 
                            TRANSP_TILE, NPP_FT_ACC, RESP_W_FT_ACC )

! Lestevens Sept2012 - Call CASA-CNP
!      if (l_casacnp) then
!      !if (l_casacnp .and. ktau_gl==kend_gl ) then
!        if (knode_gl==0 .and. ktau_gl==kend_gl) then
!        !if (knode_gl==0) then
!         print *, '  '; print *, 'CASA_log:'
!         print *, '  Calling CasaCNP - Poolout '
!         print *, '  l_casacnp = ',l_casacnp
!         print *, '  ktau_gl, kend_gl = ',ktau_gl,kend_gl
!         print *, 'End CASA_log:'; print *, '  '
!        endif
!       CALL casa_poolout_unpk(casapool,casaflux,casamet,casabal,phen,  &
!                              CPOOL_TILE,NPOOL_TILE,PPOOL_TILE, &
!                              GLAI,PHENPHASE)
!      endif
       
! DEPENDS ON: cable_hyd_driver
   call cable_hyd_driver( SNOW_TILE, LYING_SNOW, SURF_ROFF, SUB_SURF_ROFF,  &
                             TOT_TFALL )
   vfrac = pack(um1%tile_frac,um1%l_tile_pts)

   cable_runtime%um_implicit = .FALSE.

!end if ! cycleno
  
   IF(cable_user%run_diag_level == "BASIC")                                    &     
      CALL basic_diag(subr_name, "Done.") 

END SUBROUTINE cable_implicit_driver


!========================================================================= 
!========================================================================= 
!========================================================================= 
        

! rml 30/10/15 rml: rewrite carbon fluxes - previous version seemed
! unnecessarily complicated by allowing for fluxes from sea-ice when don't
! need to account for any carbon flux not from land fraction.
SUBROUTINE implicit_unpack( TSOIL, TSOIL_TILE, SMCL, SMCL_TILE,                &
                            SMVCST, STHF, STHF_TILE, STHU, STHU_TILE,          &
                            snow_tile, SNOW_RHO1L ,ISNOW_FLG3L, SNOW_DEPTH3L,  &
                            SNOW_MASS3L, SNOW_RHO3L, SNOW_TMP3L, SNOW_COND,    &
                            FTL_1, FTL_TILE, FQW_1,  FQW_TILE, TSTAR_TILE,     &
                            SURF_HT_FLUX_LAND, ECAN_TILE,        &
                            ESOIL_TILE, EI_TILE, RADNET_TILE, TOT_ALB,         &
                            SNAGE_TILE, CANOPY_TILE, GS, T1P5M_TILE,  &
                            !args increased following merge
                            !SNAGE_TILE, CANOPY_TILE, GS, GS_TILE, T1P5M_TILE,  &
                            Q1P5M_TILE, CANOPY_GB, FLAND, MELT_TILE, DIM_CS1,  &
                            DIM_CS2, NPP, NPP_FT, GPP, GPP_FT, RESP_S,         &
                            RESP_S_TOT, RESP_S_TILE, RESP_P, RESP_P_FT, G_LEAF,&
                            TRANSP_TILE, NPP_FT_ACC, RESP_W_FT_ACC )
 
   USE cable_def_types_mod, ONLY : mp
   USE cable_data_module,   ONLY : PHYS
   USE cable_um_tech_mod,   ONLY : um1 ,canopy, rad, soil, ssnow, air,         &
                                   basic_diag
   USE cable_common_module, ONLY : cable_runtime, cable_user, fudge_out,       &
                                   L_fudge
   IMPLICIT NONE
 
   character(len=*), parameter :: subr_name = "cable_implicit_unpack"
        
   !jhan:these need to be cleaned out to what is actualllly passed
   INTEGER :: DIM_CS1 ,DIM_CS2 

   REAL, DIMENSION(um1%land_pts) ::                                            &
      GS,         &  ! OUT "Stomatal" conductance to
      SMVCST,     &  ! IN Volumetric saturation point
      FLAND          ! IN Land fraction on land tiles
   
   real, dimension(um1%ROW_LENGTH,um1%ROWS) ::                                 &
      !--- Net downward heat flux at surface over land.
      !--- fraction of gridbox (W/m2).
      SURF_HT_FLUX_LAND,           &
      !--- Moisture flux between layers. (kg/m^2/sec).
      !--- FQW(,1) is total water flux from surface, 'E'.
      FQW_1,       &  
      !--- FTL(,K) =net turbulent sensible heat flux into layer K
      !--- from below; so FTL(,1) = surface sensible heat, H.(W/m2)
      FTL_1         

   REAL, DIMENSION(um1%land_pts,um1%ntiles) ::                                 &
      !___Surface FTL, FQL for land tiles
      FTL_TILE, FQW_TILE,                 &  
      !___(tiled) latent heat flux, melting, stomatatal conductance
      LE_TILE, MELT_TILE, GS_TILE,     &  
      RADNET_TILE, & ! INOUT Surface net radiation on tiles (W/m2)
      TOT_ALB,     & ! total albedo
      EI_TILE,     & ! OUT EI for land tiles.
      ECAN_TILE,   & ! OUT ECAN for snow-free land tiles
      ESOIL_TILE     ! evapotranspiration from soil moisture store (kg/m2/s) 

   !___ soil prognostics: moisture, frozen, unfrozen content, soil temp.
   !___ runoff ??
   REAL, DIMENSION(um1%land_pts,um1%sm_levels) ::                              &
      SMCL,       & !
      STHF,       &
      STHU,       &
      TSOIL       

   !___(tiled) soil prognostics: as above 
   REAL, DIMENSION(um1%land_pts,um1%ntiles,um1%sm_levels) ::                   &
      SMCL_TILE,  & 
      STHU_TILE,  &
      TSOIL_TILE, &
      STHF_TILE  

   !___flag for 3 layer snow pack
   INTEGER :: ISNOW_FLG3L(um1%LAND_PTS,um1%NTILES)
   
   !___(tiled, 3 layer) Snow depth (m), mass, density, temp., conductivity
   REAL, DIMENSION(um1%land_pts,um1%ntiles,3) ::                               &
      SNOW_DEPTH3L,  &
      SNOW_MASS3L,   &
      SNOW_RHO3L,    &
      SNOW_TMP3L,    &
      SNOW_COND 

   REAL, dimension(um1%land_pts,um1%ntiles) ::                                 &
!     FRS_TILE,   & ! Local
!     NEE_TILE,   & ! Local
!     NPP_TILE,   & ! Local
!     GPP_TILE,   & ! Local
      SURF_HTF_T_CAB
!     GLEAF_TILE, & ! Local, kdcorbin, 10/10
!     FRP_TILE

   REAL, dimension(um1%land_pts,um1%ntiles) ::                           &
      NPP_FT,     &
!     NPP_FT_old, &
      GPP_FT       
!     GPP_FT_old

   REAL, dimension(um1%land_pts) ::                                            &
      RESP_P,     & 
      NPP,        & 
      GPP

   REAL, dimension(um1%land_pts) ::                                            &
      SNOW_GRD,   &  
      CANOPY_GB

   REAL, DIMENSION(um1%land_pts,um1%ntiles) ::                                 &
      NPP_FT_ACC,    & ! sresp for CASA-CNP
      RESP_W_FT_ACC    ! presp for CASA-CNP
   
   REAL, DIMENSION(um1%land_pts,um1%ntiles) ::                                 &
      SNOW_TILE,     & !
      SNOW_RHO1L,    & ! Mean snow density
      SNAGE_TILE,    & !
      CANOPY_TILE,   & !
      T1P5M_TILE,    &
      Q1P5M_TILE,    &
      TSTAR_TILE,    &
      RESP_S_TILE,   & 
      RESP_P_FT,     &
!     RESP_P_FT_old, &
      G_LEAF,        &
      TRANSP_TILE

   REAL ::                                                                     &
      RESP_S(um1%LAND_PTS,DIM_CS1),    & !
      RESP_S_old(um1%LAND_PTS,DIM_CS1),& !
      RESP_S_TOT(DIM_CS2)                !
  
   REAL, DIMENSION(mp) ::                                                                     &
      fe_dlh,    & !
      fes_dlh,   & !
      fev_dlh      !

   !--- Local vars
   INTEGER :: i,j,l,k,n

   REAL, DIMENSION(um1%land_pts,um1%ntiles) ::                                 &
         !--- Local buffer surface FTL, FQL @ prev dt
         FTL_TILE_old, FQW_TILE_old

   INTEGER:: i_miss = 0
   REAL :: miss = 0.0
   
   REAL, POINTER :: TFRZ
   
!  LOGICAL, SAVE :: first_call = .TRUE.
   
   IF(cable_user%run_diag_level == "BASIC")                                    &     
      CALL basic_diag(subr_name, "Called.") 

      TFRZ => PHYS%TFRZ
  
      !--- set UM vars to zero
      SMCL_TILE = 0.; STHF_TILE = 0.; STHU_TILE = 0.
      TSOIL_TILE = 0.

      DO j = 1,um1%SM_LEVELS
         TSOIL_TILE(:,:,j)= UNPACK(ssnow%tgg(:,j), um1%L_TILE_PTS, miss)
         SMCL_TILE(:,:,j)= UNPACK(REAL(ssnow%wb(:,j)), um1%L_TILE_PTS, miss)
         SMCL_TILE(:,:,j)=SMCL_TILE(:,:,j)*soil%zse(j)*um1%RHO_WATER
         STHF_TILE(:,:,j)= UNPACK(REAL(ssnow%wbice(:,j)), um1%L_TILE_PTS, miss)
         SMCL(:,j) = SUM(um1%TILE_FRAC * SMCL_TILE(:,:,j),2)
         TSOIL(:,j) = SUM(um1%TILE_FRAC * TSOIL_TILE(:,:,j),2)
         
         DO N=1,um1%NTILES
            DO K=1,um1%TILE_PTS(N)
               I = um1%TILE_INDEX(K,N)
               IF ( SMVCST(I) > 0. ) THEN ! Exclude permanent ice - mrd
                  STHF_TILE(I,N,J)= STHF_TILE(I,N,J)/SMVCST(I)
                  STHU_TILE(I,N,J)= MAX( 0., SMCL_TILE(I,N,J) -                &
                                    STHF_TILE(I,N,J) * SMVCST(I) * soil%zse(J) &
                                    * um1%RHO_WATER ) / ( soil%zse(J) *        &
                                    um1%RHO_WATER * SMVCST(I) )
               ENDIF
            ENDDO
         ENDDO

         STHF(:,J) = SUM(um1%TILE_FRAC * STHF_TILE(:,:,J),2)
         STHU(:,J) = SUM(um1%TILE_FRAC * STHU_TILE(:,:,J),2)
      ENDDO

      !--- unpack snow vars 
      SNOW_RHO1L  = UNPACK(ssnow%ssdnn, um1%L_TILE_PTS, miss)
      ISNOW_FLG3L = UNPACK(ssnow%isflag, um1%L_TILE_PTS, i_miss)
      MELT_TILE   = UNPACK(ssnow%smelt, um1%L_TILE_PTS, miss)
      SNOW_TILE= UNPACK(ssnow%snowd, um1%L_TILE_PTS, miss)
      SNOW_GRD=  SUM(um1%TILE_FRAC * SNOW_TILE,2)  ! gridbox snow mass & snow below canopy 

      !--- unpack layered snow vars 
      do k = 1,3
        SNOW_TMP3L(:,:,k) = UNPACK(ssnow%tggsn(:,k), um1%L_TILE_PTS, miss)
        SNOW_MASS3L(:,:,k)= UNPACK(ssnow%smass(:,k), um1%L_TILE_PTS, miss)
        SNOW_RHO3L(:,:,k) = UNPACK(ssnow%ssdn(:,k), um1%L_TILE_PTS, miss)
        SNOW_COND(:,:,k)  = UNPACK(ssnow%sconds(:,k),um1%L_TILE_PTS,miss)
        SNOW_DEPTH3L(:,:,k)  = UNPACK(ssnow%sdepth(:,k),um1%L_TILE_PTS,miss)
      enddo

      
      canopy%gswx_T = canopy%gswx_T/air%cmolar
      GS_TILE = UNPACK(canopy%gswx_T,um1%L_TILE_PTS,miss)
      GS =  SUM(um1%TILE_FRAC * GS_TILE,2)

      !---preserve fluxes from the previous time step for the coastal grids
      FTL_TILE_old = FTL_TILE
      FQW_TILE_old = FQW_TILE
      !___return fluxes
      FTL_TILE = UNPACK(canopy%fh,  um1%l_tile_pts, miss)
!      fe_dlh = canopy%fe/(air%rlam*ssnow%cls)
      fes_dlh = canopy%fes/(air%rlam*ssnow%cls)
      fev_dlh = canopy%fev/air%rlam
      fe_dlh =  fev_dlh + fes_dlh

      !---update fluxes 
      FQW_TILE = UNPACK(fe_dlh, um1%l_tile_pts, miss)


      !___return temp and roughness
      TSTAR_TILE = UNPACK(rad%trad, um1%l_tile_pts, miss)

      !___return miscelaneous 
      RADNET_TILE = unpack( canopy%rnet , um1%l_tile_pts, miss)

     SURF_HTF_T_CAB = UNPACK(canopy%ga,um1%L_TILE_PTS,miss)

     TOT_ALB=UNPACK(rad%albedo_T,um1%L_TILE_PTS, miss) 
     ESOIL_TILE = UNPACK(fes_dlh, um1%L_TILE_PTS, miss)
     ECAN_TILE = UNPACK(fev_dlh,  um1%L_TILE_PTS, miss)
     EI_TILE = 0.
     SNAGE_TILE = UNPACK(ssnow%snage, um1%L_TILE_PTS, miss) 
     TRANSP_TILE = UNPACK(canopy%fevc, um1%L_TILE_PTS, miss) 

     !unpack screen level (1.5m) variables
     !Convert back to K 
     t1p5m_tile     = UNPACK(canopy%tscrn+tfrz, um1%L_TILE_PTS, miss)
     q1p5m_tile     = UNPACK(canopy%qscrn, um1%L_TILE_PTS, miss)
     CANOPY_TILE    = UNPACK(canopy%cansto, um1%L_TILE_PTS, miss)
     CANOPY_GB      = SUM(um1%TILE_FRAC * CANOPY_TILE,2)

     !initialse full land grids and retain coastal grid fluxes
      DO N=1,um1%NTILES
         DO K=1,um1%TILE_PTS(N)
           L = um1%TILE_INDEX(K,N)
           J=(um1%LAND_INDEX(L)-1)/um1%ROW_LENGTH + 1
           I = um1%LAND_INDEX(L) - (J-1)*um1%ROW_LENGTH
           IF( FLAND(L) == 1.0) THEN 
             FTL_1(I,J) =  0.0
             FQW_1(I,J) =  0.0
           ELSE
             !retain sea/ice contribution and remove land contribution
             FTL_1(I,J) = FTL_1(I,J) - FLAND(L) * um1%TILE_FRAC(L,N) *         &
                          FTL_TILE_old(L,N)
             FQW_1(I,J) = FQW_1(I,J) - FLAND(L) * um1%TILE_FRAC(L,N) *         &
                          FQW_TILE_old(L,N)
           ENDIF
           SURF_HT_FLUX_LAND(I,J) = 0.
         ENDDO
     ENDDO

      DO N=1,um1%NTILES
         DO K=1,um1%TILE_PTS(N)
           L = um1%TILE_INDEX(K,N)
           J=(um1%LAND_INDEX(L)-1)/um1%ROW_LENGTH + 1
           I = um1%LAND_INDEX(L) - (J-1)*um1%ROW_LENGTH
           FTL_1(I,J) = FTL_1(I,J)+FLAND(L)*um1%TILE_FRAC(L,N)*FTL_TILE(L,N)
           FQW_1(I,J) = FQW_1(I,J)+FLAND(L)*um1%TILE_FRAC(L,N)*FQW_TILE(L,N)
           SURF_HT_FLUX_LAND(I,J) = SURF_HT_FLUX_LAND(I,J) +                   &
                                    FLAND(L)*um1%TILE_FRAC(L,N) *              &
                                    SURF_HTF_T_CAB(L,N)
         ENDDO
      ENDDO

!    rml 30/10/15 Initialise grid-cell carbon fields that are accumulated over tiles
     RESP_P = 0.
     NPP = 0.
     GPP = 0.
     RESP_S = 0.

     ! Lestevens - Passing CO2 from CABLE to bl_trmix_dd.F90
!    rml 30/10/15 unpack cable carbon variables straight into um variables
!    (previously to local)
     RESP_S_TILE    = UNPACK(canopy%frs, um1%L_TILE_PTS, miss)
!    NEE_TILE       = UNPACK(canopy%fnee, um1%L_TILE_PTS, miss)
     NPP_FT         = UNPACK(canopy%fnpp, um1%L_TILE_PTS, miss)
     G_LEAF         = UNPACK(canopy%frday,um1%L_TILE_PTS, miss)
     RESP_P_FT      = UNPACK(canopy%frp, um1%L_TILE_PTS, miss)

!    rml - probably want to retire this switch and always have GPP as real GPP
!    (off case)
     IF( cable_user%leaf_respiration == 'on' .OR.                             &
           cable_user%leaf_respiration == 'ON') THEN
        GPP_FT = UNPACK(canopy%fnpp+canopy%frp, um1%L_TILE_PTS, miss)
     ELSE
        GPP_FT = UNPACK(canopy%fnpp+canopy%frp+canopy%frday,  &
                            um1%L_TILE_PTS, miss)
     ENDIF

!    convert from CABLE units (gC/m2/s) to UM units (kgC/m2/s)
     RESP_S_TILE = RESP_S_TILE*1.e-3
     G_LEAF      = G_LEAF*1.e-3
     NPP_FT      = NPP_FT*1.e-3
     GPP_FT      = GPP_FT*1.e-3
     RESP_P_FT   = RESP_P_FT*1.e-3

     ! Lestevens 23apr13 - possible miss match ntiles<-->npft
!    If CASA-CNP used, plant and soil resp need to be passed into 
!    variables that are dumped to restart, because CASA-CNP only run daily
     NPP_FT_ACC = RESP_S_TILE
     RESP_W_FT_ACC = RESP_P_FT

      DO N=1,um1%NTILES 
         DO K=1,um1%TILE_PTS(N)
            L = um1%TILE_INDEX(K,N)
            NPP(L)=NPP(L)+FLAND(L)*um1%TILE_FRAC(L,N)*NPP_FT(L,N)
            GPP(L)=GPP(L)+FLAND(L)*um1%TILE_FRAC(L,N)*GPP_FT(L,N)
            RESP_P(L)=RESP_P(L)+FLAND(L)*um1%TILE_FRAC(L,N)*RESP_P_FT(L,N)

            !loop for soil resp. although DIM_CS1=1 (not 1 for triffid)
            DO I=1,DIM_CS1
               RESP_S(L,I) = RESP_S(L,I) + &
                             FLAND(L)*um1%TILE_FRAC(L,N)*RESP_S_TILE(L,N)
            ENDDO
            RESP_S_TOT(L)=sum(RESP_S(L,:))
         ENDDO
      ENDDO

   IF(cable_user%run_diag_level == "BASIC")                                    &     
      CALL basic_diag(subr_name, "Call fudge_out on Logic switch") 

      if(L_fudge) then
         call fudge_out( 1,1,1,  tsoil_tile,    'tsoil_tile',  .TRUE., 278. )
         call fudge_out( 1,1,1,  smcl_tile,     'smcl_tile',   .TRUE., 0. )
         call fudge_out( 1,1,    tsoil,         'tsoil',       .TRUE., 278. )
         call fudge_out( 1,1,1,  sthf_tile,     'sthf_tile',   .TRUE., 0. )
         call fudge_out( 1,1,1,  sthu_tile,     'sthu_tile',   .TRUE., 0. )
         call fudge_out( 1,1,    sthf,          'sthf',        .TRUE., 0. )
         call fudge_out( 1,1,    sthu,          'sthu',        .TRUE., 0. )
         call fudge_out( 1,1,    snow_rho1l,    'snow_rho1l',  .TRUE., 0. )
         call fudge_out( 1,1,    ISNOW_FLG3L,   'ISNOW_FLG3L', .TRUE., 0 )     
         call fudge_out( 1,1,    MELT_TILE,     'MELT_TILE',   .TRUE., 0. )   
         call fudge_out( 1,1,    snow_TILE,     'snow_TILE',   .TRUE., 0. )   
         call fudge_out( 1,      snow_grd,      'snow_grd',    .TRUE., 0. )   
         call fudge_out( 1,1,1,  snow_Tmp3L,    'snow_Tmp3L',  .TRUE., 0. )      
         call fudge_out( 1,1,1,  snow_mass3L,   'snow_mass3L', .TRUE., 0. )      
         call fudge_out( 1,1,1,  snow_rho3L,    'snow_rho3L',  .TRUE., 0. )      
         call fudge_out( 1,1,1,  snow_depth3L,  'snow_depth3L',.TRUE., 0. )         
         call fudge_out( 1,1,1,  snow_cond,     'snow_cond',   .TRUE., 0. )         
         call fudge_out( 1,1,    GS_TILE,       'GS_TILE',     .TRUE., 0. )      
         call fudge_out( 1,      GS,            'GS',          .TRUE., 0. )            
         call fudge_out( 1,1,    FTL_TILE,      'FTL_TILE',    .TRUE., 0. )                        
         call fudge_out( 1,1,    Fqw_TILE,      'Fqw_TILE',    .TRUE., 0. )                                 
         call fudge_out( 1,1,    tstar_TILE,    'tstar_TILE',  .TRUE., 280. )                              
         call fudge_out( 1,1,    radnet_TILE,   'radnet_TILE', .TRUE., 0. )                                    
         call fudge_out( 1,1,    TOT_ALB,       'TOT_ALB',     .TRUE., 0. )                                 
         call fudge_out( 1,1,    ESOIL_TILE,    'ESOIL_TILE',  .TRUE., 0. )                                 
         call fudge_out( 1,1,    ECAN_TILE,     'ECAN_TILE',   .TRUE., 0. )                                    
         call fudge_out( 1,1,    snage_TILE,    'snage_TILE',  .TRUE., 0. )                                 
         call fudge_out( 1,1,    t1p5m_TILE,    't1p5m_TILE',  .TRUE., 0. )                           
         call fudge_out( 1,1,    q1p5m_TILE,    'q1p5m_TILE',  .TRUE., 0. )                                       
         call fudge_out( 1,1,    canopy_TILE,   'canopy_TILE', .TRUE., 0. )               
         call fudge_out( 1,      canopy_GB,     'canopy_GB',   .TRUE., 0. )                           
!        call fudge_out( 1,1,    frs_TILE,      'frs_TILE',    .TRUE., 0. )                                 
!        call fudge_out( 1,1,    NEE_TILE,      'NEE_TILE',    .TRUE., 0. )                                 
!        call fudge_out( 1,1,    NPP_TILE,      'NPP_TILE',    .TRUE., 0. )                                          
!        call fudge_out( 1,1,    Gleaf_TILE,    'Gleaf_TILE',  .TRUE., 0. )                                          
!        call fudge_out( 1,1,    GPP_TILE,      'GPP_TILE',    .TRUE., 0. )                                    
!        call fudge_out( 1,1,    frp_TILE,      'frp_TILE',    .TRUE., 0. )                                                
         call fudge_out( 1,1,    fTL_1,         'fTL_1',       .TRUE., 0. )                                                
         call fudge_out( 1,1,    fqw_1,         'fqw_1',       .TRUE., 0. )                                                
         call fudge_out( 1,1,SURF_HT_FLUX_LAND, 'SURF_HT_FLUX_LAND',  .TRUE., 0. )                                         
         call fudge_out( 1,1,    RESP_S_TILE,   'RESP_S_TILE', .TRUE., 0. )
         call fudge_out( 1,1,    g_leaf,        'g_leaf',      .TRUE., 0. )
         call fudge_out( 1,1,    NPP_ft,        'NPP_ft',      .TRUE., 0. )
         call fudge_out( 1,1,    GPP_ft,        'GPP_ft',      .TRUE., 0. )
         call fudge_out( 1,1,    RESP_S,        'RESP_S',      .TRUE., 0. )                        
         call fudge_out( 1,      RESP_S_tot,    'RESP_S_tot',  .TRUE., 0. )                              
         call fudge_out( 1,1,    RESP_P_FT,     'RESP_p_FT',   .TRUE., 0. )                              
         call fudge_out( 1,      RESP_p,        'RESP_p',      .TRUE., 0. )                                       
      endif

   IF(cable_user%run_diag_level == "BASIC")                           &     
      CALL basic_diag(subr_name, "Done.") 

!print *, "jhan:cable_imUN 7"
   !write( 6,'("End of cable_implicit_UNPACK")' )
END SUBROUTINE Implicit_unpack


