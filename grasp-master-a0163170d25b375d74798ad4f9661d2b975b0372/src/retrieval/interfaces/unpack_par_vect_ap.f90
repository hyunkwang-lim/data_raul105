! **
! **  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.  
! **  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
! **

! file contains :

! subroutine unpack_parameter_vector_ap
! subroutine get_SURF_wl
! subroutine get_REFI_wl
! subroutine get_SD_normalized
! subroutine get_AVP_altitudes_retr
! subroutine unpack_SURF
! subroutine unpack_C0
! subroutine unpack_SHD
! subroutine unpack_REFI
! subroutine unpack_AVP
! subroutine unpack_lidar_calibr_coeff
! subroutine unpack_SD
! subroutine unpack_AVP_std
! subroutine normalize_vector_par_single_pixel
!
! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

module mod_forward_model_characteristics

    use mod_par_inv, only : KW,KSHAPE,KBF,KIDIM2,KIDIM3,KVERTM       
    use mod_par_OS,  only : KSD !,KSD_clouds

    implicit none

!!PL Example of particle component type definition for future development	
!    type :: forward_model_characteristics_component 
!! forward_model_parameters --  type contains all perticle caracteristics driving forward model
!      integer                       ::  par_type_SD
!      real,dimension(KIDIM3,KIDIM2) ::	RIMAG,RREAL  
!      real,dimension(KVERTM)    ::  H0
!      real           			::  C0,sigma
!      integer			        ::  NBIN
!      !real,dimension(KPARS,KSD)	   ::	RADIUS,SD
!      real,dimension(KIDIM3)	::	RADIUS,SD
!      integer				    ::  NSHAPE
!      real,dimension(KSHAPE)	::	RATIO,SHD
!      integer                   ::  NHVP_retr
!      real,dimension(KVERTM)    ::  HVP_retr_km
!      real,dimension(KW)        ::  CL 
!      real       				::  fclouds
!    end type forward_model_characteristics_component

!PL Example of particle components type definition for future development
! 	type :: forward_model_characteristics_components_all
!	 type(forward_model_characteristics_component),dimension(KSD_clouds)		:: character_components 
!	end type forward_model_characteristics_components_all


    type :: forward_model_characteristics_particles 
! forward_model_parameters --  type contains all perticle caracteristics driving forward model
      integer                       ::  par_type_SD
      real,dimension(KIDIM3,KIDIM2) ::	RIMAG,RREAL  
      real,dimension(KVERTM,KSD)    ::  H0
      real,dimension(KSD)           ::  C0,sigma
      integer,dimension(KSD)        ::  NBIN
      !real,dimension(KPARS,KSD)	   ::	RADIUS,SD
      real,dimension(KIDIM3,KSD)	  ::	RADIUS,SD
      integer,dimension(KSD)	      ::  NSHAPE
      real,dimension(KSHAPE,KSD)	  ::	RATIO,SHD
      integer                       ::  NHVP_retr
      real,dimension(KVERTM)        ::  HVP_retr_km
      real,dimension(KW)            ::  CL 
      real,dimension(KSD)           ::  fclouds
    end type forward_model_characteristics_particles

    type :: forward_model_characteristics_surface 
! forward_model_parameters --  type contains all surface caracteristics driving forward model
      real,dimension(KIDIM3,KIDIM2) ::	BRF_land,BRP_land,BRM_water
    end type forward_model_characteristics_surface

contains   

!!PL Continue example with components
!      subroutine initialize_forward_model_characteristics_components(forw_components)
!
!		type(forward_model_characteristics_component),dimension(KSD_clouds)  ::  forw_components
!
!        forw_components(:)%par_type_SD = 0      
!        forw_components(:)%NBIN     = 0 
!        forw_components(:)%RADIUS(:)= 0.0      
!        forw_components(:)%SD(:)    = 0.0
!        forw_components(:)%C0       = 0.0
!            
!        forw_components(:)%NSHAPE  	= 0       
!        forw_components(:)%RATIO(:) = 0.0
!        forw_components(:)%SHD(:)   = 0.0
!      
!        forw_components(:)%RIMAG(:) = 0.0
!        forw_components(:)%RREAL(:) = 0.0
!                              
!        forw_components(:)%sigma  = 0.0
!        forw_components(:)%H0(:)   = 0.0
!        forw_components(:)%fclouds = 0.0
        
!        forw_components(:)%NHVP_retr      = 0
!        forw_components(:)%HVP_retr_km(:) = 0.0            
! CL - lidar calibration coefficient is not a particle characteristic. 
!! There is no better solution for this moment.
!        forw_components(:)%CL(:)          = 0.0 
!
!      end subroutine initialize_forward_model_characteristics_components


      subroutine initialize_forward_model_characteristics_particles(forw_particles)

        type(forward_model_characteristics_particles), intent(inout)  ::  forw_particles

        forw_particles%par_type_SD = 0      
        forw_particles%NBIN(:)     = 0 
        forw_particles%RADIUS(:,:) = 0.0      
        forw_particles%SD(:,:)     = 0.0
        forw_particles%C0(:)       = 0.0
            
        forw_particles%NSHAPE(:)  = 0       
        forw_particles%RATIO(:,:) = 0.0
        forw_particles%SHD(:,:)   = 0.0
      
        forw_particles%RIMAG(:,:) = 0.0
        forw_particles%RREAL(:,:) = 0.0
                              
        forw_particles%sigma(:)   = 0.0
        forw_particles%H0(:,:)    = 0.0
        forw_particles%fclouds(:) = 0.0
        
        forw_particles%NHVP_retr      = 0
        forw_particles%HVP_retr_km(:) = 0.0            
! CL - lidar calibration coefficient is not a particle characteristic. 
! There is no better solution for the moment.
        forw_particles%CL(:)          = 1.0

      end subroutine initialize_forward_model_characteristics_particles

      subroutine initialize_forward_model_characteristics_surface(forw_surface)

        type(forward_model_characteristics_surface), intent(inout)  ::  forw_surface

        forw_surface%BRF_land(:,:)  = 0.0
        forw_surface%BRP_land(:,:)  = 0.0
        forw_surface%BRM_water(:,:) = 0.0 
                                                                                                
      end subroutine initialize_forward_model_characteristics_surface      

! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      subroutine unpack_parameter_vector_ap(  RIN, APSING,                              &
                                              forw_aerosol, forw_clouds, forw_surface,  &
                                              pixel_fit, HGR_km, NHVP_meas, HVP_meas_km )      

      use mod_sdata_derived_type
      use mod_retr_settings_derived_type
      use mod_retr_general_output_derived_type
      
      implicit none
!	------------------------------------------------------------------------------------------------------   
      !integer,                            intent(in)     ::  iu_iP_main_output
      real,dimension(KPARS),              intent(in)     ::  APSING
      type(retr_input_settings),          intent(in)     ::  RIN
      type(forward_model_characteristics_particles),intent(inout)  ::  forw_aerosol ! Forward model chracteristics
      type(forward_model_characteristics_particles),intent(inout)  ::  forw_clouds  ! Forward model characteristics
      type(forward_model_characteristics_surface  ),intent(inout)  ::  forw_surface ! Forward model characteristics
      type(pixel),               optional,intent(in)     ::  pixel_fit
      real,                      optional,intent(in)     ::  HGR_km    ! site height above sea level 
      integer,                   optional,intent(in)     ::  NHVP_meas ! number of heights for lidar measurements
      real,dimension(KVERTM),    optional,intent(in)     ::  HVP_meas_km  ! heights for lidar measurements

!	------------------------------------------------------------------------------------------------------         
      integer  ::  IDIM1, par_type
      integer  ::  nbins1, nbins2, nbins
      logical  ::  icv ! presence of concentration characteristic
!	------------------------------------------------------------------------------------------------------
      forw_aerosol%par_type_SD = -999
      icv = .false.

      call initialize_forward_model_characteristics_surface(forw_surface)
      call initialize_forward_model_characteristics_particles(forw_aerosol)
      !if(RIN%NSD_clouds .gt. 0)  & 
      !call initialize_forward_model_characteristics_particles(forw_clouds)

! Unpack parameter vector ap (APSING) into characteristics driving forward model            
      do IDIM1=1,RIN%NDIM%n1
        par_type = RIN%NDIM%par_type(IDIM1)        
! Aerosol
        if(par_type .gt. par_type_SD_beg .and. par_type .lt. par_type_SD_end) then
            call unpack_SD ( RIN, par_type, IDIM1, APSING, forw_aerosol%NBIN, forw_aerosol%RADIUS, forw_aerosol%SD )   ! Size Distribution
            forw_aerosol%par_type_SD = par_type
        elseif(par_type .gt. par_type_RERI_beg .and. par_type .lt. par_type_RERI_end) then             
            call unpack_REFI ( RIN, par_type, IDIM1, APSING, forw_aerosol%RREAL ) ! Real part of refr.index n(wl) or Chemistry
        elseif(par_type .gt. par_type_IMRI_beg .and. par_type .lt. par_type_IMRI_end) then
            call unpack_REFI ( RIN, par_type, IDIM1, APSING, forw_aerosol%RIMAG ) ! Imaginary part of refr.index k(wl)
        elseif(par_type .gt. par_type_SHD_beg  .and. par_type .lt. par_type_SHD_end) then
            call unpack_SHD ( RIN, par_type, IDIM1, APSING, forw_aerosol%NSHAPE, forw_aerosol%RATIO, forw_aerosol%SHD ) ! Sphericity or Shape distribution
        elseif(par_type .gt. par_type_AVP_beg  .and. par_type .lt. par_type_AVP_end) then
            if(present(HGR_km) .and. present(NHVP_meas) .and. present(HVP_meas_km)) then
            forw_aerosol%NHVP_retr = RIN%NDIM%n3(1,IDIM1)
            call unpack_AVP ( RIN, par_type, IDIM1, APSING, RIN%NSD, forw_aerosol%H0 ) ! Aerosol Vertical Profile
            call get_AVP_altitudes_retr ( HGR_km, NHVP_meas, HVP_meas_km, &
                                          forw_aerosol%NHVP_retr, &
                                          forw_aerosol%HVP_retr_km )
            endif
        elseif(par_type .gt. par_type_Cv_beg .and. par_type .lt. par_type_Cv_end) then
            call unpack_C0 ( RIN, par_type, IDIM1, APSING, forw_aerosol%C0 )   ! Aerosol Concentration
            icv = .true.
        elseif(par_type .gt. par_type_ADD_beg  .and. par_type .lt. par_type_ADD_end) then
            if(par_type .eq. par_type_CL .and. present(pixel_fit))  &  ! Calibration coefficient
            call unpack_lidar_calibr_coeff ( RIN, pixel_fit, IDIM1, APSING, forw_aerosol%CL )
            if(par_type .eq. par_type_AVP_par_std)  &  ! Standard deviation of vertical profile
            call unpack_AVP_std ( RIN, par_type, IDIM1, APSING, forw_aerosol%sigma )
! Surface
        elseif(par_type .gt. par_type_SURF1_land_beg .and. par_type .lt. par_type_SURF1_land_end) then
            call unpack_SURF ( RIN, IDIM1, APSING, forw_surface%BRF_land )  ! Land surface ( RPV, Ross_Li, Litvinov)
        elseif(par_type .gt. par_type_SURF2_land_beg .and. par_type .lt. par_type_SURF2_land_end) then
            call unpack_SURF ( RIN, IDIM1, APSING, forw_surface%BRP_land )  ! Land surface ( Maignan, Litvinov_BPDF, Nadal-Breon)
        elseif(par_type .gt. par_type_SURF_water_beg .and. par_type .lt. par_type_SURF_water_end) then
            call unpack_SURF ( RIN, IDIM1, APSING, forw_surface%BRM_water )  ! Water surface ( Cox_Munk_iso, Cox_Munk_ani, Litvinov)
! Clouds
!        elseif(par_type .gt. par_type_SD_clouds_beg .and. par_type .lt. par_type_SD_clouds_end) then
!            call unpack_SD ( RIN,par_type,IDIM1,APSING,forw_clouds%NBIN,forw_clouds%RADIUS,forw_clouds%SD )   ! Cloud Size Distribution
!            forw_clouds%par_type_SD = par_type
!        elseif(par_type .gt. par_type_RERI_clouds_beg .and. par_type .lt. par_type_RERI_clouds_end) then
!            call unpack_REFI ( RIN,par_type,IDIM1,APSING,forw_clouds%RREAL ) ! Real part of refr.index n(wl) or Chemistry
!        elseif(par_type .gt. par_type_IMRI_clouds_beg .and. par_type .lt. par_type_IMRI_clouds_end) then
!            call unpack_REFI ( RIN,par_type,IDIM1,APSING,forw_clouds%RIMAG ) ! Imaginary part of refr.index k(wl)
!        elseif(par_type .gt. par_type_SHD_clouds_beg  .and. par_type .lt. par_type_SHD_clouds_end) then
!            call unpack_SHD ( RIN,par_type,IDIM1,APSING,forw_clouds%NSHAPE,forw_clouds%RATIO,forw_clouds%SHD ) ! Sphericity or Shape distribution
!        elseif(par_type .gt. par_type_CVP_beg  .and. par_type .lt. par_type_CVP_end) then
!            select case(par_type)
!            case(par_type_CVP_par_height)
!              forw_clouds%NHVP_retr = 1
!            case(par_type_CVP_prof)
!              write(*,*) 'IDIM1=',IDIM1,'  par_type =',par_type,'  - value is not supported for clouds'
!              stop 'stop in forw_IMAGE_I'
!            end select
!            call unpack_AVP ( RIN,par_type,IDIM1,APSING,RIN%NSD_clouds,forw_clouds%H0 ) ! Cloud Vertical Profile
!        elseif(par_type .gt. par_type_Cv_clouds_beg .and. par_type .lt. par_type_Cv_clouds_end) then
!            call unpack_C0 ( RIN,par_type,IDIM1,APSING,forw_clouds%C0 )   ! Cloud particle Concentration
!        elseif(par_type .gt. par_type_ADD_clouds_beg  .and. par_type .lt. par_type_ADD_clouds_end) then
!            if(par_type .eq. par_type_CVP_par_std)  &  ! Standard deviation of vertical profile
!            call unpack_AVP_std ( RIN,par_type,IDIM1,APSING,forw_clouds%sigma )
!        elseif(par_type .gt. par_type_FRACT_clouds_beg  .and. par_type .lt. par_type_FRACT_clouds_end) then
!            if(par_type .eq. par_type_fract_clouds)  &
!            call unpack_C0 ( RIN,par_type,IDIM1,APSING,forw_clouds%fclouds ) ! Cloud Fraction
        else
            write(*,*) 'IDIM1=',IDIM1,'  par_type =',par_type,'  - value is not valid'
            stop 'stop in unpack_parameter_vector_ap'
        endif ! par_type .gt. 
      enddo ! IDIM1

! Size Distribution normalization or concentration as a parameter of lognormal SD
      if( (forw_aerosol%par_type_SD .eq. par_type_SD_LN) .or. (icv .eqv. .true.) ) then
            call get_SD_normalized ( RIN, forw_aerosol%par_type_SD, forw_aerosol%C0, &
                                     forw_aerosol%NBIN, forw_aerosol%RADIUS, forw_aerosol%SD )
      endif
!      if( any(RIN%NDIM%par_type(1:RIN%NDIM%n1) .eq. par_type_SD_clouds_LN) .or. &
!                   any(RIN%NDIM%par_type(1:RIN%NDIM%n1) .eq. par_type_Cv_clouds) ) then
!            call get_SD_normalized ( RIN,forw_clouds%par_type_SD,forw_clouds%C0, &
!                                     forw_clouds%NBIN,forw_clouds%RADIUS,        &
!                                     forw_clouds%SD )
!      endif

! Extract precomputed lognormal bins for fine and coarse modes
      if( forw_aerosol%par_type_SD .eq. par_type_SD_LB .and. RIN%NSD .eq. 2 ) then
            nbins1 = forw_aerosol%NBIN(1) ! number of bins in fine mode
            nbins2 = forw_aerosol%NBIN(2) ! number of bins in coarse mode
            nbins =  nbins1 + nbins2
            forw_aerosol%RADIUS(nbins1+1:nbins,1) = forw_aerosol%RADIUS(1:nbins2,2)
            forw_aerosol%RADIUS(1:nbins,2)        = forw_aerosol%RADIUS(1:nbins,1)
            !forw_aerosol%SD(1:nbins1,1)   = forw_aerosol%SD(1:nbins1,1)
            forw_aerosol%SD(nbins1+1:nbins,2) = forw_aerosol%SD(1:nbins2,2)
            forw_aerosol%SD(1:nbins1,2)       = 0.0
            forw_aerosol%NBIN(1) = nbins
            forw_aerosol%NBIN(2) = nbins
!write(*,*) 'TEST LB NSD=2 in unpack_parameter_vector_ap'
!write(*,*) nbins1, nbins2, nbins,'  nbins1, nbins2, nbins'
!write(*,*) forw_aerosol%NBIN(1), forw_aerosol%NBIN(2), &
!        '  forw_aerosol%NBIN(1), forw_aerosol%NBIN(2)'
!write(*,*) 'forw_aerosol%RADIUS(1:nbins,1): ',forw_aerosol%RADIUS(1:nbins,1)
!write(*,*) 'forw_aerosol%RADIUS(1:nbins,2): ',forw_aerosol%RADIUS(1:nbins,2)
!write(*,*) 'forw_aerosol%SD(1:nbins,1): ',forw_aerosol%SD(1:nbins,1)
!write(*,*) 'forw_aerosol%SD(1:nbins,2): ',forw_aerosol%SD(1:nbins,2)
!stop 'stop test in unpack_parameter_vector_ap'
      endif

      return
      end subroutine unpack_parameter_vector_ap

! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

end module mod_forward_model_characteristics

!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! Surface parameters at IW-th wave length

      subroutine get_SURF_wl ( RIN,IW,ind_wl,BRF1,BRP1,BRM1,BRF,BRP,BRM )
      
      use mod_par_inv, only : KPARS,KIDIM2,KIDIM3,KBF
      use mod_retr_settings_derived_type 
           
      implicit none
! ----------------------------------------------------------------	  
      type(retr_input_settings),    intent(in)  ::  RIN
      integer,                      intent(in)  ::  IW
      integer,                      intent(in)  ::  ind_wl
      real,dimension(KIDIM3,KIDIM2),intent(in)  ::	BRP1,BRF1,BRM1
      real,dimension(KBF),          intent(out) ::	BRF,BRP,BRM

! ----------------------------------------------------------------	  
      integer  :: IDIM1,II,par_type
! ----------------------------------------------------------------	  
      BRP(:) = 0.0
      BRF(:) = 0.0
      BRM(:) = 0.0
      
      par_type = 0
      do IDIM1=1,RIN%NDIM%n1
        par_type = RIN%NDIM%par_type(IDIM1)
        if(par_type .gt. par_type_SURF1_land_beg .and. par_type .lt. par_type_SURF1_land_end) then             
          do II=1,RIN%NDIM%n2(IDIM1)
            BRF(II)=BRF1(ind_wl,II)
          enddo
        elseif(par_type .gt. par_type_SURF2_land_beg .and. par_type .lt. par_type_SURF2_land_end) then 
          do II=1,RIN%NDIM%n2(IDIM1)
            BRP(II)=BRP1(ind_wl,II)
          enddo
        elseif(par_type .gt. par_type_SURF_water_beg .and. par_type .lt. par_type_SURF_water_end) then 
          do II=1,RIN%NDIM%n2(IDIM1)
            BRM(II)=BRM1(ind_wl,II)
          enddo

        endif ! par_type .gt. par_type_SURF1_land_beg .and.
      enddo ! IDIM1

      if(par_type .eq. 0) then
        write(*,*) 'in get_SURF_wl: par_type=',par_type
        stop 'stop in get_SURF_wl'
      endif
          
      return
      end subroutine get_SURF_wl 

!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! Complex refractive index at IW-th wave length

      subroutine get_REFI_wl ( RIN,GOUT_chem_pixel,IW,ind_wl, &
                               RREALL,RIMAGL,RREAL,RIMAG )
      
      use mod_par_inv, only : KPARS,KIDIM2,KIDIM3,KW
      use mod_par_OS,  only : KSD
      use mod_retr_settings_derived_type
      use mod_retr_general_output_derived_type
      
      implicit none
! ----------------------------------------------------------------	  
      integer,                        intent(in)    ::  IW
      integer,                        intent(in)    ::  ind_wl
      type(retr_input_settings),      intent(in)    ::  RIN
      type(output_pixel_chemistry),   intent(inout) ::  GOUT_chem_pixel
      real,dimension(KIDIM3,KIDIM2),  intent(in)    ::  RIMAGL,RREALL      
      real,dimension(KSD),            intent(out)   ::  RREAL,RIMAG 

!	------------------------------------------------------------------------------------------------------
      real (selected_real_kind(p=15)),dimension(RIN%NW) :: RREALMM,RIMAGMM
      real (selected_real_kind(p=15)),dimension(RIN%NW) :: WAVELM
!      real (selected_real_kind(p=15))	::	rh_mxtr       ! Valid range: 0--1
!      real (selected_real_kind(p=15))	::	fract_inslbl  ! Volume fraction of insoluble inclusions
!      real (selected_real_kind(p=15))	::	fract_inslbl1 ! Volume fraction of insoluble inclusions 
!      real (selected_real_kind(p=15))	::	fract_carbon  ! Volume fraction of mixture that is composed of soot+Brc .
!      real (selected_real_kind(p=15))	::	fract_iron    ! Volume fraction of mixture that is composed of iron (aka Hematite).
!      real (selected_real_kind(p=15))	::	fract_slbl    ! Volume fraction of soluble inclusions
!      real (selected_real_kind(p=15))	::	fract_wtr     ! Volume fraction of mixture that is composed of water.
!      character(12)                   ::  instrument
!**********************************14/11/2016
      real (selected_real_kind(p=15))	::	rh_mxtr_f, rh_mxtr_c, rh_mxtr       ! Valid range: 0--1
      real (selected_real_kind(p=15))	::	fract_bc, fract_brc     ! Volume fraction of soot and brown carbon .
      real (selected_real_kind(p=15))	::	fract_inslbl_f, fract_inslbl_c  ! Volume fraction of insoluble inclusions
      real (selected_real_kind(p=15))	::	fract_quartz_f, fract_quartz_c ! Volume fraction of insoluble inclusions
      real (selected_real_kind(p=15))	::	fract_iron    ! Volume fraction of mixture that is composed of iron (aka Hematite).
      real (selected_real_kind(p=15))	::	fract_slbl_f, fract_slbl_c   ! Volume fraction of soluble inclusions
      real (selected_real_kind(p=15))	::	fract_wtr_f, fract_wtr_c    ! Volume fraction of mixture that is composed of water.
      character(12)                   ::  instrument
!**********************************14/11/2016
!	------------------------------------------------------------------------------------------------------
      integer :: par_type
!tl      real,dimension(KW+1,3)    ::	REALM, RIMAGM
      integer :: II, IDIM1, ISD, ind_wl1
      real :: AA, BB, AAT
      integer :: nfract, ifract
      real, dimension(5) :: vfract
      character(12), dimension(5) :: chem_name
! ----------------------------------------------------------------
      RREAL(:) = 0.0
      RIMAG(:) = 0.0
      par_type = 0
      do IDIM1=1,RIN%NDIM%n1
        par_type = RIN%NDIM%par_type(IDIM1)
        if(par_type .gt. par_type_RERI_beg .and. par_type .lt. par_type_RERI_end) then
            select case( par_type )
            case( par_type_RERI_spect )
              do ISD=1,RIN%NSD
                RREAL(ISD) = RREALL(ind_wl,ISD)
                RIMAG(ISD) = RIMAGL(ind_wl,ISD)
              enddo ! ISD
            case( par_type_RERI_const )
              ind_wl1 = 1
              do ISD=1,RIN%NSD
                RREAL(ISD) = RREALL(ind_wl1,ISD)
                RIMAG(ISD) = RIMAGL(ind_wl1,ISD)
              enddo ! ISD
            case( par_type_CXRI_nmix )
              WAVELM(1:RIN%NW)=RIN%WAVE(1:RIN%NW)
              ind_wl1 = 1
              do ISD=1,RIN%NSD
                nfract = RIN%NDIM%n3(ISD,IDIM1)
! compute normalized volume fractions vfract(1:nfract) of chem components 
! chem_name(1:nfract) and set them to generalized output structure
                call mixture_chemical_components( RIN%NSD,ISD,nfract,RREALL(1:nfract,ISD), &
                                                  GOUT_chem_pixel,vfract(1:nfract),chem_name(1:nfract) )
! compute spectral complex refractive index
                call complex_refr_index_practical_mixture( nfract,vfract(1:nfract),chem_name(1:nfract), &
                                                           ind_wl1,WAVELM(ind_wl), &
                                                           RREALMM(ind_wl1),RIMAGMM(ind_wl1) )
                RREAL(ISD) = RREALMM(ind_wl1)
                RIMAG(ISD) = RIMAGMM(ind_wl1)
                !write(*,'(/,f10.4,i4,f9.4,e12.4,12x,a,/)') WAVELM(ind_wl), &
                !ISD,RREAL(ISD),RIMAG(ISD),'  - wl,isd,refr,refi      in get_REFI_wl'

              enddo ! ISD
!stop 'stop: delete after testing practical mixture complex refractive index'
!            case( par_type_CXRI_chem )
!              do ISD=1,RIN%NSD
!                rh_mxtr      = RREALL(1,ISD)
!                fract_inslbl = RREALL(2,ISD)
!                fract_carbon   = RREALL(3,ISD)
!                fract_iron   = RREALL(4,ISD)
!                fract_inslbl1 = fract_inslbl*(1.0-fract_carbon-fract_iron)
!                GOUT_chem_pixel%finslbl(ISD) = fract_inslbl1
!                GOUT_chem_pixel%fsoot(ISD)   = fract_carbon
!                GOUT_chem_pixel%firon(ISD)   = fract_iron
!                GOUT_chem_pixel%rh(ISD)      = rh_mxtr
!
!                WAVELM(1:RIN%NW)=RIN%WAVE(1:RIN%NW)
!                instrument = 'user_defined'  ! aeronet, parasol, eumetsat, or user_defined.
!!                 write(*,*) rh_mxtr,RREALL(1,ISD),RREALL(2,ISD),RREALL(3,ISD),RREALL(4,ISD), fract_inslbl1
!!                 stop
!!                 ' rh_mxtr,RREALL1,RREALL2,RREALL3,RREALL4 inget_REFI_wl'
!
!! Wavelength assumptions for 'instrument':
!! user_defined:  wvln_um   can be anything, between 0.25 and 4 um.
!!     aeronet :  wvln_um = (/  .440,                .675,  .870, 1.020                      /)
!!     parasol :  wvln_um = (/  .443,  .490,  .565,  .675,  .870, 1.020                      /)
!!     eumetsat:  wvln_um = (/  .443,  .510,  .640,         .865,        1.380, 1.640, 2.130 /)
!
!!CD        WRITE(*,*) 'before mixing'
!                call MIXING ( instrument,RIN%NW,WAVELM,rh_mxtr,   & 
!                              fract_inslbl1,fract_carbon,           &
!                              fract_iron,fract_slbl,fract_wtr,    &
!                              RREALMM,RIMAGMM                     &
!                            )
!                GOUT_chem_pixel%fwtr(ISD)  = fract_wtr
!                GOUT_chem_pixel%fslbl(ISD) = fract_slbl
!
!                RREAL(ISD) = RREALMM(ind_wl)
!                RIMAG(ISD) = RIMAGMM(ind_wl)
!!               WRITE(*,*) RREAL(ISD),RIMAG(ISD),ISD,IW,' RREAL, RIMAG,ISD,IW'
!!               WRITE(*,*) 'WAVELM(1:RIN%NW): ',WAVELM(1:RIN%NW)
!


!********************************14/11/2016
            case( par_type_CXRI_chem )
                instrument = 'user_defined'  ! aeronet, parasol, eumetsat, or user_defined.

              do ISD=1,RIN%NSD
                rh_mxtr = RREALL(1,1)
              if ( ISD .eq. 1) then
!                rh_mxtr_f   = RREALL(1,ISD)
                fract_bc   = RREALL(2,ISD)
                fract_brc   = RREALL(3,ISD)
                fract_inslbl_f = RREALL(4,ISD)
                fract_quartz_f = fract_inslbl_f*(1.0-fract_bc-fract_brc)

                GOUT_chem_pixel%fquartz(ISD) = fract_quartz_f
                GOUT_chem_pixel%fsoot(ISD)   = fract_bc
                GOUT_chem_pixel%fbrc(ISD)    = fract_brc
                GOUT_chem_pixel%rh(ISD)      = rh_mxtr

                GOUT_chem_pixel%fwtr(ISD)  = fract_wtr_f
                GOUT_chem_pixel%fslbl(ISD) = fract_slbl_f

                WAVELM(1:RIN%NW)=RIN%WAVE(1:RIN%NW)

                call MIXING_f ( instrument,RIN%NW,WAVELM, ISD, rh_mxtr,   &
                              fract_quartz_f,fract_bc,           &
                              fract_brc,fract_slbl_f,fract_wtr_f,    &
                              RREALMM,RIMAGMM )

                RREAL(ISD) = RREALMM(ind_wl)
                RIMAG(ISD) = RIMAGMM(ind_wl)

               else
!                rh_mxtr_c   = RREALL(1,ISD)
                fract_iron   = RREALL(1,ISD)
                fract_inslbl_c = RREALL(2,ISD)
                fract_quartz_c = fract_inslbl_c*(1.0-fract_iron)

                GOUT_chem_pixel%fquartz(ISD) = fract_quartz_c
                GOUT_chem_pixel%firon(ISD)   = fract_iron
                GOUT_chem_pixel%rh(ISD)      = rh_mxtr
                GOUT_chem_pixel%fwtr(ISD)  = fract_wtr_c
                GOUT_chem_pixel%fslbl(ISD) = fract_slbl_c

                WAVELM(1:RIN%NW)=RIN%WAVE(1:RIN%NW)

                call MIXING_c ( instrument,RIN%NW,WAVELM, ISD, rh_mxtr,   &
                              fract_quartz_c,   &
                              fract_iron,fract_slbl_c,fract_wtr_c,    &
                              RREALMM,RIMAGMM )

                RREAL(ISD) = RREALMM(ind_wl)
                RIMAG(ISD) = RIMAGMM(ind_wl)

                endif
!********************************14/11/2016




              enddo ! ISD
            end select
        exit
        endif ! par_type .gt. par_type_RERI_beg
        if(par_type .eq. 0) then
          write(*,*) 'in get_REFI_wl: par_type=',par_type
          stop 'stop in get_REFI_wl'
        endif
      enddo ! IDIM1
                
      return
      end subroutine get_REFI_wl 

!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! Surface parameters 

      subroutine unpack_SURF ( RIN,IDIM1,APSING1,B1 )
      
      use mod_par_inv, only : KPARS,KIDIM2,KIDIM3
      use mod_retr_settings_derived_type 
           
      implicit none
! ----------------------------------------------------------------	  
      type(retr_input_settings),    intent(in)  ::  RIN
      integer,                      intent(in)  ::  IDIM1
      real,dimension(KPARS),        intent(in)  ::  APSING1
      real,dimension(KIDIM3,KIDIM2),intent(out) ::	B1

! ----------------------------------------------------------------	  
      integer  :: IDIM2,IDIM3,II
! ----------------------------------------------------------------	  
      B1(:,:) = 0.0
            
      do IDIM2=1,RIN%NDIM%n2(IDIM1) 
        do IDIM3=1,RIN%NDIM%n3(IDIM2,IDIM1)
            II = RIN%NDIM%ISTARSING(IDIM2,IDIM1)+IDIM3-1                     
            if(RIN%NDIM%n3(IDIM2,IDIM1) .eq. 1) then
              B1(1:RIN%NW,IDIM2) = APSING1(II)
            else
              B1(IDIM3,IDIM2) = APSING1(II)
            endif
        enddo ! IDIM3
      enddo ! IDIM2   
          
      return
      end subroutine unpack_SURF 

!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! Aerosol concentration

      subroutine unpack_C0 ( RIN,par_type,IDIM1,APSING1,C0 )
      
      use mod_par_inv, only : KPARS
      use mod_par_OS,  only : KSD
      use mod_retr_settings_derived_type
            
      implicit none
! ----------------------------------------------------------------	  
      type(retr_input_settings),  intent(in)    ::  RIN
      integer,                    intent(in)    ::  par_type,IDIM1
      real,dimension(KPARS),      intent(in)    ::  APSING1
      real,dimension(KSD),        intent(out)   ::  C0      
! ----------------------------------------------------------------	  
      integer              :: IDIM2,IDIM3,II
! ----------------------------------------------------------------	  
      C0(:)   = 0.0
      
      do IDIM2=1,RIN%NDIM%n2(IDIM1) 
        do IDIM3=1,RIN%NDIM%n3(IDIM2,IDIM1)
            II = RIN%NDIM%ISTARSING(IDIM2,IDIM1)+IDIM3-1                     
            C0(IDIM2) = APSING1(II)
!write(*,*) par_type,'  - par_type, II=',II,'  IDIM2=',IDIM2,'  C0(IDIM2)=',C0(IDIM2),'  in unpack_C0'      
        enddo ! IDIM3
      enddo ! IDIM2   

      return
      end subroutine unpack_C0 

!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! Altitudes for retrieved Aerosol Vertical Profile

      subroutine get_AVP_altitudes_retr ( HGR_km,                & 
                                          NHVP_meas,HVP_meas_km, &
                                          NHVP_retr,HVP_retr_km )
      
      use mod_par_inv, only : KVERTM
      use mod_retr_settings_derived_type
            
      implicit none
! ----------------------------------------------------------------	  
      integer,                    intent(in)    ::  NHVP_meas
      real,dimension(KVERTM),     intent(in)    ::  HVP_meas_km
      real,                       intent(in)    ::  HGR_km
      integer,                    intent(inout) ::  NHVP_retr
      real,dimension(KVERTM),     intent(out)   ::  HVP_retr_km
! ----------------------------------------------------------------	  
      integer  ::  i
! ----------------------------------------------------------------	  
      HVP_retr_km(:) = 0.0
      
      if (NHVP_retr .eq. 1) then
         if(NHVP_meas .ne. NHVP_retr) then
         write(*,*) 'NHVP_meas=',NHVP_meas,' .ne. ','NHVP_retr=',NHVP_retr
         stop 'stop in forw_IMAGE_I'
         endif
         NHVP_retr = NHVP_meas
         HVP_retr_km(1:NHVP_retr) = HVP_meas_km(1:NHVP_meas)      
      else
         if(NHVP_retr-NHVP_meas .gt. 1 .or. NHVP_retr .lt. NHVP_meas) then
            write(*,*) 'NHVP_retr=',NHVP_retr,'  NHVP_meas=',NHVP_meas, &
            '  - inconsistent number of heights'
            stop 'stop in get_AVP_altitudes_retr'
         endif
         HVP_retr_km(1:NHVP_meas) = HVP_meas_km(1:NHVP_meas)
         if(NHVP_retr .gt. NHVP_meas) then
            NHVP_retr              = NHVP_meas + 1
            HVP_retr_km(NHVP_retr) = HGR_km 
            if(HVP_retr_km(NHVP_retr-1) .lt. HGR_km) then
              write(*,*) 'HVP_retr_km(NHVP_retr-1)=',HVP_retr_km(NHVP_retr-1),  &
              ' can not be less than   HGR_km=',HGR_km
              stop 'stop in get_AVP_altitudes_retr'              
            endif
         endif
      endif

      !write(*,*) 'NHVP_retr=',NHVP_retr,'  NHVP_meas=', NHVP_meas
      !do i=1,NHVP_retr
      !write(*,*) i,HVP_retr_km(i),H0(i,1:2),'  - i,HVP_retr_km(i),H0(i,1:NSD)'
      !enddo
      !stop
      
      return
      end subroutine get_AVP_altitudes_retr 

!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! Spectral calibration coefficient for lidar

      subroutine unpack_lidar_calibr_coeff ( RIN,pixel_fit,IDIM1,APSING1,CL )
      
      use mod_par_inv, only : KPARS,KW,KVERTM
      use mod_par_OS,  only : KSD
      use mod_retr_settings_derived_type
      use mod_sdata
      
      implicit none
! ----------------------------------------------------------------	  
      type(retr_input_settings),  intent(in)    ::  RIN
      integer,                    intent(in)    ::  IDIM1
      real,dimension(KPARS),      intent(in)    ::  APSING1
      type(pixel),                intent(in)    ::  pixel_fit
      real,dimension(KW),         intent(inout)   ::  CL       
! ----------------------------------------------------------------	  
      real,dimension(KW)   :: CL_temp 
      integer              :: IDIM2,IDIM3,II,IW,IP,IW1
! ----------------------------------------------------------------	  
      !CL(:)   = 0.0
      
      do IDIM2=1,RIN%NDIM%n2(IDIM1) 
        if(IDIM2 .gt. 1) then
          write(*,*) 'in unpack_lidar_calibr_coeff: IDIM2=',IDIM2,' can not be > 1'
          stop 'stop in unpack_lidar_calibr_coeff'
        endif
        do IDIM3=1,RIN%NDIM%n3(IDIM2,IDIM1)
            II = RIN%NDIM%ISTARSING(IDIM2,IDIM1)+IDIM3-1                     
            CL(IDIM3) = APSING1(II) ! lidar calibration parameters
!write(*,*) 'in unpack_lidar_calibr_coeff: II=',II,'IDIM3=',IDIM3,'  IDIM2=',IDIM2,'  CL(IDIM3)=',CL(IDIM3)
!write(*,*) 'in unpack_lidar_calibr_coeff: II=',II, 'APSING=', APSING1(II)
        enddo ! IDIM3
      enddo ! IDIM2   

! assingning CL to the correct wavelength in a full set of WL 
! write(*,*) 'in unpack_lidar_calibr_coeff: CL: ',CL(1:RIN%NW)      
      IW1=1
      CL_temp(:)=CL(:)
      do IW=1,RIN%NW
! write(*,*) 'in unpack_lidar_calibr_coeff: NIP=',pixel_fit%meas(IW)%NIP      
        do IP=1,pixel_fit%meas(IW)%NIP
          if(pixel_fit%meas(IW)%meas_type(IP) .ge. meas_type_LS .and. &
            pixel_fit%meas(IW)%meas_type(IP) .le. meas_type_RL ) then
            CL(IW)=CL_temp(IW1)
            IW1=IW1+1
          else
            CL(IW)=1.0
          endif ! meas_type .ge. meas_type_LS .and.
        enddo ! IP
      enddo ! IW                            
!write(*,'(a,10e14.6)') 'in unpack_lidar_calibr_coeff: CL: ',CL(1:RIN%NW)      
!stop          
      
      return
      end subroutine unpack_lidar_calibr_coeff 
!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! Profile standard deviation
      subroutine unpack_AVP_std ( RIN, par_type, IDIM1, APSING1, sigma )
      
      use mod_par_inv, only : KPARS
      use mod_par_OS,  only : KSD
      use mod_retr_settings_derived_type
            
      implicit none
! ----------------------------------------------------------------	  
      type(retr_input_settings),  intent(in)    ::  RIN
      integer,                    intent(in)    ::  par_type, IDIM1
      real,dimension(KPARS),      intent(in)    ::  APSING1
      real,dimension(KSD),        intent(out)   ::  sigma      
! ----------------------------------------------------------------	  
      integer              :: IDIM2,IDIM3,II
! ----------------------------------------------------------------	  
      sigma(:)   = 0.0
      
      do IDIM2=1,RIN%NDIM%n2(IDIM1) 
        do IDIM3=1,RIN%NDIM%n3(IDIM2,IDIM1)
            II = RIN%NDIM%ISTARSING(IDIM2,IDIM1)+IDIM3-1                     
            if(RIN%NDIM%n2(IDIM1) .lt. RIN%NSD) then
              sigma(1:RIN%NSD) = APSING1(II)
            else
              sigma(IDIM2) = APSING1(II)
            endif ! RIN%NDIM%n2(IDIM1) .lt. RIN%NSD
!write(*,*) par_type,'  - par_type, II=',II,'  IDIM2=',IDIM2,'  sigma(IDIM2)=',sigma(IDIM2),'  in unpack_AVP_par_std'      
        enddo ! IDIM3
      enddo ! IDIM2   

      return
      end subroutine unpack_AVP_std 

!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! Aerosol or Cloud Vertical Profile 

      subroutine unpack_AVP ( RIN,par_type,IDIM1,APSING1,NSD,H0 )
      
      use mod_par_inv, only : KPARS,KVERTM
      use mod_par_OS,  only : KSD
      use mod_retr_settings_derived_type
      
      implicit none
! ----------------------------------------------------------------	  
      type(retr_input_settings),  intent(in)    ::  RIN
      integer,                    intent(in)    ::  NSD,par_type,IDIM1
      real,dimension(KPARS),      intent(in)    ::  APSING1
      real,dimension(KVERTM,KSD), intent(out)   ::  H0
! ----------------------------------------------------------------	  
      integer :: IDIM2, IDIM3, II
! ----------------------------------------------------------------
      H0(:,:) = 0.0

      do IDIM2=1,RIN%NDIM%n2(IDIM1) 
        do IDIM3=1,RIN%NDIM%n3(IDIM2,IDIM1)
            II = RIN%NDIM%ISTARSING(IDIM2,IDIM1)+IDIM3-1                     
            if(RIN%NDIM%n2(IDIM1) .lt. NSD) then
              H0(IDIM3,1:NSD) = APSING1(II)
!write(*,*) '1: II=',II,'  IDIM3',IDIM3,'  IDIM2=',IDIM2,'  H0=',APSING1(II)
            else
              H0(IDIM3,IDIM2) = APSING1(II)
!write(*,*) '2: II=',II,'  IDIM3',IDIM3,'  IDIM2=',IDIM2,'  H0=',APSING1(II)
            endif ! RIN%NDIM%n2(IDIM1) .lt. NSD
        enddo ! IDIM3
      enddo ! IDIM2   
      
      return
      end subroutine unpack_AVP 

!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! Shape of aerosol particles

      subroutine unpack_SHD ( RIN,par_type,IDIM1,APSING1,NSHAPE1,RATIO,SHD )
      
      use mod_par_inv, only : KPARS,KIDIM2,KIDIM3,KSHAPE
      use mod_par_OS,  only : KSD
      use mod_retr_settings_derived_type
      
      implicit none
! ----------------------------------------------------------------	  
      type(retr_input_settings),  intent(in)  ::  RIN
      integer,                    intent(in)  ::  par_type,IDIM1
      real,dimension(KPARS),      intent(in)  ::  APSING1
      real,dimension(KSHAPE,KSD), intent(out) ::	RATIO, SHD
      integer,dimension(KSD),     intent(out) ::	NSHAPE1      
! ----------------------------------------------------------------	  
      real,dimension(KIDIM3,KIDIM2)        :: SHDL
      integer,dimension(KSD)               :: NSHAPE
      integer                              :: IDIM2,IDIM3,II,ISD
! ----------------------------------------------------------------	  
      RATIO(:,:) = 0.0
      NSHAPE1(:) = 0.0
      SHD(:,:)   = 0.0
      
      do IDIM2=1,RIN%NDIM%n2(IDIM1) 
        do IDIM3=1,RIN%NDIM%n3(IDIM2,IDIM1)
            II = RIN%NDIM%ISTARSING(IDIM2,IDIM1)+IDIM3-1                     
            if(RIN%NDIM%n2(IDIM1) .lt. RIN%NSD) then
              NSHAPE(1:RIN%NSD) = RIN%NDIM%n3(IDIM2,IDIM1)
              SHDL(IDIM3,1:RIN%NSD) = APSING1(II)
            else
              NSHAPE(IDIM2) = RIN%NDIM%n3(IDIM2,IDIM1)
              SHDL(IDIM3,IDIM2) = APSING1(II)
            endif ! RIN%NDIM%n3(IDIM2,IDIM1) .ne. RIN%NSD
          
          if(par_type .eq. par_type_SHD_fsph) then ! retrieve % spherical particles
              NSHAPE1(IDIM2) = 2
              SHDL(2,IDIM2) = 1.0-SHDL(1,IDIM2)
              if(RIN%NDIM%n2(IDIM1) .lt. RIN%NSD) then
                NSHAPE1(1:RIN%NSD) = NSHAPE1(IDIM2)
                do ISD=1,RIN%NSD
                RATIO(1:NSHAPE1(IDIM2),ISD) = RIN%RATIO1(1:NSHAPE1(IDIM2),IDIM2)
                SHD(1:NSHAPE1(IDIM2),ISD)   = SHDL(1:NSHAPE1(IDIM2),IDIM2)
                enddo ! ISD
              else
                RATIO(1:NSHAPE1(IDIM2),IDIM2) = RIN%RATIO1(1:NSHAPE1(IDIM2),IDIM2)
                SHD(1:NSHAPE1(IDIM2),IDIM2)   = SHDL(1:NSHAPE1(IDIM2),IDIM2)
              endif
          else if(par_type .eq. par_type_SHD_distr) then ! retrieve axis ratio distribution
              if(RIN%NDIM%n2(IDIM1) .lt. RIN%NSD) then
                NSHAPE1(1:RIN%NSD) = NSHAPE(IDIM2)                          
                do ISD=1,RIN%NSD
                RATIO(1:NSHAPE1(IDIM2),ISD) = RIN%RATIO1(1:NSHAPE1(IDIM2),IDIM2)
                SHD(1:NSHAPE1(IDIM2),ISD)   = SHDL(1:NSHAPE1(IDIM3),IDIM2)
                enddo ! ISD
              else
                NSHAPE1(IDIM2) = NSHAPE(IDIM2)                          
                RATIO(IDIM3,IDIM2) = RIN%RATIO1(IDIM3,IDIM2)
                SHD(IDIM3,IDIM2)   = SHDL(IDIM3,IDIM2)              
              endif
          endif ! par_type .eq.
!           write(*,*) 'test: RATIO :',RATIO(1:NSHAPE1(IDIM2),IDIM2),' SHD: ',SHD(1:NSHAPE1(IDIM2),IDIM2)
!OD         do II=1,NSHAPE1(IDIM2)
!OD         write(*,*) II,SHD(II,IDIM2),'IB,SHD(II,IDIM2)'
!OD         enddo
        enddo ! IDIM3
      enddo ! IDIM2   
          
      return
      end subroutine unpack_SHD 

!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! Real part of complex refractive index n(wl) or Chemistry components

      subroutine unpack_REFI ( RIN,par_type,IDIM1,APSING1,REFI )
      
      use mod_par_inv, only : KPARS,KIDIM2,KIDIM3
      use mod_par_OS,  only : KSD
      use mod_retr_settings_derived_type
      
      implicit none
! ----------------------------------------------------------------	  
      type(retr_input_settings),  intent(in)  ::  RIN
      integer,                    intent(in)  ::  par_type,IDIM1
      real,dimension(KPARS),      intent(in)  ::  APSING1
      real,dimension(KIDIM3,KIDIM2),intent(out)  ::	REFI
! ----------------------------------------------------------------	  
      integer                              :: IDIM2,IDIM3,II
! ----------------------------------------------------------------	  
      REFI(:,:) = 0.0
     
      do IDIM2=1,RIN%NDIM%n2(IDIM1) 
        do IDIM3=1,RIN%NDIM%n3(IDIM2,IDIM1)
          II = RIN%NDIM%ISTARSING(IDIM2,IDIM1)+IDIM3-1                     
          if(RIN%NDIM%n2(IDIM1) .lt. RIN%NSD) then
            REFI(IDIM3,1:RIN%NSD) = APSING1(II)
          else
            REFI(IDIM3,IDIM2) = APSING1(II)
          endif ! RIN%NDIM%n2(IDIM1) .lt. RIN%NSD
        enddo ! IDIM3
      enddo ! IDIM2   
          
      return
      end subroutine unpack_REFI 

!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! Notmalization of Size Distribution for forward model

      subroutine get_SD_normalized ( RIN,par_type,C0,NBIN,RADIUS,SD )
      
      use mod_par_inv, only : KIDIM3
      use mod_par_OS,  only : KSD
      use mod_retr_settings_derived_type
      
      implicit none
! ----------------------------------------------------------------	  
      type(retr_input_settings),  intent(in)    ::  RIN
      integer,                    intent(in)    ::  par_type
      integer,dimension(KSD),     intent(in)    ::	NBIN
      real,dimension(KIDIM3,KSD), intent(in)    ::	RADIUS      
      real,dimension(KSD),        intent(in)    ::	C0
      real,dimension(KIDIM3,KSD), intent(inout) ::	SD 
! ----------------------------------------------------------------	  
      real      :: AA
      integer   :: ISD
! ----------------------------------------------------------------	  
      AA = 0.0
      do ISD=1,RIN%NSD
        if(par_type .eq. par_type_SD_LN) then 
          SD(NBIN(ISD)+1,ISD) = C0(ISD)
!          write(*,*) 'in unpack_SD_normalized 1: ',ISD,C0(ISD)
!          write(*,*) 'in unpack_SD_normalized 2: ',ISD,SD(1:NBIN(ISD)+1,ISD)
        else
          AA = SUM(SD(1:NBIN(ISD),ISD))
          if(par_type .eq. par_type_SD_TB) AA = AA*( LOG(RADIUS(2,ISD))-LOG(RADIUS(1,ISD)) )          
          SD(1:NBIN(ISD),ISD) = SD(1:NBIN(ISD),ISD)/AA*C0(ISD)
!          write(*,*) 'in unpack_SD_normalized 1: ',ISD,C0(ISD),AA
!          write(*,*) 'in unpack_SD_normalized 2: ',ISD,SD(1:NBIN(ISD),ISD)
!          write(*,*) 'in unpack_SD_normalized 3: ',ISD,par_type,par_type_SD_TB
      endif
      enddo ! ISD	  
          
      return
      end subroutine get_SD_normalized

!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! Size Distribution for forward model

      subroutine unpack_SD ( RIN,par_type,IDIM1,APSING1,   &
                             NBIN,RADIUS,SD )
      
      use mod_par_inv, only : KPARS,KIDIM3 ! KIDIM2,
      use mod_par_OS,  only : KSD
      use mod_retr_settings_derived_type
      
      implicit none
! ----------------------------------------------------------------	  
      type(retr_input_settings),  intent(in)  ::  RIN
      integer,                    intent(in)  ::  par_type,IDIM1
      real,dimension(KPARS),      intent(in)  ::  APSING1
      integer,dimension(KSD),     intent(out) ::	NBIN
      real,dimension(KIDIM3,KSD), intent(out) ::	RADIUS,SD      
! ----------------------------------------------------------------	  
      real     ::  AA
      integer  ::  IDIM2,IDIM3,II
! ----------------------------------------------------------------	  
      NBIN(:)     = 0
      RADIUS(:,:) = 0.0
      SD(:,:)     = 0.0
      
      do IDIM2=1,RIN%NDIM%n2(IDIM1) 
        NBIN(IDIM2)=RIN%NDIM%n3(IDIM2,IDIM1)
        do IDIM3=1,NBIN(IDIM2)
            II = RIN%NDIM%ISTARSING(IDIM2,IDIM1)+IDIM3-1                     
            SD(IDIM3,IDIM2) = APSING1(II)
            RADIUS(IDIM3,IDIM2)=RIN%RADIUS1(IDIM3,IDIM2)
        enddo ! IDIM3
      enddo ! IDIM2  
          
      return
      end subroutine unpack_SD

!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! Normalize particular retrieved characteristics in vector ap (APSING) 
! for single pixel :
!     - size distribution
!     - shape distribution
!     - aerosol vertical profile
!     - chemistry, linear volume mixture

      subroutine normalize_vector_par_single_pixel(  RIN, segment_meas, ipix, APSING )

      use mod_par_inv, only : KVERTM
      use mod_par_OS, only : HMAX_atm
      use mod_retr_settings_derived_type
      use mod_sdata_derived_type
      use mod_sdata, only : get_HVP_lidar

      implicit none
!	------------------------------------------------------------------------------------------------------
      type(retr_input_settings), intent(in) :: RIN
      type(segment_data), intent(in) :: segment_meas
      integer, intent(in) :: ipix
      real, dimension(KPARS), intent(inout) :: APSING
!	------------------------------------------------------------------------------------------------------
      integer ::  IDIM1, IDIM2, par_type
      integer ::  NBIN, ibeg, iend
      real    ::  xnorm, dlnr, HGR_km
      integer :: NHVP_meas, NH
      real, dimension(KVERTM) :: HVP_meas_m, HVP_retr_km, H_km, prof_temp
      character(len=20) :: distr_type
!	------------------------------------------------------------------------------------------------------
      do IDIM1=1,RIN%NDIM%n1
        par_type = RIN%NDIM%par_type(IDIM1)        
        if(par_type .eq. par_type_SD_TB) then
        ! Size distribution (triangle bins)
            if(any(RIN%NDIM%par_type(1:RIN%NDIM%n1) .eq. par_type_Cv)) then
              do IDIM2=1,RIN%NDIM%n2(IDIM1)
                dlnr = log(RIN%RADIUS1(2,IDIM2))-log(RIN%RADIUS1(1,IDIM2))
                NBIN = RIN%NDIM%n3(IDIM2,IDIM1)
                ibeg = RIN%NDIM%ISTARSING(IDIM2,IDIM1)
                iend = ibeg + NBIN - 1
                xnorm = sum(APSING(ibeg:iend))
                APSING(ibeg:iend) = APSING(ibeg:iend) / (xnorm*dlnr)
              enddo ! IDIM2
            endif
        elseif(par_type .eq. par_type_SD_LB) then
        ! Size distribution (precomputed lognormal bins)
            if(any(RIN%NDIM%par_type(1:RIN%NDIM%n1) .eq. par_type_Cv)) then
              do IDIM2=1,RIN%NDIM%n2(IDIM1)
                NBIN = RIN%NDIM%n3(IDIM2,IDIM1)
                ibeg = RIN%NDIM%ISTARSING(IDIM2,IDIM1)
                iend = ibeg + NBIN - 1
                xnorm = sum(APSING(ibeg:iend))
                APSING(ibeg:iend) = APSING(ibeg:iend) / xnorm
              enddo ! IDIM2
            endif
        elseif(par_type .eq. par_type_SHD_distr) then
        ! Shape distribution
              do IDIM2=1,RIN%NDIM%n2(IDIM1)
                NBIN = RIN%NDIM%n3(IDIM2,IDIM1)
                ibeg = RIN%NDIM%ISTARSING(IDIM2,IDIM1)
                iend = ibeg + NBIN - 1
                xnorm = sum(APSING(ibeg:iend))
                APSING(ibeg:iend) = APSING(ibeg:iend) / xnorm
              enddo ! IDIM2
        elseif(par_type .eq. par_type_AVP_prof) then
        ! Aerosol Vertical Profile
              HGR_km = segment_meas%pixels(ipix)%MASL*0.001
              call get_HVP_lidar ( segment_meas, NHVP_meas, HVP_meas_m )
              distr_type = 'lut'
              do IDIM2=1,RIN%NDIM%n2(IDIM1)
                NBIN = RIN%NDIM%n3(IDIM2,IDIM1)
                ibeg = RIN%NDIM%ISTARSING(IDIM2,IDIM1)
                iend = ibeg + NBIN - 1
                call get_AVP_altitudes_retr ( HGR_km,                      &
                                              NHVP_meas, HVP_meas_m*0.001, &
                                              NBIN, HVP_retr_km )
                ! Altitudes for aerosol profile normalization
                call grid_altitudes_LUT ( HGR_km, HMAX_atm*0.001,     & ! IN
                                          NBIN, HVP_retr_km(1:NBIN),  &
                                          NH, H_km(1:NBIN+2)          & ! OUT
                                        )
                ! Aerosol concentration profile normalization
                call discrvd_single_atm_component ( distr_type, 0.0, 0.0,       &
                                                    NBIN, HVP_retr_km(1:NBIN),  &
                                                    APSING(ibeg:iend),          &
                                                    NH, H_km(1:NH),             &
                                                    xnorm, prof_temp(1:NH)      & ! OUT
                                                  )
                APSING(ibeg:iend) = APSING(ibeg:iend) / xnorm
              enddo ! IDIM2
        elseif(par_type .eq. par_type_CXRI_nmix) then
        ! Chemistry, linear volume mixture
              do IDIM2=1,RIN%NDIM%n2(IDIM1)
                NBIN = RIN%NDIM%n3(IDIM2,IDIM1)
                ibeg = RIN%NDIM%ISTARSING(IDIM2,IDIM1)
                iend = ibeg + NBIN - 1
                xnorm = sum(APSING(ibeg:iend))
                APSING(ibeg:iend) = APSING(ibeg:iend) / xnorm
              enddo ! IDIM2
        endif ! par_type .eq.
      enddo ! IDIM1

      return
      end subroutine normalize_vector_par_single_pixel

!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
