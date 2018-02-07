! **
! **  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.  
! **  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
! **

!	sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine retr_input_initialization ( RIN )

      use mod_retr_settings_derived_type

      implicit none
! -----------------------------------------------------------------------------------------
	  
      type(retr_input_settings), intent(out) :: RIN
      integer   :: i	  
! -----------------------------------------------------------------------------------------
  
      RIN%KNSING  = 0 
      RIN%KNSINGF = 0 
      RIN%KL      = 0
      RIN%ISTOP   = .false. 
      RIN%IPRI_additional_info   = .false. 
      RIN%IPRI_iter_result = .false. 
      RIN%IPRI_verbose = .true.  

      RIN%NSD     = 0
      RIN%NLYRS   = 0 
      RIN%ipplane = 0 
! Polarization fitting option (and principle plane linear polarization measurements)
! Stokes parameter vector
  ! =1 Q,U; =2 Q/I,U/I; =3 sqrt(Q*Q+U*U); 
  ! =4 sqrt(Q*Q+U*U)/I; =5 P/I meas in sdata
! Phase matrix
  ! =1 p12; =2 -p12/I
! Lidar
  ! =1 LPER,LPAR; =2 LPER/LPAR; =3 LPER/(LPER+LPAR)
      RIN%iPOBS = 0
! Radiance measurement option
  ! =1 I; =2 I/sum(I(:))
      RIN%iIOBS = 1
      RIN%isurf_land(1:2) = -999
      RIN%isurf_water     = -999
      RIN%Aexp_iwl(1:2) = 0
      RIN%ndvi_iwl(1:2) = 0
              
      RIN%SHIFT   = 0.0 
      RIN%NW      = 0
      RIN%WAVE(:) = 0.0
      RIN%IBIN    = 0
      RIN%IMQ     = 0 
      RIN%IPSTOP  = 0
      RIN%INPUT     = .false.

      RIN%MAXP = 999
      RIN%EPSP = 0.0
      RIN%EPSQ = 0.0
      RIN%DL   = 0.0

      RIN%DLSF%IWL        = 0
      RIN%DLSF%key        = 0
      RIN%DLSF%keyEL      = 0       
      RIN%DLSF%distname_O = ' '
      RIN%DLSF%distname_N = ' '            
      RIN%DLSF%internal_file_path = ' '
      RIN%DLSF%external_file_path = ' '
         
      RIN%NDIM%n1               = 0
      RIN%NDIM%n2(:)            = 0
      RIN%NDIM%n3(:,:)          = 0
      RIN%NDIM%ISTARSING(:,:)   = 0
      RIN%NDIM%par_type(:)      = 0
      RIN%NDIM%par_retr(:)      = .false.
         
      RIN%OSHF%IMSC = 0
      RIN%OSHF%NG   = 0
      RIN%OSHF%NN   = 0
      RIN%OSHF%NF   = 0
      RIN%OSHD%IMSC = 0
      RIN%OSHD%NG   = 0
      RIN%OSHD%NN   = 0
      RIN%OSHD%NF   = 0         
      RIN%eps_err   = 0.0001
               
      RIN%SPCA%IO(:,:,:)  = 0
      RIN%SPCA%GSM(:,:,:) = 0.0
      RIN%SPCS%IO(:,:)    = 0
      RIN%SPCS%GSM(:,:)   = 0.0
         
      RIN%MPCS%IOT(:,:)  = 0
      RIN%MPCS%IOX(:,:)  = 0
      RIN%MPCS%IOY(:,:)  = 0
      RIN%MPCS%GSMT(:,:) = 0.0
      RIN%MPCS%GSMX(:,:) = 0.0
      RIN%MPCS%GSMY(:,:) = 0.0
               
      RIN%NOISE%INOISE  = 0
      RIN%NOISE%SGMS(:) = 0.0
      RIN%NOISE%INN(:)  = 0
      RIN%NOISE%DNN(:)  = 0.0		   
      RIN%NOISE%NMT(:)      = 0   
      RIN%NOISE%MT(:,:)     = 0       
      RIN%NOISE%NWLP(:,:)   = 0   
      RIN%NOISE%IWLP(:,:,:) = 0  

      RIN%IPFP%INVSING = 0
      RIN%IPFP%TXY_group = .false.
      RIN%IPFP%TTEST   = 0.0 
      RIN%IPFP%XTEST   = 0.0 
      RIN%IPFP%YTEST   = 0.0

      RIN%APSING(:) = 0.0
      RIN%APSMIN(:) = 0.0
      RIN%APSMAX(:) = 0.0
      RIN%RATIO1(:,:)  = 0.0 
      RIN%RADIUS1(:,:) = 0.0
      RIN%IWW_SINGL(:) = 0
      RIN%NBIN(:)      = 0

      RIN%use_models = .false.

      RIN%edges%nx = 0
      RIN%edges%ny = 0 
      RIN%edges%nt = 0

      RIN%products%retrieval%res = .false.
      RIN%products%retrieval%par = .false.
      RIN%products%retrieval%fit = .false.

      RIN%products%surface%surf = .false.

      RIN%products%aerosol%sd2m_mph = .false.
      RIN%products%aerosol%sd2m_ext = .false.
      RIN%products%aerosol%opt      = .false.
      RIN%products%aerosol%rind     = .false.
      RIN%products%aerosol%phmx     = .false.
      RIN%products%aerosol%lidar    = .false.
      RIN%products%aerosol%chem     = .false.
      RIN%products%aerosol%pm       = .false.
      RIN%products%aerosol%types    = .false.

      RIN%products%clouds%sd2m_mph = .false.
      RIN%products%clouds%sd2m_ext = .false.
      RIN%products%clouds%opt      = .false.
      RIN%products%clouds%rind     = .false.
      RIN%products%clouds%phmx     = .false.
      RIN%products%clouds%lidar    = .false.
      RIN%products%clouds%chem     = .false.
      RIN%products%clouds%pm       = .false.
      RIN%products%clouds%types    = .false.

      RIN%products%errest%par           = .false.
      RIN%products%errest%aerosol%opt   = .false.
      RIN%products%errest%aerosol%lidar = .false.
      RIN%products%errest%clouds%opt    = .false.
      RIN%products%errest%clouds%lidar  = .false.
                  
      RIN%products%forcing%bbflux   = .false.
      RIN%products%forcing%forcing  = .false.
            	   	  
      return
      end subroutine retr_input_initialization

!	sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine set_RIN_fields ( iu_main_output,RIN )

      use mod_retr_settings_derived_type

      implicit none
! -----------------------------------------------------------------------------------------
      type(retr_input_settings),  intent(inout) ::  RIN	  
      integer,                    intent(in)    ::  iu_main_output

! -----------------------------------------------------------------------------------------
! Set number of parameters driving forward model KNSING and 
!     number of retrieved parameters KNSINGF
      call set_RIN_retr_par_number ( RIN,iu_main_output )

! Set difference order for single pixel a priori estimate constraints
      call set_RIN_diff_order_apriori ( RIN,iu_main_output )

! Set radii
      call set_RIN_radii ( RIN,iu_main_output )

! Set surface model flags for RT_OSH 
      call set_RIN_RT_OSH_flags_surf ( RIN,iu_main_output )

! Set flag for edges 
      !call set_RIN_edges_flag ( RIN,iu_main_output )

      return
      end subroutine set_RIN_fields 

!	sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine set_RIN_retr_par_number ( RIN,iu_main_output )

      use mod_retr_settings_derived_type

      implicit none
! -----------------------------------------------------------------------------------------
      type(retr_input_settings),  intent(inout) ::  RIN	  
      integer,                    intent(in)    ::  iu_main_output

      integer   :: IDIM1,IDIM2	  
! -----------------------------------------------------------------------------------------

! Calculate KNSING and KNSINGF      
      RIN%KNSING = 0
      do IDIM1=1,RIN%NDIM%n1 
          do IDIM2=1,RIN%NDIM%n2(IDIM1) 
            RIN%KNSING = RIN%KNSING+RIN%NDIM%n3(IDIM2,IDIM1)
          enddo ! IDIM2
      enddo ! IDIM1
	  
      RIN%KNSINGF=0
      do IDIM1=1,RIN%NDIM%n1 
        if(RIN%NDIM%par_retr(IDIM1)) then
          do IDIM2=1,RIN%NDIM%n2(IDIM1) 
            RIN%KNSINGF = RIN%KNSINGF+RIN%NDIM%n3(IDIM2,IDIM1)
          enddo ! IDIM2
        endif ! RIN%NDIM%pat_retr(IDIM1)
      enddo ! IDIM1      

      if(RIN%IPRI_verbose) then
        write (iu_main_output,'(a)') 'in set_RIN_retr_par_number:'
        write (iu_main_output,'(4x,a,i0,a)') 'KNSING  = ',RIN%KNSING,' - number of parameters driving forward model for each pixel'
        write (iu_main_output,'(4x,a,i0,a)') 'KNSINGF = ',RIN%KNSINGF,' - number of retrieved parameters for each pixel'
      endif
      	   	  
      return
      end subroutine set_RIN_retr_par_number
            
!	sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine set_RIN_diff_order_apriori ( RIN,iu_main_output )

      use mod_retr_settings_derived_type

      implicit none
! -----------------------------------------------------------------------------------------
      type(retr_input_settings),  intent(inout) ::  RIN	  
      integer,                    intent(in)    ::  iu_main_output
! -----------------------------------------------------------------------------------------
! Difference order for single pixel a priori estimate constraints can be only '0'
      RIN%SPCA%IO(:,:,:) = 0

      return
      end subroutine set_RIN_diff_order_apriori
            
!	sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine set_RIN_nsph_package_flags_SD ( RIN,iu_main_output )

      use mod_retr_settings_derived_type

      implicit none
! -----------------------------------------------------------------------------------------
      type(retr_input_settings),  intent(inout) ::  RIN	  
      integer,                    intent(in)    ::  iu_main_output

      integer   :: IDIM1,par_type	  
! -----------------------------------------------------------------------------------------
! The parameters for DLS_bin code
! DLSF (DLS Flags) 
!   IWL=0, key=0 - use original Kernels 22/16 bins or lognormal SD   
!   IWL=1, key=2 - use spectral Kernels for precalculated lognormal bins
!   IWL=0, key=3 - retrieve parameters of lognormal SD
! -----------------------------------------------------------------------------------------

! Spheroid package options
      do IDIM1=1,RIN%NDIM%n1
        if(RIN%NDIM%par_type(IDIM1) .gt. par_type_SD_beg .and. RIN%NDIM%par_type(IDIM1) .lt. par_type_SD_end) then
          par_type = RIN%NDIM%par_type(IDIM1)        
          if(par_type     .eq. par_type_SD_TB) then  ! triangle bins
            RIN%DLSF%IWL = 0
            RIN%DLSF%key = 0        
          elseif(par_type .eq. par_type_SD_LB) then  ! precompued lognormal bins
            RIN%DLSF%IWL = 1
            RIN%DLSF%key = 2                
          elseif(par_type .eq. par_type_SD_LN) then  ! parameters of lognormal SD
            RIN%DLSF%IWL = 0
            RIN%DLSF%key = 3
          endif
        endif ! RIN%NDIM%par_type(IDIM1) .gt. par_type_SD_beg .and.
      enddo ! IDIM1
      if(RIN%IPRI_verbose)  then
        write(iu_main_output,*) 'IWL, key, keyEL :'
        write(iu_main_output,'(3i5,a)') RIN%DLSF%IWL,RIN%DLSF%key,RIN%DLSF%keyEL, &
        '  in set_RIN_nsph_package_flags_SD'
      endif
  
      if(((RIN%DLSF%key .eq. 0  .or.   RIN%DLSF%key .eq. 3) .and. (RIN%DLSF%IWL .ne. 0)) .or.  &
	       ((RIN%DLSF%key .eq. 2) .and. (RIN%DLSF%IWL .ne. 1))) then
        write(*,*) 'key=',RIN%DLSF%key,' IWL=',RIN%DLSF%IWL
        write(*,*) '((key.eq.0/3.and.IWL.ne.0).or.(key.eq.2.and.IWL.ne.1))'
        write(*,*) ' SD options for spheroid package are not correct'
        stop 'stop in set_RIN_nsph_package_flags_SD'
      endif ! key and RIN%DLSF%IWL
      	   	  
      return
      end subroutine set_RIN_nsph_package_flags_SD
      
!	sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine set_RIN_RT_OSH_flags_surf ( RIN,iu_main_output )

      use mod_retr_settings_derived_type

      implicit none
! -----------------------------------------------------------------------------------------
      type(retr_input_settings),  intent(inout) ::  RIN	  
      integer,                    intent(in)    ::  iu_main_output

      integer   :: IDIM1	  
! -----------------------------------------------------------------------------------------
! The parameters for radiative transfer module 
!   iBRDF    
!   iBPDF 
! -----------------------------------------------------------------------------------------
! Surface model options

      RIN%isurf_land(1:2) = -999
      RIN%isurf_water     = -999
           
      do IDIM1=1,RIN%NDIM%n1
        if(RIN%NDIM%par_type(IDIM1) .gt. par_type_SURF1_land_beg .and.     &
           RIN%NDIM%par_type(IDIM1) .lt. par_type_SURF1_land_end) then
                      RIN%isurf_land(1) = RIN%NDIM%par_type(IDIM1)        
        elseif(RIN%NDIM%par_type(IDIM1) .gt. par_type_SURF2_land_beg .and. &
               RIN%NDIM%par_type(IDIM1) .lt. par_type_SURF2_land_end) then
                      RIN%isurf_land(2) = RIN%NDIM%par_type(IDIM1)        
        elseif(RIN%NDIM%par_type(IDIM1) .gt. par_type_SURF_water_beg .and. &
               RIN%NDIM%par_type(IDIM1) .lt. par_type_SURF_water_end) then
                      RIN%isurf_water   = RIN%NDIM%par_type(IDIM1)        
        endif ! RIN%NDIM%par_type(IDIM1) .gt. par_type_SURF1_land_beg .and.
      enddo ! IDIM1

      if(RIN%IPRI_verbose) then
        write(iu_main_output,'(a)') 'in set_RIN_RT_OSH_flags_surf:'
        write(iu_main_output,'(4x,a)') 'isurf_land(1), isurf_land(2), isurf_water:'
        write(iu_main_output,'(1x,3i8,a)') RIN%isurf_land(1),RIN%isurf_land(2),RIN%isurf_water
      endif
      	   	  
      return
      end subroutine set_RIN_RT_OSH_flags_surf

!	sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss      
! Wavelength indices for NDVI

!      subroutine set_RIN_wavelength_index_NDVI ( iu_main_output, RIN )

!      use mod_retr_settings_derived_type

!      implicit none
!	-------------------------------------------------------------------------
!      integer,                      intent(in)     :: iu_main_output
!      type(retr_input_settings),    intent(inout)  :: RIN
!	-------------------------------------------------------------------------      
!      integer           :: iw,idim1,par_type
!      real,parameter    :: wl1_min = 0.650, wl1_max = 0.675,  &
!                           wl2_min = 0.850, wl2_max = 0.875
!	-------------------------------------------------------------------------
!	-------------------------------------------------------------------------
!      RIN%ndvi_iwl(1) = 0
!      RIN%ndvi_iwl(2) = 0
                  
!loop_par_type:  do idim1=1,RIN%NDIM%n1

!      par_type = RIN%NDIM%par_type(idim1)
!      if(par_type .gt. par_type_SURF1_land_beg .and. par_type .lt. par_type_SURF_water_end) then
!        do iw=1,RIN%nw
!          if(RIN%WAVE(iw) .ge. wl1_min .and. RIN%WAVE(iw) .le. wl1_max) then
!            RIN%ndvi_iwl(1) = iw
!          endif
!          if(RIN%WAVE(iw) .ge. wl2_min .and. RIN%WAVE(iw) .le. wl2_max) then
!            RIN%ndvi_iwl(2) = iw
!          endif
!        enddo ! iw
!        if(any(RIN%ndvi_iwl(1:2) .eq. 0)) then
!          RIN%ndvi_iwl(1:2) = 1
!          return
!        elseif(RIN%ndvi_iwl(1) .gt. RIN%ndvi_iwl(2)) then
!          write(iu_main_output,*) 'ndvi_iwl(1)=',RIN%ndvi_iwl(1),  &
!                               ' .gt. ndvi_iwl(2)=',RIN%ndvi_iwl(2)
!          stop 'stop in wavelength_index_NDVI'
!        endif

!      exit
!      endif ! par_type .gt.
!enddo loop_par_type

!      if(any(RIN%ndvi_iwl(1:2) .eq. 0)) then
!          RIN%ndvi_iwl(1:2) = 1
!      endif

!      return
!      end subroutine set_RIN_wavelength_index_NDVI

!	sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine set_RIN_radii ( RIN,iu_main_output )

      use mod_retr_settings_derived_type

      implicit none
! -----------------------------------------------------------------------------------------
      type(retr_input_settings),  intent(inout) ::  RIN	  
      integer,                    intent(in)    ::  iu_main_output

      integer ::  IDIM1,IDIM2,IDIM3,par_type
      real    ::  AD
! -----------------------------------------------------------------------------------------
      do IDIM1=1,RIN%NDIM%n1 
        par_type = RIN%NDIM%par_type(IDIM1)
        if(par_type .eq. par_type_SD_TB .or. par_type .eq. par_type_SD_LN) then
          do IDIM2=1,RIN%NDIM%n2(IDIM1)
            if(RIN%IBIN .eq. -1) then                                  
              AD = (LOG(RIN%RMAX(IDIM2))-LOG(RIN%RMIN(IDIM2)))/(RIN%NDIM%n3(IDIM2,IDIM1)-1)
            else ! RIN%IBIN .EQ. 1                                    
              AD = (RIN%RMAX(IDIM2)-RIN%RMIN(IDIM2))/(RIN%NDIM%n3(IDIM2,1)-1)
            endif
            do IDIM3=1,RIN%NDIM%n3(IDIM2,IDIM1)
              if(RIN%IBIN .eq. -1) then                                     
                RIN%RADIUS1(IDIM3,IDIM2) = EXP(LOG(RIN%RMIN(IDIM2))+AD*(IDIM3-1))
              else ! RIN%IBIN .EQ. 1                                        
                RIN%RADIUS1(IDIM3,IDIM2) = RIN%RMIN(IDIM2)+AD*(IDIM3-1)
              endif
              if(RIN%IPRI_additional_info) then
                write(iu_main_output,*) RIN%RADIUS1(IDIM3,IDIM2),IDIM2,IDIM3,' - RADIUS in set_RIN_radii'
              endif ! RIN%IPRI_additional_info
            enddo ! IDIM3
          enddo ! IDIM2
        exit 
        endif ! par_type .eq. par_type_SD_TB
      enddo ! IDIM1

      return
      end subroutine set_RIN_radii

!	sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss 

      subroutine set_input_settings ( iu_main_output, segment_meas, RIN )
      
      use mod_sdata_derived_type
      use mod_sdata, only : set_RIN_wave_length_array
      use mod_retr_settings_derived_type
      use mod_stop_report

      implicit none
! -----------------------------------------------------------------------------------------
      integer,                    intent(in)    ::  iu_main_output
      type(segment_data),         intent(in)    ::  segment_meas
      type(retr_input_settings),  intent(inout) ::  RIN
! -----------------------------------------------------------------------------------------
! Set RIN fields (radii, surface model flags for RT_OSH, number of retrieved parameters )
      call set_RIN_fields ( iu_main_output, RIN )
      
! Set RIN with wave length array for retrieved parameters (combined array of wavelengths in case  
!     different pixels contain niether diff numbers or diff values of wavelengths)
      call set_RIN_wave_length_array ( iu_main_output, segment_meas, RIN )    

! Search for NDVI wavelength indices
! Wavelength indices for ndvi are set from settings file. 
! Subroutine set_RIN_wavelength_index_NDVI and the call can be deleted
!      call set_RIN_wavelength_index_NDVI ( iu_main_output, RIN )

      return
      end subroutine set_input_settings

!	sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss



