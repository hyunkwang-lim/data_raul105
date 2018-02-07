! **
! **  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.  
! **  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
! **

!> @file inversion.f90
! ********************************************************************************************
! **  The GRASP - Generalized Retrieval of Aerosol and Surface Properties                  ***
! **                                                                                       ***
! **  This program implements the retrieval detailed properties of aerosol and surface     ***
! **  from diverse types of remote sensing (passive and active) and  laboratory            ***
! **  (e.g. scattering phase matrix) observations.                                         ***
! **  The introductory information can be found in the following publications:             ***
! **  - Dubovik, O., T. Lapyonok, P. Litvinov, M. Herman, D. Fuertes, F. Ducos,            *** 
! **    A. Lopatin, A. Chaikovsky, B. Torres, Y. Derimian, X. Huang, M. Aspetsberger,      *** 
! **    and C. Federspiel, “GRASP: a versatile algorithm for characterizing the atmosphere”,**
! **    SPIE: Newsroom, DOI:10.1117/2.1201408.005558, Published Online: September 19, 2014. ** 
! **    http://spie.org/x109993.xml                                                           **
! **                                                                                       ***
! **  - Dubovik, O., M. Herman, A. Holdak, T. Lapyonok, D. Tanré, J. L. Deuzé, F. Ducos,   ***
! **    A. Sinyuk, and A. Lopatin, “Statistically optimized inversion algorithm for        ***
! **    enhanced retrieval of aerosol properties from spectral multi-angle polarimetric    *** 
! **    satellite observations”, Atmos. Meas. Tech., 4, 975-1018, 2011.                    ***
! **                                                                                        **
! *** NOTE: This program inherits many methodological aspects from the retrieval code      *** 
! **        previously developed for AERONET processing. Nonetheless, the program is       *** 
! **        entirely rewritten and redesigned and includes several conceptually new        *** 
! **        features and options:
! **        - the program is highly modular and flexible in a sense that it can be applied *** 
! **          to a variety of measurements and can be set to retrieve different parameters ***
! **        - the inversion can be implemented for each "observed pixel" independently and *** 
! **          also in multi-pixel regime, i.e. when the retrieval is implemented for       *** 
! **          a segment of data (a set of pixels) simultaneously. The multi-pixel retrieval***
! **          allows for improved constraining of retrieved aerosol properties by a priori *** 
! **          limitations of space and time variability of aerosol and surface properties. ***
! **                                                                                       ***
! ********************************************************************************************
! **       This program implements general scheme of non-linear Multi-Term LSM Fitting     ***
! **                    (see Dubovik 2004, Dubovik et al. 2011, etc.)                      ***
! **                                                                                       ***
! **       The general scheme of p and q iterations is realized                            ***
! ********************************************************************************************
#include "../constants_set/mod_globals.inc"
      module mod_grasp

      use mod_sdata
      use mod_retr_settings_derived_type
      use mod_retr_general_output_derived_type
      use mod_index_cloud
      use mod_alloc_arrays
      use mod_edges
      use inversion_subsystem
      use mod_molecular_scattering, only: rayleia
      use mod_stop_report
      use mod_globals

!#if defined(GRASP_SETTINGS) 
!      use mo_grasp_settings
!#endif
      implicit none
    
      contains

! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
        !> main routine for calling the GRASP algorithm 
        !>
        !> @author Tatsiana Lapionak
        !> 
        !> @param[in] inversion_input_settings_file TO_BE_DOCUMENTED
        !> @param[in] sdata_file TO_BE_DOCUMENTED
        !> @param[in] main_output_file TO_BE_DOCUMENTED
        !> @param[in] plotting_output_file TO_BE_DOCUMENTED
        !> @param[in] sdata_sim_file TO_BE_DOCUMENTED
      subroutine grasp_main ( inversion_input_settings_file, sdata_file, main_output_file,   &
                              plotting_output_file,sdata_sim_file                            &
                            )
      use mod_par_inv, only : KPARS,KIMAGE
      use mod_sdata, only : read_input_sdata

      character(*), intent(in)  ::  inversion_input_settings_file
      character(*), intent(in)  ::  sdata_file ! the satellite input data
      character(*), intent(in)  ::  main_output_file
      character(*), intent(in)  ::  plotting_output_file
      character(*), intent(in)  ::  sdata_sim_file ! simulated satellite input data
!	------------------------------------------------------------------------------------------
      type(retr_input_settings)       ::  input_settings
      type(segment_data)              ::  segment_meas
      type(output_segment_general)    ::  output_general
      type(kernels_triangle_bin)      ::  KERNELS1
      type(kernels_lognormal_bin)     ::  KERNELS2
      real,dimension(:,:,:),allocatable ::  US
      real,dimension(KPARS,KIMAGE)    ::  iguess
      type(segment_edges)             ::  edges
!	------------------------------------------------------------------------------------------
!      real :: start1,finish1          
!	------------------------------------------------------------------------------------------
      integer                         ::  irun
      integer                         ::  iu_main_output
!   ------------------------------------------------------------------------------------------

!	------------------------------------------------------------------------------------------
!call CPU_TIME(start1)
! ***  READ Retrieval setting input            *******************************	   

!#if defined(GRASP_SETTINGS)    
!   integer (kind = C_INT)          :: read_status !status of readYamlInput function
!   type(grasp_settings), pointer   :: settings    ! Declare RIN var     
!
!   read_status = grasp_settings_read( settings, inversion_input_settings_file)
!
!   if (read_status .ge. 0) then 
!       if (read_status .eq. 0) then
!            write(*, '(a)') 'Config file read successfully'
!       else
!           write(*, '(a,i0,a)') "Input file read succesfully but it contains ", read_status, " error(s)." 
!       endif
!   else
!       write(*, '(a)') 'ERROR: Config file does not read successfully'
!   endif
!
!   input_settings=grasp_settings_rin(settings);
!   
!   input_settings%main_output_file     = main_output_file
!   input_settings%plotting_output_file = plotting_output_file
!   input_settings%sdata_sim_file       = sdata_sim_file   
!   
!    if (input_settings%main_output_file .eq. '-') then
!         iu_main_output = 6
!    else   
!         open (newunit=iu_main_output, FILE=trim(input_settings%main_output_file), STATUS='replace',recl=20000)         
!    endif        
!#elif defined(CLASSIC_SETTINGS)
    
    call read_input_settings (  inversion_input_settings_file, & ! IN
                                main_output_file,              & ! IN
                                plotting_output_file,          &
                                sdata_sim_file,                &  
                                input_settings,iu_main_output  & ! OUT
                             )
!#else
!#error "Settings module undefined or not exist"
!#endif
            
      if(input_settings%IPRI_additional_info)   &
      write(iu_main_output,*) "after read_input_settings"

! ***  READ Measurements from data file    ******************************* 

      call read_input_sdata ( iu_main_output, &
                              sdata_file,     &
                              segment_meas 	  &
                            )

      if(input_settings%IPRI_additional_info)   &
      write(iu_main_output,*) "after read_input_sdata"

! Edge test
!     Read input file for edges  
        call edges_initialization(edges)
!      if(ledges) then
        call read_edges(iu_main_output,input_settings,edges)
        if(input_settings%IPRI_additional_info)   &
        write(iu_main_output,*) "after read_edges"
!      endif

!     Allocate arrays for spheroid package (DLS_bin) and matrix of Jacobians
      call alloc_arrays(input_settings,KERNELS1,KERNELS2,US)
!     Allocate general output arrays 
      !call alloc_output_general(input_settings,output_general)
!	------------------------------------------------------------------------------------------

      call inversion_init(input_settings)

! ***  CALL multipixel aerosol/surface inversion     ******************************
      iguess(:,:) = -999.0
      do irun=1,1
         call initialize_stop_report()
         call set_segment_for_stop_report(segment_meas)
         call prepare_segment_settings ( iu_main_output, segment_meas, input_settings )
         if ( stop_report%status ) return

!         call print_segment('segment_meas', segment_meas) ! print for debugging 
         call inversion ( input_settings,iu_main_output, & ! IN
                          segment_meas,iguess,edges,     &
                          output_general,	               & ! OUT
                          KERNELS1,KERNELS2,US           & ! INOUT
                        )
         if ( stop_report%status ) return

! Print output results
         if(input_settings%IPRI_verbose) then
           call print_output_results ( iu_main_output, input_settings,  &
                                       segment_meas, output_general     &
                                     )
           if ( stop_report%status ) return
         endif
! Print retrieval results into plotting output file
         if (trim(input_settings%plotting_output_file) .ne. '') &
         call print_classic_plot ( input_settings, segment_meas, output_general )
      enddo ! irun

      call inversion_finalize()
     
!     Allocate arrays for spheroid package (DLS_bin) and matrix of Jacobians
      call dealloc_arrays(input_settings,KERNELS1,KERNELS2,US)
!     Allocate general output arrays 
      !call dealloc_output_general(input_settings,output_general)

!call CPU_TIME(finish1)
!write(*,*) 'in main: CPU_time(sec)=',finish1-start1

      stop 
      end subroutine grasp_main

! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine inversion (  RIN,                  &
                              iu_main_output,       &
                              segment_meas,         &
                              iguess,edges,         &
                              GOUT,		              &
                              KERNELS1,KERNELS2,US  &
                           )
      use mod_index_cloud
      use mod_par_inv
      use mod_par_OS
      use mod_os
      use mod_forward_model
      use mod_inversion_utils
      use mod_par_DLS_bin, only : NRR,NRC
      use mod_par_DLS,     only : KMpar
      use mod_fisher_matrix_ccs
      use mod_retr_settings_derived_type
      use mod_retr_general_output_derived_type
      use mod_molecular_scattering
      use mod_edges
      use mod_covariance_matrix
!XH   for broad band flux
      use mod_bbflux, only : bbflux_pixel
      use mod_sdata_meas_type
      use mod_sdata, only : write_sdata_pixels, print_segment, get_HVP_lidar, set_segment_meas_vector_fs, &
      assign_noise_index, add_rnoise_segment, set_index_clouds, read_input_sdata

      implicit none

!	--------------------------------------------------------------------------------------------                
      type(retr_input_settings),     intent(inout)  ::  RIN
      type(segment_data),            intent(inout)  ::  segment_meas
      real,dimension(KPARS,KIMAGE),  intent(in)     ::  iguess
      type(segment_edges),           intent(in)     ::  edges
      integer,                       intent(in)     ::  iu_main_output 
      type(output_segment_general),  intent(inout)  ::  GOUT
      type(kernels_triangle_bin),    intent(inout)  ::  KERNELS1
      type(kernels_lognormal_bin),   intent(inout)  ::  KERNELS2
      real,dimension(KMESS,KPARS,KIMAGE),intent(inout) ::  US
!	--------------------------------------------------------------------------------------------
      type(ind_clouds)               :: index_clouds
      real,dimension(KKNOISE)        :: resa,resr,resat,resrt
      character(20)                  :: INOISEC
      character(155)                 :: CFMT
      type(nonzero)                  :: UFNZ
      integer                        :: nnz,nnz_err
      real                           :: solver_timer
!	------------------------------------------------------------------------------------------
      integer, DIMENSION(KITIME)     :: IT_m
      integer                        :: npixels,IT
      integer                        :: ipix,ipixstart,ipixstop,KN,KNF,IW
      integer                        :: I,IS1,J,J1,         & 
                                        IAC,IKQ,IIQ,        &
                                        IKS,IKS1,NISM,IKSIM
!	--------------------------------------------------------------------------------------------
      real                           :: ccor_max,ccor,EPSINT,  &  
                                        AQBEFORE,AQAFTER,AAQ,  & 
                                        AGEN,AGENP,AGENQ,      & 
                                        TCINT,EPSPP,EPSQ_test, & 
                                        AFS0, AB_edge

      integer                        :: IP2,LP   
      logical                        :: apriori_estim,apriori_smooth
      integer                        :: nangles
      real,dimension(KMpar)          :: sca_angles	  	  
      integer                        :: KMIMAGE	  	  													           
      integer,dimension(KKNOISE,KIMAGE)       :: IKI ! number of meas for noise (:,:) 
      logical,dimension(KKNOISE,KIMAGE)       :: IKI_shift ! presence of SHIFT
      integer,dimension(KMESS,KKNOISE,KIMAGE) :: KNOISEI
      real,dimension(KPARS)          :: APMIN,APMAX
      type(pixel_vector)             :: pixel_vec_fit_deriv(KIMAGE)
      real                           :: DL1
      REAL, DIMENSION(KPARS,KIMAGE) :: AP,APQ,AP0
      REAL, DIMENSION(KMESS,KIMAGE) :: CS
      REAL, DIMENSION(KPARS,KIMAGE) :: FFMSQ,FFMS0
      REAL, DIMENSION(KPARS,KIMAGE) :: FFMS,QS,AAI, FF_edge

      REAL, DIMENSION(KIMAGE)       :: AMEAS,ASMO,AEST,  &   
                                       ALS,ALSP,ALSQ 
								 
      REAL, DIMENSION(KPARS,KPARS)                :: SMSING,SMSING0     
      REAL, DIMENSION(KIMAGE,KPARS,KIMAGE)        :: SMMULTI

      REAL, DIMENSION(KPARS,KPARS,KIMAGE) :: UFS      

      REAL, DIMENSION(:,:), allocatable :: UF !,UF1

      INTEGER                                :: alloc_stat

      REAL, DIMENSION(KPAR)                  :: FFMI0,QIM 

      INTEGER,DIMENSION(KMPSM)               :: KSM
      INTEGER,DIMENSION(KMPSM,KIMAGE)        :: IMSM

      integer,dimension(:),   allocatable    :: NIPPN  
      integer,dimension(:,:), allocatable    :: IPPN           
      logical		                             :: status_funct 
      !integer,dimension(:,:,:),allocatable   :: MNOISEI   
      integer,dimension(KIP,KWM,KIMAGE)      :: MNOISEI
!	------------------------------------------------------------------------------------------
! Modified normal system
! for the case when some retrieval parameters are assumed the same
      type(nonzero)                  :: UFNZ1       
!	------------------------------------------------------------------------------------------
      type(pixel_vector),dimension(KIMAGE) :: segment_vec_meas
      type(pixel_vector),dimension(KIMAGE) :: segment_vec_fit
!	------------------------------------------------------------------------------------------
! *** deep_random_switch -- this switch enables to set more random noise generations
      logical, parameter				  ::	deep_random_switch = .true.
!	------------------------------------------------------------------------------------------
      integer                     :: INVSING,NW,IWB,IWE
      real, dimension(KW,KIMAGE)  :: tau_mol
      real*8                      :: wl1,h1,lat1 
      integer                     :: NHVP_meas ! number of heights for vertical profile
      real, dimension(KVERTM)     :: HVP_meas ! heights for vertical profile
!	------------------------------------------------------------------------------------------
      integer,save                :: inversion_run_count = 1

      integer, dimension(KIMAGE)         :: KM1_pix
      integer                            :: KM_segm,KM1_segm,KM1_segm_temp
      logical, dimension(KIMAGE)         :: errest_pix
      logical                            :: errest_segm                            

      real,dimension(KPARS,KIMAGE)       :: FFS0,FFS  
!      real,dimension(KMESS,KPARS,KIMAGE) :: US ! allocated before call inversion routine
      real                               :: ALM
      integer                            :: IERR
      real,dimension(KIMAGE)             :: ARES2
!	------------------------------------------------------------------------------------------
! PhF
      real*8, dimension(KPAR)     :: b	
!	------------------------------------------------------------------------------------------ 
      character(len=20)  ::  sparse_solver_name
      logical            ::  test_UFNZ
!	------------------------------------------------------------------------------------------
      logical            ::  IPRI_additional_info
      logical, parameter ::  lcmatrix = .false.!! !.true.: with covarience matrix
      real               ::  temp
! lresult=.true. only single scattering properties are calculated for parameters retrieved 
! at all wavelengths RIN%WAVE(1:RIN%NW)
      logical            ::  lresult = .false., ledges = .false.
      logical, dimension(KPARS,KIMAGE) :: dermask
      logical            ::  edges_present, SvR_meas_present
!	------------------------------------------------------------------------------------------
      logical            ::  molprofskipped = .false. ! AL for output when molecular profile provided but not used
!  ------------------------------------------------------------------------------------------
      integer            :: ipm
!  ------------------------------------------------------------------------------------------

#if defined(SUPERLU_MT)
         sparse_solver_name = 'SUPERLU_MT'
#elif defined(SUPERLU_43)
         sparse_solver_name = 'SUPERLU_43'
#elif defined(VIENNA_CL)
         sparse_solver_name = 'VIENNA_CL'
#elif defined(MUMPS)
         sparse_solver_name = 'MUMPS'
#else
#error No SPARSE_SOLVER configured!
#endif
! ------------------------------------------------------------------------------------------
!#warning "inversion_run_count - TEMP"
      !write(*,'(a,i4,a)') '***** start inversion_run_count =', inversion_run_count,' **********'

!      write(*,*) ".out_x =",edges%group_y(1)%out_x(1,1,1)
!      write(*,*) ".icloud=",edges%group_y(1)%icloud(1,1,1)

      solver_timer = 0.0
      test_UFNZ = .false. ! for comparison of UF,UF1 and UFNZ,UFNZ1	 		                            

      npixels = segment_meas%npixels
      INVSING = RIN%IPFP%INVSING
      KN      = npixels*RIN%KNSING  
      KNF     = npixels*RIN%KNSINGF
      NW      = RIN%NW
      IPRI_additional_info = RIN%IPRI_additional_info

      if(RIN%IPRI_verbose) then
        write(iu_main_output,'(a)') 'in inversion:'
        write(iu_main_output,'(4x,a,i0)') 'npixels = ',npixels
        write(iu_main_output,'(4x,a,i0)') 'nsd = ',RIN%NSD
        if( SvR_meas_present( segment_meas ) ) then
          write(iu_main_output,'(4x,a,f8.5)') 'eps_err =',RIN%eps_err
          write(iu_main_output,'(4x,a,i0,2x,a)') 'aerosol profile   type = ', &
                        RIN%aer_prof_type,'(0 - exponential, 1 - gaussian)'
          write(iu_main_output,'(4x,a,i0,2x,a)') 'molecular profile type = ', &
                        RIN%mol_prof_type,'(0 - exponential, 1 - stdatm)'
        endif
      endif
!	------------------------------------------------------------------------------------------
      if(IPRI_additional_info) then
        call print_segment('segment_meas', segment_meas) ! print for debugging 
        call print_input_settings ( RIN ) ! print classic input settings file for debugging
      endif
! Validate constants
      call validator_constants ( RIN )
      if ( error_present() ) return
! Validate Radiative Transfer (OSH) parameters
      call validator_SOS_RT_parameters ( RIN )
      if ( stop_report%status ) return
! Validate Inversion parameters
      call validator_inversion_settings ( RIN )
      if ( stop_report%status ) return
! Validate presence of surf. and vert.distr. characteristics in settings
      call validator_settings_sdata ( RIN, segment_meas )
      if ( stop_report%status ) return
! Validate sdata (measurements)
      call validator_sdata ( RIN, segment_meas )
      if ( stop_report%status ) return
!	------------------------------------------------------------------------------------------
! Assign General output product flags
      call set_GOUT_product_flags (RIN,GOUT%products)
!	------------------------------------------------------------------------------------------
! Set mask for retrieved surface parameters (land or water only) in order
! to avoid calling forward model for surface derivatives equal to 0.
      call set_surf_derivative_mask ( RIN, segment_meas, dermask )
!	------------------------------------------------------------------------------------------
      if(RIN%IPRI_verbose) write(*,'(4x,a,20f10.5)') 'wavelengths(um): ',RIN%WAVE(1:RIN%NW)
!	------------------------------------------------------------------------------------------
! Molecular scattering
      if( SvR_meas_present( segment_meas ) ) then
        do ipix=1,npixels
          h1   = segment_meas%pixels(ipix)%MASL
          lat1 = segment_meas%pixels(ipix)%y
          do IW=1,NW
          wl1  = RIN%WAVE(IW)
          call rayleia(wl1,h1,lat1,tau_mol(IW,ipix))
          enddo ! IW
        enddo ! ipix
      endif
!	------------------------------------------------------------------------------------------
!AL Checking if any molecular profiles given and putting a warning
      do ipix=1,npixels
         do IW=1,NW
            if (any(segment_meas%pixels(ipix)%meas(IW)%IFMP(:) .eq. 1)) then
               wl1=segment_meas%pixels(ipix)%meas(IW)%wl
               molprofskipped = .true.
               write(iu_main_output,'(A,F5.3,A)') 'WARNING: Molecular profile at ', wl1, 'um ignored!'
            endif 
         enddo !IW
         if (molprofskipped) then
            write(iu_main_output,*) 'Unsupported feature, standard atmosphere will be used instead.'
         endif
      enddo !ipix
!  ------------------------------------------------------------------------------------------
! Define Retrieval OUTput structure with general for segment variables 
      call set_segm_retr_output ( RIN,                & ! IN
                                  segment_meas,       &
                                  GOUT                & ! OUT
                                )
!	------------------------------------------------------------------------------------------
! Check if aerosol profile is a part of retrieved parameters (presence of lidar measurements)
      call get_HVP_lidar ( segment_meas,       & ! IN
                           NHVP_meas,HVP_meas  & ! INOUT
                         )
!      write(*,*) 'AFTER get_HVP_lidar: NHVP_meas=',NHVP_meas
!	------------------------------------------------------------------------------------------
! Initialize and fill in segment measurement vector FS  
      call set_segment_meas_vector_FS ( RIN,INVSING,       & ! IN
                                        segment_meas,      & 
                                        segment_vec_meas   & ! OUT
                                      ) 

! Assign noise indices to measurement types
      call assign_noise_index(RIN,segment_meas,MNOISEI)

! If aplicable, add random noise to segment measurement vector
      if(any(RIN%NOISE%SGMS(1:RIN%NOISE%INOISE) .gt. 0.0)) &
      call add_rnoise_segment ( RIN,deep_random_switch,    & ! IN
                                segment_meas,              &  
                                segment_vec_meas,          & ! INOUT
                                MNOISEI                    & ! IN 
                              )
! CS - calculate covariance matrix
! IKI(I,:)       - total number of measurements of i-th source
! IKI_shift(I,:) - presence of shift for measurements of i-th source
! KNOISE(1,K)  - specific numbers affected by i-th source        
! ARES - parameters adjusting the weight of a priori constraints during iterations
! if RESIDIAL > ARES, the weight of a priori constraints is increased, once
! RESIDIAL < ARES the weight of a priori constraints is fixed
      call covariance_matrix_segment (RIN,INVSING,           & ! IN
                                      segment_meas,          & 
                                      segment_vec_meas,      &
                                      MNOISEI(:,:,:),        & 
                                      IKI(:,:),              &
                                      IKI_shift(:,:),        &
                                      KNOISEI(:,:,:),        & ! OUT
                                      CS(:,:),ARES2(:)       &
                                     ) 	  

!      write(*,'(a)') 'ARES2 :'
!      write(*,'(10e12.4)') ARES2(1:npixels)
!	------------------------------------------------------------------------------------------

      segment_vec_fit = segment_vec_meas

!	------------------------------------------------------------------------------------------
! Definition of a priori and smoothness matrices 

! Initial guess and
! Definition of Vector of a priori astimates 

!      if(inversion_run_count .eq. 2) then
! delete or comment after testing edges
!        RIN%INPUT = .true.
!      endif

      call iguess_apriori (inversion_run_count,  & ! delete
                           iu_main_output,       & ! IN
                           RIN,                  &
                           npixels,iguess,       &
                           AP,AP0,APMIN,APMAX    &	! OUT
	                       )

      !if(inversion_run_count .eq. 2) then
        !write(iu_main_output,*) 'New test'
        !AP(21,1) = log(0.16819E+00)
        !AP(22,1) = log(0.17209E+00)
        !AP(23,1) = log(0.18540E+00)
        !AP(24,1) = log(0.22048E+00)
        !AP(25,1) = log(0.37086E+00)
        !AP(26,1) = log(0.45338E+00)
      !endif

         !if(RIN%IPRI_verbose)  write(iu_main_output,*) 'after iguess_apriori, AP,AP0:'
         !DO I=1,RIN%KNSING
         !if(RIN%IPRI_verbose)  write(iu_main_output,'(i5,1000e14.5)')     &
             !I,(exp(AP(I,ipix)),exp(AP0(I,ipix)),RIN%APSING(I),ipix=1,1)!npixels)
         !ENDDO

      if(IPRI_additional_info) write(iu_main_output,*) 'Before smoothterm_single_pixel'

! Single pixel constrains: a priori smoothness constraints for single-pixel
      call smoothterm_single_pixel_smoothness ( RIN,RIN%SPCS, & ! IN
                                                IKS,SMSING,apriori_smooth  & ! OUT
                                              )
      if ( stop_report%status ) return
! Single pixel constrains: a priori estimates
      call smoothterm_single_pixel_apriori ( RIN,RIN%SPCA, & ! IN
                                             IKS1,SMSING0,apriori_estim  & ! OUT
                                           )
      IP2 = 0
      if(INVSING .gt. 0) then

! Multi pixel constrains: a priori smoothness multi-pixel constraints
         if(IPRI_additional_info)  write(iu_main_output,*) 'Before smoothterm_multi_pixel'
         call set_index_clouds ( RIN, segment_meas, index_clouds )
         call smoothterm_multi_pixel ( iu_main_output,RIN, & ! IN   
                                       index_clouds,       &
                                       NISM,KSM,IMSM,      & ! OUT
                                       IKSIM,SMMULTI       &
                                     )
!write(*,*) 'after smoothterm_multi_pixel: NISM, IKSIM - ',NISM,IKSIM
        if(IPRI_additional_info)  write(iu_main_output,*) 'After smoothterm_multi_pixel'
!*** Multi pixel constrains: a priori estimates of retrieved parameters near
!***                         the edges of the inverteded segment
        FF_edge  (:,:) = 0.0
        AB_edge        = 0.0

        ledges =  &
        edges_present(iu_main_output,segment_meas%NX,segment_meas%NY,segment_meas%NT,RIN,edges)
        if(IPRI_additional_info)  write(iu_main_output,*) 'Before smoothterm_mutli_pixel_edges, ledges=',ledges
        if(ledges)  &
        call smoothterm_mutli_pixel_edges(  iu_main_output,RIN,index_clouds,edges,  &   ! IN
                                            NISM,KSM,IKSIM,IMSM,SMMULTI,            &   ! INOUT
                                            FF_edge,AB_edge ) ! OUT
        if(IPRI_additional_info)  write(iu_main_output,*) 'After smoothterm_mutli_pixel_edges'
!write(*,*) 'after smoothterm_mutli_pixel_edges NISM,IKSIM: ',NISM,IKSIM
!write(*,*) 'in inversion:  FF_edge(:,1): '
!write(*,'(10e14.4)')  FF_edge(:,1)

!write(*,*) 'AB_edge =',AB_edge,'  ledges =',ledges ! delete or comment after edges testing
	
        ! Selection of the group of the same parameters
        ! Search for Time and Space qroups in order to reduce Jacobian size             
        if(RIN%IPFP%TXY_group) then
          allocate(NIPPN(KPAR),stat=alloc_stat)
          if (alloc_stat/=0) stop 'error while trying to allocate NIPPN'
          allocate(IPPN(KPAR,KPAR),stat=alloc_stat)
          if (alloc_stat/=0) stop 'error while trying to allocate IPPN'
          call time_space_group ( iu_main_output,      & ! IN
                                  RIN,                 & 
                                  index_clouds,        &
                                  npixels,             &                      
                                  IP2,NIPPN,IPPN       & ! OUT
                                )
           write(*,*) 'after  time_space_group IP2=',IP2      

           if(IP2 .gt. 0) then
           !Allocate Compressed Column Storage (sparse matrix)
              call allocate_sparse_matrix_storage ( UFNZ1 )
           else
              deallocate(NIPPN,stat=alloc_stat)
              if (alloc_stat /= 0) stop 'error while trying to deallocate NIPPN'
              deallocate(IPPN,stat=alloc_stat)
              if (alloc_stat /= 0) stop 'error while trying to deallocate IPPN'
           endif ! IP2 .gt. 0
        endif ! RIN%IPFP%TXY_group
      endif ! INVSING .GT. 0
	  
      if(RIN%IMQ .eq. 2 .or. test_UFNZ) then
!     Allocate arrays for UF matrix test
         allocate(UF(KPAR,KPAR),stat=alloc_stat)
            if (alloc_stat /= 0) stop 'error while trying to allocate UF'
!         allocate(UF1(KPAR,KPAR),stat=alloc_stat)
!            if (alloc_stat /= 0) stop 'error while trying to allocate UF1'
      endif ! RIN%IMQ .le. 2 .or. 
						  
!     Allocate Compressed Column Storage (sparse matrix)
      call allocate_sparse_matrix_storage ( UFNZ )

!	------------------------------------------------------------------------------------------
! Copy segment_meas structure into segment_fit structure

      GOUT%retrieval%fit%segment_fit = segment_meas ! place in the code is important
!	------------------------------------------------------------------------------------------
! Numbers of virtual degree of freedom for single pixel and multipixel to be used 
! for residual calculation and 
! for scalling of a priori constrains by measurement residual  
      errest_pix(1:npixels) = .true.
      errest_segm           = .true.

      KM_segm = SUM(segment_vec_meas(1:npixels)%KMIMAGE)
      KM1_segm_temp = 0
      do ipix=1,npixels
         KMIMAGE = segment_vec_meas(ipix)%KMIMAGE
         KM1_pix(ipix) = KMIMAGE - RIN%KNSINGF
         if(apriori_smooth) KM1_pix(ipix) = KM1_pix(ipix) + IKS
         if(apriori_estim)  KM1_pix(ipix) = KM1_pix(ipix) + IKS1
         KM1_segm_temp = KM1_segm_temp + KM1_pix(ipix)

         if(KM1_pix(ipix) .lt. RIN%KNSINGF) then
            errest_pix(ipix) = .false.
         endif
         if(KM1_pix(ipix) .le. 0) then
            KM1_pix(ipix)  = KMIMAGE
         endif
      enddo ! ipix 
      IERR = SUM(KM1_pix(1:npixels))

      if(RIN%IPRI_verbose .and. KM1_segm_temp .ne. IERR) &
      write(*,'(2(a,i8))') &
      '!!! WARNING: number of measurements for residual calculation KM1 =', &
      KM1_segm_temp,'  KM =',IERR
      if(INVSING .gt. 0)  then
        KM1_segm = KM1_segm_temp + IKSIM
        if(KM1_segm .lt. KNF) then
          errest_segm = .false.
        endif
        if(KM1_segm .le. 0) then
          KM1_segm    = KM_segm
        endif
      endif
!	------------------------------------------------------------------------------------------

      ipixstart = 1
      select case(INVSING)
      case(0)   ! single pixel
         ipixstop = 1
      case(1,2) ! multi pixel
         ipixstop = npixels  
      case default  
        write(tmp_message,'(a,i0,a)') 'INVSING = ',INVSING,' value is not valid'
        G_ERROR(trim(tmp_message))
      end select

!C******************************************************
!C***         p - iterations
!C******************************************************
!C*** LP -number of iterations
!C*** LTP-number of iterations for choosing TP
!C******************************************************
!C*** DETERMINING of delta AP()     ***
!C******************************************************
!C**  Calculating "measurements" model from AP() vector:
!C******************************************************
  
      ccor_max = 100.0
      ccor     = 0.0

      !ALSP(:)   = 0.0
      !FFMS(:,:) = 0.0
      !US(:,:,:) = 0.0
      QS(:,:)   = 0.0 ! do not remove
      
      IF(RIN%KNSING .GT. RIN%KNSINGF) THEN 
        DO ipix=1,npixels
        APQ(RIN%KNSINGF+1:RIN%KNSING,ipix) = AP(RIN%KNSINGF+1:RIN%KNSING,ipix)
        ENDDO
      ENDIF

177   CONTINUE

!do while(ipixstop .LT. npixels)

      LP = 0

! *********************************************************************************************************************
! *********************************************************************************************************************
! Forward model and Residual for INITIAL GUESS

      ALS(ipixstart:ipixstop) = 0.0
      ASMO(:) = 0.0
      AEST(:) = 0.0 

      do ipix=ipixstart,ipixstop
         segment_vec_fit(ipix)%FS(:) = 0.0
      enddo

      call inversion_forward_model(                              &
                  ipixstart,                                     &
                  ipixstop,                                      &
                  RIN,                                           &
                  RIN%OSHF,                                      &
                  0,                                             &
                  lresult,                                       &
                  tau_mol,                                       &
                  NHVP_meas,                                     &
                  HVP_meas,                                      &
                  AP,                                            &
                  GOUT%retrieval%fit%segment_fit%pixels,         &
                  segment_vec_fit,                               &
                  GOUT%aerosol,                                  &
                  GOUT%clouds,                                   &
                  GOUT%surface,                                  &
                  nangles,                                       &
                  sca_angles,                                    &
                  KERNELS1,                                      &
                  KERNELS2                                       &
            )

      do ipix=ipixstart,ipixstop
         KMIMAGE=segment_vec_meas(ipix)%KMIMAGE

!write(*,*) 'LP=0 AP:'
!write(*,'(10e14.4)') AP(1:RIN%KNSING,ipix)

         !DO I=1,KMIMAGE
         !!WRITE(iu_main_output,*)   &
         !WRITE(0,*)   &
         !ipix,I,segment_vec_fit(ipix)%FS(I),'   - ipix,I,FPS(I,ipix)) AFTER forward_model_pixel in inversion'
         !ENDDO

         if(IPRI_additional_info .and. INVSING .eq. 0) then
         do I=1,RIN%KNSING
            write(iu_main_output,'(2i5,e13.5,a)')   &
            ipix,I,EXP(AP(I,ipix)),' ipix,I,EXP(AP(I)) AFTER forward_model_pixel in inversion'
         enddo
         endif ! IPRI_additional_info .AND.

! residual: measurment term (F-FP)^T (C)^(-1)(F-FP)

         call residual_meas_term (  IPRI_additional_info,iu_main_output,                 & ! IN
                                    KMIMAGE,                              &
                                    segment_vec_meas(ipix)%FS(1:KMIMAGE), &
                                    segment_vec_fit (ipix)%FS(1:KMIMAGE), &
                                    CS(1:KMIMAGE,ipix),                & 
                                    AMEAS(ipix)                        & ! OUT
                                 )
!if (ipix .eq. 3) then
!      call print_segment('segment_meas', segment_meas) ! print for debagging 
!write(*,*) 'FS_meas:'
!write(*,'(10e14.6)') segment_vec_meas(ipix)%FS(1:KMIMAGE)
!write(*,*) 'FS_fit:'
!write(*,'(10e14.6)') segment_vec_fit (ipix)%FS(1:KMIMAGE)
!write(*,*) 'CS:'
!write(*,'(10e14.6)') CS(1:KMIMAGE,ipix)
!write(*,*) 'ipix',ipix,'  AMEAS(ipix)=',AMEAS(ipix)
!stop 'test: after residual_meas_term'
!endif

! residual: a priori smoothness term (single-pixel): (AP)^T(D_single)^T(D_single) AP
         if (apriori_smooth) &
         call residual_apriori_smooth_term ( IPRI_additional_info,iu_main_output,  & ! IN
                                             RIN%KNSINGF,                  &
                                             AP(1:RIN%KNSINGF,ipix),       &
                                             SMSING(1:RIN%KNSINGF,1:RIN%KNSINGF), &
                                             ASMO(ipix)                    & ! OUT
                                           ) 
! residual: (AP-AP0)^T C^(-1) (AP-AP0); A PRIORI ESTIMATE TERM of the residual
         if (apriori_estim) &
         call residual_apriori_estim_term ( IPRI_additional_info,iu_main_output,  & ! IN
                                            RIN%KNSINGF,                  &
                                            AP(1:RIN%KNSINGF,ipix),  &
                                            AP0(1:RIN%KNSINGF,ipix), &
                                            SMSING0(1:RIN%KNSINGF,1:RIN%KNSINGF), &
                                            AEST(ipix)               & ! OUT
                                         ) 
! residual: provide residuals: provides simple mes residuals (F-FP)^T (F-FP) for each noise:
!           the residual inlcludes both redidual for absolute and relative vlaues
         call residual_details ( RIN,                      & ! IN
                                 npixels,ipix,             & 
                                 IKI,IKI_shift,KNOISEI,    & 						 
                                 segment_vec_meas,         &
                                 segment_vec_fit,          &
                                 resa,resr,resat,resrt     & ! OUT 		
                               )
!C***      Prepare residual values required for the output
         call set_pixel_retr_output_residual (  INVSING,LP,ipix,0.0,         & ! IN
                                                RIN%NOISE%INOISE,resa,resr,  &
                                                GOUT%retrieval%res           & ! INOUT
                                             )
     
         if(RIN%IPRI_verbose .and. INVSING .ge. 1)  then
            write(INOISEC,*) RIN%NOISE%INOISE
            CFMT = '(10x,'//trim(adjustl(INOISEC))//'(i5,a,e14.5,f15.5,a),4x,a,i0,6x,a)'
            write(iu_main_output,TRIM(CFMT))                            & 
            (I,':',resa(I),resr(I)*100.0,' %',I=1,RIN%NOISE%INOISE),'  pixel # ', ipix,  &
            'Residual using INITIAL GUESS'
         endif ! RIN%IPRI_verbose

        ccor = AMEAS(ipix)/(ARES2(ipix)*KMIMAGE)
        if(ccor .gt. ccor_max) ccor = ccor_max
        ALS(ipix) = AMEAS(ipix)
        if(apriori_smooth)  ALS(ipix) = ALS(ipix) + ccor*ASMO(ipix)
        if(apriori_estim)   ALS(ipix) = ALS(ipix) + ccor*AEST(ipix) 
        ALS (ipix) = ALS (ipix)/KM1_pix(ipix)            
        ALSP(ipix) = ALS(ipix)
      enddo ! ipix=ipixstart,ipixstop
      !write(*,'(a,10e14.4)') '1: ALS:     ',ALS(1:npixels) 
      !write(*,'(a,10i14)') '1: KM1_pix: ',KM1_pix(1:npixels) 
                         
      select case(INVSING)
      case(0)   ! pixel by pixel scenario
         AGEN = ALS(ipixstop) !/KM1_pix(ipixstop)      
      case(1,2) ! multi pixel scenario
         ccor = SUM(AMEAS(1:npixels))/  &
                SUM(ARES2(1:npixels)*segment_vec_meas(1:npixels)%KMIMAGE)  
         if(ccor .gt. ccor_max) ccor = ccor_max
         AGEN = SUM(AMEAS(1:npixels))
         if(apriori_smooth)  AGEN = AGEN+ccor*SUM(ASMO(1:npixels))
         if(apriori_estim)   AGEN = AGEN+ccor*SUM(AEST(1:npixels))
!C*****************************************************************************************
!C***   inter-pixel smoothness term: (AP)^T (D_inter)^T(D_inter) AP) + edge terms (see subroutine)
!C***                                 See  Dubovik et al. [2011]
!C***      call residual_inter_pixel_term(IPRI_additional_info,iu_iP_main_output,RIN%KNSINGF,npixels,AP,SMMULTI,AFS0)
!C*****************************************************************************************

         call residual_inter_pixel_term(IPRI_additional_info,iu_main_output,RIN%KNSINGF,npixels,AP,SMMULTI, &
                                        ledges,FF_edge,AB_edge,AFS0)
         AGEN = (AGEN+ccor*AFS0)/KM1_segm

         if(IPRI_additional_info) then
            write(iu_main_output,*) SQRT(AGEN)*100.0,'  AGEN ERR %'
            write(iu_main_output,'(2e16.7,i12,a)') AFS0,ccor,KM1_segm,'  AFTER0=AFS0,ccor,KM1_segm'
         endif ! IPRI_additional_info
      end select ! INVSING 
!C*****************************************************************************************
!C***      Prepare residual values required for the output
!C*****************************************************************************************
      call set_total_retr_output_residual (  LP,AGEN,                      & ! IN
                                             RIN%NOISE%INOISE,resat,resrt, &
                                             GOUT%retrieval%res            & ! INOUT
                                          )
      if(RIN%IPRI_verbose) then 
         write(INOISEC,*) RIN%NOISE%INOISE
         if(INVSING .eq. 0) then
           CFMT = '(/,f10.5,'//trim(adjustl(INOISEC))//'(i5,a,e14.5,f15.5,a),4x,a,i0,6x,a)'
           write(iu_main_output,TRIM(CFMT)) SQRT(AGEN)*100.0,      & 
           (I,':',resa(I),resr(I)*100.0,' %',I=1,RIN%NOISE%INOISE),'  pixel # ',ipixstop,  &
           'Residual using INITIAL GUESS'
         else
           CFMT = '(/,f10.5,'//trim(adjustl(INOISEC))//'(i5,a,e14.5,f15.5,a),6x,a,/)'
           write(iu_main_output,TRIM(CFMT)) SQRT(AGEN)*100.0,      &
           (I,':',resat(I),resrt(I)*100.0,' %',I=1,RIN%NOISE%INOISE), &
           'Residual using INITIAL GUESS for TOTAL SEGMENT'
         endif ! INVSING .eq. 0
      endif ! RIN%IPRI_verbose
                  
!      write(*,'(a,10e14.4)') '0: ALS:     ',ALS(1:npixels)
!      write(*,'(a,10i14)')   '0: KM1_pix: ',KM1_pix(1:npixels)
!      write(*,'(a,10e14.4)') '0: AMEAS:   ',AMEAS(1:npixels)                                         
!      write(*,'(a,10e14.4)') '0: ASMO:    ',ASMO(1:npixels)                                         
!      write(*,'(a,10e14.4)') '0: AEST:    ',AEST(1:npixels)                                         
!      write(*,'(a,10e14.4)') '0: ARES2:   ',ARES2(1:npixels)                                         
!      write(*,'(a,10i14)')   '0: KMIMAGE: ',segment_vec_meas(1:npixels)%KMIMAGE                                         
!      write(*,'(2(a,f14.5))')'0: AGEN=',SQRT(AGEN)*100.,'  AFS0=',AFS0                                         

      if(RIN%ISTOP .or. AGEN .lt. 1e-8) then
        select case(INVSING)
        case(0)
          if(ipixstop .lt. npixels) then
            ipixstart = ipixstart + 1
            ipixstop  = ipixstop  + 1
            goto 177
          else
            goto 221
          endif
        case(1,2)
          goto 221
        end select
      endif ! RIN%ISTOP

12    continue

!C***********************************************************
!C*** DERIVATIVES CALCULATION:                            ***
!C***********************************************************

      DO ipix=ipixstart,ipixstop
        US(:,:,ipix)  = 0.0  !- matrix of Jacobians U
        UFS(:,:,ipix) = 0.0  !- Fisher matrix: (U)^T(C)^(-1) U +...
        FFMS(:,ipix)  = 0.0  !- gradient: contribution of measurements
        FFS(:,ipix)   = 0.0  !- contribution of single-pixel smoothness constraints to the gradient
        FFS0(:,ipix)  = 0.0  !- contribution of a priori estimates to the gradient
      ENDDO

      IF(RIN%OSHF%IMSC .NE. RIN%OSHD%IMSC) THEN
        !write(*,*) 'before inversion_forward_model in inversion'
        call inversion_forward_model(                              &
                    ipixstart,                                     &
                    ipixstop,                                      &
                    RIN,                                           &
                    RIN%OSHD,                                      &
                    0,                                             &
                    lresult,                                       &
                    tau_mol,                                       &
                    NHVP_meas,                                     &
                    HVP_meas,                                      &
                    AP,                                            &
                    GOUT%retrieval%fit%segment_fit%pixels,         &
                    pixel_vec_fit_deriv,                           &
                    GOUT%aerosol,                                  &
                    GOUT%clouds,                                   &
                    GOUT%surface,                                  &
                    nangles,                                       &
                    sca_angles,                                    &
                    KERNELS1,                                      &
                    KERNELS2                                       &
              )
!write(*,*) 'after inversion_forward_model in inversion
!do ipix=ipixstart,ipixstop
!   KMIMAGE=segment_vec_meas(ipixstop)%KMIMAGE
!   write(*,*) 'deriv AP:'
!   write(*,'(10e14.4)') AP(1:RIN%KNSING,ipix)
!   write(*,*) 'FS:'
!   write(*,'(10e14.4)') segment_vec_fit(ipix)%FS(1:KMIMAGE)
!   write(*,*) 'FS_deriv:'
!   write(*,'(10e14.4)') pixel_vec_fit_deriv%FS(1:KMIMAGE)
!   stop
!enddo
      ENDIF ! RIN%OSHF%IMSC .NE. RIN%OSHD%IMSC

      DO ipix=ipixstart,ipixstop
      IF (apriori_smooth) THEN
! smoothness term in gradient of minimized residual
        CALL gradient_apriori_smooth_term ( RIN%KNSINGF,                 &  ! IN
                                            AP(1:RIN%KNSINGF,ipix),      &
                                            SMSING(1:RIN%KNSINGF,1:RIN%KNSINGF), &
                                            FFS(1:RIN%KNSINGF,ipix)      & ! OUT
                                          ) 
      ENDIF ! apriori_smooth

      IF (apriori_estim) THEN
! a priori estimates term in gradient of minimized residual
        CALL gradient_apriori_estim_term ( RIN%KNSINGF,                & ! IN
                                           AP(1:RIN%KNSINGF,ipix),     &
                                           AP0(1:RIN%KNSINGF,ipix),    &
                                           SMSING0(1:RIN%KNSINGF,1:RIN%KNSINGF), &
                                           FFS0(1:RIN%KNSINGF,ipix)    & ! OUT
                                         ) 
      ENDIF ! apriori_estim
      ENDDO

      call inversion_jacobian_matrix(                             &
                  iu_main_output,                                 &
                  lp,                                             &
                  ipixstart,                                      &
                  ipixstop,                                       &
                  RIN,                                            &
                  lresult,                                        &
                  tau_mol,                                        &
                  NHVP_meas,                                      &
                  HVP_meas,                                       &
                  AP,                                             &
                  ALS,                                            &
                  GOUT%retrieval%fit%segment_fit%pixels,          &
                  pixel_vec_fit_deriv,                            &
                  GOUT%aerosol,                                   &
                  GOUT%clouds,                                    &
                  GOUT%surface,                                   &
                  nangles,                                        &
                  sca_angles,                                     &
                  KERNELS1,                                       &
                  KERNELS2,                                       &
                  APMIN,                                          &
                  APMAX,                                          &
                  dermask,                                        &
                  US                                              &
            )

      DO ipix=ipixstart,ipixstop
!******************************************************************************
! Measurement term in Fisher matrix for each pixel (Up)^T W^(-1) (Up);  UFS(:,:,ipix)
! Measurement term in gradient of minimized residual for each pixel (Up)^T W^(-1) (fp - f)); FFMS(:,ipix)
!******************************************************************************

!       KMIMAGE = segment_vec_meas(ipixstop)%KMIMAGE
       KMIMAGE = segment_vec_meas(ipix)%KMIMAGE
!write(*,*) 'ipix, KMIMAGE', ipix, KMIMAGE

       CALL FISHMX ( KMIMAGE,RIN%KNSINGF,   &
                     segment_vec_meas(ipix)%FS(:), &
                     segment_vec_fit(ipix)%FS(:),  &
                     CS(:,ipix),            &
                     US(:,:,ipix),          & 
                     UFS(:,:,ipix),         & ! OUT
                     FFMS(:,ipix)           & 
                  )
!stop
!write(*,*) 'ipix: '
!write(*,*) ipix
!write(*,*) 'FS: '
!write(*,*) segment_vec_meas(ipix)%FS(:)
!write(*,*) 'FPS: '
!write(*,*) segment_vec_fit(ipix)%FS(:)
!if (ipix .eq.1 ) then
!write(*,*) 'CS: '
!write(*,*) CS(:,ipix)
!write(*,*) CS(1:13,ipix)
!endif
!write(*,*) 'US: '
!write(*,*) US(:,:,ipix)
!if (ipix .eq. 1) then
!write(*,*) 'UFS: '
!write(*,*) UFS(:,:,ipix)

!write(*,'(10es12.4)') UFS(14,1:32,1)
!write(*,*) 'US: '
!write(*,'(10es12.4)') US(1:13,14,ipix)
!endif !ipix
!write(*,*) 'FFMS1: '  
!write(*,*) FFMS(:,ipix)

!write(*,*) 'ipix=',ipix
!write(*,*)
!write(*,*) 'FS: '
!write(*,'(10e12.4)') segment_vec_meas(ipix)%FS(1:KMIMAGE)
!write(*,*) 'FPS: '
!write(*,'(10e12.4)') segment_vec_fit(ipix)%FS(1:KMIMAGE)

!write(*,*) 'UFS:'
!write(*,'(10e12.4)') UFS(1:RIN%KNSINGF,1:RIN%KNSINGF,ipix)
!write(*,*) 'FFMS:'
!write(*,'(10e12.4)') FFMS(1:RIN%KNSINGF,ipix)

!C*****************************************************
!C*   INCLUSION of SMOOTHNESS  AND A PRIORI  ESTIMATES *
!C*   HERE the weight of a priori term is enhanced at
!C*   earlier iterations and decreasing  to assumed one
!C*   with decrease of total residual to expected
!C*   the methodology described by  Dubovik [2004]    *
!C********************************************************************************

!write(*,*) 'FFS:'
!write(*,'(10e12.4)') FFS(1:RIN%KNSING,ipix)      
!C******************************************************************************
! Adding A priori smoothness term  to Fisher matrix for each pixel (single-pixel):
!       (Up)^T W^(-1) (Up) + (D_singl)^T(D_singl) : - UFS(:,:,ipix)
!!******************************************************************************
! Adding A priori smoothnes contribution in gradient of minimized residual for each pixel (single-pixel):
!       (Up)^T W^(-1) (fp - f)) + ((D_singl)^T(D_singl)Ap-0) : - FFMS(:,ipix)
!C*****************************************************************************
      IF(apriori_smooth) THEN
!*** calculation of contribution of single-pixel smoothness constraints to
!*** the gradient of minimized residual
        CALL gradient_apriori_smooth_term ( RIN%KNSINGF,                 &  ! IN
                                            AP(1:RIN%KNSINGF,ipix),      &
                                            SMSING(1:RIN%KNSINGF,1:RIN%KNSINGF), &
                                            FFS(1:RIN%KNSINGF,ipix)      & ! OUT
                                          ) 
!*** inlcuding smooth conributions to Fisher matrix and gradient
         DO I=1, RIN%KNSINGF
         UFS(I,1:RIN%KNSINGF,ipix) = UFS(I,1:RIN%KNSINGF,ipix) + ccor*SMSING(I,1:RIN%KNSINGF)
         FFMS(I,ipix) = FFMS(I,ipix) + ccor*FFS(I,ipix)
         ENDDO ! I
      ENDIF
!C******************************************************************************
! Adding A priori estimate term  to Fisher matrix for each pixel (single-pixel)
!       (Up)^T W^(-1) (Up) + (D_singl)^T(D_singl)+ (Wa)^(-1) : - UFS(:,:,ipix)
!C******************************************************************************
! Adding A priori estimates contribution in gradient of minimized residual for each pixel (single-pixel):
!       (Up)^T W^(-1) (fp - f)) + ((D_singl)^T(D_singl)Ap - 0) + (Wa)^(-1) (Ap - A0)) : - FFMS(:,ipix)
!C****************************************************************************
      IF(apriori_estim) THEN
!*** calculation of contribution of single-pixel a priori estimates to
!*** the gradient of minimized residual
      CALL gradient_apriori_estim_term (   RIN%KNSINGF,                & ! IN
                                           AP(1:RIN%KNSINGF,ipix),     &
                                           AP0(1:RIN%KNSINGF,ipix),    &
                                           SMSING0(1:RIN%KNSINGF,1:RIN%KNSINGF), &
                                           FFS0(1:RIN%KNSINGF,ipix)    & ! OUT
                                         ) 
!*** inlcuding a priori estimation conributions to Fisher matrix and gradient
         DO I=1, RIN%KNSINGF
         UFS(I,1:RIN%KNSINGF,ipix) = UFS(I,1:RIN%KNSINGF,ipix) + ccor*SMSING0(I,1:RIN%KNSINGF)
         FFMS(I,ipix) = FFMS(I,ipix) + ccor*FFS0(I,ipix)
         ENDDO ! I
      ENDIF          
!C*******************************************************************************
!C* Lev-Mqrdt type constraints inclusion                                        *
!C*******************************************************************************
!C*  This correction is included according                                      *
!C*  the methodology described by Dubovik & King [2000]     
!C*  and Dubovik [2004]
!C*******************************************************************************
      IF(LP .LE. RIN%IPSTOP) THEN
        DO I=1,RIN%KNSINGF
          ALM=(APMAX(I)-APMIN(I))*0.5
          !UFS(I,I,ipix)=UFS(I,I,ipix)+ALS(ipix)/(ALM*ALM*RIN%KNSINGF)
          UFS(I,I,ipix)=UFS(I,I,ipix)+(ALS(ipix)*KM1_pix(ipix))/(ALM*ALM*RIN%KNSINGF)
          !write(*,*) ipix,ALS(ipix),AGEN,'  ipix,ALS(ipix),AGEN'
        ENDDO ! I
        !write(*,'(a,10e14.4)') '2: ALS:     ',ALS(1:npixels)
        !write(*,'(a,10i14)') '2: KM1_pix: ',KM1_pix(1:npixels)

      ENDIF ! LP .LE. IPSTOP         
!C*************************************************************
      ENDDO ! ipix
!write(*,*) 'FFMS2: '
!write(*,'(10e16.6)') FFMS
!write(*,*) 'FFS: '
!write(*,'(10e16.6)') FFS
!write(*,*) 'FFS0: '
!write(*,'(10e16.6)') FFS0

!stop 'stop in iP_MTG.f90 after normal_system'
!write(*,*) 'FFMS3:'
!write(*,'(10e16.6)') FFMS

!C********************************************************************************
!C**           Solving linear system of equation for
!C**                    each p-iterations
!C********************************************************************************
!   The solution is organized as q-iterations (Dubovik and King, 2000; Dubovik 2004)
!   The following linear system is solved: A d(AP) = F(AP)-F
!   the p-th sollution approximation is searched as:    AP+1=AP - d(AP)
!            (d(AP))q+1=(d(AP))q - QS, where QS is iterative solution of system: A QS= A(d(AP))q-[F(AP)-F], i.e.
!            (d(AP))q+1=(d(AP))q - H [A(d(AP))q-[F(AP)-F]];
!   *NOTE* that if H=A^(-1), then we get normal Newton-Gauss:AP+1=AP-A^(-1)[F(AP)-F]
!C********************************************************************************
SELECT CASE(INVSING)
CASE(0) ! Single Pixel inversion
      QS(:,ipixstop)  = 0.0  ! do not delete
      if(RIN%IMQ .eq. 2) then
!C********************************************************************************
!C*** The linear system A QS = F(AP)-F using SVD  solver :
!C********************************************************************************
        do ipix=ipixstart,ipixstop
!write(*,*)
!write(*,*) 'UFS:'
!write(*,'(10e12.4)') UFS(1:RIN%KNSINGF,1:RIN%KNSINGF,ipix)
!write(*,*) 'FFMS:'
!write(*,'(10e12.4)') FFMS(1:RIN%KNSINGF,ipix)
!write(*,*) 
!write(*,*) 'SVD UFS diagonal: LP=',LP
!do i=1,RIN%KNSINGF
!write(*,'(e12.4)') UFS(i,i,ipix)
!enddo
          do i=1,RIN%KNSINGF
            if(UFS(i,i,ipix) .eq. 0.0) then
              !UFS(i,i,ipix) = 1e-30
              if(RIN%IPRI_verbose)  & 
              write(*,*) 'Warning in inversion: UFS diagonal element is ZERO for',ipix,i,'  ipix,i'
            endif 
          enddo

          call ITERQ (  RIN%IMQ,RIN%KNSINGF,                    & ! IN
                        UFS(1:RIN%KNSINGF,1:RIN%KNSINGF,ipix),  &
                        FFMS(1:RIN%KNSINGF,ipix),               &
                        QS(1:RIN%KNSINGF,ipix)                  & ! OUT
                     ) 
        enddo ! ipix
      elseif(RIN%IMQ .eq. 3) then
!C********************************************************************************
!C*** The linear system A QS = F(AP)-F using sparse matrix solver :
!C********************************************************************************
        do ipix=ipixstart,ipixstop
!write(*,*) 'SuperLU UFS diagonal:: LP=',LP
!do i=1,RIN%KNSINGF
!write(*,'(e12.4)') UFS(i,i,ipix)
!enddo
          do i=1,RIN%KNSINGF
            if(UFS(i,i,ipix) .eq. 0.0) then
              !UFS(i,i,ipix) = 1e-30
              if(RIN%IPRI_verbose)  & 
              write(*,*) 'Warning in inversion: UFS diagonal element is ZERO for',ipix,i,'  ipix,i'
            endif 
          enddo
!C*** preparation of matrix A = UF as large vertor with non-zero ellements   :
          call UF_nonzero_one_pixel ( RIN%KNSINGF,UFS(:,:,ipix),  & ! IN
                                      UFNZ,nnz                    & ! OUT
                                    ) 
          if ( stop_report%status ) return
          b(1:RIN%KNSINGF) = FFMS(1:RIN%KNSINGF,ipix)
!write(*,*) 'b:'
!write(*,'(10e16.5)') b(1:KNF)
!write(*,*) 'UFNZ:'
!do i=1,nnz
!write(*,'(2i5,e16.5)') UFNZ%col(i),UFNZ%row(i),UFNZ%val(i)
!enddo 
          !if(RIN%IPRI_verbose)  write(iu_main_output,*)  &
          !'UF inversion, before sparse_solver ',trim(sparse_solver_name),'  nnz=',nnz  

          call sparse_solver ( iu_main_output,       & ! IN 
                               RIN%KNSINGF,nnz,UFNZ, &
                               b,                    & ! INOUT
                               solver_timer          &
                             )
          if ( stop_report%status ) return
          !if(RIN%IPRI_verbose)  write(iu_main_output,*)  &
          !'UF inversion, after  sparse_solver ',trim(sparse_solver_name)
         
          QS(1:RIN%KNSINGF,ipix) = b(1:RIN%KNSINGF)
        enddo ! ipix
      elseif(RIN%IMQ .eq. 1) then
        write(tmp_message,'(2(a,i0),a)') 'Pixel by pixel inversion INVSING = ',INVSING,  &
                  ' and simple linear iterations IMQ = ',RIN%IMQ,'  - not supported'
        G_ERROR(trim(tmp_message))
      endif ! RIN%IMQ .eq. 2
!**********************************************************************************************
CASE(1,2) !***  Multi Pixel inversion ****
!**********************************************************************************************
      FFMS0(:,:) = 0.0 ! ask Oleg
      DO I=1,RIN%KNSINGF
        DO IS1=1,npixels
!**********************************************************************************************
!*** ! Adding inter-pixel smothness contributions ONLY!!! to the diagonal ellement of Fisher matrix
          UFS(I,I,IS1) = UFS(I,I,IS1) + ccor*SMMULTI(IS1,I,IS1)
!**********************************************************************************************
!*** ! Adding inter-pixel smothness contributions to gradient estimation at p-th iteration
!*** ! * NOTE*: if there are a priori estimates values near adges, outsode of inverted segment
!*** !          part of their contribtuon is inlcuded in this term
          FFMS0(I,IS1) = FFMS0(I,IS1) + ccor*DOT_PRODUCT(SMMULTI(IS1,I,1:npixels),AP(I,1:npixels))
!**********************************************************************************************
!*** !   adding second part of contribtuon form a priori estimates near edges to the gradient
!od '-' replaced with'+'
          !IF(ledges) FFMS0(I,IS1) = FFMS0(I,IS1) - ccor*FF_edge(I,IS1)
          IF(ledges) FFMS0(I,IS1) = FFMS0(I,IS1) + ccor*FF_edge(I,IS1)
        ENDDO ! IS1
      ENDDO ! I
!**********************************************************************************************
      IKQ =1
      IF(RIN%IMQ .eq. 3) THEN ! Sparse matrix normal system solver for segment
          !write(*,*) 'before UF_nonzero, IP2=',IP2,  &
          !'  size(UFNZ)(bit)=',size(UFNZ%val)*8+size(UFNZ%row)*4+size(UFNZ%col)*4+size(UFNZ%nval)*4
!**********************************************************************************************
!C*** preparation of matrix UF = [UFS +(D_inter)^T (D_inter)] as large vertor with non-zero ellements:
!C *NOTE* that in "UF_nonzero" non-diagonal ellement of (D_inter)^T (D_inter) are added before
!C        generating final UF_nonzero vector
!**********************************************************************************************
          CALL UF_nonzero ( RIN%KNSINGF,npixels,UFS, & ! IN
                            ccor,SMMULTI,IP2,        &
                            UFNZ,nnz                 & ! OUT
                          )
          if ( stop_report%status ) return
          nnz_err = nnz
          !WRITE(*,*) 'nnz=',nnz
          if(test_UFNZ .and. IP2 .eq. 0) call UFNZ_symmetry_check ( KNF,nnz,UFNZ,UF )
      ELSE IF(RIN%IMQ .eq. 2 .or. test_UFNZ)  THEN ! SVD normal system solver for segment
          ! filling in UF matrix 
          call UF_matrix (  RIN,KNF,npixels,UFS,        & ! IN
                            ccor,SMMULTI,NISM,KSM,IMSM, &
                            UF,nnz                      & ! OUT
                         )
      ELSE IF(RIN%IMQ .eq. 1)  THEN ! simple linear iterations ~ Dubovik et al. [1995] ~
          call matrix_Q_iter (RIN%KNSINGF,npixels,UFS,SMMULTI,NISM,KSM,IMSM,ccor,AAI)
          IKQ = INT(1e+7)/npixels  !IKQ = 10
      ENDIF ! RIN%IMQ .eq. 3		 

      IIQ = 0
      !QS(:,:) = 0.0
! begin solving equation for multi-pixel scenario and its correction
121   CONTINUE 
!**********************************************************************************************
! inclusion of INTER-PIXEL SMOOTHNESS CONSTRAINTS See Dubovik et al.[2011, 2008])
!**********************************************************************************************
      call inter_pix_smooth_constr_incl ( RIN%KNSINGF,npixels,      &
                                          UFS,QS,FFMS,FFMS0,        &
                                          SMMULTI,NISM,KSM,IMSM,ccor, &
                                          FFMSQ                     &
                                        )
      AQAFTER = 0.0
      DO ipix=1,npixels
        CALL residual ( IPRI_additional_info,iu_main_output,     & ! IN
                        RIN%KNSINGF,                 &
                        FFMSQ(1:RIN%KNSINGF,ipix),   &
                        FFMSQ(1:RIN%KNSINGF,ipix),   &
                        AAQ                          & ! OUT
                      ) 
        AQAFTER = AQAFTER+AAQ      
      ENDDO  ! ipix
      IF(IIQ .EQ. 0) AQBEFORE = AQAFTER

      EPSQ_test = (AQAFTER-AQBEFORE)
      if(RIN%IMQ .eq. 1) EPSQ_test = ABS(EPSQ_test)
      EPSQ_test = EPSQ_test/AQBEFORE

      IF(IPRI_additional_info) THEN
        write(iu_main_output,*)     &
        SQRT(AQAFTER/KNF),SQRT(AQBEFORE/KNF),IIQ,IKQ,EPSQ_test,RIN%EPSQ,  &
        '  AQAFTER,AQBEFORE,IIQ,IKQ,EPSQ_test,EPSQ' 
      ENDIF
      
      IF(IIQ .EQ. 0 .OR. (IIQ .LT. IKQ .AND. EPSQ_test .GT. RIN%EPSQ)) THEN		
        AQBEFORE = AQAFTER
        IF(RIN%IMQ .GE. 2) THEN
          FFMI0(:) = 0.0
          DO ipix=1,npixels
            DO I=1,RIN%KNSINGF
              FFMI0((ipix-1)*RIN%KNSINGF+I) = FFMSQ(I,ipix)
            ENDDO ! I
          ENDDO ! ipix
          IF(IP2 .EQ. 0)  THEN       

          if(RIN%IPRI_verbose) &
          write(iu_main_output,'(a,i0,a)') 'KNF = ',KNF, &
          ' - number of retrieved parameters for TOTAL SEGMENT'

            IF(RIN%IMQ .eq. 2) THEN
              if(RIN%IPRI_verbose)  write(iu_main_output,*) 'UF inversion by SVD'  
              CALL ITERQ (  RIN%IMQ,KNF,           &
                            UF(1:KNF,1:KNF),       &
                            FFMI0(1:KNF),          &
                            QIM(1:KNF)             & ! OUT
                         )
            ELSEIF(RIN%IMQ .EQ. 3) THEN
              b(1:KNF) = FFMI0(1:KNF)

              !write(*,*) 'b:'
              !write(*,'(10e16.5)') b(1:KNF)
              !write(*,*) 'UFNZ:'
              !do i=1,nnz
                !write(*,'(2i5,e16.5)') UFNZ%col(i),UFNZ%row(i),UFNZ%val(i)
              !enddo 
              !if(RIN%IPRI_verbose)  write(iu_main_output,*)  &
              !'UF inversion, before sparse_solver ',trim(sparse_solver_name),'  nnz=',nnz  
              call sparse_solver (  iu_main_output, & ! IN 
                                    KNF,nnz,UFNZ,      &
                                    b,                 & ! INOUT
                                    solver_timer       &							   
                                 )
              if ( stop_report%status ) return
              !if(RIN%IPRI_verbose)  write(iu_main_output,*)  &
              !'UF inversion, after  sparse_solver ',trim(sparse_solver_name)         
              QIM(1:KNF) = b(1:KNF)
            ENDIF ! RIN%IMQ .eq. 2
          ELSE  ! IF(IP2 .EQ. 0)   

! (Appendix C in Dubovik et al. 2011) - modifying the normal system 
!  for the case when some retrieval parameters are assumed the same      
            call modified_normal_system ( iu_main_output,RIN,KNF,    & ! IN
                                          npixels,IP2,NIPPN,IPPN,       &      
                                          test_UFNZ,UFNZ,UFNZ1,UF,FFMI0,QIM,  & 
                                          solver_timer,                 &
                                          sparse_solver_name            & ! OUT 
                                        )
          ENDIF ! IP2 .EQ. 0 
        ENDIF ! RIN%IMQ .ge. 2 

! Solution correction
        SELECT CASE(RIN%IMQ)
        CASE(1) 
          DO ipix=1,npixels
            DO I=1,RIN%KNSINGF
              QS(I,ipix) = QS(I,ipix)-FFMSQ(I,ipix)/AAI(I,ipix)
            ENDDO ! I
          ENDDO ! ipix
        CASE(2,3)
          DO ipix=1,npixels
            DO I=1,RIN%KNSINGF
              QS(I,ipix) = QS(I,ipix)-QIM((ipix-1)*RIN%KNSINGF+I)
            ENDDO ! I
          ENDDO ! ipix
        CASE DEFAULT
          write(tmp_message,'(a,i0,2x,a)') 'Q-interetions: IMQ = ',RIN%IMQ,'unknown value'
          G_ERROR(trim(tmp_message))
        END SELECT ! RIN%IMQ

        IIQ = IIQ+1
        GO TO 121
      ENDIF ! IIQ .EQ. 0 .OR. (IIQ .LT. IKQ .AND. EPSQ_test .GT. RIN%EPSQ)
! end solving equation for multi-pixel scenario and its correction
      IF(IPRI_additional_info) WRITE(iu_main_output,*) SQRT(AQAFTER/KNF),SQRT(AQBEFORE/KNF),IIQ,IKQ,  &
                                              '  AQAFTER, AQBEFORE, IIQ, IKQ'
END SELECT ! INVSING  

!C*************************************************************
!C***   END of Solving linear system of equation for 
!C**                    each p-iterations
!C*************************************************************

!C*************************************************************
!C** Determining the optimum length of the correction 
!C**             DELTA AP()
!C*************************************************************     
      LP     = LP+1
      TCINT  = 1.0
      IAC    = 0

! FD hardcoded for debugging
!      IF(RIN%IPRI_verbose) write(6,'(a,i3)') 'LP iteration # ',LP

LOOP_IAC: DO WHILE (IAC .LT. 5)

      APQ(1:RIN%KNSINGF,ipixstart:ipixstop) =                  &
                      AP(1:RIN%KNSINGF,ipixstart:ipixstop) -   &
                      TCINT*QS(1:RIN%KNSINGF,ipixstart:ipixstop) 		 
!      DO ipix=ipixstart,ipixstop 
!        IF(RIN%KL .EQ. 1) THEN
!          WHERE(APQ(1:RIN%KNSINGF,ipix) .LT. APMIN(1:RIN%KNSINGF))  &
!          APQ(1:RIN%KNSINGF,ipix)=APMIN(1:RIN%KNSINGF)
!          WHERE(APQ(1:RIN%KNSINGF,ipix) .GT. APMAX(1:RIN%KNSINGF))  &
!          APQ(1:RIN%KNSINGF,ipix)=APMAX(1:RIN%KNSINGF)
!        ELSE
!          WHERE(APQ(1:RIN%KNSINGF,ipix) .LT. RIN%APSMIN(1:RIN%KNSINGF)) &
!          APQ(1:RIN%KNSINGF,ipix)=RIN%APSMIN(1:RIN%KNSINGF)
!          WHERE(APQ(1:RIN%KNSINGF,ipix) .GT. RIN%APSMAX(1:RIN%KNSINGF)) &
!          APQ(1:RIN%KNSINGF,ipix)=RIN%APSMAX(1:RIN%KNSINGF)
!        ENDIF ! RIN%KL .EQ. 1
!      ENDDO ! ipix

!write(iu_main_output,*) IAC,TCINT,'  IAC, TCINT'
!write(*,*) "QS:"
!write(),'(10e12.4)') QS(1:RIN%KNSINGF,ipixstart:ipixstop)

      ALSQ(:) = 0.0
      ASMO(:) = 0.0
      AEST(:) = 0.0

      DO ipix=ipixstart,ipixstop
         segment_vec_fit(ipix)%FS(:) = 0.0
         call MIN_AP_MAX(1,RIN%KNSINGF,APMIN,APMAX,APQ(:,ipix))
      ENDDO

      call inversion_forward_model(                              &
                  ipixstart,                                     &
                  ipixstop,                                      &
                  RIN,                                           &
                  RIN%OSHF,                                      &
                  0,                                             &
                  lresult,                                       &
                  tau_mol,                                       &
                  NHVP_meas,                                     &
                  HVP_meas,                                      &
                  APQ,                                           &
                  GOUT%retrieval%fit%segment_fit%pixels,         &
                  segment_vec_fit,                               &
                  GOUT%aerosol,                                  &
                  GOUT%clouds,                                   &
                  GOUT%surface,                                  &
                  nangles,                                       &
                  sca_angles,                                    &
                  KERNELS1,                                      &
                  KERNELS2                                       &
            )

      DO ipix=ipixstart,ipixstop
        KMIMAGE=segment_vec_meas(ipix)%KMIMAGE
       !DO I=1,KMIMAGE
        !!WRITE(iu_main_output,*)   &
        !WRITE(*,*)   &
        !ipix,I,segment_vec_fit(ipix)%FS(I),'   - ipix,I,FPS(I,ipix)) AFTER forward_model_pixel'
       !ENDDO

! residual (F-FP)^T C (F-FP)                  
         CALL residual_meas_term (  IPRI_additional_info,iu_main_output,  & ! IN
                                    KMIMAGE,                  &
                                    segment_vec_meas(ipix)%FS(1:KMIMAGE), &
                                    segment_vec_fit (ipix)%FS(1:KMIMAGE), &
                                    CS(1:KMIMAGE,ipix),       & 
                                    AMEAS(ipix)               & ! OUT
                                 )

!write(*,*) 'AMEAS =',AMEAS(ipix)
!write(*,*) 'FS_meas:'
!write(*,'(10e14.6)') segment_vec_meas(ipix)%FS(1:KMIMAGE)
!write(*,*) 'FS_fit:'
!write(*,'(10e14.6)') segment_vec_fit (ipix)%FS(1:KMIMAGE)
!write(*,*) 'CS:'
!write(*,'(10e14.6)') CS(1:KMIMAGE,ipix)
!stop 'test in inversion'

! residual AP^T SMSING AP
         IF (apriori_smooth) &
         CALL residual_apriori_smooth_term ( IPRI_additional_info,iu_main_output,  & ! IN
                                             RIN%KNSINGF,              &
                                             APQ(1:RIN%KNSINGF,ipix),  &
                                             SMSING(1:RIN%KNSINGF,1:RIN%KNSINGF), &
                                             ASMO(ipix)                & ! OUT
                                           ) 

! residual (AP-AP0)^T SMSING0 (AP-AP0); A PRIORI ESTIMATE TERM of the residual
         IF (apriori_estim) &
         CALL residual_apriori_estim_term ( IPRI_additional_info,iu_main_output,  & ! IN
                                            RIN%KNSINGF,              &
                                            APQ(1:RIN%KNSINGF,ipix),  &
                                            AP0(1:RIN%KNSINGF,ipix),  &
                                            SMSING0(1:RIN%KNSINGF,1:RIN%KNSINGF), &
                                            AEST(ipix)                & ! OUT
                                          ) 
           
!write(*,*) 'ipix=',ipix,'  ASMO(ipix)=',ASMO(ipix),'  AEST(ipix)=',AEST(ipix)
!******* END of A PRIORI TERMS *******************************

!tl      write(*,'(a,5i5)') 'ipix,KMIMAGE,KNSING,IKS,IKS1: ',ipix,KMIMAGE,RIN%KNSING,IKS,IKS1
!tl      write(*,'(a,3e13.5)') 'AMEAS,ASMO,AESTQ: ',AMEAS(ipix),ASMO(ipix),AEST(ipix)
	  	  	                   
        ccor = AMEAS(ipix)/(ARES2(ipix)*KMIMAGE)
        if(ccor .gt. ccor_max) ccor = ccor_max
        ALSQ(ipix) = AMEAS(ipix)
        if(apriori_smooth)  ALSQ(ipix) = ALSQ(ipix) + ccor*ASMO(ipix)
        if(apriori_estim)   ALSQ(ipix) = ALSQ(ipix) + ccor*AEST(ipix) 
        ALSQ(ipix) = ALSQ(ipix)/KM1_pix(ipix)
        ENDDO ! ipix=ipixstart,ipixstop
!      write(*,'(a,10e14.4)') '3: ALSQ:    ',ALSQ(1:npixels)
!      write(*,'(a,10i14)')   '3: KM1_pix: ',KM1_pix(1:npixels)
!      write(*,'(a,10e14.4)') '3: AMEAS:   ',AMEAS(1:npixels)                                         
!      write(*,'(a,10e14.4)') '3: ASMO:    ',ASMO(1:npixels)                                         
!      write(*,'(a,10e14.4)') '3: AEST :   ',AEST(1:npixels)                                         
!      write(*,'(a,10e14.4)') '3: ARES2:   ',ARES2(1:npixels)                                         
!      write(*,'(a,10i14)')   '3: KMIMAGE: ',segment_vec_meas(1:npixels)%KMIMAGE                                         
                  
      select case(INVSING)
      case(0)    ! pixel by pixel scenario
         AGENQ = ALSQ(ipixstop) !/KM1_pix(ipixstop)      
      case(1,2)  ! multi pixel scenario
         ccor = SUM(AMEAS(1:npixels))/  &
                SUM(ARES2(1:npixels)*segment_vec_meas(1:npixels)%KMIMAGE)  
         if(ccor .gt. ccor_max) ccor = ccor_max
         AGENQ = SUM(AMEAS(1:npixels))
         if(apriori_smooth)  AGENQ = AGENQ+ccor*SUM(ASMO(1:npixels))
         if(apriori_estim)   AGENQ = AGENQ+ccor*SUM(AEST(1:npixels)) 
! inter-pixel smoothnes constraint term of residual, See  Dubovik et al. [2011]
         call residual_inter_pixel_term(IPRI_additional_info,iu_main_output,RIN%KNSINGF,npixels,APQ,SMMULTI,  &
                                                            ledges,FF_edge,AB_edge,AFS0)
         AGENQ = (AGENQ+ccor*AFS0)/KM1_segm
         if(IPRI_additional_info) then
            write(iu_main_output,*) SQRT(AGENQ)*100.0,'  AGENQ (ERR %)'
            write(iu_main_output,'(2e16.7,i12,a)') AFS0,ccor,KM1_segm,'  AFTER0=AFS0,ccor,KM1_segm'
         endif ! IPRI_additional_info
      end select ! INVSING 

      EPSINT = (AGEN-AGENQ)/AGEN
      if(IPRI_additional_info) then
        write(iu_main_output,*) EPSINT,'  EPSINT'
        write(iu_main_output,*) SQRT(AGENQ)*100.,SQRT(AGEN)*100.,'  AGENQ,AGEN (ERR%)'
      endif ! IF(IPRI_additional_info)

      IF(IPRI_additional_info)  &
      WRITE(iu_main_output,*) IAC,TCINT,EPSINT,'  IAC, TCINT,EPSINT'
      
      IF(EPSINT .LE. 0.) THEN
        IAC = IAC+1
        TCINT = TCINT*0.5
      ELSE
        EXIT LOOP_IAC
      ENDIF ! EPSINT .LE. 0.

ENDDO LOOP_IAC

      IF(EPSINT .LE. 0.) THEN
        IAC   = 0
          call inversion_forward_model(                              &
                      ipixstart,                                     &
                      ipixstop,                                      &
                      RIN,                                           &
                      RIN%OSHF,                                      &
                      0,                                             &
                      lresult,                                       &
                      tau_mol,                                       &
                      NHVP_meas,                                     &
                      HVP_meas,                                      &
                      AP,                                            &
                      GOUT%retrieval%fit%segment_fit%pixels,         &
                      segment_vec_fit,                               &
                      GOUT%aerosol,                                  &
                      GOUT%clouds,                                   &
                      GOUT%surface,                                  &
                      nangles,                                       &
                      sca_angles,                                    &
                      KERNELS1,                                      &
                      KERNELS2                                       &
                  )

          ALSP(ipixstart:ipixstop) = ALS(ipixstart:ipixstop)
        AGENP = AGEN    
      ELSE            
        do ipix=ipixstart,ipixstop
! Vector of retrieved parameters after LP-th ineration
          AP(1:RIN%KNSINGF,ipix) = APQ(1:RIN%KNSINGF,ipix)
          ALSP(ipix) = ALSQ(ipix)
        enddo  ! ipix
        AGENP = AGENQ
      ENDIF ! EPSINT .LE. 0.

!	------------------------------------------------------------------------------------------------------
!	------------------------------------------------------------------------------------------------------

      SELECT CASE(INVSING)
      CASE(0)
         DO ipix=ipixstart,ipixstop
            EPSPP=(SQRT(ALS(ipix))-SQRT(ALSP(ipix)))/SQRT(ALS(ipix)) 
! TL 02-12-2014
            if(EPSPP .eq. 1. .and. ALSP(ipix) .eq. 0.0 ) EPSPP = RIN%EPSP
            IF(IPRI_additional_info)  THEN
               if(RIN%IPRI_verbose) WRITE(iu_main_output,*)  &
               EPSPP,SQRT(ALSP(ipix))*100.0, SQRT(ALS(ipix))*100.0,ipix, &
               '  EPSPP,ALSP(ipix),ALS(ipix),ipix'
            ENDIF ! IPRI_additional_info
         ENDDO ! ipix
      CASE(1,2)
          EPSPP=(SQRT(AGEN)-SQRT(AGENP))/SQRT(AGEN)
! TL 02-12-2014
          if(EPSPP .eq. 1. .and. AGENP .eq. 0.0 ) EPSPP = RIN%EPSP 
!         WRITE(iu_main_output,*) EPSPP,'EPSPP'
         IF(IPRI_additional_info)  THEN
            if(RIN%IPRI_verbose)  write(iu_main_output,*) &
            SQRT(AGEN)*100.0,SQRT(AGENP)*100.0,'  AGEN,AGENP ERR%'
         ENDIF
      END SELECT

      DO ipix=ipixstart,ipixstop
        call residual_details ( RIN,                    & ! IN
                                npixels,ipix,           & 
                                IKI,IKI_shift,KNOISEI,  & 					 	 
                                segment_vec_meas,       & 
                                segment_vec_fit,        &
                                resa,resr,resat,resrt   & ! OUT 		
                              )

          call set_pixel_retr_output_residual ( INVSING,LP,ipix,ALSP(ipix),  & ! IN
                                                RIN%NOISE%INOISE,resa,resr,  &
                                                GOUT%retrieval%res           & ! INOUT
                                              )

          if(RIN%IPRI_verbose .and. INVSING .ge. 1)  then
              write(INOISEC,*) RIN%NOISE%INOISE  
              CFMT = '(10x,'//trim(adjustl(INOISEC))//'(i5,a,e14.5,f15.5,a),2(6x,a,i0))'
              write(iu_main_output,TRIM(CFMT))     &
              (I,':',resa(I),resr(I)*100.0,' %',I=1,RIN%NOISE%INOISE),'pixel # ', ipix,  &
              'Residual after iteration # ', LP
          endif ! RIN%IPRI_verbose
      ENDDO ! ipix=ipixstart,ipixstop
 
      call set_total_retr_output_residual (  LP,AGENP,                     & ! IN
                                             RIN%NOISE%INOISE,resat,resrt, &
                                             GOUT%retrieval%res            & ! INOUT
                                          )
      if(RIN%IPRI_verbose) then 
        write(INOISEC,*) RIN%NOISE%INOISE
        if(INVSING .eq. 0) then
          CFMT = '(f10.5,'//trim(adjustl(INOISEC))//'(i5,a,e14.5,f15.5,a),2(6x,a,i0))'
          write(iu_main_output,TRIM(CFMT)) SQRT(AGENP)*100.0, &
          (I,':',resa(I),resr(I)*100.0,' %',I=1,RIN%NOISE%INOISE),'pixel # ', ipixstop, &
          'Residual after iteration # ', LP
        else
          CFMT = '(/,f10.5,'//trim(adjustl(INOISEC))//'(i5,a,e14.5,f15.5,a),6x,a,i0,a/)'
          write(iu_main_output,TRIM(CFMT)) SQRT(AGENP)*100.0, &
          (I,':',resat(I),resrt(I)*100.0,' %',I=1,RIN%NOISE%INOISE), &
          'Residual after iteration # ', LP,'  for TOTAL SEGMENT'
        endif ! INVSING .eq. 0
      endif ! RIN%IPRI_verbose	   
               
!      write(*,'(a,10e14.4)') 'last: ALSQ:    ',ALSQ(1:npixels)
!      write(*,'(a,10i14)')   'last: KM1_pix: ',KM1_pix(1:npixels)
!      write(*,'(a,10e14.4)') 'last: AMEAS:   ',AMEAS(1:npixels)                                         
!      write(*,'(a,10e14.4)') 'ladt: ASMO:    ',ASMO(1:npixels)                                         
!      write(*,'(a,10e14.4)') 'last: AEST:    ',AEST(1:npixels)                                         
!      write(*,'(a,10e14.4)') 'last: ARES2:   ',ARES2(1:npixels)                                         
!      write(*,'(a,10i14)')   'last: KMIMAGE: ',segment_vec_meas(1:npixels)%KMIMAGE                                         
!      write(*,'(3(a,f14.5))')'last: AGEN=',SQRT(AGEN)*100.,'  AGENP=',SQRT(AGENP)*100., &
!                                  '  AFS0=',AFS0
!      do ipix=1,npixels
!      write(*,*) &
!      ipix,ALSP(ipix),AMEAS(ipix),ASMO(ipix),AEST(ipix),ARES2(ipix),segment_vec_meas(ipix)%KMIMAGE,KM1_pix(ipix), &
!      '  - 3:  ipix,ALSP,AMEAS,ASMO,AEST,ARES2,KMIMAGE,KM1_pix'
!      enddo
!      write(*,*) 'AGENP=',AGENP

! Update residuals
      DO ipix=ipixstart,ipixstop
        ALS(ipix) = ALSP(ipix)
      ENDDO ! ipix      
      AGEN = AGENP

 221  CONTINUE

      IF(ipixstop .EQ. npixels .OR. RIN%ISTOP) THEN

      IF(INVSING .GT. 0 .OR. (EPSPP .LE. RIN%EPSP .OR. LP .EQ. RIN%MAXP) .OR. RIN%ISTOP) THEN
 !       WRITE(iu_main_output,*)
 !       DO I=1,RIN%KNSING
 !       WRITE(iu_main_output,'(i5,1000e14.5)')  &
 !        I,EXP(AP(I,1:npixels))
 !       ENDDO ! I

! **  Calculate optical characteristics (ext,ssa,lr) at RIN%WAVE wavelengths for retrieved parameters
         lresult = .true.
         call inversion_forward_model(                              &
                     1,                                             &
                     npixels,                                       &
                     RIN,                                           &
                     RIN%OSHF,                                      &
                     0,                                             &
                     lresult,                                       &
                     tau_mol,                                       &
                     NHVP_meas,                                     &
                     HVP_meas,                                      &
                     AP,                                            &
                     GOUT%retrieval%fit%segment_fit%pixels,         &
                     segment_vec_fit,                               &
                     GOUT%aerosol,                                  &
                     GOUT%clouds,                                   &
                     GOUT%surface,                                  &
                     nangles,                                       &
                     sca_angles,                                    &
                     KERNELS1,                                      &
                     KERNELS2                                       &
                  )
         lresult = .false.

         do ipix=1,npixels
            call set_gout_pixel(  iu_main_output, RIN, segment_meas, & ! IN
                                  ipix, nangles, sca_angles, AP,     &
                                  GOUT,                      & ! INOUT
                                  KERNELS1,KERNELS2          &
                               )
         enddo ! ipix

! **  Write retrieval results after each iteration
        if(RIN%IPRI_iter_result)  then
          call print_main_output ( iu_main_output, RIN, segment_meas, GOUT )
        endif

      ENDIF ! INVSING .GT. 0 .OR. EPSPP .LE. EPSP 
      ENDIF ! ipixstop .EQ. npixels .OR.

!C*****************************************************
!C*** CHECKING if more iterations are necessary***
!C*** INVSING=0 - pixel-by-pixel inversion;
!C*** INVSING>0 - multi-pixel inversion
!C*** The condition below is set that 
!C    multi-pixel inversion is followed by
!C    pixel-by-pixel inversion using solution of 
!C    multi-pixel inversion and then it stops
!C*****************************************************

      IF(.NOT. RIN%ISTOP) THEN
        IF(EPSPP .GT. RIN%EPSP .AND. LP .LT. RIN%MAXP) THEN
          if(INVSING .gt. 0) then ! multipixel scenario  
            if(lcmatrix) then
! recalculate covariance matrix CS
!do ipix=1,npixels
!write(*,*) 'ipix=',ipix
!write(*,*) 'CM1: MNOISEI  :'
!write(*,'(10i8)') MNOISEI(1:3,1:NW,ipix) 
!write(*,*) 'CM1: IKI      :'
!write(*,'(10i8)') IKI(1:2,ipix)
!write(*,*) 'CM1: IKI_shift:'
!write(*,'(10l8)') IKI_shift(1:2,ipix)
!write(*,*) 'CM1: KNOISEI  :'
!write(*,'(10i8)') KNOISEI(1:3,1:NW,ipix)
!enddo 
              call covariance_matrix_segment (RIN,INVSING,           & ! IN
                                              segment_meas,          & 
                                              segment_vec_meas,      &
                                              MNOISEI(:,:,:),        & 
                                              IKI(:,:),              &
                                              IKI_shift(:,:),        &
                                              KNOISEI(:,:,:),        & ! OUT
                                              CS(:,:),ARES2(:),      &
                                              GOUT%retrieval%res     &
                                             )
!write(*,*) '**********************************'
!do ipix=1,npixels
!write(*,*) 'ipix=',ipix
!write(*,*) 'CM2: MNOISEI  :'
!write(*,'(10i8)') MNOISEI(1:3,1:NW,ipix) 
!write(*,*) 'CM2: IKI      :'
!write(*,'(10i8)') IKI(1:2,ipix)
!write(*,*) 'CM2: IKI_shift:'
!write(*,'(10l8)') IKI_shift(1:2,ipix)
!write(*,*) 'CM2: KNOISEI  :'
!write(*,'(10i8)') KNOISEI(1:3,1:NW,ipix)
!enddo

! residual (F-FP)^T C (F-FP)                  
              do ipix=1,npixels
                KMIMAGE=segment_vec_fit(ipix)%KMIMAGE
                call residual_meas_term ( IPRI_additional_info,iu_main_output,  & ! IN
                                          KMIMAGE,                  &
                                          segment_vec_meas(ipix)%FS(1:KMIMAGE),  &
                                          segment_vec_fit(ipix)%FS(1:KMIMAGE),   &
                                          CS(1:KMIMAGE,ipix),       & 
                                          AMEAS(ipix)               & ! OUT
                                        )
! residual AP^T SMSING AP
                if (apriori_smooth) &
                call residual_apriori_smooth_term ( IPRI_additional_info,iu_main_output, & ! IN 
                                                    RIN%KNSINGF,             &
                                                    AP(1:RIN%KNSINGF,ipix),  &
                                                    SMSING(1:RIN%KNSINGF,1:RIN%KNSINGF), &
                                                    ASMO(ipix)               & ! OUT
                                                  ) 

! residual (AP-AP0)^T SMSING0 (AP-AP0); A PRIORI ESTIMATE TERM of the residual
                if (apriori_estim) &
                call residual_apriori_estim_term (  IPRI_additional_info,iu_main_output,  & ! IN
                                                    RIN%KNSINGF,              &
                                                    AP(1:RIN%KNSINGF,ipix),   &
                                                    AP0(1:RIN%KNSINGF,ipix),  &
                                                    SMSING0(1:RIN%KNSINGF,1:RIN%KNSINGF), &
                                                    AEST(ipix)                & ! OUT
                                                 ) 

                ccor = AMEAS(ipix)/(ARES2(ipix)*KMIMAGE)
                if(ccor .gt. ccor_max) ccor = ccor_max
                ALSP(ipix) = AMEAS(ipix)
                if(apriori_smooth)  ALSP(ipix) = ALSP(ipix) + ccor*ASMO(ipix)
                if(apriori_estim)   ALSP(ipix) = ALSP(ipix) + ccor*AEST(ipix) 
                ALSP(ipix)  = ALSP(ipix)/KM1_pix(ipix)
                enddo ! ipix
!      write(*,'(a,10e14.4)') '4: ALSP:    ',ALSP(1:npixels)
!      write(*,'(a,10i14)')   '4: KM1_pix: ',KM1_pix(1:npixels)
!      write(*,'(a,10e14.4)') '4: AMEAS:   ',AMEAS(1:npixels)                                         
!      write(*,'(a,10e14.4)') '4: ASMO:    ',ASMO(1:npixels)                                         
!      write(*,'(a,10e14.4)') '4: AEST :   ',AEST(1:npixels)                                         
!      write(*,'(a,10e14.4)') '4: ARES2:   ',ARES2(1:npixels)                                         
!      write(*,'(a,10i14)')   '4: KMIMAGE: ',segment_vec_meas(1:npixels)%KMIMAGE                                         
              ccor = SUM(AMEAS(1:npixels))/  &
                     SUM(ARES2(1:npixels)*segment_vec_meas(1:npixels)%KMIMAGE)  
              if(ccor .gt. ccor_max) ccor = ccor_max
              AGENP = SUM(AMEAS(ipixstart:ipixstop))
              if(apriori_smooth)  AGENP = AGENP+ccor*SUM(ASMO(1:npixels))
              if(apriori_estim)   AGENP = AGENP+ccor*SUM(AEST(1:npixels)) 
! inter-pixel smoothnes constraint term of residual, See  Dubovik et al. [2011]
              call residual_inter_pixel_term(IPRI_additional_info,iu_main_output,RIN%KNSINGF,npixels,AP,SMMULTI,  & 
                                                                    ledges,FF_edge,AB_edge,AFS0)
              AGENP = (AGENP+ccor*AFS0)/KM1_segm
              call set_total_retr_output_residual ( LP,AGENP,                      & ! IN
                                                    RIN%NOISE%INOISE,resat,resrt,  &
                                                    GOUT%retrieval%res             & ! INOUT
                                                  )      
!      do ipix=1,npixels
!      write(*,*) &
!      ipix,ALSP(ipix),AMEAS(ipix),ASMO(ipix),AEST(ipix),ARES2(ipix),segment_vec_meas(ipix)%KMIMAGE,KM1_pix(ipix), &
!      '  - 4:  ipix,ALSP,AMEAS,ASMO,AEST,ARES2,KMIMAGE,KM1_pix'
!      enddo
!      write(*,*) 'AGENP=',AGENP

! Redefine residiuals for next iteration
              ALS(1:npixels) = ALSP(1:npixels)
              AGEN = AGENP 
            endif ! lcmatrix .and.
          endif ! INVSING .gt. 0                
          GOTO 12
        ELSE
          IF(INVSING .NE. 2) THEN
            IF(INVSING .EQ. 0 .AND. ipixstop .LT. npixels) THEN
              ipixstart = ipixstart+1
              ipixstop  = ipixstop +1
              GOTO 177         
            ELSEIF(INVSING .EQ. 1) THEN
              INVSING   = 0
              ipixstart = 1
              ipixstop  = 1
              if(any(RIN%NOISE%SGMS(1:RIN%NOISE%INOISE) .gt. 0.0)) then
                write(tmp_message,'(a)') 'Can not add twice random noise to meas vector.'
                G_ERROR(trim(tmp_message))
              endif
              ! Assign noise indices to measurement types
              call assign_noise_index(RIN,segment_meas,MNOISEI)
              call covariance_matrix_segment (RIN,INVSING,          & ! IN
                                              segment_meas,         & 
                                              segment_vec_meas,     & 
                                              MNOISEI,              & 
                                              IKI(:,:),             & ! OUT
                                              IKI_shift(:,:),       &
                                              KNOISEI(:,:,:),       &
                                              CS(:,:),ARES2         &
                                             )
              GOTO 177         
            ENDIF ! INVSING .EQ. 0 .AND.
          ENDIF ! INVSING .NE. 2 
        ENDIF ! EPSPP .GT. RIN%EPSP 
      ENDIF ! .NOT. RIN%ISTOP
      
!     Deallocate Compressed Column Storage if allocated (sparse matrix with time_space groups)
      call deallocate_sparse_matrix_storage ( UFNZ1 )

!     Deallocate time space group arrays
      if(allocated(NIPPN)) then
        deallocate(NIPPN,stat=alloc_stat)
        if (alloc_stat /= 0) stop 'error while trying to deallocate NIPPN'
      endif
      if(allocated(IPPN))  then
        deallocate(IPPN,stat=alloc_stat)
        if (alloc_stat /= 0) stop 'error while trying to deallocate IPPN'
      endif

      GOUT%aerosol%phmx%nangle   = nangles
      GOUT%aerosol%phmx%angle(:) = sca_angles(:)

!XH   testing broadband flux calculation
      if (RIN%products%forcing%bbflux .or. RIN%products%forcing%forcing) then
         do ipix=1,npixels
!XH         define the list of levels
!            GOUT%forcing%bbflux%pixel(ipix)%NHLV =KNT
!            do i=1, KNT
!               GOUT%forcing%bbflux%pixel(ipix)%HLV(i) =(KNT-i)*100.0/(KNT-1)
!            end do
!XH         level list consitent with Yevgeny's
            GOUT%forcing%bbflux%pixel(ipix)%NHLV = 10
            GOUT%forcing%bbflux%pixel(ipix)%HLV(1:GOUT%forcing%bbflux%pixel(ipix)%NHLV)=  &
                              (/120.00,15.00,10.00,6.00,5.00,4.00,3.00,2.00,1.00,0.06/)
!XH         longitude and latitude of the pixel
            h1   = segment_meas%pixels(ipix)%MASL
            lat1 = segment_meas%pixels(ipix)%y
            call bbflux_pixel (                                                &
                                 RIN,ipix,h1,lat1,                             & ! IN
                                 NHVP_meas,HVP_meas,                           &
                                 RIN%NW,RIN%WAVE,AP(:,ipix),                   &
                                 GOUT%retrieval%fit%segment_fit%pixels(ipix),  & ! INOUT
                                 GOUT%aerosol,                                 &
                                 !GOUT%clouds,                                  &
                                 !GOUT%surface,                                 &
                                 GOUT%forcing%bbflux%pixel(ipix),              &
                                 GOUT%forcing%forcing%pixel(ipix),             &
                                 KERNELS1,KERNELS2                             &
                              )
!!XH         temporary output of forcing product
!            print *, 'Pixel at latitude: ',lat1,' longitude: ',h1,' pixel # ',ipix
!            print *, '     HLV(km)   Flux0Up(w/m^2) Flux0Down(w/m^2)   Flux Up(w/m^2) Flux Down(w/m^2) NetForcing(w/m^2)'
!            do I = GOUT%forcing%bbflux%pixel(ipix)%NHLV, 1, -1
!               print *, GOUT%forcing%bbflux%pixel(ipix)%HLV(I), GOUT%forcing%bbflux%pixel(ipix)%BBUFX0(I),  &
!                        GOUT%forcing%bbflux%pixel(ipix)%BBDFX0(I), GOUT%forcing%bbflux%pixel(ipix)%BBUFXA(I),  &
!                        GOUT%forcing%bbflux%pixel(ipix)%BBDFXA(I), GOUT%forcing%forcing%pixel(ipix)%NETFORC(I)
!            end do
         end do ! ipix=1,npixels
         GOUT%products%forcing%bbflux  = .true.
         GOUT%products%forcing%forcing = .true.
      end if
!XH   end of broadband flux calculation
      
      if(RIN%IPRI_verbose .eqv. .true.) then
      if(RIN%IMQ .eq. 3) then
        write(iu_main_output,'(3a,f16.5,/)')  &
        'Solver ',trim(sparse_solver_name),'  CPU_time(sec) =',solver_timer
      endif
      endif

!AL ************ Calculation of aerosol particulate matter ********************
      if (RIN%products%aerosol%PM) then
!        write(*,*) 'before PM!'
        do ipm=1,2
          if(RIN%PM_diam(ipm) .gt. 0.0) then
            call get_PM(RIN%PM_diam(ipm),RIN,GOUT,segment_meas,GOUT%aerosol%pm%pixel(:)%PM(ipm))
          endif ! PM_diam gt 0.0
        enddo !ipm
        GOUT%products%aerosol%pm = .true.
!        stop
!        write(*,*) 'PM asserted!'
      endif
!AL ************ Aerosol classification ********************
      if (RIN%products%aerosol%types) then
      if (GOUT%products%aerosol%sd2m_mph .or. GOUT%products%aerosol%sd2m_mph) then
          call get_aerosol_types(RIN,GOUT,segment_meas%npixels,GOUT%aerosol%types%pixel(:)%index)
          GOUT%products%aerosol%types = .true.
!      write(*,*) 'types asserted!'
      endif
      endif

!***********************************************************
!***  Calculation of Error and Bias Estimates     ************************************
!***********************************************************
      if( .not. RIN%ISTOP .and. (RIN%products%errest%par .or. &
          RIN%products%errest%aerosol%opt .or. RIN%products%errest%aerosol%lidar) ) then
      if((INVSING .gt. 0 .and. errest_segm) .or.   &
         (INVSING .eq. 0 .and. any(errest_pix(1:npixels))) ) then 
         if(RIN%IMQ .eq. 1) then
            CALL UF_nonzero ( RIN%KNSINGF,npixels,UFS,    & ! IN
                              ccor,SMMULTI,IP2,           &
                              UFNZ,nnz                    & ! OUT
                            )
            if ( stop_report%status ) return
            nnz_err = nnz
            FFMI0(:) = 0.0
            DO ipix=1,npixels
              DO I=1,RIN%KNSINGF
                FFMI0((ipix-1)*RIN%KNSINGF+I) = FFMSQ(I,ipix)
              ENDDO ! I
            ENDDO ! ipix
            b(1:KNF) = FFMI0(1:KNF)
            call sparse_solver (  iu_main_output, & ! IN 
                                  KNF,nnz,UFNZ,      &
                                  b,                 & ! INOUT
                                  solver_timer       &							   
                               )
            if ( stop_report%status ) return
            QIM(1:KNF) = b(1:KNF)
         endif ! RIN%IMQ .eq. 1 

!j=1
!write(*,*) j,'  j  UFS:'
!do i=1,RIN%KNSING
!write(*,*) 'i=',i
!write(*,'(10f12.4)') UFS(1:RIN%KNSING,i,j)
!enddo 
!write(*,*)

         !write(*,*) 'before PAR_ERR_estimates'
         call output_err_estim_initialization ( GOUT%errest )
         
         if(RIN%products%errest%par) then
            solver_timer = 0.0
            call PAR_ERR_estimates (iu_main_output,          & ! IN
                                    sparse_solver_name,      &
                                    IP2,                     &
                                    UF,UFS,nnz_err,          &
                                    QS,QIM,AGEN,ALS,         & !ALSP,
                                    RIN,INVSING,segment_meas, &
                                    errest_pix,              &
                                    UFNZ,GOUT%errest%par,    & ! INOUT
                                    solver_timer             &
                                   )
            if ( stop_report%status ) return
            GOUT%products%errest%par = .true.
            !write(*,*) 'after  PAR_ERR_estimates'
            if(RIN%IPRI_verbose .eqv. .true.) then
            if(INVSING .eq. -2) then
            if(RIN%IMQ .eq. 3) then
              write(iu_main_output,'(3a,f16.5,/)') 'PAR_ERR_estimates: solver ', &
              trim(sparse_solver_name),'  CPU_time(sec)=',solver_timer
            endif
            endif
            endif
         endif
         if(RIN%products%errest%aerosol%opt .or. RIN%products%errest%aerosol%lidar) then
            solver_timer = 0.0
!            write(*,*) 'before OPT_ERR_estimates, AP:'
!            write(*,'(10e14.5)') AP(1:RIN%KNSING,1)
            call OPT_ERR_estimates (iu_main_output,         & ! IN
                                  sparse_solver_name,     &
                                  UF, UFS, nnz_err, UFNZ, &
                                  QS, QIM, AGEN, ALS, KM1_pix,  &
                                  APMIN, APMAX,                 &
                                  RIN, INVSING, segment_meas,   &
                                  errest_pix, AP,               &
                                  tau_mol, NHVP_meas, HVP_meas, &
                                  GOUT,                   &
                                  solver_timer,           &
                                  KERNELS1,KERNELS2       &
                                  )
            if ( stop_report%status ) return
            !if(RIN%products%errest%aerosol%opt) &
            GOUT%products%errest%aerosol%opt = .true.
            if(RIN%products%errest%aerosol%lidar .and. (RIN%DLSF%keyEL .gt. 0)) &
            GOUT%products%errest%aerosol%lidar = .true.
            !write(*,*) 'after  OPT_ERR_estimates'
            if(RIN%IPRI_verbose .eqv. .true.) then
            if(INVSING .eq. -2) then
            if(RIN%IMQ .eq. 3) then
              write(iu_main_output,'(3a,f16.5,/)') 'OPT_ERR_estimates: solver ', &
              trim(sparse_solver_name),'  CPU_time(sec)=',solver_timer
            endif
            endif
            endif
         endif
      else
         if(RIN%IPRI_verbose) then
            write(*,*) 'Warning in inversion !!!' 
            write(*,*) 'Do not calculate Error and Bias Estimates because of '
            write(*,*) 'number of measurements is less than number of retrieved parameters.'
            write(*,'(x,a,i3,a)') 'INVSING=',INVSING,'  errest_segm = .F. or errest_pix(1:npixels) = .F.'
         endif !
      endif ! INVSING .eq. 2 .and. errest_segm) .or.
      endif ! .not. RIN%ISTOP .and. (RIN%product%errest_mph .or.

      if( inversion_run_count .eq. -1 ) then
! can be deleted after testing edges
        if( .not. RIN%INPUT ) then
        write(*,*) 'input file : ',trim(RIN%DLSF%internal_file_path)//"input_iguess_1.dat"
        open (33, FILE=trim(RIN%DLSF%internal_file_path)//"input_iguess_1.dat",status='unknown')
          do I=1,RIN%KNSING
            write(33,*) I,( GOUT%retrieval%par%pixel(1:npixels)%par(I) )
          enddo
        close (33)
        endif
      elseif(inversion_run_count .eq. -2) then
! can be deleted after testing edges
        if( .not. RIN%INPUT ) then
        open (33, FILE=trim(RIN%DLSF%internal_file_path)//"input_iguess_2.dat",status='unknown')
          do I=1,RIN%KNSING
            write(33,*) I,( GOUT%retrieval%par%pixel(1:npixels)%par(I) )
          enddo
        close (33)
        endif
      endif

!***************************************************************************
!*** Print simulated data into file         ********************************
      if(RIN%ISTOP .and. RIN%sdata_sim_file .ne. '') then
         !call write_simul_pixels (  iu_main_output,  & ! IN
         !                           RIN%sdata_sim_file, & 
         !                           npixels,            &
         !                           GOUT%retrieval%fit%segment_fit%pixels, &
         !                           index_clouds        &
         !                        )
         status_funct = &
				 write_sdata_pixels ( RIN%sdata_sim_file, npixels, GOUT%retrieval%fit%segment_fit )
      endif ! RIN%ISTOP .and. RIN%sdata_sim_file

!*** FITTING (FS and FPS (measurements and modeled measurements calculated for retrieved parameters)) ***
      !if(RIN%IPRI_verbose) then
        !if(RIN%products%retrieval%fit)  then
            !call print_fitting_FS ( iu_main_output, RIN, segment_meas, &
            !                        segment_vec_meas, segment_vec_fit )
        !endif ! RIN%products%retrieval%fit
      !endif ! RIN%IPRI_verbose

!*** FITTING (FS and FPS (measurements and modeled measurements calculated for retrieved parameters)) ***
      if(RIN%IPRI_verbose) then
        !if(.not. RIN%products%retrieval%fit)  then
! If aplicable, add random noise to segment measurement vector
          if(any(RIN%NOISE%SGMS(1:RIN%NOISE%INOISE) .gt. 0.0)) then
            call print_fitting_FS ( iu_main_output, RIN, segment_meas, &
                                    segment_vec_meas, segment_vec_fit )
          endif ! any(RIN%NOISE%SGMS
        !endif ! RIN%products%retrieval%fit
      endif ! RIN%IPRI_verbose
!***************************************************************************
!     Deallocate arrays UF,UF1

      if(allocated(UF)) then
         deallocate(UF,stat=alloc_stat)
         if (alloc_stat /= 0) stop 'error while trying to deallocate UF'
      endif
      !if(allocated(UF1)) then
         !deallocate(UF1,stat=alloc_stat)
         !if (alloc_stat /= 0) stop 'error while trying to deallocate UF1'
      !endif
!	-------------------------------------------------------------------------
      if(RIN%main_output_file .ne. '-') close(iu_main_output) 
!	-------------------------------------------------------------------------

!     Deallocate Compressed Column Storage (sparse matrix)
      call deallocate_sparse_matrix_storage ( UFNZ ) 
!	-------------------------------------------------------------------------
!	-------------------------------------------------------------------------
      !write(iu_main_output,'(a,i4,a)') '*****   end inversion_run_count =', inversion_run_count,' **********'
      !inversion_run_count = inversion_run_count +1

      return
      end subroutine inversion

! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

end module mod_grasp 


