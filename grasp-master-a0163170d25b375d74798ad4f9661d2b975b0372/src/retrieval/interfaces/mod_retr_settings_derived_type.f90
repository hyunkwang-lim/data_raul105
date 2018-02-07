! **
! **  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.  
! **  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
! **

! mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm

#ifdef WARN_DRY
#warning "__RETRIEVAL_SETTINGS_DEFINITION__ binded"
#endif                  
                
module mod_retr_settings_derived_type
      use iso_c_binding
      use mod_par_inv, only   : KIDIM1,KIDIM2,KIDIM3,   &
                                KKNOISE,KPARS,KW,KWM,KIP            
      use mod_par_OS,   only  : KSD
      use mod_globals, only  : GBL_FILE_PATH_LEN
      use mod_par_type_aerosol
      use mod_par_type_surface
      use mod_par_type_cloud 
                 	  
      implicit none

!	----------------------------------------------------------------------------------------
! DLS_bin code parameters: 
! keyEL         - number of scattering matrix elements  
! IWL=0, key=0  - use original Kernels 
! IWL=1, key=2  - use spectral Kernels for precalculated lognormal bins     
     type, bind(C) :: iP_flags_for_DLS
        integer(kind=C_INT)             ::  IWL
        integer(kind=C_INT)             ::  key
        integer(kind=C_INT)             ::  keyEL
        character(kind=C_CHAR,len=GBL_FILE_PATH_LEN)  ::  distname_O
        character(kind=C_CHAR,len=GBL_FILE_PATH_LEN)  ::  distname_N
        character(kind=C_CHAR,len=GBL_FILE_PATH_LEN)  ::  internal_file_path
        character(kind=C_CHAR,len=GBL_FILE_PATH_LEN)  ::  external_file_path
     end type iP_flags_for_DLS 
!	----------------------------------------------------------------------------------------
!  The the atmospheric characteristics that are retrieved for EACH PIXES    

! NDIM1 - number of different retrieved characteristics
! NDIM2 - number of atmospheric layers
! NDIM3 - number of aerosol components at each layer
! NDIM4 - number of parameters used for size distribution  (e.g. number of size bins)
! ISTARSING - definition of the retrieved parameter places of in the segment of the vector of unknown AP for single pixel 
! par_type - aerosol parameter type
! par_retr - define if aerosol parameter is retrieved
     type, bind(C) :: retr_par_number_NDIM
         integer(kind=C_INT)          :: n1
         integer(kind=C_INT)          :: n2(KIDIM1)
         integer(kind=C_INT)          :: n3(KIDIM2,KIDIM1)
         integer(kind=C_INT)          :: ISTARSING(KIDIM2,KIDIM1)
         integer(kind=C_INT)          :: par_type(KIDIM1)
         logical(kind=C_BOOL)         :: par_retr(KIDIM1)
     end type retr_par_number_NDIM      
!	----------------------------------------------------------------------------------------
! The parameters define LIGHT SCATTERING REGIMEs for OS 
! of BOTH modeling and retrieving      !!!

! IMSC = 0 - multiple scattering regime for forward calculations
!      = 1 - single scattering regime 
!      = 2 - multiple scattering regime for derivative calculations
! NG1
! NG2
      type, bind(C) :: OSH_par
         integer(kind=C_INT) :: IMSC
         integer(kind=C_INT) :: NG
         integer(kind=C_INT) :: NN
         integer(kind=C_INT) :: NF
      end type OSH_par
!	----------------------------------------------------------------------------------------
! single_pixel_constraints
!      type, bind(C) :: single_pixel_constraints
!         integer(kind=C_INT) :: IO(KIDIM2,KIDIM1)
!         real(kind=C_FLOAT)  :: GSM(KIDIM2,KIDIM1)
!      end type single_pixel_constraints
!	----------------------------------------------------------------------------------------
! single_pixel_apriori_estimates
      type, bind(C) :: single_pixel_constraints_apriori
         integer(kind=C_INT) :: IO(KIDIM3,KIDIM2,KIDIM1)
         real(kind=C_FLOAT)  :: GSM(KIDIM3,KIDIM2,KIDIM1)
      end type single_pixel_constraints_apriori
!	----------------------------------------------------------------------------------------
! single_pixel_apriori_smoothness
      type, bind(C) :: single_pixel_constraints_smoothness
         integer(kind=C_INT) :: IO(KIDIM2,KIDIM1)
         real(kind=C_FLOAT)  :: GSM(KIDIM2,KIDIM1)
      end type single_pixel_constraints_smoothness
!	----------------------------------------------------------------------------------------
! multi_pixel_constraints
      type, bind(C) :: multi_pixel_constraints	  
         integer(kind=C_INT) :: IOT(KIDIM2,KIDIM1)
         integer(kind=C_INT) :: IOX(KIDIM2,KIDIM1)
         integer(kind=C_INT) :: IOY(KIDIM2,KIDIM1)
         real(kind=C_FLOAT)  :: GSMT(KIDIM2,KIDIM1)
         real(kind=C_FLOAT)  :: GSMX(KIDIM2,KIDIM1)
         real(kind=C_FLOAT)  :: GSMY(KIDIM2,KIDIM1)
      end type multi_pixel_constraints
!	----------------------------------------------------------------------------------------
! The parameters define NOISE
! INOISE  - the number of different noise sources              
! SGMS    - std of noise in i -th source                      
! INN     - EQ.1.THEN error is absolute with                  
!         - EQ.0 THEN error assumed relative
! DNN     - variation of the noise of the i-th source
! NMT     - number of meas types with particular noise
! MT      - measurement types
! NWLP    - number of wavelength in pixel for meas type
! IWLP    - index of wavelength in pixel
!	----------------------------------------------------------------------------------------
      type, bind(C) :: NOISE_par
         integer(kind=C_INT)  :: INOISE
         real(kind=C_FLOAT)   :: SGMS(KKNOISE)
         integer(kind=C_INT)  :: INN(KKNOISE)
         real(kind=C_FLOAT)   :: DNN(KKNOISE)
         integer(kind=C_INT)  :: NMT(KKNOISE)          ! number of meas types with particular noise
         integer(kind=C_INT)  :: MT(KIP,KKNOISE)       ! measurement types
         integer(kind=C_INT)  :: NWLP(KIP,KKNOISE)     ! number of wavelength in pixel for meas type  
         integer(kind=C_INT)  :: IWLP(KWM,KIP,KKNOISE) ! index of wavelength in pixel 
      end type NOISE_par
!	----------------------------------------------------------------------------------------
!  Inter-pixel  fitting
      type, bind(C) :: inter_pixel_fitting
         integer(kind=C_INT) :: INVSING
         logical(kind=C_BOOL):: TXY_group           
         real(kind=C_FLOAT)  :: TTEST
         real(kind=C_FLOAT)  :: XTEST
         real(kind=C_FLOAT)  :: YTEST
      end type inter_pixel_fitting
!	----------------------------------------------------------------------------------------
!  Edge sizes
      type, bind(C) :: edges_size
         integer(kind=C_INT) :: nx
         integer(kind=C_INT) :: ny        
         integer(kind=C_INT) :: nt         
      end type edges_size
!	----------------------------------------------------------------------------------------
!  Output products
#ifdef WARN_DRY
#warning "__RETRIEVAL_PRODUCTS_DEFINITION__ binded"
#endif      
                 
      type, bind(C) :: output_segment_products_retrieval
         logical(kind=C_BOOL) ::  res
         logical(kind=C_BOOL) ::  par
         logical(kind=C_BOOL) ::  fit
      end type output_segment_products_retrieval

      type, bind(C) :: output_segment_products_particles
         logical(kind=C_BOOL) ::  opt    ! ext ssa Aexp
         logical(kind=C_BOOL) ::  rind
         logical(kind=C_BOOL) ::  chem
         logical(kind=C_BOOL) ::  phmx   
         logical(kind=C_BOOL) ::  lidar  ! lr ldpr
         logical(kind=C_BOOL) ::  sd2m_mph
         logical(kind=C_BOOL) ::  sd2m_ext
         logical(kind=C_BOOL) ::  PM ! particulate matter
         logical(kind=C_BOOL) ::  types ! aerosol type
      end type output_segment_products_particles

      type, bind(C) :: output_segment_products_surface
         logical(kind=C_BOOL) ::  surf
      end type output_segment_products_surface

      type, bind(C) :: output_segment_products_forcing
         logical(kind=C_BOOL) ::  bbflux
         logical(kind=C_BOOL) ::  forcing
      end type output_segment_products_forcing

      type, bind(C) :: output_segment_products_errest_particles
         logical(kind=C_BOOL) ::  opt        ! ext ssa
         logical(kind=C_BOOL) ::  lidar      ! lr     
      end type output_segment_products_errest_particles

      type, bind(C) :: output_segment_products_errest
         logical(kind=C_BOOL)                           ::  par
         type(output_segment_products_errest_particles) ::  aerosol  
         type(output_segment_products_errest_particles) ::  clouds   
      end type output_segment_products_errest

      type, bind(C) :: output_segment_products
         type(output_segment_products_retrieval)  ::  retrieval
         type(output_segment_products_particles)  ::  aerosol
         type(output_segment_products_particles)  ::  clouds
         type(output_segment_products_surface)    ::  surface
         type(output_segment_products_forcing)    ::  forcing          
         type(output_segment_products_errest)     ::  errest
      end type output_segment_products
!	----------------------------------------------------------------------------------------
!  Retrieval INput structure (RIN)
      type, bind(C) :: retr_input_settings
        integer(kind=C_INT)            ::  KNSING 
        integer(kind=C_INT)            ::  KNSINGF                                             
        integer(kind=C_INT)            ::  KL                                             
        logical(kind=C_BOOL)           ::  ISTOP                                             
        logical(kind=C_BOOL)           ::  IPRI_additional_info                                             
        logical(kind=C_BOOL)           ::  IPRI_iter_result                                             
        logical(kind=C_BOOL)           ::  IPRI_verbose                                             
        integer(kind=C_INT)            ::  NSD                                             
        integer(kind=C_INT)            ::  NLYRS                                             
        integer(kind=C_INT)            ::  ipplane                                             
        integer(kind=C_INT)            ::  iPOBS                                             
        integer(kind=C_INT)            ::  iIOBS
        integer(kind=C_INT)            ::  isurf_land(2)
        integer(kind=C_INT)            ::  isurf_water
        integer(kind=C_INT)            ::  Aexp_iwl(2)
        integer(kind=C_INT)            ::  ndvi_iwl(2)
        real(kind=C_FLOAT)             ::  SHIFT        
        integer(kind=C_INT)            ::  NW       
        real(kind=C_FLOAT)             ::  WAVE(KW) 
        integer(kind=C_INT)            ::  IBIN
        integer(kind=C_INT)            ::  IMQ
        integer(kind=C_INT)            ::  IPSTOP
        logical(kind=C_BOOL)           ::  INPUT
        logical(kind=C_BOOL)           ::  ITRONC
        integer(kind=C_INT)            ::  MAXP
        real(kind=C_FLOAT)             ::  EPSP
        real(kind=C_FLOAT)             ::  EPSQ
        real(kind=C_FLOAT)             ::  DL
        integer(kind=C_INT)            ::  mol_prof_type ! type of the molecular vertical profile used
        integer(kind=C_INT)            ::  aer_prof_type ! type of the aerosol vertical profile used
        real(kind=C_FLOAT)             ::  PM_diam(2)    ! diameters at wich PM is calculated
        integer(kind=C_INT)            ::  nPM_diam
        logical(kind=C_BOOL)           ::  use_models  ! use modeled phase matrices instead of LB Kernels

        type(iP_flags_for_DLS)         ::  DLSF
        type(retr_par_number_NDIM)     ::  NDIM
        type(OSH_par)                  ::  OSHF
        type(OSH_par)                  ::  OSHD
        type(single_pixel_constraints_apriori)    ::  SPCA
        type(single_pixel_constraints_smoothness) ::  SPCS
        type(multi_pixel_constraints)  ::  MPCS        
        type(NOISE_par)                ::  NOISE
        type(inter_pixel_fitting)      ::  IPFP   ! INVSING,TTEST,XTEST,YTEST
! **
        real(kind=C_FLOAT)             ::  APSING(KPARS)
        real(kind=C_FLOAT)             ::  APSMIN(KPARS)
        real(kind=C_FLOAT)             ::  APSMAX(KPARS) 
        
        real(kind=C_FLOAT)             ::  RMIN(KSD),RMAX(KSD)
        real(kind=C_FLOAT)             ::  RATIO1(KIDIM3,KSD)
        real(kind=C_FLOAT)             ::  RADIUS1(KIDIM3,KSD)
        integer(kind=C_INT)            ::  IWW_SINGL(KPARS)
        integer(kind=C_INT)            ::  NBIN(KSD)

        type(edges_size)               ::  edges
                    
        character(kind=C_CHAR,len=GBL_FILE_PATH_LEN) ::  plotting_output_file
        character(kind=C_CHAR,len=GBL_FILE_PATH_LEN) ::  main_output_file
        character(kind=C_CHAR,len=GBL_FILE_PATH_LEN) ::  sdata_sim_file

! eps_err - absolute value of truncation threshold of Fourier and order-of-scattering 
!           series expansions in radiative transfer calculations
        real(kind=C_FLOAT)             ::  eps_err   
        
        type(output_segment_products)  :: products                    

        !character(kind=C_CHAR,len=GBL_FILE_PATH_LEN) ::  internal_file_path
        !character(kind=C_CHAR,len=GBL_FILE_PATH_LEN) ::  external_file_path         
	  end type retr_input_settings
!	----------------------------------------------------------------------------------------
!	----------------------------------------------------------------------------------------

      contains
      
      subroutine for_avoiding_warnings_in_an_empty_module_settings_derived_type()
      end subroutine for_avoiding_warnings_in_an_empty_module_settings_derived_type                

end module mod_retr_settings_derived_type

