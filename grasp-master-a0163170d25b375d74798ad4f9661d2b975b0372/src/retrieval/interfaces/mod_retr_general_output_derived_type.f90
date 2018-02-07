! **
! **  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.  
! **  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
! **

! mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm

      module mod_retr_general_output_derived_type
      
      use mod_time_utils, only : KIND_TIME 
      use mod_index_cloud
      use mod_retr_settings_derived_type
      use mod_par_inv,     only : KPARS,KIDIM3,KPAR,KVERTM,KW,KKNOISE
      use mod_par_DLS,     only : KMpar,KNpar
      use mod_par_OS,      only : KSD,KNT
      use mod_par_DLS_bin, only : NRR,NRC
      use mod_sdata_derived_type

      implicit none

#ifdef WARN_DRY
#warning "__RETRIEVAL_OUTPUT_DEFINITION__ binded"
#endif   
              
      !integer,parameter :: KF = 3 ! number of functions for retrieved parameters: 
                                  ! 1 - ext, 2 - ssa, 3 - lr
! Error estimates: indexes for UW0 matrix 
      integer,parameter	::	index_ext = 1
      integer,parameter	::	index_ssa = 2
      integer,parameter	::	index_lr  = 3
! -----------------------------------------------------------------------------------------
! -----------------------------------------------------------------------------------------
! Retrieval output
! -----------------------------------------------------------------------------------------
! Structure contains retrieval output data :
! res       - single pixel total residual (meas + smothness constrain )
! resa,resr - detailed absolut and relative measurement residuals for single pixel 
! niter       - number of iterations
! rest        - total residual for multi-pixel retrieval
! resat,resrt - detailed absolut and relative measurement residuals for segment
! par         - retrieved aerosol and surface reflectance parameters
! ngrid       - nimber of grid radii for SD (saved to print output)
! radius      - grid radii (saved to print output)
! SDL         - if applicable, precalculated lognormal bins SD (saved to print output)
      type, bind(C) :: output_pixel_residual 	
         integer(kind=C_INT) ::  niter
         real(kind=C_FLOAT)  ::  res
         real(kind=C_FLOAT)  ::  resa(KKNOISE)
         real(kind=C_FLOAT)  ::  resr(KKNOISE)
      end type output_pixel_residual
      type, bind(C) :: output_segment_residual 
         integer(kind=C_INT)      ::  niter
         real(kind=C_FLOAT)       ::  rest  
         real(kind=C_FLOAT)       ::  resat(KKNOISE) 
         real(kind=C_FLOAT)       ::  resrt(KKNOISE)
         type(output_pixel_residual)  ::  pixel(KIMAGE)  
      end type output_segment_residual

      type, bind(C) :: output_pixel_retr_par 	
         real(kind=C_FLOAT)  ::  par(KPARS)
      end type output_pixel_retr_par
      type, bind(C) :: output_segment_retr_par 
! ngrid - number of grid radii for precalculated lognormal bins
! radius(NRR) - grid radii
! SDL(NRR,NRC) - values of SD for each precalculated lognormal bin
         integer(kind=C_INT)          ::  ngrid
         real(kind=C_FLOAT)           ::  radius(NRR)
         real(kind=C_FLOAT)           ::  SDL(NRR,NRC)
         type(output_pixel_retr_par)  ::  pixel(KIMAGE)  
      end type output_segment_retr_par

! -----------------------------------------------------------------------------------------
! Fitting
! -----------------------------------------------------------------------------------------
      type, bind(C) :: output_segment_fitting
         type(segment_data) :: segment_fit
      end type output_segment_fitting
! -----------------------------------------------------------------------------------------

      type, bind(C) :: output_segment_retrieval 
         type(output_segment_residual)  ::  res  
         type(output_segment_retr_par)  ::  par  
         type(output_segment_fitting)   ::  fit
      end type output_segment_retrieval
! -----------------------------------------------------------------------------------------
! -----------------------------------------------------------------------------------------
! Optical characteristics      
! -----------------------------------------------------------------------------------------
! ext,ssa,aext   - spectral ext, ssa and aext for each aerosol component
! Aexp      - Angstrom exponent (wl(4)/wl(5))
      type, bind(C) :: output_pixel_opt_wl 	
         real(kind=C_FLOAT)  ::  extt 
         real(kind=C_FLOAT)  ::  ssat 
         real(kind=C_FLOAT)  ::  aextt
         real(kind=C_FLOAT)  ::  ext(KSD)
         real(kind=C_FLOAT)  ::  ssa(KSD)
         real(kind=C_FLOAT)  ::  aext(KSD)
      end type output_pixel_opt_wl
      type, bind(C) :: output_pixel_opt 	
         real(kind=C_FLOAT)         ::  Aexp
         type(output_pixel_opt_wl)  ::  wl(KW)
      end type output_pixel_opt
      type, bind(C) :: output_segment_opt 
         type(output_pixel_opt)  ::  pixel(KIMAGE)
      end type output_segment_opt
! -----------------------------------------------------------------------------------------
! Refractive index
      type, bind(C) :: output_pixel_rindex_wl 	
         real(kind=C_FLOAT)  ::  mreal(KSD)
         real(kind=C_FLOAT)  ::  mimag(KSD)
      end type output_pixel_rindex_wl
      type, bind(C) :: output_pixel_rindex 	
         type(output_pixel_rindex_wl)  ::  wl(KW)
      end type output_pixel_rindex
      type, bind(C) :: output_segment_rindex 
         type(output_pixel_rindex)  ::  pixel(KIMAGE)
      end type output_segment_rindex
! -----------------------------------------------------------------------------------------
! Phase matrix :
      type, bind(C) :: output_pixel_ph_matrix_wl      
         real(kind=C_FLOAT),dimension(KMpar,KSD) ::  ph11,ph12,ph22,ph33,ph34,ph44
         real(kind=C_FLOAT),dimension(KMpar)     ::  pht11,pht12,pht22,pht33,pht34,pht44
      end type output_pixel_ph_matrix_wl
      type, bind(C) :: output_pixel_ph_matrix      
         type(output_pixel_ph_matrix_wl)  ::  wl(KW)
      end type output_pixel_ph_matrix
      type, bind(C) :: output_segment_ph_matrix 
         integer(kind=C_INT)                 ::  nangle
         real(kind=C_FLOAT),dimension(KMpar) ::  angle
         type(output_pixel_ph_matrix)        ::  pixel(KIMAGE)
      end type output_segment_ph_matrix
! Lidar and depolarization ratios :
      type, bind(C) :: output_pixel_lidar_ratio_wl      
         real(kind=C_FLOAT),dimension(KSD)   ::  lr,ldpr
         real(kind=C_FLOAT)                  ::  lrt,ldprt 
      end type output_pixel_lidar_ratio_wl
      type, bind(C) :: output_pixel_lidar_ratio      
         type(output_pixel_lidar_ratio_wl)  ::  wl(KW)
      end type output_pixel_lidar_ratio
      type, bind(C) :: output_segment_lidar_ratio 
         type(output_pixel_lidar_ratio)     ::  pixel(KIMAGE)
      end type output_segment_lidar_ratio
! -----------------------------------------------------------------------------------------
! Chemistry parameters
      type, bind(C) :: output_pixel_chemistry      
!         real(kind=C_FLOAT)  ::  rh(KSD)
!         real(kind=C_FLOAT)  ::  fwtr(KSD)
!         real(kind=C_FLOAT)  ::  fslbl(KSD)
!         real(kind=C_FLOAT)  ::  finslbl(KSD)
!         real(kind=C_FLOAT)  ::  fsoot(KSD)
!         real(kind=C_FLOAT)  ::  firon(KSD)

!***********14/11/2016
         real(kind=C_FLOAT)  ::  rh(KSD)
         real(kind=C_FLOAT)  ::  fwtr(KSD)
         real(kind=C_FLOAT)  ::  fslbl(KSD)
         real(kind=C_FLOAT)  ::  fquartz(KSD)
         real(kind=C_FLOAT)  ::  fsoot(KSD)
         real(kind=C_FLOAT)  ::  firon(KSD)
         real(kind=C_FLOAT)  ::  fbrc(KSD)
!***********14/11/2016


      end type output_pixel_chemistry
      type, bind(C) :: output_segment_chemistry 
         type(output_pixel_chemistry)  :: pixel(KIMAGE)
      end type output_segment_chemistry
! -----------------------------------------------------------------------------------------
! Two mode aerosol characteristics
! Structure contains output data (microphysical parameters):
! 0 - total, 1 - fine mode, 2 - coarse mode 
! reff        - volume median radius
! std         - standard deviation
! cv          - concentration 
! rm          - median radius 
! ext         - ext each aerosol component
      type, bind(C) :: output_pixel_sd2m_mph 	
         real(kind=C_FLOAT),dimension(0:2)     ::  cv,std,rm,reff
      end type output_pixel_sd2m_mph 
      type, bind(C) :: output_segment_sd2m_mph 
         type(output_pixel_sd2m_mph)  ::  pixel(KIMAGE)
      end type output_segment_sd2m_mph

      type, bind(C) :: output_pixel_sd2m_opt 	
         real(kind=C_FLOAT),dimension(KW,1:2)  ::  ext
      end type output_pixel_sd2m_opt 
      type, bind(C) :: output_segment_sd2m_opt 
         type(output_pixel_sd2m_opt)  ::  pixel(KIMAGE)
      end type output_segment_sd2m_opt

      type, bind(C) :: output_segment_sd2m 
         type(output_segment_sd2m_mph)  ::  mph
         type(output_segment_sd2m_opt)  ::  opt
      end type output_segment_sd2m      
! -----------------------------------------------------------------------------------------
! -----------------------------------------------------------------------------------------
! Particulate Matter
      type, bind(C) :: output_pixel_PM
         real(kind=C_FLOAT)  ::  pm(2)
      end type output_pixel_PM
      type, bind(C) :: output_segment_PM
         type(output_pixel_PM)  :: pixel(KIMAGE)
      end type output_segment_PM
! -----------------------------------------------------------------------------------------
! Typing (aerosol and presumably clouds)
! 0 – Complex mixture
! 1 – Background Aerosol
! 2 – Water/Maritime
! 3 — Urban Polluted
! 4 – Mixed aerosol
! 5 – Urban Clean
! 6 – Smoke Smoldering
! 7 – Smoke flaming
! 8 – Mineral dust
      type, bind(C) :: output_pixel_types
         integer(kind=C_INT)  ::  index  ! describes aerosol type
      end type output_pixel_types
      type, bind(C) :: output_segment_types
         type(output_pixel_types)  :: pixel(KIMAGE)
      end type output_segment_types
! -----------------------------------------------------------------------------------------

! -----------------------------------------------------------------------------------------
! Surface      
! -----------------------------------------------------------------------------------------
! Surface albedo 
      type, bind(C) :: output_pixel_surface_wl      
         real(kind=C_FLOAT)  ::  salbedo
      end type output_pixel_surface_wl
      type, bind(C) :: output_pixel_surface 
         real(kind=C_FLOAT)             ::  NDVI           
         type(output_pixel_surface_wl)  ::  wl(KW)
      end type output_pixel_surface
      type, bind(C) :: output_segment_surface 
         type(output_pixel_surface)  :: pixel(KIMAGE)
      end type output_segment_surface
! -----------------------------------------------------------------------------------------
! -----------------------------------------------------------------------------------------
! Radiative forcing
! -----------------------------------------------------------------------------------------
! Radiative broad band flux and forcing
      type, bind(C) :: output_pixel_bbflux
         integer(kind=C_INT)                 ::  NHLV
         real(kind=C_FLOAT), dimension(KNT)  ::  BBUFX0,BBDFX0,BBUFXA,BBDFXA,HLV
      end type output_pixel_bbflux
      type, bind(C) :: output_segment_bbflux
         type(output_pixel_bbflux)  :: pixel(KIMAGE)
      end type output_segment_bbflux

      type, bind(C) :: output_pixel_forcing
         integer(kind=C_INT)                 ::  NHLV
         real(kind=C_FLOAT), dimension(KNT)  ::  NETFORC,FORCEFF,HLV
      end type output_pixel_forcing
      type, bind(C) :: output_segment_forcing
         type(output_pixel_forcing)  :: pixel(KIMAGE)
      end type output_segment_forcing

      type, bind(C) :: output_segment_rad_forcing
         type(output_segment_bbflux)   :: bbflux
         type(output_segment_forcing)  :: forcing
      end type output_segment_rad_forcing
! -----------------------------------------------------------------------------------------
! -----------------------------------------------------------------------------------------
! Error estimations
! -----------------------------------------------------------------------------------------
! ERRP - Standard deviations of retrieved parameter logarithms (~relative errors)
! BIASP - Standard deviation of systematic errors of retrieved parameter logarithms
! ERR_ - Standard deviations of retrieved optical characteristic logarithms (~relative errors) 
! BIAS_ - Standard deviations of systematic errors of retrieved optical characteristic logarithms
! structure par      contains BIAS & ERR for all retrieved (of aerosol, clouds and surface)
! structure aerosol1 contains BIAS & ERR for ext & ssa - optical thickness and single scattering albedo of aerosol
! structure aerosol2 contains BIAS & ERR for lr        - lidar ratio of aerosol
! structure cloud1   contains BIAS & ERR for ext & ssa - optical thickness and single scattering albedo of clouds
! structure cloud2   contains BIAS & ERR for lr        - lidar ratio of clouds
      type, bind(C) :: output_pixel_err_estim_par  
         real(kind=C_FLOAT) ::  ERRP(KPARS)
         real(kind=C_FLOAT) ::  BIASP(KPARS)
      end type output_pixel_err_estim_par 
      type, bind(C) :: output_segment_err_estim_par 
         type(output_pixel_err_estim_par)  ::  pixel(KIMAGE)
      end type output_segment_err_estim_par

      type, bind(C) :: output_pixel_err_estim_particles_opt  
         real(kind=C_FLOAT) ::  ERR_ext(KW,KSD)
         real(kind=C_FLOAT) ::  BIAS_ext(KW,KSD)
         real(kind=C_FLOAT) ::  ERR_extt(KW)
         real(kind=C_FLOAT) ::  BIAS_extt(KW)
         real(kind=C_FLOAT) ::  ERR_ssa(KW,KSD)
         real(kind=C_FLOAT) ::  BIAS_ssa(KW,KSD)
         real(kind=C_FLOAT) ::  ERR_ssat(KW)
         real(kind=C_FLOAT) ::  BIAS_ssat(KW)
      end type output_pixel_err_estim_particles_opt 
      type, bind(C) :: output_segment_err_estim_particles_opt 
         type(output_pixel_err_estim_particles_opt)  ::  pixel(KIMAGE)
      end type output_segment_err_estim_particles_opt
      
      type, bind(C) :: output_pixel_err_estim_particles_lidar  
         real(kind=C_FLOAT) ::  ERR_lr(KW,KSD)
         real(kind=C_FLOAT) ::  BIAS_lr(KW,KSD)       
         real(kind=C_FLOAT) ::  ERR_lrt(KW)
         real(kind=C_FLOAT) ::  BIAS_lrt(KW)
      end type output_pixel_err_estim_particles_lidar
      type, bind(C) :: output_segment_err_estim_particles_lidar 
         type(output_pixel_err_estim_particles_lidar)  ::  pixel(KIMAGE)
      end type output_segment_err_estim_particles_lidar

      type, bind(C) :: output_segment_err_estim_particles 
         type(output_segment_err_estim_particles_opt)    ::  opt
         type(output_segment_err_estim_particles_lidar)  ::  lidar
      end type output_segment_err_estim_particles

      type, bind(C) :: output_segment_err_estim 
         type(output_segment_err_estim_par)          ::  par
         type(output_segment_err_estim_particles)    ::  aerosol   
         type(output_segment_err_estim_particles)    ::  clouds     
      end type output_segment_err_estim
! -----------------------------------------------------------------------------------------
! -----------------------------------------------------------------------------------------
! Aerosol
! -----------------------------------------------------------------------------------------
!      type, bind(C) :: output_segment_aerosol 
      type, bind(C) :: output_segment_particles 
         type(output_segment_opt)         ::  opt
         type(output_segment_rindex)      ::  rind
         type(output_segment_ph_matrix)   ::  phmx
         type(output_segment_lidar_ratio) ::  lidar
         type(output_segment_sd2m)        ::  sd2m
         type(output_segment_chemistry)   ::  chem
         type(output_segment_PM)          ::  pm
         type(output_segment_types)       ::  types
      end type output_segment_particles
!      end type output_segment_aerosol
! -----------------------------------------------------------------------------------------
! -----------------------------------------------------------------------------------------
! Clouds
! -----------------------------------------------------------------------------------------
      !type, bind(C) :: output_segment_clouds 
         !type(output_segment_opt)         ::  opt
         !type(output_segment_rindex)      ::  rind
         !type(output_segment_ph_matrix)   ::  phmx
         !type(output_segment_lidar_ratio) ::  lidar
         !type(output_segment_sd2m)        ::  sd2m
         !type(output_segment_chemistry)   ::  chem
      !end type output_segment_clouds
! -----------------------------------------------------------------------------------------
! -----------------------------------------------------------------------------------------
! General output
! -----------------------------------------------------------------------------------------
! FRAMEWORK INTERFACE: The bound structure called output_segment_general is placed in output/grasp_output.h
      type, bind(C) :: output_segment_general
          type(output_segment_products)      ::  products
          type(output_segment_retrieval)     ::  retrieval
          type(output_segment_particles)     ::  aerosol
          type(output_segment_particles)     ::  clouds
          type(output_segment_surface)       ::  surface
          type(output_segment_err_estim)     ::  errest
          type(output_segment_rad_forcing)   ::  forcing
      endtype output_segment_general
! -----------------------------------------------------------------------------------------

      contains
      
      subroutine for_avoiding_warnings_in_an_empty_module_general_output_derived()
      end subroutine for_avoiding_warnings_in_an_empty_module_general_output_derived                

end module mod_retr_general_output_derived_type

! mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm


