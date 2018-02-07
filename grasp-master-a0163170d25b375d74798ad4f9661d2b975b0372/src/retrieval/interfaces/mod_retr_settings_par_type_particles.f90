! **
! **  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.  
! **  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
! **

! mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm

module mod_par_type_aerosol

  implicit none
  
#ifdef WARN_DRY
#warning "__CHARACTERISTIC_TYPE__ duplicated"
#endif
        
  integer,parameter ::	par_type_aerosol_beg = 10000
! Size Distribution 
      integer,parameter ::	par_type_SD_beg = 10100
        integer,parameter ::	par_type_SD_TB = 10101   ! Size distribution dV/dlnr at "Triangle" Bins
        integer,parameter ::	par_type_SD_LB = 10102   ! Size distribution dV/dlnr for precompued Lognormal Bins
        integer,parameter ::	par_type_SD_LN = 10103   ! Parameters of LogNormal size distribution dV/dlnr
      integer,parameter ::	par_type_SD_end = 10200
! Real part of complex Refractive Index or complex refractive index for Chemistry and mixture
      integer,parameter ::	par_type_RERI_beg = 10200
        integer,parameter ::	par_type_RERI_spect = 10201  ! Spectral dependent real part of complex refractive index
        integer,parameter ::	par_type_RERI_const = 10202  ! Real part of complex refractive index is constant
        integer,parameter ::	par_type_CXRI_nmix  = 10203  ! Complex refractive index is mixture.
        integer,parameter ::	par_type_CXRI_chem  = 10204  ! Complex refractive index. Chemistry: rh,finslbl,fsoot,firon
      integer,parameter ::	par_type_RERI_end = 10300
! Imaginary part of complex Refractive Index
      integer,parameter ::	par_type_IMRI_beg = 10300
        integer,parameter ::	par_type_IMRI_spect = 10301  ! Spectral dependent imaginary part of complex refractive index
        integer,parameter ::	par_type_IMRI_const = 10302  ! Imaginary part of complex refractive index is constant
      integer,parameter ::	par_type_IMRI_end = 10400
! Particles shape (nonsphericity)
      integer,parameter ::	par_type_SHD_beg = 10400
        integer,parameter ::	par_type_SHD_fsph  = 10401     ! Fraction of spherical particles
        integer,parameter ::	par_type_SHD_distr = 10402     ! Axis Ratio Distribution
        !integer,parameter ::	par_type_SHD_prof  = 10403     ! Sphericity profile
      integer,parameter ::	par_type_SHD_end = 10500
! Aerosol profile
      integer,parameter ::	par_type_AVP_beg = 10500
        integer,parameter ::	par_type_AVP_par_height = 10501  ! Parameter of Aerosol Vertical Profile (...)
        integer,parameter ::	par_type_AVP_prof       = 10502  ! Normalized Aerosol Vertical Profile
      integer,parameter ::	par_type_AVP_end = 10600
! Aerosol concentration
      integer,parameter ::	par_type_Cv_beg = 10600
        integer,parameter ::	par_type_Cv = 10601     ! Aerosol concentration
      integer,parameter ::	par_type_Cv_end = 10700
! Additional parameters
      integer,parameter ::	par_type_ADD_beg = 10700
        integer,parameter ::	par_type_CL          = 10701  ! Calibration coefficient for lidar 
        integer,parameter ::	par_type_AVP_par_std = 10702  ! Standard deviation for vertical profile 
      integer,parameter ::	par_type_ADD_end = 10800
  integer,parameter ::	par_type_aerosol_end = 20000

end module mod_par_type_aerosol

! mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm

module mod_par_type_cloud

  implicit none
  
  integer,parameter ::	par_type_cloud_beg = 30000
! Size Distribution 
      integer,parameter ::	par_type_SD_cloud_beg = 30100
        integer,parameter ::	par_type_SD_cloud_par = 30101   ! Parameters of Cloud Size Distribution 
      integer,parameter ::	par_type_SD_cloud_end = 30200
! Cloud profile
      integer,parameter ::	par_type_CVP_beg = 30400
        integer,parameter ::	par_type_CVP_par    = 30401  ! Parameters of Cloud Vertical Profile (height)
      integer,parameter ::	par_type_CVP_end = 30500
! Cloud fraction
      integer,parameter ::	par_type_fcloud_beg = 30500
        integer,parameter ::	par_type_fcloud_par    = 30501  ! Pixel Cloud Fraction
      integer,parameter ::	par_type_fcloud_end = 30600
  integer,parameter ::	par_type_cloud_end = 40000

      contains
      
      subroutine dummy_subroutine_for_avoiding_warnings_in_an_empty_module()
      end subroutine dummy_subroutine_for_avoiding_warnings_in_an_empty_module                

end module mod_par_type_cloud

! mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
