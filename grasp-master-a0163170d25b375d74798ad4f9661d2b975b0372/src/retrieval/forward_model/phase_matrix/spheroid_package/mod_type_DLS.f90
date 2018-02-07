! **
! **  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.  
! **  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
! **

MODULE mod_type_DLS

      use mod_par_DLS,     only : KN1par,KM1par,KR1par,  &
	                                KNpar,KMpar,KRpar
      use mod_par_DLS_bin, only : KWLpar,KCpar
	  
      implicit none
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
! DLS_kernel_location: 
  
! distname_O    - directory name for original Kernels
! distname_N    - directory name for Spectral Kernels with lognormal bins
! NRATN         - number of axis ratios for reading Kernels
! comm_name(1:KR1par)  - common part of Kernel names (Axis ratio dependent) 
!
      type DLS_kernel_location

		   character(LEN=255)                    ::  internal_file_path
		   character(LEN=255)                    ::  external_file_path
		   character(len=255)                    ::  distname_O  
		   character(len=255)                    ::  distname_N  
		   character(len=255),dimension(KR1par)  ::  comm_name

      end type DLS_kernel_location      

! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
! DLS_kernel_grid1: 
  
! KN1    - number of grid radii in original Kernels
! KM1    - number of grid angles in original Kernels
! NRATN         - number of axis ratios for reading Kernels
! XWL    - wavelength of original Kernels
! grid1(1:KN1par)  - grid radii 
! ANGLE1(1:KM1par) - grid scattering angles
!
      type DLS_kernel_grid1

        integer                 :: KN1,KM1,NRATN
        real                    :: XWL
        real, dimension(KN1par) :: grid1
        real, dimension(KM1par) :: ANGLE1
        real, dimension(KR1par) :: RATIO

      end type DLS_kernel_grid1      

! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
! DLS_input_config : 
  
      type DLS_input_config

        integer                         ::  key,keyEL,keySUB,keyLS,key_RD1,key_RD, &
                                            KN,KM,KR
        real                            ::  xnmin,xnmax,xkmin,xkmax,  &		                                      
                                            pomin,pomax
        real,   dimension (KNpar)       ::  grid              
        real*4, dimension (KRpar)       ::  R
        real,   dimension (KRpar)       ::  RD
        real,   dimension (KMpar)       ::  ANGLE

        integer                         ::  KC,NWL,LB,LE                        																		
        real,   dimension(KCpar)        ::  RC2  
        real,   dimension(KMpar,KWLpar) ::  f11_bin,f12_bin,f22_bin,   &
                                            f33_bin,f34_bin,f44_bin
        real,   dimension(KWLpar)       ::  WL_bin,RN_bin,RK_bin,      & 
                                            xext_bin,xabs_bin,xsca_bin,albedo_bin
        real,   dimension(KCpar,KNpar)  ::  SD_bin
        integer,dimension(KWLpar)       ::  IWL1 

      end type DLS_input_config      

! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
! DLS_config : 
  
      type DLS_CONFIG

       type(DLS_input_config)    :: DLSIN
		   type(DLS_kernel_location) :: KERNEL
		   type(DLS_kernel_grid1)    :: XGRID1
		   
      end type DLS_CONFIG     



! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
      contains
      
        subroutine dummy_subroutine_for_avoiding_warnings_in_an_empty_module()
   	end subroutine dummy_subroutine_for_avoiding_warnings_in_an_empty_module
        
END MODULE mod_type_DLS


