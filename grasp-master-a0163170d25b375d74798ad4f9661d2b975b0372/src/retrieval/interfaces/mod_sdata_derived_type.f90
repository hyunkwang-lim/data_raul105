! **
! **  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.  
! **  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
! **

! mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
module mod_sdata_derived_type

      use mod_time_utils
      use iso_c_binding
      use mod_par_OS,  only : NBVM
      use mod_par_inv, only : KWM,KW,KVERTM,KNBVM,KIP,KSURF, &
                              KITIME,KIX,KIY,KIMAGE,KMESS,KKNOISE

      implicit none 
      ! KWM	 ! number of wavelengths for measurements
      ! NBVM	 ! number of directions 
      ! KIP	 ! number of measurement types
      ! KSURF ! 1
      ! KVERTM ! maximum of number of measured profile vert. heights and retrieved profile vert. heights
      ! KNBVM  ! max number of measurements for single wavelength
!	------------------------------------------------------------------------------------------------------
#ifdef WARN_DRY
#warning "WARN_DRY: __CHARACTERISTIC_TYPE__ duplicated (binded)"
#endif        
	type, bind(C) :: data_wl 
! *** data_wl --  type contains all wavelength dependent data

      integer(kind=C_INT)         ::  meas_type(KIP)
      real(kind=C_FLOAT)  	      ::  wl
      integer(kind=C_INT)         ::  ind_wl   ! index of wl in wave length array for retrieved parameters
      real(kind=C_FLOAT)          ::  sza      ! solar zenith angle
      real(kind=C_FLOAT)          ::  thetav(NBVM,KIP)   ! vis_IM(ITIME,IX,IY,IV)=thetav(IV)
      real(kind=C_FLOAT)          ::  phi(NBVM,KIP)      ! fis_IM(ITIME,IX,IY,IV)=phi(IV)
      integer(kind=C_INT)         ::  Nsurf    ! Number of surface components ??? 
      real(kind=C_FLOAT)          ::  groundpar(KSURF)
      real(kind=C_FLOAT)          ::  gaspar
 
      integer(kind=C_INT)         ::  NBVM(KIP)     ! number of valid measurements
      integer(kind=C_INT)         ::  NIP

      real(kind=C_FLOAT)          ::  tau(NBVM)  ! optical thickness;            meas_type = 11/12
      real(kind=C_FLOAT)          ::  p11(NBVM)  ! p11 phase matrix element;     meas_type = 21
      real(kind=C_FLOAT)          ::  p12(NBVM)  ! p11 phase matrix element;     meas_type = 22
      real(kind=C_FLOAT)          ::  p22(NBVM)  ! p11 phase matrix element;     meas_type = 23
      real(kind=C_FLOAT)          ::  p33(NBVM)  ! p11 phase matrix element;     meas_type = 24
      real(kind=C_FLOAT)          ::  p34(NBVM)  ! p11 phase matrix element;     meas_type = 25
      real(kind=C_FLOAT)          ::  p44(NBVM)  ! p11 phase matrix element;     meas_type = 26
          
      real(kind=C_FLOAT)          ::  LS(KVERTM) ! lidar signal;                 meas_type = 31
      real(kind=C_FLOAT)          ::  DP(KVERTM) ! depolarization ratio;         meas_type = 32
      real(kind=C_FLOAT)          ::  RL(KVERTM) ! Raman lidar signal;           meas_type = 33
          
      real(kind=C_FLOAT)          ::  I(NBVM)    ! I Stokes parameter;           meas_type = 41
      real(kind=C_FLOAT)          ::  Q(NBVM)    ! Q Stokes parameter;           meas_type = 42
      real(kind=C_FLOAT)          ::  U(NBVM)    ! U Stokes parameter;           meas_type = 43
      real(kind=C_FLOAT)          ::  P(NBVM)    ! linear polarization sqrt(Q*Q+U*U);   meas_type = 44
                      
      real(kind=C_FLOAT)          ::  CMTRX(KNBVM,KIP)   ! diagonal of covariance matrix (also known as OMEGA)
      real(kind=C_FLOAT)          ::  MPROF(KVERTM,KIP)  ! vertical profile of Rayleigh backscatter (beta_m) 
      integer(kind=C_INT)         ::  IFCOV(KIP)    ! 0/1 presence of covar. matrix in input data file
      integer(kind=C_INT)         ::  IFMP(KIP)
   
   end type data_wl

!	------------------------------------------------------------------------------------------------------
!	------------------------------------------------------------------------------------------------------

	type, bind(C) :: pixel 

    real(kind=C_FLOAT)      ::  HOBS ! height of observation
		integer(kind=C_INT)     ::	nwl
		!C*** ICLOUD(KITIME,KX,KY) - cloud index for each pixes
		!C***        =0  - cloudy  		
		!C***        =1  - cloud free
		integer(kind=C_INT)     ::	cloudy
		real(kind=C_FLOAT)      ::	x   ! position in input grid
		real(kind=C_FLOAT)      ::	y
		integer(KIND_TIME)      ::	t 

		integer(kind=C_INT)  	  ::	ix   ! position in segment
		integer(kind=C_INT)  	  ::	iy
		integer(kind=C_INT)  	  ::	it
		integer(kind=C_INT)  	  ::  out_x; ! position in output grid
		integer(kind=C_INT)  	  ::  out_y;
		integer(kind=C_INT)  	  ::  out_t;                
		real(kind=C_FLOAT)      ::  MASL ! Metres Above Sea Level
		real(kind=C_FLOAT)  	  ::  land_percent
		integer(kind=C_INT)  	  ::	irow
		integer(kind=C_INT)  	  ::	icol
		integer(kind=C_INT)  	  ::	file_index
                        
    integer(kind=C_INT)  	  ::  IFGAS      

		type(data_wl)           ::  meas(KWM)
    real(kind=C_FLOAT)      ::  HVP(KVERTM)  ! height for vertical profile     

	end type pixel 

!	------------------------------------------------------------------------------------------------------
!	------------------------------------------------------------------------------------------------------
	type, bind(C) :: segment_data 

      integer(kind=C_INT)  	  ::	npixels ! number of pixels in meas data segment
! NX*NY*NT - meas data segment dimension
      integer(kind=C_INT)  	  ::	NT  ! segment time dimension
      integer(kind=C_INT)  	  ::	NX  ! cell lon dimension (number of pixels)
      integer(kind=C_INT)  	  ::	NY  ! cell lat dimension (number of pixels)
      integer(kind=C_INT)     ::  id_inversion  !  identification of current segment inside the tile
      type(pixel)             ::  pixels(KIMAGE)

	end type segment_data 

!	------------------------------------------------------------------------------------------------------
!	------------------------------------------------------------------------------------------------------
	type :: pixel_vector 
      !real                      ::	FS(KNBVM*KWM*KIP)
      real                      ::	FS(KMESS)
      integer                   ::  nFS(KWM)
      integer                   ::  KMIMAGE
	end type pixel_vector 

!	------------------------------------------------------------------------------------------------------
!	------------------------------------------------------------------------------------------------------

      contains
      
      subroutine dummy_subroutine_for_avoiding_warnings_in_an_empty_module()
      end subroutine dummy_subroutine_for_avoiding_warnings_in_an_empty_module                

end module mod_sdata_derived_type

! mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
