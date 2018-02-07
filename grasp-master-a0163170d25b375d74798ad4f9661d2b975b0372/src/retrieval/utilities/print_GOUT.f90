! **
! **  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.  
! **  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
! **

#include "../constants_set/mod_globals.inc"
!	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
module mod_print_array
  use mod_stop_report

contains

   subroutine iprint_array_int(iu, array, element_format, from_element, to_element)

     integer,      intent(in) :: iu
     integer,      intent(in) :: array(:)
     character(*), intent(in) :: element_format
     integer,      intent(in) :: from_element
     integer,      intent(in) :: to_element
     integer :: i
     integer :: nelements
     
     nelements = to_element - from_element + 1
     if(nelements .gt. size(array)) then
        write(tmp_message,'(a,i0,3x,a,i8)') &
        'nelements = ',nelements,'size(array) = ',size(array)
        G_ERROR(trim(tmp_message))
     endif

     do i = from_element, to_element
        write(iu,fmt=element_format, advance='no') array(i)
     end do
     write(iu,'("")')
     
   end subroutine iprint_array_int

   subroutine rprint_array_real(iu, array, element_format, from_element, to_element)

     integer,      intent(in) :: iu
     real,         intent(in) :: array(:)
     character(*), intent(in) :: element_format
     integer,      intent(in) :: from_element
     integer,      intent(in) :: to_element
     integer :: i
     integer :: nelements
     
     nelements = to_element - from_element + 1
     if(nelements .gt. size(array)) then
        write(tmp_message,'(a,i0,3x,a,i8)') &
        'nelements = ',nelements,'size(array)=',size(array)
        G_ERROR(trim(tmp_message))
     endif

     do i = from_element, to_element
        write(iu,fmt=element_format, advance='no') array(i)
     end do
     write(iu,'("")')
     
   end subroutine rprint_array_real

end module mod_print_array

!	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm

! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      logical function print_status (iu, RIN_flag, GOUT_flag, product_name)

!      use mod_retr_settings_derived_type
      use iso_c_binding

      implicit none
! ----------------------------------------------------------------	  
      integer, intent(in) :: iu
      logical(kind=C_BOOL), intent(in) :: RIN_flag
      logical(kind=C_BOOL), intent(in) :: GOUT_flag
      character(*), intent(in) :: product_name
! ----------------------------------------------------------------
      character(len=50) :: a1, a2
! ----------------------------------------------------------------	 
      a1 = '!!! WARNING!!! products.'
      a2 = ' is not available to be printed'

      print_status = .false.
      if(RIN_flag .eqv. .true.) then
        if(GOUT_flag .eqv. .true.) then
          print_status = .true.
        else
          write(iu,'(a)') trim(a1)//trim(product_name)//trim(a2)
        endif
      endif

      return
      end function print_status

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      subroutine set_GOUT_print_flags (iu, RIN_products, GOUT_products, PRN_products)

      use mod_retr_settings_derived_type
            
      implicit none
! ----------------------------------------------------------------	  
      integer, intent(in) ::  iu
      type(output_segment_products), intent(in)  ::  RIN_products
      type(output_segment_products), intent(in)  ::  GOUT_products
      type(output_segment_products), intent(out) ::  PRN_products
! ----------------------------------------------------------------
      logical :: print_status
      character(len=50) :: product_name
! ----------------------------------------------------------------
! retrieval
      ! final inversion residual
        product_name = 'retrieval.res'
      PRN_products%retrieval%res = &
      print_status(iu, RIN_products%retrieval%res, GOUT_products%retrieval%res, product_name)
      ! retrieved parameters
        product_name = 'retrieval.par'
      PRN_products%retrieval%par = &
      print_status(iu, RIN_products%retrieval%par, GOUT_products%retrieval%par, product_name)
      ! real and fitted measurements
        product_name = 'retrieval.fit'
      PRN_products%retrieval%fit = &
      print_status(iu, RIN_products%retrieval%fit, GOUT_products%retrieval%fit, product_name)

! aerosol
      ! aerosol optical depth and absorption optical depth
        product_name = 'aerosol.opt'
      PRN_products%aerosol%opt = &
      print_status(iu, RIN_products%aerosol%opt, GOUT_products%aerosol%opt, product_name)
      ! complex refractive index
        product_name = 'aerosol.rind'
      PRN_products%aerosol%rind = &
      print_status(iu, RIN_products%aerosol%rind, GOUT_products%aerosol%rind, product_name)
      ! particle chemical fractions if retrieved
        product_name = 'aerosol.chem'
      PRN_products%aerosol%chem = &
      print_status(iu, RIN_products%aerosol%chem, GOUT_products%aerosol%chem, product_name)
      ! lognormal size distribution parameters for theoretical fine and coarse aerosol particles
        product_name = 'aerosol.sd2m_mph'
      PRN_products%aerosol%sd2m_mph = &
      print_status(iu, RIN_products%aerosol%sd2m_mph, GOUT_products%aerosol%sd2m_mph, product_name)
      ! aod for theoretical fine and coarse aerosol particles
        product_name = 'aerosol.sd2m_ext'
      PRN_products%aerosol%sd2m_ext = &
      print_status(iu, RIN_products%aerosol%sd2m_ext, GOUT_products%aerosol%sd2m_ext, product_name)
      ! elements of scattering matrix
        product_name = 'aerosol.phmx'
      PRN_products%aerosol%phmx = &
      print_status(iu, RIN_products%aerosol%phmx, GOUT_products%aerosol%phmx, product_name)
      ! lidar ratio
        product_name = 'aerosol.lidar'
      PRN_products%aerosol%lidar = &
      print_status(iu, RIN_products%aerosol%lidar, GOUT_products%aerosol%lidar, product_name)
      ! particulate matter
        product_name = 'aerosol.pm'
      PRN_products%aerosol%pm = &
      print_status(iu, RIN_products%aerosol%pm, GOUT_products%aerosol%pm, product_name)
      ! aerosol type
        product_name = 'aerosol.types'
      PRN_products%aerosol%types = &
      print_status(iu, RIN_products%aerosol%types, GOUT_products%aerosol%types, product_name)

! clouds; future development
      ! cloud optical depth and absorption optical depth
        product_name = 'clouds.opt'
      PRN_products%clouds%opt = &
      print_status(iu, RIN_products%clouds%opt, GOUT_products%clouds%opt, product_name)
      ! complex refractive index
        product_name = 'clouds.rind'
      PRN_products%clouds%rind = &
      print_status(iu, RIN_products%clouds%rind, GOUT_products%clouds%rind, product_name)
      ! particle chemical fractions if retrieved
        product_name = 'clouds.chem'
      PRN_products%clouds%chem = &
      print_status(iu, RIN_products%clouds%chem, GOUT_products%clouds%chem, product_name)
      ! lognormal size distribution parameters for theoretical fine and coarse aerosol particles
        product_name = 'clouds.sd2m_mph'
      PRN_products%clouds%sd2m_mph = &
      print_status(iu, RIN_products%clouds%sd2m_mph, GOUT_products%clouds%sd2m_mph, product_name)
      ! aod for theoretical fine and coarse aerosol particles
        product_name = 'clouds.sd2m_ext'
      PRN_products%clouds%sd2m_ext = &
      print_status(iu, RIN_products%clouds%sd2m_ext, GOUT_products%clouds%sd2m_ext, product_name)
      ! elements of scattering matrix
        product_name = 'clouds.phmx'
      PRN_products%clouds%phmx = &
      print_status(iu, RIN_products%clouds%phmx, GOUT_products%clouds%phmx, product_name)
      ! lidar ratio
        product_name = 'clouds.lidar'
      PRN_products%clouds%lidar = &
      print_status(iu, RIN_products%clouds%lidar, GOUT_products%clouds%lidar, product_name)
      ! particulate matter
        product_name = 'clouds.pm'
      PRN_products%clouds%pm = &
      print_status(iu, RIN_products%clouds%pm, GOUT_products%clouds%pm, product_name)
      ! aerosol type
        product_name = 'clouds.types'
      PRN_products%clouds%types = &
      print_status(iu, RIN_products%clouds%types, GOUT_products%clouds%types, product_name)

! surface
        product_name = 'surface.surf'
      PRN_products%surface%surf = &
      print_status(iu, RIN_products%surface%surf, GOUT_products%surface%surf, product_name)

! errest
      ! for retrived parameters
        product_name = 'errest.par'
      PRN_products%errest%par = &
      print_status(iu, RIN_products%errest%par, GOUT_products%errest%par, product_name)
    ! aerosol
      ! for aerosol optical properties
        product_name = 'errest.aerosol.opt'
      PRN_products%errest%aerosol%opt = &
      print_status(iu, RIN_products%errest%aerosol%opt, GOUT_products%errest%aerosol%opt, product_name)
      ! for aerosol lidar ratio
        product_name = 'errest.aerosol.lidar'
      PRN_products%errest%aerosol%lidar = &
      print_status(iu, RIN_products%errest%aerosol%lidar, GOUT_products%errest%aerosol%lidar, product_name)
    ! clouds future development
      ! for cloud optical properties
        product_name = 'errest.clouds.opt'
      PRN_products%errest%clouds%opt = &
      print_status(iu, RIN_products%errest%clouds%opt, GOUT_products%errest%clouds%opt, product_name)
      ! for cloud lidar ratio
        product_name = 'errest.clouds.lidar'
      PRN_products%errest%clouds%lidar = &
      print_status(iu, RIN_products%errest%clouds%lidar, GOUT_products%errest%clouds%lidar, product_name)

! forcing
        product_name = 'forcing.bbflux'
      PRN_products%forcing%bbflux = &
      print_status(iu, RIN_products%forcing%bbflux, GOUT_products%forcing%bbflux, product_name)
        product_name = 'forcing.forcing'
      PRN_products%forcing%forcing = &
      print_status(iu, RIN_products%forcing%forcing, GOUT_products%forcing%forcing, product_name)

      return
      end subroutine set_GOUT_print_flags

! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! Print Output results
      
      subroutine print_output_results( iu_main_output, RIN, segment_meas, GOUT )

      use mod_retr_settings_derived_type
      use mod_retr_general_output_derived_type
      use mod_sdata_derived_type
      use mod_stop_report

      implicit none  
!----------------------------------------------------------------------------------------
! IN :
      integer, intent(in) :: iu_main_output
      type(retr_input_settings), intent(in) :: RIN
      type(segment_data), intent(inout) :: segment_meas
      type(output_segment_general), intent(in) :: GOUT
!----------------------------------------------------------------------------------------
      type(output_segment_products):: PRN_products
!----------------------------------------------------------------------------------------
! Set print GOUT flags
      call set_GOUT_print_flags ( iu_main_output, RIN%products, GOUT%products, PRN_products )

! Print detailed residual for final iteration
      if(PRN_products%retrieval%res) &
      call print_final_iteration_residuals ( iu_main_output, RIN, segment_meas, GOUT )

! Print final retrieval output
      call print_main_output ( iu_main_output, RIN, segment_meas, GOUT )

! Print general output
      ! Theoretical 2 modes properties
      if(PRN_products%aerosol%sd2m_mph .or. PRN_products%aerosol%sd2m_ext) &
      call print_sd2m ( iu_main_output, RIN, segment_meas, GOUT )
      if ( stop_report%status ) return

      ! Print particulate matter (PM)
      if(PRN_products%aerosol%PM) &
      call print_aerosol_PM ( iu_main_output, RIN, segment_meas, GOUT )

      ! Print aerosol types
      if(PRN_products%aerosol%types) &
      call print_aerosol_types ( iu_main_output, RIN, segment_meas, GOUT )

      ! Print phase matrix
      if(PRN_products%aerosol%phmx) &
      call print_phmx ( iu_main_output, RIN, segment_meas, GOUT )
      !call print_phmx_file ( RIN, segment_meas, GOUT )

      ! Print fitting (FS and FPS (measurements and modeled measurements calculated for retrieved parameters))
      if(PRN_products%retrieval%fit) &
      call print_fitting ( iu_main_output, RIN, segment_meas, GOUT%retrieval%fit%segment_fit )

      ! Print error estimates
      if(PRN_products%errest%par .or. PRN_products%errest%aerosol%opt .or. PRN_products%errest%aerosol%lidar) &
      call print_error_estimates ( iu_main_output, RIN, segment_meas, GOUT )
      if ( stop_report%status ) return

      ! Print broadband flux
      if(PRN_products%forcing%bbflux) &
      call print_forcing_bbflux ( iu_main_output, RIN, segment_meas, GOUT )

      ! Print net forcing
      if(PRN_products%forcing%forcing) &
      call print_forcing_forcing ( iu_main_output, RIN, segment_meas, GOUT )

      return
      end subroutine print_output_results 

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine print_main_output (  iu,RIN,           & 
                                      segment_meas,     &
                                      GOUT              &
                                   )

      use mod_par_inv,     only : KIDIM2,KIDIM3,KIMAGE,KVERTM
      use mod_par_DLS_bin, only : NRR,NRC
      use mod_par_DLS,     only : KNpar
      use mod_par_OS,      only : KSD	  
      use mod_retr_settings_derived_type
      use mod_retr_general_output_derived_type
      use mod_time_utils
      use mod_sdata_derived_type
      use mod_stop_report

      implicit none  
!----------------------------------------------------------------------------------------
! IN :
      integer,                     intent(in)  ::  iu
      type(retr_input_settings),   intent(in)  ::  RIN
      type(segment_data),          intent(in)  ::  segment_meas
      type(output_segment_general),intent(in)  ::  GOUT
!----------------------------------------------------------------------------------------
! LOCAL :
      integer                              :: i,i2,KNSING,IDIM1,IDIM2,NDIM3
      real,dimension(KIDIM3,KIDIM2,KIMAGE) :: SD,CHEM,RRE,RIM,sph,h01,CV,BRF,BRP,CL
!----------------------------------------------------------------------------------------
      integer                              :: ISD,ILR,IW,ipix,npixels
      real,dimension(NRR,KIDIM2,KIMAGE)    :: SDout 
      character(LEN=20)                 	 :: NBINC
      character(LEN=12),dimension(KIMAGE)  :: pdate,ptime
      real                                 :: DLNR,RSD,HGR_km
      integer                              :: par_type,par_type_SD,par_type_CV1
      integer                              :: IDIM1_CV,IDIM1_SD
      integer                              :: NHVP_meas,NHVP_retr
      real,dimension(KVERTM)               :: HVP_meas,HVP_retr_km
      logical                              :: lclouds = .false.
      type(output_segment_products)        :: PRN_products
!----------------------------------------------------------------------------------------
! Set print GOUT flags
      call set_GOUT_print_flags(iu, RIN%products, GOUT%products, PRN_products)

      KNSING = RIN%KNSING
      npixels=segment_meas%npixels

! Print output for clouds
      if(PRN_products%clouds%opt) lclouds = .true.
      if(PRN_products%clouds%rind) lclouds = .true.
      if(PRN_products%clouds%chem) lclouds = .true.
      if(PRN_products%clouds%sd2m_mph) lclouds = .true.
      if(PRN_products%clouds%sd2m_ext) lclouds = .true.
      if(PRN_products%clouds%phmx) lclouds = .true.
      if(PRN_products%clouds%lidar) lclouds = .true.
      if(PRN_products%clouds%pm) lclouds = .true.
      if(PRN_products%clouds%types) lclouds = .true.
      if(lclouds) then
        write(iu,'(a,/,a)') 'There is no cloud product yet to be printed.', &
        'Cloud properties retrieval is under development. '
      endif

! Print a table of retrieved parameters
      call print_date_time(iu,segment_meas)
      call print_coordinates(iu,segment_meas)      
      if(PRN_products%retrieval%par) &
      call print_retrieved_parameters(iu,RIN,segment_meas,GOUT)

! Print detailed parameters            
      write(iu,'(/,a)') '*** DETAILED PARAMETERS ***'
      call print_date_time(iu,segment_meas)
      call print_coordinates(iu,segment_meas)

      ! Size Distribution
      call print_size_distribution(iu,RIN,segment_meas,GOUT)
      if ( stop_report%status ) return

      ! Aerosol Concentration 
      call print_aerosol_concentration(iu,RIN,segment_meas,GOUT)

      ! Sphere fraction or Shape distribution
      call print_shape_distribution(iu,RIN,segment_meas,GOUT)

      ! Aerosol Vertical Profile
      call print_aerosol_profile(iu,RIN,segment_meas,GOUT) 

      ! Calibration coefficient
      call print_lidar_calibration_coeffitient(iu,RIN,segment_meas,GOUT)

      ! Angstrom exponent
      call print_angstrom_exponent(iu,RIN,segment_meas,GOUT)

      ! Aerosol optical depth
      if(PRN_products%aerosol%opt) &
      call print_optical_thickness(iu,RIN,segment_meas,GOUT)

      ! Single Scattering Albedo
      if(PRN_products%aerosol%opt) &
      call print_single_scattering_albedo(iu,RIN,segment_meas,GOUT)

      ! Absorption optical depth
      if(PRN_products%aerosol%opt) &
      call print_absorption(iu,RIN,segment_meas,GOUT)

      ! Real part of refr.index n(wl)
      if(PRN_products%aerosol%rind) &
      call print_refractive_index_real(iu,RIN,segment_meas,GOUT)

      ! Imaginary part of refr.index k(wl)
      if(PRN_products%aerosol%rind) &
      call print_refractive_index_imaginary(iu,RIN,segment_meas,GOUT)

      ! Particle chemical component fractions
      if(PRN_products%aerosol%chem) &
      call print_chemistry(iu,RIN,segment_meas,GOUT)

      ! Lidar Ratio  	  
      if(PRN_products%aerosol%lidar) &
      call print_lidar_ratio(iu,RIN,segment_meas,GOUT)

      ! Land percent
      if(PRN_products%surface%surf) &
      call print_surface_land_percent(iu,RIN,segment_meas,GOUT)

      ! Surface ndvi
      if(PRN_products%surface%surf) &
      call print_surface_ndvi(iu,RIN,segment_meas,GOUT)

      ! Surface albedo 
      if(PRN_products%surface%surf) &
      call print_surface_albedo(iu,RIN,segment_meas,GOUT)

      ! Surface parameters 
      if(PRN_products%surface%surf) &
      call print_surface_parameters(iu,RIN,segment_meas,GOUT)

      flush(iu)
      return
      end subroutine print_main_output
!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! Unpack retrieved parameters from retrieved parameter vector

      subroutine unpack_retr_param_to_print ( RIN,GOUT_retrieval_par,IDIM1,npixels,par )
      
      use mod_par_inv, only : KIDIM3,KIDIM2,KIMAGE
      use mod_retr_settings_derived_type 
      use mod_retr_general_output_derived_type
                 
      implicit none
! ----------------------------------------------------------------	  
      integer,                      intent(in)  ::  IDIM1,npixels
      type(retr_input_settings),    intent(in)  ::  RIN
      type(output_segment_retr_par),intent(in)  ::  GOUT_retrieval_par
      real,dimension(KIDIM3,KIDIM2,KIMAGE),intent(out) ::	par
! ----------------------------------------------------------------	  
      integer  :: NDIM3,IDIM2,ipix,i1,i2
! ----------------------------------------------------------------	  
      par(:,:,:)=0.0
      do ipix=1,npixels
        do IDIM2=1,RIN%NDIM%n2(IDIM1) 
            NDIM3=RIN%NDIM%n3(IDIM2,IDIM1)
            i1=RIN%NDIM%ISTARSING(IDIM2,IDIM1)
            i2=i1+NDIM3-1
            par(1:NDIM3,IDIM2,ipix) = GOUT_retrieval_par%pixel(ipix)%par(i1:i2)
          enddo ! IDIM2
        enddo ! ipix
          
      return
      end subroutine unpack_retr_param_to_print 

!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! OUTPUT print: date&time

      subroutine print_date_time ( iu,segment_meas )
      
      use mod_par_inv, only  : KIMAGE
      use mod_sdata_derived_type
            
      implicit none
! ----------------------------------------------------------------	  
      integer,             intent(in)  ::  iu
      type(segment_data),  intent(in)  ::  segment_meas
! ----------------------------------------------------------------	  
      character(len=12),dimension(KIMAGE)  ::  pdate,ptime
      integer  :: ipix,npixels
! ----------------------------------------------------------------	  
      npixels = segment_meas%npixels 
      call yyyymmdd_hhmmss ( segment_meas,pdate,ptime )      
      write(iu,'(a5,10x,20000(3x,A10,x))')  "Date:",(pdate(ipix),ipix=1,npixels)
      write(iu,'(a5,12x,20000(3x,A8,3x))')  "Time:",(ptime(ipix),ipix=1,npixels)
      
      return
      end subroutine print_date_time 
      
! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! OUTPUT print: coordinates

      subroutine print_coordinates ( iu,segment_meas )
      
      use mod_sdata_derived_type
      use mod_print_array  
          
      implicit none
! ----------------------------------------------------------------	  
      integer,                       intent(in)  ::  iu
      type(segment_data),            intent(in)  ::  segment_meas
! ----------------------------------------------------------------	  
      integer  :: npixels
! ----------------------------------------------------------------	  
      npixels = segment_meas%npixels 
      write(iu,'(a10,4x)', advance='no') "Longitude:"
      call rprint_array_real(iu,segment_meas%pixels(:)%x,"(f14.4)", 1, npixels)
      write(iu,'(a10,4x)', advance='no') "Latitude :"
      call rprint_array_real(iu,segment_meas%pixels(:)%y,"(f14.4)", 1, npixels)
      
      return
      end subroutine print_coordinates 
      
! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! OUTPUT print: sd2m structure print

      subroutine print_sd2m ( iu,RIN,segment_meas,GOUT )
      
      use mod_retr_settings_derived_type 
      use mod_retr_general_output_derived_type
      use mod_sdata_derived_type
      use mod_stop_report

      implicit none
! ----------------------------------------------------------------	  
      integer,                       intent(in)  ::  iu
      type(retr_input_settings),     intent(in)  ::  RIN
      type(segment_data),            intent(in)  ::  segment_meas
      type(output_segment_general),  intent(in)  ::  GOUT
! ----------------------------------------------------------------	  
      integer  :: ipix,npixels
      integer,parameter  :: par_total  = 0
      integer,parameter  :: par_fine   = 1
      integer,parameter  :: par_coarse = 2      
! ----------------------------------------------------------------	  
      npixels = segment_meas%npixels 
! ***  Date and time for npixels           *************************
      if(RIN%products%aerosol%sd2m_mph .or. RIN%products%aerosol%sd2m_ext) then
        write(iu,*)
        call print_date_time(iu,segment_meas)
        call print_coordinates(iu,segment_meas)
      endif
      if(RIN%products%aerosol%sd2m_mph) then
        if(RIN%NSD .eq. 1)  & 
        call print_sd2m_mph(iu,RIN,segment_meas,GOUT,par_total )      
        call print_sd2m_mph(iu,RIN,segment_meas,GOUT,par_fine  )      
        call print_sd2m_mph(iu,RIN,segment_meas,GOUT,par_coarse) 
      endif
      if(RIN%products%aerosol%sd2m_ext .and. RIN%NSD .eq. 1) then
        call print_sd2m_ext(iu,RIN,segment_meas,GOUT,par_fine  )       
        call print_sd2m_ext(iu,RIN,segment_meas,GOUT,par_coarse)       
        if ( stop_report%status ) return
      endif
      
      return
      end subroutine print_sd2m 
      
! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine print_sd2m_mph (iu,RIN,segment_meas,GOUT,index)

      use mod_retr_general_output_derived_type
      use mod_print_array
      use mod_sdata_derived_type
            
      implicit none
! ----------------------------------------------------------------	  
      integer,                     intent(in)  ::  iu,index
      type(retr_input_settings),   intent(in)  ::  RIN
      type(segment_data),          intent(in)  ::  segment_meas
      type(output_segment_general),intent(in)  ::  GOUT
! ----------------------------------------------------------------	  
      integer                 ::  i1,i2,ipix,npixels
      character(len=6)        ::  aaa      
! ----------------------------------------------------------------	 
      npixels = segment_meas%npixels

      select case(index)
      case(0) 
        aaa = 'total'
      case(1) 
        aaa = 'fine'
      case(2) 
        aaa = 'coarse'
      end select 
      
      if(RIN%NSD .eq. 1) then
        write(iu,'(a14)', advance='no') 'cv   '//aaa
        call rprint_array_real(iu,GOUT%aerosol%sd2m%mph%pixel(:)%cv(index),"(e14.5)", 1, npixels)
        write(iu,'(a14)', advance='no') 'rv   '//aaa
        call rprint_array_real(iu,GOUT%aerosol%sd2m%mph%pixel(:)%rm(index),"(e14.5)", 1, npixels)
        write(iu,'(a14)', advance='no') 'std  '//aaa
        call rprint_array_real(iu,GOUT%aerosol%sd2m%mph%pixel(:)%std(index),"(e14.5)", 1, npixels)
      endif ! RIN%NSD .eq. 1
        write(iu,'(a14)', advance='no') 'reff '//aaa
        call rprint_array_real(iu,GOUT%aerosol%sd2m%mph%pixel(:)%reff(index),"(e14.5)", 1, npixels)
      return
      end subroutine print_sd2m_mph 

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      
      subroutine print_sd2m_ext (iu,RIN,segment_meas,GOUT,index)

      use mod_retr_general_output_derived_type
      use mod_retr_settings_derived_type
      use mod_print_array
      use mod_sdata_derived_type
      use mod_stop_report

      implicit none
! ----------------------------------------------------------------	  
      integer,                     intent(in)  ::  iu,index
      type(retr_input_settings),   intent(in)  ::  RIN
      type(segment_data),          intent(in)  ::  segment_meas
      type(output_segment_general),intent(in)  ::  GOUT
! ----------------------------------------------------------------	  
      integer           ::  i1,i2,ipix,iw,iv,npixels
      character(len=6)  ::  aaa
! ----------------------------------------------------------------	 
      npixels = segment_meas%npixels
      select case(index)
      case(1) 
        aaa = 'fine'
      case(2) 
        aaa = 'coarse'
      case default
        write(tmp_message,'(a,i0,a)') 'index = ',index,'  - unknown'
        G_ERROR(trim(tmp_message))
      end select
      
      write(iu,'(a)') 'Wavelength (um),  aod_'//aaa
      do iw=1,RIN%NW
        write(iu,'(f14.5)', advance='no') RIN%WAVE(iw)
        call rprint_array_real(iu,GOUT%aerosol%sd2m%opt%pixel(:)%ext(iw,index),"(e14.5)", 1, npixels)
      enddo ! iw
      
      return
      end subroutine print_sd2m_ext 

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! Print detailed residual for final iteration

      subroutine print_final_iteration_residuals (iu, RIN, segment_meas, GOUT)

      use mod_retr_general_output_derived_type
      use mod_print_array
      use mod_sdata_derived_type
            
      implicit none  
!----------------------------------------------------------------------------------------
      integer,                      intent(in)  ::  iu
      type(retr_input_settings),    intent(in)  ::  RIN
      type(output_segment_general), intent(in)  ::  GOUT
      type(segment_data),           intent(in)  ::  segment_meas
!----------------------------------------------------------------------------------------
      integer  ::  i, ipix, npixels, INOISE
      character(LEN=20) :: KNC, INOISEC
      character(LEN=150) :: CFMT
!----------------------------------------------------------------------------------------
      npixels  = segment_meas%npixels
      INOISE = RIN%NOISE%INOISE
      write(INOISEC,*) INOISE

      write(iu,'(/,a)') 'Detailed residuals after final iteration'
      CFMT = '(8x,'//trim(adjustl(INOISEC))//'(4x,a))'
      write(iu,trim(CFMT)) ('noise          abs            rel', i=1,INOISE)

      if(RIN%IPFP%INVSING .ge. 2) then

        loop_pixel1: do ipix=1,npixels
            CFMT = '(10x,'//trim(adjustl(INOISEC))//'(i5,a,e14.5,f15.5,a),2(6x,a,i0))'
            write(iu,trim(CFMT)) &
            (i,':',GOUT%retrieval%res%pixel(ipix)%resa(i), &
            GOUT%retrieval%res%pixel(ipix)%resr(i)*100.0,' %', &
            i=1,INOISE),'pixel # ', ipix,  &
            'Residual after iteration # ',GOUT%retrieval%res%niter
        enddo loop_pixel1
        CFMT = '(/,f10.5,'//trim(adjustl(INOISEC))//'(i5,a,e14.5,f15.5,a),6x,a,i0,a/)'
        write(iu,trim(CFMT)) GOUT%retrieval%res%rest*100.0, &
        (i,':',GOUT%retrieval%res%resat(i), &
        GOUT%retrieval%res%resrt(i)*100.0,' %', &
        i=1,INOISE),'Residual after iteration # ', GOUT%retrieval%res%niter,'  for TOTAL SEGMENT'

      elseif(RIN%IPFP%INVSING .eq. 0 .or. RIN%IPFP%INVSING .eq. 1) then

        loop_pixel2: do ipix=1,npixels
            CFMT = '(f10.5,'//trim(adjustl(INOISEC))//'(i5,a,e14.5,f15.5,a),2(6x,a,i0))'
            write(iu,trim(CFMT)) GOUT%retrieval%res%pixel(ipix)%res*100.0, &
            (i,':',GOUT%retrieval%res%pixel(ipix)%resa(i), &
            GOUT%retrieval%res%pixel(ipix)%resr(i)*100.0,' %', &
            i=1,INOISE),'pixel # ', ipix, &
            'Residual after iteration # ', GOUT%retrieval%res%pixel(ipix)%niter
        enddo loop_pixel2

      endif ! INVSING .ge. 2
      write(iu,*)

      return
      end subroutine print_final_iteration_residuals

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! Print vector of Retreived parameters for segment  	  

      subroutine print_retrieved_parameters (iu,RIN,segment_meas,GOUT)

      use mod_retr_general_output_derived_type
      use mod_print_array
      use mod_sdata_derived_type
            
      implicit none  
!----------------------------------------------------------------------------------------
      integer,                      intent(in)  ::  iu
      type(retr_input_settings),    intent(in)  ::  RIN
      type(output_segment_general), intent(in)  ::  GOUT
      type(segment_data),           intent(in)  ::  segment_meas
      integer  ::  i,KNSING,npixels
!----------------------------------------------------------------------------------------
      npixels  = segment_meas%npixels
      KNSING   = RIN%KNSING

      !write(iu,*)
      write(iu,'(a)') 'Parameter #, Vector of retrieved parameters'
      do i=1,KNSING
        write(iu,'(i14)', advance='no') i
        call rprint_array_real (iu,GOUT%retrieval%par%pixel(:)%par(i), '(e14.5)', 1, npixels)
      enddo ! i
            
      return
      end subroutine print_retrieved_parameters 

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! Print Size Distribution   	  

      subroutine print_size_distribution (iu,RIN,segment_meas,GOUT)

      use mod_retr_settings_derived_type
      use mod_retr_general_output_derived_type
      use mod_print_array
      use mod_sdata_meas_type
      use mod_sdata_derived_type
      use mod_stop_report

      implicit none  
!----------------------------------------------------------------------------------------
! IN :
      integer,                      intent(in)  ::  iu
      type(retr_input_settings),    intent(in)  ::  RIN
      type(output_segment_general), intent(in)  ::  GOUT
      type(segment_data),           intent(in)  ::  segment_meas
!----------------------------------------------------------------------------------------
! LOCAL :
      integer          ::  IDIM1,par_type
!----------------------------------------------------------------------------------------

      do IDIM1=1,RIN%NDIM%n1
        par_type = RIN%NDIM%par_type(IDIM1)        
        if(par_type .gt. par_type_SD_beg .and. par_type .lt. par_type_SD_end) then
          select case(par_type)
          case(par_type_SD_TB)
            call print_size_distribution_TB (iu,RIN,segment_meas,GOUT)
          case(par_type_SD_LB)
            call print_size_distribution_LB (iu,RIN,segment_meas,GOUT)
          case(par_type_SD_LN)
            call print_size_distribution_LN (iu,RIN,segment_meas,GOUT)
            if ( stop_report%status ) return
          end select
        exit
        endif ! par_type .gt. par_type_SD_beg .and.
      enddo ! IDIM1

      return
      end subroutine print_size_distribution 

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! Print Triangle Bin (TB) Size Distribution   	  

      subroutine print_size_distribution_TB (iu,RIN,segment_meas,GOUT)

      use mod_par_inv,     only : KIDIM2,KIDIM3,KIMAGE
      use mod_retr_settings_derived_type
      use mod_retr_general_output_derived_type
      use mod_print_array
      use mod_sdata_meas_type
      use mod_sdata_derived_type

      implicit none  
!----------------------------------------------------------------------------------------
! IN :
      integer,                      intent(in)  ::  iu
      type(retr_input_settings),    intent(in)  ::  RIN
      type(output_segment_general), intent(in)  ::  GOUT
      type(segment_data),           intent(in)  ::  segment_meas
!----------------------------------------------------------------------------------------
! LOCAL :
      integer                              :: IDIM1,IDIM2,NDIM3
      real,dimension(KIDIM3,KIDIM2,KIMAGE) :: SD,CV
      integer                              :: i,ipix,npixels
      real                                 :: DLNR
      real,dimension(KIMAGE)               :: temp
      integer                              :: par_type
      logical                              :: icv = .false.
!----------------------------------------------------------------------------------------
      temp(:) = 0.0
      npixels  = segment_meas%npixels
      
      !do i=1,RIN%NDIM%n1
      !  if(RIN%NDIM%par_type(i) .gt. par_type_Cv_beg .and. &
      !     RIN%NDIM%par_type(i) .lt. par_type_Cv_end) then
      !    ! if Aerosol Concentration retreived 
      !    icv = .true.
      !  exit
      !  endif ! RIN%NDIM%par_type(i) .gt. par_type_SD_beg .and. 
      !enddo ! i

      do IDIM1=1,RIN%NDIM%n1
        par_type = RIN%NDIM%par_type(IDIM1)        
        if(par_type .eq. par_type_SD_TB) then
          call unpack_retr_param_to_print ( RIN,GOUT%retrieval%par,IDIM1,npixels,SD )
          do IDIM2=1,RIN%NDIM%n2(IDIM1)    ! aerosol component loop
            NDIM3 = RIN%NDIM%n3(IDIM2,IDIM1)                   
            write(iu,'(a,1x,i0,a)') 'Size Distribution dV/dlnr (normalized to 1) for', &
                                     IDIM2,' - fraction'
            !if(icv) then
                DLNR=LOG(RIN%radius1(2,IDIM2))-LOG(RIN%radius1(1,IDIM2))
                do i=1,NDIM3
                write(iu,'(e14.5)',advance='no') RIN%radius1(i,IDIM2)
                do ipix=1,npixels
                temp(ipix) = SD(i,IDIM2,ipix)/(SUM(SD(1:NDIM3,IDIM2,ipix))*DLNR)
                enddo ! ipix
                call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
                enddo ! i
            !else
            !    do i=1,NDIM3
            !    write(iu,'(e14.5)',advance='no') RIN%radius1(i,IDIM2)
            !    call rprint_array_real (iu,SD(i,IDIM2,:), '(e14.5)', 1, npixels)
            !    enddo ! i                  
            !endif ! icv
          enddo  ! IDIM2
        exit
        endif ! par_type .eq. par_type_SD_TB
      enddo ! IDIM1

      return
      end subroutine print_size_distribution_TB 

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! Print precalculated Lognormal Bin (LB) Size Distribution   	  

      subroutine print_size_distribution_LB (iu,RIN,segment_meas,GOUT)

      use mod_par_inv,     only : KIDIM2,KIDIM3,KIMAGE
      use mod_par_DLS_bin, only : NRR
      use mod_retr_settings_derived_type
      use mod_retr_general_output_derived_type
      use mod_print_array
      use mod_sdata_derived_type

      implicit none  
!----------------------------------------------------------------------------------------
! IN :
      integer,                      intent(in)  ::  iu
      type(retr_input_settings),    intent(in)  ::  RIN
      type(output_segment_general), intent(in)  ::  GOUT
      type(segment_data),           intent(in)  ::  segment_meas
!----------------------------------------------------------------------------------------
! LOCAL :
      integer                              :: IDIN1,IDIM2,NDIM3
      real,dimension(KIDIM3,KIDIM2,KIMAGE) :: SD
      integer                              :: i,ipix,npixels
      real,dimension(NRR,KIDIM2,KIMAGE)    :: SDout 
      integer                              :: par_type,IDIM1
      real,dimension(KIDIM2,KIMAGE)        :: temp

!----------------------------------------------------------------------------------------
      npixels  = segment_meas%npixels
      
      do IDIM1=1,RIN%NDIM%n1
        par_type = RIN%NDIM%par_type(IDIM1)        
        if(par_type .eq. par_type_SD_LB) then
          call unpack_retr_param_to_print ( RIN,GOUT%retrieval%par,IDIM1,npixels,SD )
          ! precalculated lognormal bins
            call SD_precalcul_ln_bins (  RIN,GOUT%retrieval%par,npixels,IDIM1,SD,SDout )
            do ipix=1,npixels
              do IDIM2=1,RIN%NDIM%n2(IDIM1)
                NDIM3 = RIN%NDIM%n3(IDIM2,IDIM1)
                temp(IDIM2,ipix) = sum(SD(1:NDIM3,IDIM2,ipix))              
              enddo ! IDIM2
            enddo ! ipix
            do IDIM2=1,RIN%NDIM%n2(IDIM1)    ! aerosol component loop
              if(.not. RIN%use_models) then
                write(iu,'(a,1x,i0,a)') 'Radius (um), Size Distribution dV/dlnr (normalized to 1) for', &
                                        IDIM2,' - fraction'
                do i=1,GOUT%retrieval%par%ngrid
                  write(iu,'(e14.5)',advance='no') GOUT%retrieval%par%radius(i)
                  SDout(i,IDIM2,1:npixels) = SDout(i,IDIM2,1:npixels)/temp(IDIM2,1:npixels)
                  call rprint_array_real (iu,SDout(i,IDIM2,:), '(e14.5)', 1, npixels)
                enddo ! i
                write(iu,'(2a,1x,i0,a)') 'rv (um), Volume concentration (um^3/um^2 or um^3/um^3) ', &
                'of precomputed lognormal bins for',IDIM2,' - fraction'
              else
                write(iu,'(a,1x,i0,a)') '      #, Model fraction in total concentration for', &
                                                    IDIM2,' - fraction'
              endif

              NDIM3 = RIN%NDIM%n3(IDIM2,IDIM1)
              do i=1,NDIM3
                write(iu,'(e14.5)',advance='no') RIN%radius1(i,IDIM2)
                call rprint_array_real (iu,SD(i,IDIM2,:), '(e14.5)', 1, npixels)
              enddo ! i
            enddo  ! IDIM2
        exit
        endif ! par_type .eq. par_type_SD_LB
      enddo ! IDIM1

      return
      end subroutine print_size_distribution_LB 

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! Print Lognormal Bi-Modal Size Distribution   	  

      subroutine print_size_distribution_LN (iu,RIN,segment_meas,GOUT)

      use mod_par_inv,     only : KIDIM2,KIDIM3,KIMAGE
      use mod_retr_settings_derived_type
      use mod_retr_general_output_derived_type
      use mod_print_array
      use mod_sdata_derived_type
      use mod_stop_report

      implicit none  
!----------------------------------------------------------------------------------------
! IN :
      integer,                      intent(in)  ::  iu
      type(retr_input_settings),    intent(in)  ::  RIN
      type(output_segment_general), intent(in)  ::  GOUT
      type(segment_data),           intent(in)  ::  segment_meas
!----------------------------------------------------------------------------------------
! LOCAL :
      integer                              ::  IDIM1,IDIM2,NDIM3
      real,dimension(KIDIM3,KIDIM2,KIMAGE) ::  SD,CV
      integer                              ::  i,npixels
      real,dimension(KIDIM3,KIDIM2,KIMAGE) ::  SDout 
      real,dimension(KIMAGE)               ::  temp
      integer                              ::  par_type
      logical                              ::  icv = .false.
!----------------------------------------------------------------------------------------
      temp(:) = 0.0
      npixels = segment_meas%npixels
      
      do IDIM1=1,RIN%NDIM%n1
        par_type = RIN%NDIM%par_type(IDIM1)        
        if(par_type .gt. par_type_Cv_beg .and. par_type .lt. par_type_Cv_end) then
          ! if Aerosol Concentration retreived 
          call unpack_retr_param_to_print ( RIN,GOUT%retrieval%par,IDIM1,npixels,CV ) 
          icv = .true.
        exit
        endif ! RIN%NDIM%par_type(i) .gt. par_type_SD_beg .and. 
      enddo ! i
      if(.not. icv) then
        write(tmp_message,'(2a)') &
        'Can not print Bi-modal SD. Aerosol Concentration is not', &
        ' a parameter of forward model.'
        G_ERROR(trim(tmp_message))
      endif ! .not. icv

      do IDIM1=1,RIN%NDIM%n1
        par_type = RIN%NDIM%par_type(IDIM1)        
        if(par_type .gt. par_type_SD_beg .and. par_type .lt. par_type_SD_end) then
          call unpack_retr_param_to_print ( RIN,GOUT%retrieval%par,IDIM1,npixels,SD )
          call SD_lognormal ( RIN,GOUT%retrieval%par,npixels,IDIM1,SD,CV,SDout )
          do IDIM2=1,RIN%NDIM%n2(IDIM1)    ! aerosol component loop
              NDIM3 = RIN%NDIM%n3(IDIM2,IDIM1)    
              write(iu,'(a,1x,i0,a)') 'Radius (um), Size Distribution dV/dlnr (normalized to 1) for', &
                                    IDIM2,' - fraction'
              do i=1,GOUT%retrieval%par%ngrid
              write(iu,'(f14.5)',advance='no') GOUT%retrieval%par%radius(i)
              SDout(i,IDIM2,:) = SDout(i,IDIM2,:)/CV(1,IDIM2,:)
              call rprint_array_real (iu,SDout(i,IDIM2,:), '(e14.5)', 1, npixels)
              enddo ! i                                
              write(iu,'(a,1x,i0,a)') 'Parameters of lognormal SD for',IDIM2,' - fraction'
              do i=1,NDIM3
              if(i .eq. 1) write(iu,'(a14)',advance='no') 'rv (um):'
              if(i .eq. 2) write(iu,'(a14)',advance='no') 'ln(sigma):'
              call rprint_array_real (iu,SD(i,IDIM2,:), '(e14.5)', 1, npixels)
              enddo ! i
          enddo  ! IDIM2
        exit
        endif ! par_type .gt. par_type_SD_beg .and.
      enddo ! IDIM1

      return
      end subroutine print_size_distribution_LN

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! Print Aerosol concentrations             	  

      subroutine print_aerosol_concentration (iu,RIN,segment_meas,GOUT)

      use mod_par_inv,     only : KIDIM2,KIDIM3,KIMAGE
      use mod_retr_settings_derived_type
      use mod_retr_general_output_derived_type
      use mod_print_array
      use mod_sdata_derived_type
      
      implicit none  
!----------------------------------------------------------------------------------------
! IN :
      integer,                      intent(in)  ::  iu
      type(retr_input_settings),    intent(in)  ::  RIN
      type(output_segment_general), intent(in)  ::  GOUT
      type(segment_data),           intent(in)  ::  segment_meas
!----------------------------------------------------------------------------------------
! LOCAL :
      integer                              :: IDIM1,IDIM2,NDIM3
      real,dimension(KIDIM3,KIDIM2,KIMAGE) :: CV,SD
      integer                              :: i,ipix,npixels,par_type,par_type_SD
      real,dimension(KIMAGE)               :: temp
      logical                              :: icv = .false.
!----------------------------------------------------------------------------------------
      temp(:) = 0.0
      npixels = segment_meas%npixels
      
      icv = .false.
      do IDIM1=1,RIN%NDIM%n1
        par_type = RIN%NDIM%par_type(IDIM1)        
        if(par_type .gt. par_type_Cv_beg .and. par_type .lt. par_type_Cv_end) then
          icv = .true.
        exit
        endif ! par_type .gt. par_type_Cv_beg .and. 
      enddo ! IDIM1

      !write(iu,*) 

! volume of particles in air column per unit area um^3/um^2 or volume of particles per unit air volume um^3/um^3
      write(iu,'(a)')  'Aerosol volume concentration (um^3/um^2 or um^3/um^3)'

      if(icv) then
        call unpack_retr_param_to_print ( RIN,GOUT%retrieval%par,IDIM1,npixels,CV ) 
      else
        do IDIM1=1,RIN%NDIM%n1
          par_type = RIN%NDIM%par_type(IDIM1)              
          if(par_type .gt. par_type_SD_beg .and. par_type .lt. par_type_SD_end) then
            call unpack_retr_param_to_print ( RIN,GOUT%retrieval%par,IDIM1,npixels,SD )
            par_type_SD = par_type
          exit
          endif ! par_type .gt. par_type_Cv_beg .and. 
        enddo ! IDIM1
      endif ! icv

      do IDIM2=1,RIN%NDIM%n2(IDIM1)    ! aerosol component loop
        NDIM3 = RIN%NDIM%n3(IDIM2,IDIM1)    
        write(iu,'(i14)',advance='no') IDIM2
        if(icv) then
          temp(1:npixels) = CV(1,IDIM2,1:npixels)
        else
          do ipix=1,npixels
            temp(ipix) = sum(SD(1:NDIM3,IDIM2,ipix))
            if(par_type_SD .eq. par_type_SD_TB)  &
            temp(ipix) = temp(ipix)*( LOG(RIN%radius1(2,IDIM2))-LOG(RIN%radius1(1,IDIM2)) )
          enddo ! ipix
        endif ! icv        
        call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
      enddo ! IDIM2        

!!************add by lei on 20160914 to output aerosol chemical elements mass concentration *******
!
!      write(iu,'(a)')  'Aerosol chemical elements mass concentration (mg/m^2)'
!
!      do IDIM1=1,RIN%NDIM%n1
!        par_type = RIN%NDIM%par_type(IDIM1)        
!        if(par_type .gt. par_type_RERI_beg .and. par_type .lt. par_type_RERI_end) then
!
!          if(par_type .eq. par_type_CXRI_nmix) then
!
!      do IDIM2=1,RIN%NDIM%n2(IDIM1)       ! aerosol component loop for chemical elements mass concentration
!        NDIM3 = RIN%NDIM%n3(IDIM2,IDIM1)    
!          if (IDIM2 .eq. 1) then
!        write(iu,'(a)') 'Mass concentration of BC:'
!        write(iu,'(i14)',advance='no') IDIM2
!          do ipix=1,npixels
!            temp(ipix) = sum(SD(1:NDIM3,IDIM2,ipix))
!            temp(ipix) = GOUT%aerosol%chem%pixel(ipix)%fsoot(IDIM2)*1.8*1000*sum(SD(1:NDIM3,IDIM2,ipix))
!            if(par_type_SD .eq. par_type_SD_TB)  &
!            temp(ipix) = temp(ipix)*( LOG(RIN%radius1(2,IDIM2))-LOG(RIN%radius1(1,IDIM2)) )
!          enddo ! ipi
!        call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
!
!        write(iu,'(a)') 'Mass concentration of BrC:'
!        write(iu,'(i14)',advance='no') IDIM2
!          do ipix=1,npixels
!            temp(ipix) = sum(SD(1:NDIM3,IDIM2,ipix))
!            temp(ipix) = GOUT%aerosol%chem%pixel(ipix)%fslbl(IDIM2)*1.2*1000*sum(SD(1:NDIM3,IDIM2,ipix))
!            if(par_type_SD .eq. par_type_SD_TB)  &
!            temp(ipix) = temp(ipix)*( LOG(RIN%radius1(2,IDIM2))-LOG(RIN%radius1(1,IDIM2)) )
!          enddo ! ipi
!        call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
!
!        write(iu,'(a)') 'Mass concentration of Quartz:'
!        write(iu,'(i14)',advance='no') IDIM2
!          do ipix=1,npixels
!            temp(ipix) = sum(SD(1:NDIM3,IDIM2,ipix))
!            temp(ipix) = GOUT%aerosol%chem%pixel(ipix)%finslbl(IDIM2)*2.5*1000*sum(SD(1:NDIM3,IDIM2,ipix))
!            if(par_type_SD .eq. par_type_SD_TB)  &
!            temp(ipix) = temp(ipix)*( LOG(RIN%radius1(2,IDIM2))-LOG(RIN%radius1(1,IDIM2)) )
!          enddo ! ipi
!        call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
!
!        write(iu,'(a)') 'Mass concentration of Water:'
!        write(iu,'(i14)',advance='no') IDIM2
!          do ipix=1,npixels
!            temp(ipix) = sum(SD(1:NDIM3,IDIM2,ipix))
!            temp(ipix) = GOUT%aerosol%chem%pixel(ipix)%fwtr(IDIM2)*1*1000*sum(SD(1:NDIM3,IDIM2,ipix))
!            if(par_type_SD .eq. par_type_SD_TB)  &
!            temp(ipix) = temp(ipix)*( LOG(RIN%radius1(2,IDIM2))-LOG(RIN%radius1(1,IDIM2)) )
!          enddo ! ipi
!        call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
!
!         else
!!        write(iu,'(a)') 'Mass concentration of Carbonaceous:'
!!        write(iu,'(i14)',advance='no') IDIM2
!!          do ipix=1,npixels
!!            temp(ipix) = sum(SD(1:NDIM3,IDIM2,ipix))
!!            temp(ipix) = GOUT%aerosol%chem%pixel(ipix)%fcarbon(IDIM2)*1.2*1000*sum(SD(1:NDIM3,IDIM2,ipix))
!!            if(par_type_SD .eq. par_type_SD_TB)  &
!!            temp(ipix) = temp(ipix)*( LOG(RIN%radius1(2,IDIM2))-LOG(RIN%radius1(1,IDIM2)) )
!!          enddo ! ipi
!!        call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
!!        write(iu,'(a)') 'Mass concentration of BC:'
!!        write(iu,'(i14)',advance='no') IDIM2
!!          do ipix=1,npixels
!!            temp(ipix) = sum(SD(1:NDIM3,IDIM2,ipix))
!!            temp(ipix) = GOUT%aerosol%chem%pixel(ipix)%fsoot(IDIM2)*1.8*1000*sum(SD(1:NDIM3,IDIM2,ipix))
!!            if(par_type_SD .eq. par_type_SD_TB)  &
!!            temp(ipix) = temp(ipix)*( LOG(RIN%radius1(2,IDIM2))-LOG(RIN%radius1(1,IDIM2)) )
!!          enddo ! ipi
!!        call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
!!
!!        write(iu,'(a)') 'Mass concentration of BrC:'
!!        write(iu,'(i14)',advance='no') IDIM2
!!          do ipix=1,npixels
!!            temp(ipix) = sum(SD(1:NDIM3,IDIM2,ipix))
!!            temp(ipix) = GOUT%aerosol%chem%pixel(ipix)%fslbl(IDIM2)*1.2*1000*sum(SD(1:NDIM3,IDIM2,ipix))
!!            if(par_type_SD .eq. par_type_SD_TB)  &
!!            temp(ipix) = temp(ipix)*( LOG(RIN%radius1(2,IDIM2))-LOG(RIN%radius1(1,IDIM2)) )
!!          enddo ! ipi
!!        call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
!
!        write(iu,'(a)') 'Mass concentration of Quartz:'
!        write(iu,'(i14)',advance='no') IDIM2
!          do ipix=1,npixels
!            temp(ipix) = sum(SD(1:NDIM3,IDIM2,ipix))
!            temp(ipix) = GOUT%aerosol%chem%pixel(ipix)%finslbl(IDIM2)*2.5*1000*sum(SD(1:NDIM3,IDIM2,ipix))
!            if(par_type_SD .eq. par_type_SD_TB)  &
!            temp(ipix) = temp(ipix)*( LOG(RIN%radius1(2,IDIM2))-LOG(RIN%radius1(1,IDIM2)) )
!          enddo ! ipi
!        call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
!
!        write(iu,'(a)') 'Mass concentration of Iron:'
!        write(iu,'(i14)',advance='no') IDIM2
!          do ipix=1,npixels
!            temp(ipix) = sum(SD(1:NDIM3,IDIM2,ipix))
!            temp(ipix) = GOUT%aerosol%chem%pixel(ipix)%firon(IDIM2)*4.8*1000*sum(SD(1:NDIM3,IDIM2,ipix))
!            if(par_type_SD .eq. par_type_SD_TB)  &
!            temp(ipix) = temp(ipix)*( LOG(RIN%radius1(2,IDIM2))-LOG(RIN%radius1(1,IDIM2)) )
!          enddo ! ipi
!        call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
!
!        write(iu,'(a)') 'Mass concentration of Water:'
!        write(iu,'(i14)',advance='no') IDIM2
!          do ipix=1,npixels
!            temp(ipix) = sum(SD(1:NDIM3,IDIM2,ipix))
!            temp(ipix) = GOUT%aerosol%chem%pixel(ipix)%fwtr(IDIM2)*1*1000*sum(SD(1:NDIM3,IDIM2,ipix))
!            if(par_type_SD .eq. par_type_SD_TB)  &
!            temp(ipix) = temp(ipix)*( LOG(RIN%radius1(2,IDIM2))-LOG(RIN%radius1(1,IDIM2)) )
!          enddo ! ipi
!        call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
!
!           endif
!
!       enddo ! IDIM2
!
!!   elseif(par_type .eq. par_type_CXRI_chem) then
!!
!!      do IDIM2=1,RIN%NDIM%n2(IDIM1)       ! aerosol component loop for chemical elements mass concentration
!!        NDIM3 = RIN%NDIM%n3(IDIM2,IDIM1)    
!!     if (IDIM2 .eq. 1) then
!!        write(iu,'(a)') 'Mass concentration of BC:'
!!        write(iu,'(i14)',advance='no') IDIM2
!!          do ipix=1,npixels
!!            temp(ipix) = sum(SD(1:NDIM3,IDIM2,ipix))
!!            temp(ipix) = GOUT%aerosol%chem%pixel(ipix)%fsoot(IDIM2)*1.8*1000*sum(SD(1:NDIM3,IDIM2,ipix))
!!            if(par_type_SD .eq. par_type_SD_TB)  &
!!            temp(ipix) = temp(ipix)*( LOG(RIN%radius1(2,IDIM2))-LOG(RIN%radius1(1,IDIM2)) )
!!          enddo ! ipi
!!        call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
!!
!!        write(iu,'(a)') 'Mass concentration of BrC:'
!!        write(iu,'(i14)',advance='no') IDIM2
!!          do ipix=1,npixels
!!            temp(ipix) = sum(SD(1:NDIM3,IDIM2,ipix))
!!            temp(ipix) = GOUT%aerosol%chem%pixel(ipix)%fbrc(IDIM2)*1.2*1000*sum(SD(1:NDIM3,IDIM2,ipix))
!!            if(par_type_SD .eq. par_type_SD_TB)  &
!!            temp(ipix) = temp(ipix)*( LOG(RIN%radius1(2,IDIM2))-LOG(RIN%radius1(1,IDIM2)) )
!!          enddo ! ipi
!!        call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
!!
!!        write(iu,'(a)') 'Mass concentration of Quartz:'
!!        write(iu,'(i14)',advance='no') IDIM2
!!          do ipix=1,npixels
!!            temp(ipix) = sum(SD(1:NDIM3,IDIM2,ipix))
!!            temp(ipix) = GOUT%aerosol%chem%pixel(ipix)%finslbl(IDIM2)*2.5*1000*sum(SD(1:NDIM3,IDIM2,ipix))
!!            if(par_type_SD .eq. par_type_SD_TB)  &
!!            temp(ipix) = temp(ipix)*( LOG(RIN%radius1(2,IDIM2))-LOG(RIN%radius1(1,IDIM2)) )
!!          enddo ! ipi
!!        call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
!!
!!        write(iu,'(a)') 'Mass concentration of Soluble:'
!!        write(iu,'(i14)',advance='no') IDIM2
!!          do ipix=1,npixels
!!            temp(ipix) = sum(SD(1:NDIM3,IDIM2,ipix))
!!            temp(ipix) = GOUT%aerosol%chem%pixel(ipix)%fslbl(IDIM2)*1.73*1000*sum(SD(1:NDIM3,IDIM2,ipix))  !! NH4NO3
!!            if(par_type_SD .eq. par_type_SD_TB)  &
!!            temp(ipix) = temp(ipix)*( LOG(RIN%radius1(2,IDIM2))-LOG(RIN%radius1(1,IDIM2)) )
!!          enddo ! ipi
!!        call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
!!
!!
!!        write(iu,'(a)') 'Mass concentration of Water:'
!!        write(iu,'(i14)',advance='no') IDIM2
!!          do ipix=1,npixels
!!            temp(ipix) = sum(SD(1:NDIM3,IDIM2,ipix))
!!            temp(ipix) = GOUT%aerosol%chem%pixel(ipix)%fwtr(IDIM2)*1*1000*sum(SD(1:NDIM3,IDIM2,ipix))
!!            if(par_type_SD .eq. par_type_SD_TB)  &
!!            temp(ipix) = temp(ipix)*( LOG(RIN%radius1(2,IDIM2))-LOG(RIN%radius1(1,IDIM2)) )
!!          enddo ! ipi
!!        call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
!!
!!         else
!!        write(iu,'(a)') 'Mass concentration of Carbonaceous:'
!!        write(iu,'(i14)',advance='no') IDIM2
!!          do ipix=1,npixels
!!            temp(ipix) = sum(SD(1:NDIM3,IDIM2,ipix))
!!            temp(ipix) = GOUT%aerosol%chem%pixel(ipix)%fcarbon(IDIM2)*1.2*1000*sum(SD(1:NDIM3,IDIM2,ipix))
!!            if(par_type_SD .eq. par_type_SD_TB)  &
!!            temp(ipix) = temp(ipix)*( LOG(RIN%radius1(2,IDIM2))-LOG(RIN%radius1(1,IDIM2)) )
!!          enddo ! ipi
!!        call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
!!
!!        write(iu,'(a)') 'Mass concentration of Quartz:'
!!        write(iu,'(i14)',advance='no') IDIM2
!!          do ipix=1,npixels
!!            temp(ipix) = sum(SD(1:NDIM3,IDIM2,ipix))
!!            temp(ipix) = GOUT%aerosol%chem%pixel(ipix)%finslbl(IDIM2)*2.5*1000*sum(SD(1:NDIM3,IDIM2,ipix))
!!            if(par_type_SD .eq. par_type_SD_TB)  &
!!            temp(ipix) = temp(ipix)*( LOG(RIN%radius1(2,IDIM2))-LOG(RIN%radius1(1,IDIM2)) )
!!          enddo ! ipi
!!        call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
!!
!!        write(iu,'(a)') 'Mass concentration of Iron:'
!!        write(iu,'(i14)',advance='no') IDIM2
!!          do ipix=1,npixels
!!            temp(ipix) = sum(SD(1:NDIM3,IDIM2,ipix))
!!            temp(ipix) = GOUT%aerosol%chem%pixel(ipix)%firon(IDIM2)*4.8*1000*sum(SD(1:NDIM3,IDIM2,ipix))
!!            if(par_type_SD .eq. par_type_SD_TB)  &
!!            temp(ipix) = temp(ipix)*( LOG(RIN%radius1(2,IDIM2))-LOG(RIN%radius1(1,IDIM2)) )
!!          enddo ! ipi
!!        call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
!!
!!        write(iu,'(a)') 'Mass concentration of Soluble:'
!!        write(iu,'(i14)',advance='no') IDIM2
!!          do ipix=1,npixels
!!            temp(ipix) = sum(SD(1:NDIM3,IDIM2,ipix))
!!            temp(ipix) = GOUT%aerosol%chem%pixel(ipix)%fslbl(IDIM2)*2.17*1000*sum(SD(1:NDIM3,IDIM2,ipix))   !!sea salt
!!            if(par_type_SD .eq. par_type_SD_TB)  &
!!            temp(ipix) = temp(ipix)*( LOG(RIN%radius1(2,IDIM2))-LOG(RIN%radius1(1,IDIM2)) )
!!          enddo ! ipi
!!        call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
!!
!!
!!        write(iu,'(a)') 'Mass concentration of Water:'
!!        write(iu,'(i14)',advance='no') IDIM2
!!          do ipix=1,npixels
!!            temp(ipix) = sum(SD(1:NDIM3,IDIM2,ipix))
!!            temp(ipix) = GOUT%aerosol%chem%pixel(ipix)%fwtr(IDIM2)*1*1000*sum(SD(1:NDIM3,IDIM2,ipix))
!!            if(par_type_SD .eq. par_type_SD_TB)  &
!!            temp(ipix) = temp(ipix)*( LOG(RIN%radius1(2,IDIM2))-LOG(RIN%radius1(1,IDIM2)) )
!!          enddo ! ipi
!!        call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
!!
!!           endif
!!
!!        enddo ! IDIM2
!!
!      endif
!
!     endif
!
!   enddo  ! IDIM1

!!************add by lei on 20160914 to output aerosol chemical elements mass concentration *******

      return
      end subroutine print_aerosol_concentration 

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! Print Sphericity or Shape Distribution               	  

      subroutine print_shape_distribution (iu,RIN,segment_meas,GOUT)

      use mod_par_inv,     only : KIDIM2,KIDIM3,KIMAGE
      use mod_retr_settings_derived_type
      use mod_retr_general_output_derived_type
      use mod_print_array
      use mod_sdata_derived_type
      
      implicit none  
!----------------------------------------------------------------------------------------
! IN :
      integer,                   intent(in)  ::  iu
      type(retr_input_settings), intent(in)  ::  RIN
      type(output_segment_general), intent(in)  ::  GOUT
      type(segment_data),        intent(in)  ::  segment_meas
!----------------------------------------------------------------------------------------
! LOCAL :
      integer                              :: IDIM1,IDIM2,NDIM3
      real,dimension(KIDIM3,KIDIM2,KIMAGE) :: sph
      real,dimension(KIMAGE)               :: temp
      integer                              :: i,npixels,par_type
!----------------------------------------------------------------------------------------
      temp(:) = 0.0
      npixels = segment_meas%npixels

      do IDIM1=1,RIN%NDIM%n1
        par_type = RIN%NDIM%par_type(IDIM1)        
        if(par_type .gt. par_type_SHD_beg  .and. par_type .lt. par_type_SHD_end) then 
          call unpack_retr_param_to_print ( RIN,GOUT%retrieval%par,IDIM1,npixels,sph )      
          if(RIN%NDIM%n2(IDIM1) .eq. 1 .and. RIN%NDIM%n3(1,IDIM1) .eq. 2) & 
          sph(2,:,:) = sph(1,:,:)     	  
          do IDIM2=1,RIN%NDIM%n2(IDIM1)          
            write(iu,'(a,1x,i0,a)') '% of spherical particles',IDIM2,' - fraction'
            NDIM3 = RIN%NDIM%n3(IDIM2,IDIM1)
            do i=1,NDIM3
              write(iu,'(i14)',advance='no') i
              temp(1:npixels) = sph(i,IDIM2,1:npixels)*100.
              call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
            enddo ! i
          enddo ! IDIM2
        exit
        endif ! par_type .gt. par_type_SHD_beg  .and. 
      enddo ! IDIM1

      return
      end subroutine print_shape_distribution 

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! Print Aerosol Vertical Profile (if lidar measurements are present) or 
! mean height of the profole                          	  

      subroutine print_aerosol_profile (iu,RIN,segment_meas,GOUT)

      use mod_par_inv, only : KIDIM2, KIDIM3, KIMAGE, KVERTM !,max_NH
      use mod_par_OS,  only : KSD, NBVM, HMAX_atm
      use mod_retr_settings_derived_type
      use mod_retr_general_output_derived_type
      use mod_print_array
      use mod_sdata_derived_type
      use mod_sdata, only : get_HVP_lidar
      use mod_os
      
      implicit none  
!----------------------------------------------------------------------------------------
! IN :
      integer,                     intent(in)  ::  iu
      type(retr_input_settings),   intent(in)  ::  RIN
      type(segment_data),          intent(in)  ::  segment_meas
      type(output_segment_general),intent(in)  ::  GOUT
!----------------------------------------------------------------------------------------
! LOCAL :
      integer                              :: IDIM1,IDIM2,NDIM3,npixels
      real,dimension(KIDIM3,KIDIM2,KIMAGE) :: h01,std
      integer                              :: NHVP_meas,NHVP_retr
      real                                 :: HGR_km
      real,dimension(KVERTM)               :: HVP_meas,HVP_retr_km
      integer                              :: i,ipix,par_type
      real                                 :: xnorm
      integer                              :: NH
      real,dimension(KVERTM)               :: H_km
      character(len=20) :: distr_type
      real, dimension(KVERTM) :: prof_temp
      real,dimension(KIDIM2,KIMAGE) :: h01_HGR
!----------------------------------------------------------------------------------------
      npixels = segment_meas%npixels

      do IDIM1=1,RIN%NDIM%n1
        par_type = RIN%NDIM%par_type(IDIM1)        
        if(par_type .gt. par_type_AVP_beg  .and. par_type .lt. par_type_AVP_end) then 
          call unpack_retr_param_to_print ( RIN,GOUT%retrieval%par,IDIM1,npixels,h01 ) 	  
          if(par_type .eq. par_type_AVP_prof) then
            call get_HVP_lidar (  segment_meas,       & ! IN
                                  NHVP_meas, HVP_meas & ! INOUT
                               )
            distr_type = 'lut'
            do ipix=1,npixels
              HGR_km = segment_meas%pixels(ipix)%MASL*0.001
              do IDIM2=1,RIN%NDIM%n2(IDIM1)
                NHVP_retr = RIN%NDIM%n3(IDIM2,IDIM1)
                call get_AVP_altitudes_retr ( HGR_km,                    &
                                              NHVP_meas, HVP_meas*0.001, &
                                              NHVP_retr, HVP_retr_km )
                ! Altitudes for aerosol concentration profile normalization
                call grid_altitudes_LUT ( HGR_km, HMAX_atm*0.001,               & ! IN
                                          NHVP_retr, HVP_retr_km(1:NHVP_retr),  &
                                          NH, H_km(1:NHVP_retr+2)               & ! OUT
                                        )
                ! Aerosol concentration profile normalization
                call discrvd_single_atm_component ( distr_type, 0.0, 0.0,                 &
                                                    NHVP_retr, HVP_retr_km(1:NHVP_retr),  &
                                                    h01(1:NHVP_retr,IDIM2,ipix),          &
                                                    NH, H_km(1:NH),                       &
                                                    xnorm, prof_temp(1:NH)                & ! OUT
                                                  )
                h01(1:NHVP_retr,IDIM2,ipix) = h01(1:NHVP_retr,IDIM2,ipix)/xnorm
                h01(1:NHVP_retr,IDIM2,ipix) = h01(1:NHVP_retr,IDIM2,ipix)*0.001 ! AL 0.001 to make aerosol vertical profile in 1/m
                if(HVP_retr_km(NHVP_retr) .ne. HGR_km) h01_HGR(IDIM2,ipix) = prof_temp(NH)/xnorm*0.001 ! vertical profile at ground level (assumed)
              enddo ! IDIM2
            enddo ! ipix
          endif ! par_type .eq. par_type_AVP_prof 
          do IDIM2=1,RIN%NDIM%n2(IDIM1)    ! aerosol component loop
            NDIM3 = RIN%NDIM%n3(IDIM2,IDIM1)    
            if(par_type .eq. par_type_AVP_par_height) then
              write(iu,'(a)') 'Aerosol profile mean height (m)'
              do i=1,NDIM3
                write(iu,'(i14)',advance='no')  IDIM2
                call rprint_array_real (iu,h01(i,IDIM2,:), '(e14.5)', 1, npixels)
              enddo ! i
            elseif(par_type .eq. par_type_AVP_prof) then
              write(iu,'(a,1x,i0,a)') 'Aerosol vertical profile [1/m]', IDIM2,' - fraction'
              do i=1,NDIM3
                write(iu,'(i4,f10.2)',advance='no') i,HVP_retr_km(i)*1000.
                call rprint_array_real (iu,h01(i,IDIM2,:), '(e14.5)', 1, npixels)
              enddo ! I
              if(HVP_retr_km(NHVP_retr) .ne. HGR_km) then
                ! print vertical profile at ground level (assumed)
                write(iu,'(a,i3,f10.2)',advance='no') '*',NHVP_retr+1,HGR_km*1000.
                call rprint_array_real (iu,h01_HGR(IDIM2,:), '(e14.5)', 1, npixels)
              endif
            endif ! par_type .eq. par_type_AVP_par_height
          enddo ! NDIM2
        elseif(par_type .eq. par_type_AVP_par_std) then
          call unpack_retr_param_to_print ( RIN,GOUT%retrieval%par,IDIM1,npixels,std ) 	  
          write(iu,'(a)') 'Aerosol profile standard deviation (m)'
          do IDIM2=1,RIN%NDIM%n2(IDIM1)    ! aerosol component loop
            NDIM3 = RIN%NDIM%n3(IDIM2,IDIM1)    
            do i=1,NDIM3
              write(iu,'(i14)',advance='no') IDIM2
              call rprint_array_real (iu,std(i,IDIM2,:), '(e14.5)', 1, npixels)
            enddo ! i
          enddo ! NDIM2
        endif ! par_type .gt. par_type_AVP_beg  .and. 
      enddo ! IDIM1
      
      return
      end subroutine print_aerosol_profile 

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! Print Lidar Calibration coefficient               	  

      subroutine print_lidar_calibration_coeffitient (iu,RIN,segment_meas,GOUT)

      use mod_par_inv,     only : KIDIM2,KIDIM3,KIMAGE
      use mod_retr_settings_derived_type
      use mod_retr_general_output_derived_type
      use mod_print_array
      use mod_sdata_derived_type
      
      implicit none  
!----------------------------------------------------------------------------------------
! IN :
      integer,                     intent(in)  ::  iu
      type(retr_input_settings),   intent(in)  ::  RIN
      type(segment_data),          intent(in)  ::  segment_meas
      type(output_segment_general),intent(in)  ::  GOUT
!----------------------------------------------------------------------------------------
! LOCAL :
      integer                              :: IDIM1,IDIM2,NDIM3
      real,dimension(KIDIM3,KIDIM2,KIMAGE) :: CL
      real,dimension(KIMAGE)               :: temp
      integer                              :: i,npixels,par_type
!----------------------------------------------------------------------------------------
      temp(:) = 0.0
      npixels = segment_meas%npixels

      do IDIM1=1,RIN%NDIM%n1
        par_type = RIN%NDIM%par_type(IDIM1)        
        if(par_type .gt. par_type_ADD_beg  .and. par_type .lt. par_type_ADD_end) then    
          if(par_type .eq. par_type_CL) then
            call unpack_retr_param_to_print ( RIN,GOUT%retrieval%par,IDIM1,npixels,CL ) 	  
            do IDIM2=1,RIN%NDIM%n2(IDIM1)          
              write(iu,'(a)') 'Lidar calibration coefficient'
              NDIM3 = RIN%NDIM%n3(IDIM2,IDIM1)
              do i=1,NDIM3
                write(iu,'(i14)',advance='no') i
                call rprint_array_real (iu,CL(i,IDIM2,:), '(e14.5)', 1, npixels)
              enddo ! i
            enddo ! IDIM2
          exit
          endif 
        endif ! par_type .gt. par_type_ADD_beg  .and. 
      enddo ! IDIM1

      return
      end subroutine print_lidar_calibration_coeffitient 

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! Print Sutface Parameters               	  

      subroutine print_surface_parameters (iu,RIN,segment_meas,GOUT)

      use mod_par_inv,     only : KIDIM2,KIDIM3,KIMAGE
      use mod_retr_settings_derived_type
      use mod_retr_general_output_derived_type
      use mod_print_array
      use mod_sdata_derived_type
      
      implicit none  
!----------------------------------------------------------------------------------------
! IN :
      integer,                     intent(in)  ::  iu
      type(retr_input_settings),   intent(in)  ::  RIN
      type(segment_data),          intent(in)  ::  segment_meas
      type(output_segment_general),intent(in)  ::  GOUT
!----------------------------------------------------------------------------------------
! LOCAL :
      integer                              :: IDIM1,IDIM2,IDIM3,NDIM3,ns3
      real,dimension(KIDIM3,KIDIM2,KIMAGE) :: SURF
      real,dimension(KIMAGE)               :: temp
      integer                              :: npixels,par_type
!----------------------------------------------------------------------------------------
      temp(:) = 0.0
      npixels = segment_meas%npixels

      do IDIM1=1,RIN%NDIM%n1
        par_type = RIN%NDIM%par_type(IDIM1)        
        if(par_type .gt. par_type_SURF1_land_beg .and. par_type .lt. par_type_SURF1_land_end) then      
          call unpack_retr_param_to_print ( RIN,GOUT%retrieval%par,IDIM1,npixels,SURF )
          write(iu,'(a)') 'Wavelength (um), BRDF parameters'
          if(par_type .eq. par_type_SURF1_land_RPV_BRDF) then  
            ns3 = 3
            SURF(1:RIN%NDIM%n3(ns3,IDIM1),ns3,1:npixels) =         & 
            SURF(1:RIN%NDIM%n3(ns3,IDIM1),ns3,1:npixels) - 1.0
          endif
          do IDIM2=1,RIN%NDIM%n2(IDIM1)
            write(iu,'(i4)') IDIM2
            NDIM3 = RIN%NDIM%n3(IDIM2,IDIM1)
            do IDIM3=1,NDIM3
              write(iu,'(f14.5)',advance='no') RIN%wave(IDIM3)
              call rprint_array_real (iu,SURF(IDIM3,IDIM2,:), '(e14.5)', 1, npixels)
            enddo ! IDIM3
          enddo ! IDIM2
        elseif(par_type .gt. par_type_SURF2_land_beg .and. par_type .lt. par_type_SURF2_land_end) then
          call unpack_retr_param_to_print ( RIN,GOUT%retrieval%par,IDIM1,npixels,SURF )
          write(iu,'(a)') 'Wavelength (um), BPDF parameters'
          do IDIM2=1,RIN%NDIM%n2(IDIM1)
            write(iu,'(i4)') IDIM2
            NDIM3 = RIN%NDIM%n3(IDIM2,IDIM1)
            do IDIM3=1,NDIM3
              write(iu,'(f14.5)',advance='no') RIN%wave(IDIM3)
              call rprint_array_real (iu,SURF(IDIM3,IDIM2,:), '(e14.5)', 1, npixels)
            enddo ! IDIM3
          enddo ! IDIM2
        elseif(par_type .gt. par_type_SURF_water_beg .and. par_type .lt. par_type_SURF_water_end) then
          call unpack_retr_param_to_print ( RIN,GOUT%retrieval%par,IDIM1,npixels,SURF )
          write(iu,'(a)') 'Wavelength (um), Water surface parameters'
          do IDIM2=1,RIN%NDIM%n2(IDIM1)
            write(iu,'(i4)') IDIM2
            NDIM3 = RIN%NDIM%n3(IDIM2,IDIM1)
            do IDIM3=1,NDIM3
              write(iu,'(f14.5)',advance='no') RIN%wave(IDIM3)
              call rprint_array_real (iu,SURF(IDIM3,IDIM2,:), '(e14.5)', 1, npixels)
            enddo ! IDIM3
          enddo ! IDIM2
        endif ! par_type .gt. par_type_SURF1_land_beg .and. 
      enddo ! IDIM1
           
      return
      end subroutine print_surface_parameters 

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! Print Real part of Refractiv Index               	  

      subroutine print_refractive_index_real (iu,RIN,segment_meas,GOUT)

      use mod_par_inv,     only : KIDIM2,KIDIM3,KIMAGE
      use mod_retr_settings_derived_type
      use mod_retr_general_output_derived_type
      use mod_print_array
      use mod_sdata_derived_type
      
      implicit none  
!----------------------------------------------------------------------------------------
! IN :
      integer,                     intent(in)  ::  iu
      type(retr_input_settings),   intent(in)  ::  RIN
      type(segment_data),          intent(in)  ::  segment_meas
      type(output_segment_general),intent(in)  ::  GOUT
!----------------------------------------------------------------------------------------
! LOCAL :
      integer                              :: IDIM1,IDIM2,IDIM3,NDIM3
      real,dimension(KIDIM3,KIDIM2,KIMAGE) :: RERI
      real,dimension(KIMAGE)               :: temp
      integer                              :: iw,npixels,par_type
!----------------------------------------------------------------------------------------
      !temp(:) = 0.0
      npixels = segment_meas%npixels
      
      do IDIM1=1,RIN%NDIM%n1
        par_type = RIN%NDIM%par_type(IDIM1)        
        if(par_type .gt. par_type_RERI_beg .and. par_type .lt. par_type_RERI_end) then
          !if(par_type .ne. par_type_CXRI_chem .and. par_type .ne. par_type_CXRI_nmix) then
            !call unpack_retr_param_to_print ( RIN,GOUT%retrieval%par,IDIM1,npixels,RERI )
            !write(iu,'(a)') 'Wavelength (um), REAL Ref. Index :'
            !do IDIM2=1,RIN%NDIM%n2(IDIM1)
            !  write(iu,*) IDIM2
            !  NDIM3 = RIN%NDIM%n3(IDIM2,IDIM1)
            !  do IDIM3=1,NDIM3
            !    write(iu,'(e14.5)',advance='no') RIN%wave(IDIM3)
            !    call rprint_array_real (iu,RERI(IDIM3,IDIM2,:), '(e14.5)', 1, npixels)
            !  enddo ! IDIM3
            !enddo ! IDIM2
          !exit
          !endif ! par_type .ne. par_type_CXRI_chem
          write(iu,'(a)') 'Wavelength (um), REAL Ref. Index'
          do IDIM2=1,RIN%NDIM%n2(IDIM1)    ! aerosol component loop
            write(iu,*) IDIM2
            !NDIM3 = RIN%NDIM%n3(IDIM2,IDIM1)
            NDIM3 = RIN%NW
            do IDIM3=1,NDIM3
              write(iu,'(f14.5)',advance='no') RIN%wave(IDIM3)
              temp(1:npixels) = GOUT%aerosol%rind%pixel(1:npixels)%wl(IDIM3)%mreal(IDIM2)
              call rprint_array_real (iu,temp(:), '(f14.5)', 1, npixels)
            enddo ! IDIM3
          enddo ! IDIM2
        exit
        endif ! par_type .gt. par_type_RERI_beg .and.
      enddo ! IDIM1

      return
      end subroutine print_refractive_index_real 


! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! Print Imaginary part of Refractiv Index               	  

      subroutine print_refractive_index_imaginary (iu,RIN,segment_meas,GOUT)

      use mod_par_inv,     only : KIDIM2,KIDIM3,KIMAGE
      use mod_retr_settings_derived_type
      use mod_retr_general_output_derived_type
      use mod_print_array
      use mod_sdata_derived_type
      
      implicit none  
!----------------------------------------------------------------------------------------
! IN :
      integer,                     intent(in)  ::  iu
      type(retr_input_settings),   intent(in)  ::  RIN
      type(segment_data),          intent(in)  ::  segment_meas
      type(output_segment_general),intent(in)  ::  GOUT
!----------------------------------------------------------------------------------------
! LOCAL :
      integer                              :: IDIM1,IDIM2,IDIM3,NDIM3
      real,dimension(KIDIM3,KIDIM2,KIMAGE) :: IMRI
      real,dimension(KIMAGE)               :: temp
      integer                              :: iw,npixels,par_type
!----------------------------------------------------------------------------------------
      !temp(:) = 0.0
      npixels = segment_meas%npixels
      
      do IDIM1=1,RIN%NDIM%n1
        par_type = RIN%NDIM%par_type(IDIM1)        
        if((par_type .gt. par_type_IMRI_beg .and. par_type .lt. par_type_IMRI_end) &
            .or. &
           (par_type .eq. par_type_CXRI_chem .or. par_type .eq. par_type_CXRI_nmix)) then
          !call unpack_retr_param_to_print ( RIN,GOUT%retrieval%par,IDIM1,npixels,IMRI )
          !write(iu,'(a)') 'Wavelength (um), IMAG Ref. Index :'
          !do IDIM2=1,RIN%NDIM%n2(IDIM1)
            !write(iu,*) IDIM2
            !NDIM3 = RIN%NDIM%n3(IDIM2,IDIM1)
            !do IDIM3=1,NDIM3
              !write(iu,'(e14.5)',advance='no') RIN%wave(IDIM3)
              !call rprint_array_real (iu,IMRI(IDIM3,IDIM2,:), '(e14.5)', 1, npixels)
            !enddo ! IDIM3
          !enddo ! IDIM2
          !exit
          write(iu,'(a)') 'Wavelength (um), IMAG Ref. Index'
          do IDIM2=1,RIN%NDIM%n2(IDIM1)    ! aerosol component loop
            write(iu,*) IDIM2
            !NDIM3 = RIN%NDIM%n3(IDIM2,IDIM1)
            NDIM3 = RIN%NW
            do IDIM3=1,NDIM3
              write(iu,'(f14.5)',advance='no') RIN%wave(IDIM3)
              temp(1:npixels) = GOUT%aerosol%rind%pixel(1:npixels)%wl(IDIM3)%mimag(IDIM2)
              call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
            enddo ! IDIM3
          enddo ! IDIM2
        exit
        endif ! par_type .gt. par_type_RERI_beg .and.
      enddo ! IDIM1

      return
      end subroutine print_refractive_index_imaginary 

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! Print chemistry              	  

      subroutine print_chemistry (iu,RIN,segment_meas,GOUT)

      use mod_par_inv,     only : KIMAGE
      use mod_retr_settings_derived_type
      use mod_retr_general_output_derived_type
      use mod_print_array
      use mod_sdata_derived_type
      
      implicit none  
!----------------------------------------------------------------------------------------
! IN :
      integer,                     intent(in)  ::  iu
      type(retr_input_settings),   intent(in)  ::  RIN
      type(segment_data),          intent(in)  ::  segment_meas
      type(output_segment_general),intent(in)  ::  GOUT
!----------------------------------------------------------------------------------------
! LOCAL :
      integer                              :: IDIM1, IDIM2, IDIM3, NDIM3
      real,dimension(KIMAGE)               :: temp, aa
      integer                              :: iw, npixels, par_type
!----------------------------------------------------------------------------------------
      temp(:) = 0.0
      npixels = segment_meas%npixels
      
      do IDIM1=1,RIN%NDIM%n1
        par_type = RIN%NDIM%par_type(IDIM1)        
        if(par_type .gt. par_type_RERI_beg .and. par_type .lt. par_type_RERI_end) then
          if(par_type .eq. par_type_CXRI_chem) then
!            do IDIM2=1,RIN%NDIM%n2(IDIM1)    ! aerosol component loop
!              write(iu,'(a)') 'Relative humidity'
!                write(iu,'(i14)',advance='no') IDIM2
!                temp(1:npixels) = GOUT%aerosol%chem%pixel(1:npixels)%rh(IDIM2)
!                call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
!              write(iu,'(a)') 'Volume Fraction of Quartz+Clay (Insoluble):'   ! change 'Quartz' to 'Quartz+Clay' by lei
!                write(iu,'(i14)',advance='no') IDIM2
!                temp(1:npixels) = GOUT%aerosol%chem%pixel(1:npixels)%fquartz(IDIM2)
!                call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
!              write(iu,'(a)') 'Volume Fraction of Carbonaceous:'              ! change 'soot' to 'Carbonaceous' (Soot+BrC) by lei
!                write(iu,'(i14)',advance='no') IDIM2
!                temp(1:npixels) = GOUT%aerosol%chem%pixel(1:npixels)%fsoot(IDIM2)
!                call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
!              write(iu,'(a)') 'Volume Fraction of Iron'
!                write(iu,'(i14)',advance='no') IDIM2
!                temp(1:npixels) = GOUT%aerosol%chem%pixel(1:npixels)%firon(IDIM2)
!                call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
!              write(iu,'(a)') 'Volume Fraction of Soluble'
!                write(iu,'(i14)',advance='no') IDIM2
!                temp(1:npixels) = GOUT%aerosol%chem%pixel(1:npixels)%fslbl(IDIM2)
!                call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
!              write(iu,'(a)') 'Volume Fraction of Water'
!                write(iu,'(i14)',advance='no') IDIM2
!                temp(1:npixels) = GOUT%aerosol%chem%pixel(1:npixels)%fwtr(IDIM2)
!                call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)        
!            enddo ! IDIM2
!

!************************14/11/2016
            select case( RIN%NSD )
            case(1) ! one component aerosol
              IDIM2 = 1
              write(iu,'(a)') 'Relative humidity'
                write(iu,'(i14)',advance='no') IDIM2
                temp(1:npixels) = GOUT%aerosol%chem%pixel(1:npixels)%rh(IDIM2)
                call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
              write(iu,'(a)') 'Volume Fraction of Quartz+Clay (Insoluble):'   ! change 'Quartz' to 'Quartz+Clay' by lei
                write(iu,'(i14)',advance='no') IDIM2
                temp(1:npixels) = GOUT%aerosol%chem%pixel(1:npixels)%fquartz(IDIM2)
                call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
              write(iu,'(a)') 'Volume Fraction of Carbonaceous:'              ! change 'soot' to 'Carbonaceous' (Soot+BrC) by lei
                write(iu,'(i14)',advance='no') IDIM2
                temp(1:npixels) = GOUT%aerosol%chem%pixel(1:npixels)%fsoot(IDIM2)
                call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
              write(iu,'(a)') 'Volume Fraction of Iron'
                write(iu,'(i14)',advance='no') IDIM2
                temp(1:npixels) = GOUT%aerosol%chem%pixel(1:npixels)%firon(IDIM2)
                call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
              write(iu,'(a)') 'Volume Fraction of Soluble'
                write(iu,'(i14)',advance='no') IDIM2
                temp(1:npixels) = GOUT%aerosol%chem%pixel(1:npixels)%fslbl(IDIM2)
                call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
              write(iu,'(a)') 'Volume Fraction of Water'
                write(iu,'(i14)',advance='no') IDIM2
                temp(1:npixels) = GOUT%aerosol%chem%pixel(1:npixels)%fwtr(IDIM2)
                call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)        
            case(2) ! 2 component aerosol
              do IDIM2=1,RIN%NDIM%n2(IDIM1)
              if( IDIM2 .eq. 1 ) then
              write(iu,'(a)') 'Relative humidity'
                write(iu,'(i14)',advance='no') IDIM2
                temp(1:npixels) = GOUT%aerosol%chem%pixel(1:npixels)%rh(IDIM2)
                call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
              write(iu,'(a)') 'Volume Fraction of BC:'
                write(iu,'(i14)',advance='no') IDIM2
                temp(1:npixels) = GOUT%aerosol%chem%pixel(1:npixels)%fsoot(IDIM2)
                call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
              write(iu,'(a)') 'Volume Fraction of BrC:'
                write(iu,'(i14)',advance='no') IDIM2
                temp(1:npixels) = GOUT%aerosol%chem%pixel(1:npixels)%fbrc(IDIM2)
                call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
              write(iu,'(a)') 'Volume Fraction of Quartz'
                write(iu,'(i14)',advance='no') IDIM2
                temp(1:npixels) = GOUT%aerosol%chem%pixel(1:npixels)%fquartz(IDIM2)
                call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
              write(iu,'(a)') 'Volume Fraction of Soluble'
                write(iu,'(i14)',advance='no') IDIM2
                temp(1:npixels) = GOUT%aerosol%chem%pixel(1:npixels)%fslbl(IDIM2)
                call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
              write(iu,'(a)') 'Volume Fraction of Water'
                write(iu,'(i14)',advance='no') IDIM2
                temp(1:npixels) = GOUT%aerosol%chem%pixel(1:npixels)%fwtr(IDIM2)
                call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)        
              else
              write(iu,'(a)') 'Relative humidity'
                write(iu,'(i14)',advance='no') IDIM2
                temp(1:npixels) = GOUT%aerosol%chem%pixel(1:npixels)%rh(IDIM2)
                call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
              write(iu,'(a)') 'Volume Fraction of Quartz:'
                write(iu,'(i14)',advance='no') IDIM2
                temp(1:npixels) = GOUT%aerosol%chem%pixel(1:npixels)%fquartz(IDIM2)
                call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
              write(iu,'(a)') 'Volume Fraction of Iron'
                write(iu,'(i14)',advance='no') IDIM2
                temp(1:npixels) = GOUT%aerosol%chem%pixel(1:npixels)%firon(IDIM2)
                call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
              write(iu,'(a)') 'Volume Fraction of Soluble'
                write(iu,'(i14)',advance='no') IDIM2
                temp(1:npixels) = GOUT%aerosol%chem%pixel(1:npixels)%fslbl(IDIM2)
                call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
              write(iu,'(a)') 'Volume Fraction of Water'
                write(iu,'(i14)',advance='no') IDIM2
                temp(1:npixels) = GOUT%aerosol%chem%pixel(1:npixels)%fwtr(IDIM2)
                call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)        

              endif ! IDIM2 .eq. 1
              enddo ! IDIM2
            end select
!***************************14/11/2016

          elseif(par_type .eq. par_type_CXRI_nmix) then
            select case( RIN%NSD )
            case(1) ! one component aerosol
              IDIM2 = 1
              write(iu,'(a)') 'Fraction of Carbonaceous:'      ! change 'soot' to 'Carbonaceous' (Soot+BrC) by lei
                write(iu,'(i14)',advance='no') IDIM2
                temp(1:npixels) = GOUT%aerosol%chem%pixel(1:npixels)%fsoot(IDIM2)
                call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
!! remove BrC element in volume mixture model by lei
!              write(iu,'(a)') 'Fraction of BrC:'
!                write(iu,'(i14)',advance='no') IDIM2
!                temp(1:npixels) = GOUT%aerosol%chem%pixel(1:npixels)%fslbl(IDIM2)
!                call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
              write(iu,'(a)') 'Fraction of Quartz+Clay:'       ! change 'Quartz' to 'Quartz+Clay' by lei
                write(iu,'(i14)',advance='no') IDIM2
                temp(1:npixels) = GOUT%aerosol%chem%pixel(1:npixels)%fquartz(IDIM2)
                call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
              write(iu,'(a)') 'Fraction of Iron:'
                write(iu,'(i14)',advance='no') IDIM2
                temp(1:npixels) = GOUT%aerosol%chem%pixel(1:npixels)%firon(IDIM2)
                call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
              write(iu,'(a)') 'Fraction of Water:'
                write(iu,'(i14)',advance='no') IDIM2
                temp(1:npixels) = GOUT%aerosol%chem%pixel(1:npixels)%fwtr(IDIM2)
                call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)

!!*********improve two modes by lei on 20160608**start*************
            case(2) ! 2 component aerosol
              do IDIM2=1,RIN%NDIM%n2(IDIM1)
              if( IDIM2 .eq. 1 ) then
              write(iu,'(a)') 'Fraction of BC:'
                write(iu,'(i14)',advance='no') IDIM2
                temp(1:npixels) = GOUT%aerosol%chem%pixel(1:npixels)%fsoot(IDIM2)
                call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)

              write(iu,'(a)') 'Fraction of BrC:'
                write(iu,'(i14)',advance='no') IDIM2
                temp(1:npixels) = GOUT%aerosol%chem%pixel(1:npixels)%fslbl(IDIM2)
                call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
              write(iu,'(a)') 'Fraction of Quartz:'
                write(iu,'(i14)',advance='no') IDIM2
                temp(1:npixels) = GOUT%aerosol%chem%pixel(1:npixels)%fquartz(IDIM2)
                call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
!              write(iu,'(a)') 'Fraction of Iron:'    ! remove iron in fine mode by lei
!                write(iu,'(i14)',advance='no') IDIM2
!                temp(1:npixels) = GOUT%aerosol%chem%pixel(1:npixels)%firon(IDIM2)
!                call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
              write(iu,'(a)') 'Fraction of Water:'
                write(iu,'(i14)',advance='no') IDIM2
                temp(1:npixels) = GOUT%aerosol%chem%pixel(1:npixels)%fwtr(IDIM2)
                call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
              else
!              write(iu,'(a)') 'Fraction of Carbonaceous:'      ! change 'soot' to 'Carbonaceous' (Soot+BrC) by lei
!                write(iu,'(i14)',advance='no') IDIM2
!                temp(1:npixels) = GOUT%aerosol%chem%pixel(1:npixels)%fsoot(IDIM2)
!                call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
!! remove BrC element in volume mixture model by lei on 20160330
!
!              write(iu,'(a)') 'Fraction of BC:'
!                write(iu,'(i14)',advance='no') IDIM2
!                temp(1:npixels) = GOUT%aerosol%chem%pixel(1:npixels)%fsoot(IDIM2)
!                call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
!              write(iu,'(a)') 'Fraction of BrC:'
!                write(iu,'(i14)',advance='no') IDIM2
!                temp(1:npixels) = GOUT%aerosol%chem%pixel(1:npixels)%fslbl(IDIM2)
!                call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
              write(iu,'(a)') 'Fraction of Quartz:'
                write(iu,'(i14)',advance='no') IDIM2
                temp(1:npixels) = GOUT%aerosol%chem%pixel(1:npixels)%fquartz(IDIM2)
                call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
              write(iu,'(a)') 'Fraction of Iron:'
                write(iu,'(i14)',advance='no') IDIM2
                temp(1:npixels) = GOUT%aerosol%chem%pixel(1:npixels)%firon(IDIM2)
                call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
              write(iu,'(a)') 'Fraction of Water:'
                write(iu,'(i14)',advance='no') IDIM2
                temp(1:npixels) = GOUT%aerosol%chem%pixel(1:npixels)%fwtr(IDIM2)
                call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
              endif ! IDIM2 .eq. 1
              enddo ! IDIM2
            end select
          endif
!!*********improve two modes by lei on 20160608**end**************

!          if(par_type .eq. par_type_CXRI_chem .or. par_type .eq. par_type_CXRI_nmix) then
!            !write(iu,*)
!            write(iu,'(a)') 'Wavelength (um), REAL Ref. Index'
!            do IDIM2=1,RIN%NDIM%n2(IDIM1)    ! aerosol component loop
!              write(iu,*) IDIM2
!!              NDIM3 = RIN%NDIM%n3(IDIM2,IDIM1)
!              NDIM3 = RIN%NW
!              do IDIM3=1,NDIM3
!                write(iu,'(f14.5)',advance='no') RIN%wave(IDIM3)
!                temp(1:npixels) = GOUT%aerosol%rind%pixel(1:npixels)%wl(IDIM3)%mreal(IDIM2)
!                call rprint_array_real (iu,temp(:), '(f14.5)', 1, npixels)
!              enddo ! IDIM3
!            enddo ! IDIM2
!            !write(iu,*)
!            write(iu,'(a)') 'Wavelength (um), IMAG Ref. Index'
!            do IDIM2=1,RIN%NDIM%n2(IDIM1)    ! aerosol component loop
!              write(iu,*) IDIM2
!!              NDIM3 = RIN%NDIM%n3(IDIM2,IDIM1)
!              NDIM3 = RIN%NW
!              do IDIM3=1,NDIM3
!                write(iu,'(f14.5)',advance='no') RIN%wave(IDIM3)
!                temp(1:npixels) = GOUT%aerosol%rind%pixel(1:npixels)%wl(IDIM3)%mimag(IDIM2)
!                call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
!              enddo ! IDIM3
!            enddo ! IDIM2
!          endif
        exit
        endif ! par_type .gt. par_type_RERI_beg .and.
      enddo ! IDIM1
                        
      return
      end subroutine print_chemistry 

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! Print Optical Thickness              	  
      
      subroutine print_optical_thickness(iu,RIN,segment_meas,GOUT)

      use mod_par_inv,     only : KIMAGE
      use mod_retr_settings_derived_type
      use mod_retr_general_output_derived_type
      use mod_print_array
      use mod_sdata_derived_type
      
      implicit none  
!----------------------------------------------------------------------------------------
! IN :
      integer,                     intent(in)  ::  iu
      type(retr_input_settings),   intent(in)  ::  RIN
      type(segment_data),          intent(in)  ::  segment_meas
      type(output_segment_general),intent(in)  ::  GOUT
!----------------------------------------------------------------------------------------
! LOCAL :
      real,dimension(KIMAGE)               :: temp
      integer                              :: iw,isd,npixels
!----------------------------------------------------------------------------------------
      temp(:) = 0.0
      npixels = segment_meas%npixels
                      
      write(iu,'(a)') 'Wavelength (um), Total_AOD (unitless or 1/um)'
      do iw=1,RIN%nw
        write(iu,'(f14.5)',advance='no') RIN%wave(iw)
        temp(1:npixels) = GOUT%aerosol%opt%pixel(1:npixels)%wl(iw)%extt
        call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
      enddo ! iw

      if(RIN%NSD .gt. 1) then
        do isd=1,RIN%NSD    ! aerosol component loop
          if(isd .eq. 1) write(iu,'(a)') 'Wavelength (um), Fine_AOD (unitless or 1/um)'
          if(isd .eq. 2) write(iu,'(a)') 'Wavelength (um), Coarse_AOD (unitless or 1/um)'
          do iw=1,RIN%nw
            write(iu,'(f14.5)',advance='no') RIN%wave(iw)
            temp(1:npixels) = GOUT%aerosol%opt%pixel(1:npixels)%wl(iw)%ext(isd)
            call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
          enddo ! iw
        enddo ! isd
      endif ! RIN%NSD .gt. 1
  
      return
      end subroutine print_optical_thickness 

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! Print Single Scattering Abledo              	  
      
      subroutine print_single_scattering_albedo(iu,RIN,segment_meas,GOUT)

      use mod_par_inv,     only : KIMAGE
      use mod_retr_settings_derived_type
      use mod_retr_general_output_derived_type
      use mod_print_array
      use mod_sdata_derived_type
      
      implicit none  
!----------------------------------------------------------------------------------------
! IN :
      integer,                     intent(in)  ::  iu
      type(retr_input_settings),   intent(in)  ::  RIN
      type(segment_data),          intent(in)  ::  segment_meas
      type(output_segment_general),intent(in)  ::  GOUT
!----------------------------------------------------------------------------------------
! LOCAL :
      real,dimension(KIMAGE)               :: temp
      integer                              :: iw,isd,npixels
!----------------------------------------------------------------------------------------
      temp(:) = 0.0
      npixels = segment_meas%npixels
                      
      write(iu,'(a)') 'Wavelength (um), Total_SSA'
      do iw=1,RIN%nw
        write(iu,'(f14.5)',advance='no') RIN%wave(iw)
        temp(1:npixels) = GOUT%aerosol%opt%pixel(1:npixels)%wl(iw)%ssat
        call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
      enddo ! iw

      if(RIN%NSD .gt. 1) then
        do isd=1,RIN%NSD    ! aerosol component loop
          if(isd .eq. 1) write(iu,'(a)') 'Wavelength (um), Fine_SSA'
          if(isd .eq. 2) write(iu,'(a)') 'Wavelength (um), Coarse_SSA'
          do iw=1,RIN%nw
            write(iu,'(f14.5)',advance='no') RIN%wave(iw)
            temp(1:npixels) = GOUT%aerosol%opt%pixel(1:npixels)%wl(iw)%ssa(isd)
            call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
          enddo ! iw
        enddo ! isd
      endif ! RIN%NSD .gt. 1

      return
      end subroutine print_single_scattering_albedo 
! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! Print Absorption
      
      subroutine print_absorption(iu,RIN,segment_meas,GOUT)

      use mod_par_inv,     only : KIMAGE
      use mod_retr_settings_derived_type
      use mod_retr_general_output_derived_type
      use mod_print_array
      use mod_sdata_derived_type
      
      implicit none  
!----------------------------------------------------------------------------------------
! IN :
      integer,                     intent(in)  ::  iu
      type(retr_input_settings),   intent(in)  ::  RIN
      type(segment_data),          intent(in)  ::  segment_meas
      type(output_segment_general),intent(in)  ::  GOUT
!----------------------------------------------------------------------------------------
! LOCAL :
      real,dimension(KIMAGE)               :: temp
      integer                              :: iw,isd,npixels
!----------------------------------------------------------------------------------------
      temp(:) = 0.0
      npixels = segment_meas%npixels
                      
      write(iu,'(a)') 'Wavelength (um), Total_AAOD (unitless or 1/um)'
      do iw=1,RIN%nw
        write(iu,'(f14.5)',advance='no') RIN%wave(iw)
        temp(1:npixels) = GOUT%aerosol%opt%pixel(1:npixels)%wl(iw)%aextt
        call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
      enddo ! iw

      if(RIN%NSD .gt. 1) then
        do isd=1,RIN%NSD    ! aerosol component loop
          if(isd .eq. 1) write(iu,'(a)') 'Wavelength (um), Fine_AAOD (unitless or 1/um)'
          if(isd .eq. 2) write(iu,'(a)') 'Wavelength (um), Coarse_AAOD (unitless or 1/um)'
          do iw=1,RIN%nw
            write(iu,'(f14.5)',advance='no') RIN%wave(iw)
            temp(1:npixels) = GOUT%aerosol%opt%pixel(1:npixels)%wl(iw)%aext(isd)
            call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
          enddo ! iw
        enddo ! isd
      endif ! RIN%NSD .gt. 1
  
      return
      end subroutine print_absorption

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! Print Angstrom exponent              	  
      
      subroutine print_angstrom_exponent(iu,RIN,segment_meas,GOUT)

      use mod_par_inv,     only : KIMAGE
      use mod_retr_settings_derived_type
      use mod_retr_general_output_derived_type
      use mod_print_array
      use mod_sdata_derived_type
      
      implicit none  
!----------------------------------------------------------------------------------------
! IN :
      integer,                     intent(in)  ::  iu
      type(retr_input_settings),   intent(in)  ::  RIN
      type(segment_data),          intent(in)  ::  segment_meas
      type(output_segment_general),intent(in)  ::  GOUT
!----------------------------------------------------------------------------------------
! LOCAL :
      integer                              :: IDIM1,par_type
      real,dimension(KIMAGE)               :: temp
      integer                              :: npixels
!----------------------------------------------------------------------------------------
      temp(:) = 0.0
      npixels = segment_meas%npixels

      do IDIM1=1,RIN%NDIM%n1
        par_type = RIN%NDIM%par_type(IDIM1)        
        if(par_type .gt. par_type_SURF1_land_beg .and. par_type .lt. par_type_SURF_water_end) then      
          if(any(GOUT%aerosol%opt%pixel(1:npixels)%Aexp .ne. 0.0)) then
            write(iu,'(a,f6.3,a,f6.3,a)') 'Angstrom exp (',RIN%WAVE(RIN%Aexp_iwl(1)),  &
                                                  ' /',RIN%WAVE(RIN%Aexp_iwl(2)),' )'
            write(iu,'(14x)',advance='no') 
            temp(1:npixels) = GOUT%aerosol%opt%pixel(1:npixels)%Aexp
            call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)                
          endif
        exit
        endif        
      enddo ! IDIM1 

      return
      end subroutine print_angstrom_exponent 

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! Print Lidar Ratio              	  
      
      subroutine print_lidar_ratio(iu,RIN,segment_meas,GOUT)

      use mod_par_inv,     only : KIMAGE
      use mod_retr_settings_derived_type
      use mod_retr_general_output_derived_type
      use mod_print_array
      use mod_sdata_derived_type
      
      implicit none  
!----------------------------------------------------------------------------------------
! IN :
      integer,                     intent(in)  ::  iu
      type(retr_input_settings),   intent(in)  ::  RIN
      type(segment_data),          intent(in)  ::  segment_meas
      type(output_segment_general),intent(in)  ::  GOUT
!----------------------------------------------------------------------------------------
! LOCAL :
      integer                      ::  IDIM1,par_type
      real,dimension(KIMAGE)       ::  temp
      integer                      ::  iw,isd,npixels
!----------------------------------------------------------------------------------------
      temp(:) = 0.0
      npixels = segment_meas%npixels

      do IDIM1=1,RIN%NDIM%n1
        par_type = RIN%NDIM%par_type(IDIM1)        
        if(par_type .eq. par_type_AVP_prof) then      
          write(iu,'(a)') 'Wavelength (um),  Lidar Ratio (Total)'
          do iw=1,RIN%nw
            write(iu,'(f14.5)',advance='no') RIN%wave(iw)
            temp(1:npixels) = GOUT%aerosol%lidar%pixel(1:npixels)%wl(iw)%lrt
            call rprint_array_real (iu,temp(:), '(f14.5)', 1, npixels)
          enddo ! iw
          if(RIN%NSD .gt. 1) then
            do isd=1,RIN%NSD    ! aerosol component loop
              if(isd .eq. 1) write(iu,'(a)') 'Wavelength (um), Lidar Ratio (Fine)'
              if(isd .eq. 2) write(iu,'(a)') 'Wavelength (um), Lidar Ratio (Coarse)'
              do iw=1,RIN%nw
                write(iu,'(f14.5)',advance='no') RIN%wave(iw)
                temp(1:npixels) = GOUT%aerosol%lidar%pixel(1:npixels)%wl(iw)%lr(isd)
                call rprint_array_real (iu,temp(:), '(f14.5)', 1, npixels)
              enddo ! iw
            enddo ! isd
          endif ! RIN%NSD .gt. 1
        exit
        endif ! par_type .eq. par_type_CL
      enddo ! IDIM1
  
      return
      end subroutine print_lidar_ratio
! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! Print land percent              	  
      
      subroutine print_surface_land_percent(iu,RIN,segment_meas,GOUT)

      use mod_par_inv,     only : KIMAGE
      use mod_retr_settings_derived_type
      use mod_retr_general_output_derived_type
      use mod_print_array
      use mod_sdata_derived_type
      
      implicit none  
!----------------------------------------------------------------------------------------
! IN :
      integer,                     intent(in)  ::  iu
      type(retr_input_settings),   intent(in)  ::  RIN
      type(segment_data),          intent(in)  ::  segment_meas
      type(output_segment_general),intent(in)  ::  GOUT
!----------------------------------------------------------------------------------------
! LOCAL :
      integer                              :: IDIM1,par_type
      real,dimension(KIMAGE)               :: temp
      integer                              :: npixels
!----------------------------------------------------------------------------------------
      temp(:) = 0.0
      npixels = segment_meas%npixels

      do IDIM1=1,RIN%NDIM%n1
        par_type = RIN%NDIM%par_type(IDIM1)        
        if(par_type .gt. par_type_SURF1_land_beg .and. par_type .lt. par_type_SURF_water_end) then      
          write(iu,'(/,a)') '*** SURFACE ***'
          call print_date_time(iu,segment_meas)
          call print_coordinates(iu,segment_meas)
          write(iu,*)
          write(iu,'(a)',advance='no') 'Land percent:'
          temp(1:npixels) = segment_meas%pixels(1:npixels)%land_percent
          call rprint_array_real (iu,temp(:), '(f14.5)', 1, npixels)          
        exit
        endif        
      enddo ! IDIM1 

      return
      end subroutine print_surface_land_percent 

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! Print PM
      
      subroutine print_aerosol_PM(iu,RIN,segment_meas,GOUT)

      use mod_par_inv,     only : KIMAGE
      use mod_retr_settings_derived_type
      use mod_retr_general_output_derived_type
      use mod_print_array
      use mod_sdata_derived_type
      
      implicit none  
!----------------------------------------------------------------------------------------
! IN :
      integer,                     intent(in)  ::  iu
      type(retr_input_settings),   intent(in)  ::  RIN
      type(segment_data),          intent(in)  ::  segment_meas
      type(output_segment_general),intent(in)  ::  GOUT
!----------------------------------------------------------------------------------------
! LOCAL :
      real,dimension(KIMAGE)               :: temp
      integer                              :: npixels, ipm
!----------------------------------------------------------------------------------------
      temp(:) = 0.0
      npixels = segment_meas%npixels
! AL to do: get Pm from GOUT structure not calculate once again!!!!
      do ipm=1,RIN%nPM_diam
        if(RIN%PM_diam(ipm) .gt. 0.0) then
          write(iu,'(a,f4.1)',advance='no') 'PM ',RIN%PM_diam(ipm)
          write(iu,*)
!          call get_PM(RIN%PM_diam(ipm),RIN,GOUT,segment_meas,temp(1:npixels))
          temp(1:npixels) = GOUT%aerosol%pm%pixel(1:npixels)%PM(ipm)
          call rprint_array_real (iu,temp(:), '(f14.5)', 1, npixels)
        endif ! PM_diam(ipm) .gt. 0.0)
      enddo !ipm
!      write(iu,'(a)',advance='no') 'PM10: '
!      write(iu,*) ''
!      call get_PM(10.0,RIN,GOUT,segment_meas,temp(1:npixels))
!      call rprint_array_real (iu,temp(:), '(f14.5)', 1, npixels)
      
      return
      end subroutine print_aerosol_PM


! Print AEROSOL TYPES
      
      subroutine print_aerosol_types(iu,RIN,segment_meas,GOUT)

      use mod_par_inv,     only : KIMAGE
      use mod_retr_settings_derived_type
      use mod_retr_general_output_derived_type
      use mod_print_array
      use mod_sdata_derived_type
      
      implicit none  
!----------------------------------------------------------------------------------------
! IN :
      integer,                     intent(in)  ::  iu
      type(retr_input_settings),   intent(in)  ::  RIN
      type(segment_data),          intent(in)  ::  segment_meas
      type(output_segment_general),intent(in)  ::  GOUT
!----------------------------------------------------------------------------------------
! LOCAL :
      integer,dimension(KIMAGE)            :: temp ! will store indicies
      integer                              :: npixels, ipix
      character(len=12),dimension(0:8)     :: type_name ! strings containing names of aerosol types
!----------------------------------------------------------------------------------------
! initialize dictionaty
      type_name(0)='Complex_mix '
      type_name(1)='Background  '
      type_name(2)='Maritime    '
      type_name(3)='Urban:Poltd '
      type_name(4)='Mixed       '
      type_name(5)='Urban:Clean '
      type_name(6)='Smoke:smold '
      type_name(7)='Smoke:flame '
      type_name(8)='Miner._dust '

      temp(:) = 0
      npixels = segment_meas%npixels
      write(iu,'(a)',advance='no') 'Aerosol type: '
!      write(iu,*)
      do ipix=1,npixels
         write(iu,'(a14)', advance='no') type_name(GOUT%aerosol%types%pixel(ipix)%index)
      enddo
      write(iu,*)
!      write(iu,'(a)',advance='no') 'Aerosol type index:'
!      write(iu,*) ''
!      temp(1:npixels) = GOUT%aerosol%types%pixel(1:npixels)%index
!      call iprint_array_int (iu, temp(:), '(I2)', 1, npixels)
      
      return
      end subroutine print_aerosol_types
! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! Print Surface Abledo              	  
      
      subroutine print_surface_albedo(iu,RIN,segment_meas,GOUT)

      use mod_par_inv,     only : KIMAGE
      use mod_retr_settings_derived_type
      use mod_retr_general_output_derived_type
      use mod_print_array
      use mod_sdata_derived_type
      
      implicit none  
!----------------------------------------------------------------------------------------
! IN :
      integer,                     intent(in)  ::  iu
      type(retr_input_settings),   intent(in)  ::  RIN
      type(segment_data),          intent(in)  ::  segment_meas
      type(output_segment_general),intent(in)  ::  GOUT
!----------------------------------------------------------------------------------------
! LOCAL :
      integer                              :: IDIM1,par_type
      real,dimension(KIMAGE)               :: temp
      integer                              :: iw,npixels
!----------------------------------------------------------------------------------------
      temp(:) = 0.0
      npixels = segment_meas%npixels

      do IDIM1=1,RIN%NDIM%n1
        par_type = RIN%NDIM%par_type(IDIM1)        
        if(par_type .gt. par_type_SURF1_land_beg .and. par_type .lt. par_type_SURF_water_end) then                
          write(iu,'(a)') 'Wavelength (um), Surface ALBEDO'
          do iw=1,RIN%nw
            write(iu,'(f14.5)',advance='no') RIN%wave(iw)
            temp(1:npixels) = GOUT%surface%pixel(1:npixels)%wl(iw)%salbedo
            call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)
          enddo ! IW
        exit
        endif        
      enddo ! IDIM1 
      
      return
      end subroutine print_surface_albedo 

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! Print Surface NDVI              	  
      
      subroutine print_surface_ndvi(iu,RIN,segment_meas,GOUT)

      use mod_par_inv,     only : KIMAGE
      use mod_retr_settings_derived_type
      use mod_retr_general_output_derived_type
      use mod_print_array
      use mod_sdata_derived_type
      
      implicit none  
!----------------------------------------------------------------------------------------
! IN :
      integer,                     intent(in)  ::  iu
      type(retr_input_settings),   intent(in)  ::  RIN
      type(segment_data),          intent(in)  ::  segment_meas
      type(output_segment_general),intent(in)  ::  GOUT
!----------------------------------------------------------------------------------------
! LOCAL :
      integer                              :: IDIM1,par_type
      real,dimension(KIMAGE)               :: temp
      integer                              :: npixels
!----------------------------------------------------------------------------------------
      temp(:) = 0.0
      npixels = segment_meas%npixels

      do IDIM1=1,RIN%NDIM%n1
        par_type = RIN%NDIM%par_type(IDIM1)        
        if(par_type .gt. par_type_SURF1_land_beg .and. par_type .lt. par_type_SURF_water_end) then      
          if(any(GOUT%surface%pixel(1:npixels)%ndvi .ne. 0.0)) then
            write(iu,'(a,f6.3,a,f6.3,a)') 'Surface NDVI (',RIN%WAVE(RIN%ndvi_iwl(1)),  &
                                                  ' /',RIN%WAVE(RIN%ndvi_iwl(2)),' )'

            write(iu,'(14x)',advance='no') 
            temp(1:npixels) = GOUT%surface%pixel(1:npixels)%ndvi
            call rprint_array_real (iu,temp(:), '(e14.5)', 1, npixels)                
          endif
        exit
        endif        
      enddo ! IDIM1 
      
      return
      end subroutine print_surface_ndvi 

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine print_fitting ( iu_main_output, RIN, segment_meas, segment_fit )
	  
      use mod_par_inv, only : KW, KNBVM, NBVM
      use mod_sdata_meas_type
      use mod_sdata_derived_type
      use mod_inversion_utils
      use mod_retr_settings_derived_type
      use mod_sdata, only : get_pixel_wl, get_vert_prof_h
            										      
      implicit none  

!---------------------------------------------------------------------------
! IN :
      integer, intent(in) :: iu_main_output
      type(retr_input_settings), intent(in) :: RIN
      type(segment_data), intent(inout) :: segment_meas
      type(segment_data), intent(in) :: segment_fit

!---------------------------------------------------------------------------
! LOCAL :
      integer ::  nw, npixels
      real, dimension(KW) :: wave
      integer :: ipix, iw, ip, iv
      real, dimension(KNBVM) :: temp_meas, temp_fit, &
                                I_meas, I_fit, Q_meas, Q_fit, &
                                U_meas, U_fit, P11_meas, P11_fit     
      character(len=8) :: meas_type_char
      real :: sza
      real, dimension(NBVM) :: vis,fis                     
      real, dimension(KNBVM) :: arg
      integer :: nmeas_type, meas_type, nvalid_meas 
      logical :: status_funct 
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
! iPOBS = 
!           1    I, Q, U                      
!           2    I, Q/I, U/I                  
!           3    I, P    or sqrt(Q*Q+U*U)
!           4    I, P/I  or sqrt(Q*Q+U*U)/I
!
!---------------------------------------------------------------------------

      npixels = segment_meas%npixels
      write(iu_main_output,*)
      write(iu_main_output,'(a)') '*** FITTING ***'
loop_pixel: do ipix=1,npixels
      call get_pixel_wl ( segment_meas%pixels(ipix), nw, wave )      
loop_wl: do iw=1,nw
      write(iu_main_output,'(a)') '------------------------------------------------------------------------'
      write(iu_main_output,'(a,i4,a,i3,f10.3,a)') 'pixel # ',ipix,'   wavelength # ',iw,wave(iw),' (um)'
      write(iu_main_output,'(a)') '------------------------------------------------------------------------'
      nmeas_type = segment_meas%pixels(ipix)%meas(iw)%NIP
loop_meas_type: do ip=1,nmeas_type 
      meas_type_char(:) = ' '
      meas_type   = segment_meas%pixels(ipix)%meas(iw)%meas_type(ip)
      if(RIN%iPOBS .ge. 3 .and. meas_type .eq. meas_type_U) &
      exit loop_meas_type
      nvalid_meas = segment_meas%pixels(ipix)%meas(iw)%NBVM(ip)
      if(meas_type .ge. meas_type_LS .and. meas_type .le. meas_type_RL) then
        call get_vert_prof_h ( iw, ip, segment_meas%pixels(ipix), &
                                    segment_meas%pixels(ipix)%meas(IW)%NBVM(ip), &
                                    arg(:) )
      else
        do iv=1,nvalid_meas
          sza = segment_meas%pixels(ipix)%meas(iw)%sza
          vis(iv) = segment_meas%pixels(ipix)%meas(iw)%thetav(iv,ip)
          fis(iv) = segment_meas%pixels(ipix)%meas(iw)%phi(iv,ip)
          status_funct = geom2scat_angl ( sza,  vis(iv), fis(iv), arg(iv) )
        enddo
      endif
      select case(meas_type)
! Phase matrix
      case(meas_type_tod,meas_type_aod)
        if(meas_type .eq. meas_type_tod) meas_type_char = 'tod'
        if(meas_type .eq. meas_type_aod) meas_type_char = 'aod'
        temp_meas(1:nvalid_meas) = segment_meas%pixels(ipix)%meas(iw)%tau(1:nvalid_meas)
        temp_fit(1:nvalid_meas) = segment_fit%pixels(ipix)%meas(iw)%tau(1:nvalid_meas)
        call print_meas_fit ( iu_main_output, meas_type, meas_type_char, nvalid_meas, &
                              sza, vis(1:nvalid_meas), fis(1:nvalid_meas), arg(1:nvalid_meas), &
                              temp_meas(1:nvalid_meas), temp_fit(1:nvalid_meas) )
      case(meas_type_p11)
        meas_type_char = 'p11'
        temp_meas(1:nvalid_meas) = segment_meas%pixels(ipix)%meas(iw)%p11(1:nvalid_meas)
        temp_fit(1:nvalid_meas) = segment_fit%pixels(ipix)%meas(iw)%p11(1:nvalid_meas)
        call print_meas_fit ( iu_main_output, meas_type, meas_type_char, nvalid_meas, &
                              sza, vis(1:nvalid_meas), fis(1:nvalid_meas), arg(1:nvalid_meas), &
                              temp_meas(1:nvalid_meas), temp_fit(1:nvalid_meas) )
      case(meas_type_p12)  
        temp_meas(1:nvalid_meas) = segment_meas%pixels(ipix)%meas(iw)%p12(1:nvalid_meas)
        temp_fit(1:nvalid_meas) = segment_fit%pixels(ipix)%meas(iw)%p12(1:nvalid_meas)
        if(RIN%iPOBS .eq. 1 ) then
          meas_type_char = 'p12'
        elseif(RIN%iPOBS .eq. 2) then
          meas_type_char = '-p12/p11'
          P11_meas(1:nvalid_meas) = segment_meas%pixels(ipix)%meas(iw)%p11(1:nvalid_meas)
          P11_fit(1:nvalid_meas) = segment_fit%pixels(ipix)%meas(iw)%p11(1:nvalid_meas)
          temp_meas(1:nvalid_meas) = - temp_meas(1:nvalid_meas) / P11_meas(1:nvalid_meas)
          temp_fit(1:nvalid_meas) = - temp_fit(1:nvalid_meas) / P11_fit(1:nvalid_meas)
        elseif(RIN%iPOBS .eq. 5 ) then
          meas_type_char = '-p12/p11'
        endif
        call print_meas_fit ( iu_main_output, meas_type, meas_type_char, nvalid_meas, &
                              sza, vis(1:nvalid_meas), fis(1:nvalid_meas), arg(1:nvalid_meas), &
                              temp_meas(1:nvalid_meas), temp_fit(1:nvalid_meas) )
      case(meas_type_p22)
        temp_meas(1:nvalid_meas) = segment_meas%pixels(ipix)%meas(iw)%p22(1:nvalid_meas)
        temp_fit(1:nvalid_meas) = segment_fit%pixels(ipix)%meas(iw)%p22(1:nvalid_meas)
        if(RIN%iPOBS .eq. 1) then
          meas_type_char = 'p22'
        elseif(RIN%iPOBS .eq. 2) then
          meas_type_char = 'p22/p11'
          P11_meas(1:nvalid_meas) = segment_meas%pixels(ipix)%meas(iw)%p11(1:nvalid_meas)
          P11_fit(1:nvalid_meas) = segment_fit%pixels(ipix)%meas(iw)%p11(1:nvalid_meas)
          temp_meas(1:nvalid_meas) = temp_meas(1:nvalid_meas) / P11_meas(1:nvalid_meas)
          temp_fit(1:nvalid_meas) = temp_fit(1:nvalid_meas) / P11_fit(1:nvalid_meas)
        endif
        call print_meas_fit ( iu_main_output, meas_type, meas_type_char, nvalid_meas, &
                              sza, vis(1:nvalid_meas), fis(1:nvalid_meas), arg(1:nvalid_meas), &
                              temp_meas(1:nvalid_meas), temp_fit(1:nvalid_meas) )
      case(meas_type_p33)
        temp_meas(1:nvalid_meas) = segment_meas%pixels(ipix)%meas(iw)%p33(1:nvalid_meas)
        temp_fit(1:nvalid_meas) = segment_fit%pixels(ipix)%meas(iw)%p33(1:nvalid_meas)
        if(RIN%iPOBS .eq. 1) then
          meas_type_char = 'p33'
        elseif(RIN%iPOBS .eq. 2) then
          meas_type_char = 'p33/p11'
          P11_meas(1:nvalid_meas) = segment_meas%pixels(ipix)%meas(iw)%p11(1:nvalid_meas)
          P11_fit(1:nvalid_meas) = segment_fit%pixels(ipix)%meas(iw)%p11(1:nvalid_meas)
          temp_meas(1:nvalid_meas) = temp_meas(1:nvalid_meas) / P11_meas(1:nvalid_meas)
          temp_fit(1:nvalid_meas) = temp_fit(1:nvalid_meas) / P11_fit(1:nvalid_meas)
        endif
        call print_meas_fit ( iu_main_output, meas_type, meas_type_char, nvalid_meas, &
                              sza, vis(1:nvalid_meas), fis(1:nvalid_meas), arg(1:nvalid_meas), &
                              temp_meas(1:nvalid_meas), temp_fit(1:nvalid_meas) )
      case(meas_type_p34)
        temp_meas(1:nvalid_meas) = segment_meas%pixels(ipix)%meas(iw)%p34(1:nvalid_meas)
        temp_fit(1:nvalid_meas) = segment_fit%pixels(ipix)%meas(iw)%p34(1:nvalid_meas)
        if(RIN%iPOBS .eq. 1) then
          meas_type_char = 'p34'
        elseif(RIN%iPOBS .eq. 2) then
          meas_type_char = 'p34/p11'
          P11_meas(1:nvalid_meas) = segment_meas%pixels(ipix)%meas(iw)%p11(1:nvalid_meas)
          P11_fit(1:nvalid_meas) = segment_fit%pixels(ipix)%meas(iw)%p11(1:nvalid_meas)
          temp_meas(1:nvalid_meas) = temp_meas(1:nvalid_meas) / P11_meas(1:nvalid_meas)
          temp_fit(1:nvalid_meas) = temp_fit(1:nvalid_meas) / P11_fit(1:nvalid_meas)
        endif
        call print_meas_fit ( iu_main_output, meas_type, meas_type_char, nvalid_meas, &
                              sza, vis(1:nvalid_meas), fis(1:nvalid_meas), arg(1:nvalid_meas), &
                              temp_meas(1:nvalid_meas), temp_fit(1:nvalid_meas) )
      case(meas_type_p44)
        temp_meas(1:nvalid_meas) = segment_meas%pixels(ipix)%meas(iw)%p44(1:nvalid_meas)
        temp_fit(1:nvalid_meas) = segment_fit%pixels(ipix)%meas(iw)%p44(1:nvalid_meas)
        if(RIN%iPOBS .eq. 1) then
          meas_type_char = 'p44'
        elseif(RIN%iPOBS .eq. 2) then
          meas_type_char = 'p44/p11'
          P11_meas(1:nvalid_meas) = segment_meas%pixels(ipix)%meas(iw)%p11(1:nvalid_meas)
          P11_fit(1:nvalid_meas) = segment_fit%pixels(ipix)%meas(iw)%p11(1:nvalid_meas)
          temp_meas(1:nvalid_meas) = temp_meas(1:nvalid_meas) / P11_meas(1:nvalid_meas)
          temp_fit(1:nvalid_meas) = temp_fit(1:nvalid_meas) / P11_fit(1:nvalid_meas)
        endif
        call print_meas_fit ( iu_main_output, meas_type, meas_type_char, nvalid_meas, &
                              sza, vis(1:nvalid_meas), fis(1:nvalid_meas), arg(1:nvalid_meas), &
                              temp_meas(1:nvalid_meas), temp_fit(1:nvalid_meas) )
! Lidar
      case(meas_type_LS)  
        meas_type_char = 'LS'
        temp_meas(1:nvalid_meas) = segment_meas%pixels(ipix)%meas(iw)%LS(1:nvalid_meas)
        temp_fit(1:nvalid_meas) = segment_fit%pixels(ipix)%meas(iw)%LS(1:nvalid_meas)
        call print_meas_fit ( iu_main_output, meas_type, meas_type_char, nvalid_meas, &
                              sza, vis(1:nvalid_meas), fis(1:nvalid_meas), arg(1:nvalid_meas), &
                              temp_meas(1:nvalid_meas), temp_fit(1:nvalid_meas) )
      case(meas_type_DP)
        meas_type_char = 'DP'
      case(meas_type_RL)
        meas_type_char = 'RL'
! Stokes vector
      case(meas_type_I)
        temp_meas(1:nvalid_meas) = segment_meas%pixels(ipix)%meas(iw)%I(1:nvalid_meas)
        temp_fit(1:nvalid_meas) = segment_fit%pixels(ipix)%meas(iw)%I(1:nvalid_meas)
        if(RIN%iIOBS .eq. 1) then
          meas_type_char = 'I'
        elseif(RIN%iIOBS .eq. 2) then
          meas_type_char = 'I/sum'
          temp_meas(1:nvalid_meas) = temp_meas(1:nvalid_meas) / sum(temp_meas(1:nvalid_meas))
          temp_fit(1:nvalid_meas) = temp_fit(1:nvalid_meas) / sum(temp_fit(1:nvalid_meas))
        endif
        call print_meas_fit ( iu_main_output, meas_type, meas_type_char, nvalid_meas, &
                              sza, vis(1:nvalid_meas), fis(1:nvalid_meas), arg(1:nvalid_meas), &
                              temp_meas(1:nvalid_meas), temp_fit(1:nvalid_meas) )
      case(meas_type_Q)
        select case(RIN%iPOBS)
        case(1)
          meas_type_char = 'Q'
          temp_meas(1:nvalid_meas) = segment_meas%pixels(ipix)%meas(iw)%Q(1:nvalid_meas)
          temp_fit(1:nvalid_meas) = segment_fit%pixels(ipix)%meas(iw)%Q(1:nvalid_meas)
        case(2)
          meas_type_char = 'Q/I'
          I_meas(1:nvalid_meas) = segment_meas%pixels(ipix)%meas(iw)%I(1:nvalid_meas)
          I_fit(1:nvalid_meas) = segment_fit%pixels(ipix)%meas(iw)%I(1:nvalid_meas)
          temp_meas(1:nvalid_meas) = segment_meas%pixels(ipix)%meas(iw)%Q(1:nvalid_meas) / I_meas(1:nvalid_meas)
          temp_fit(1:nvalid_meas) = segment_fit%pixels(ipix)%meas(iw)%Q(1:nvalid_meas) / I_fit(1:nvalid_meas)
        case(3)
          meas_type_char = 'P'
          Q_meas(1:nvalid_meas) = segment_meas%pixels(ipix)%meas(iw)%Q(1:nvalid_meas)
          Q_fit(1:nvalid_meas) = segment_fit%pixels(ipix)%meas(iw)%Q(1:nvalid_meas)
          U_meas(1:nvalid_meas) = segment_meas%pixels(ipix)%meas(iw)%U(1:nvalid_meas)
          U_fit(1:nvalid_meas) = segment_fit%pixels(ipix)%meas(iw)%U(1:nvalid_meas)
          temp_meas(1:nvalid_meas) = sqrt( Q_meas(1:nvalid_meas)*Q_meas(1:nvalid_meas) + &
                                           U_meas(1:nvalid_meas)*U_meas(1:nvalid_meas) )         
          temp_fit(1:nvalid_meas) = sqrt( Q_fit(1:nvalid_meas)*Q_fit(1:nvalid_meas) + &
                                          U_fit(1:nvalid_meas)*U_fit(1:nvalid_meas) )         
        case(4)
          meas_type_char = 'P/I'
          I_meas(1:nvalid_meas) = segment_meas%pixels(ipix)%meas(iw)%I(1:nvalid_meas)
          I_fit(1:nvalid_meas) = segment_fit%pixels(ipix)%meas(iw)%I(1:nvalid_meas)
          Q_meas(1:nvalid_meas) = segment_meas%pixels(ipix)%meas(iw)%Q(1:nvalid_meas)
          Q_fit(1:nvalid_meas) = segment_fit%pixels(ipix)%meas(iw)%Q(1:nvalid_meas)
          U_meas(1:nvalid_meas) = segment_meas%pixels(ipix)%meas(iw)%U(1:nvalid_meas)
          U_fit(1:nvalid_meas) = segment_fit%pixels(ipix)%meas(iw)%U(1:nvalid_meas)
          temp_meas(1:nvalid_meas) = sqrt( Q_meas(1:nvalid_meas)*Q_meas(1:nvalid_meas) + &
                                           U_meas(1:nvalid_meas)*U_meas(1:nvalid_meas) ) / I_meas(1:nvalid_meas)         
          temp_fit(1:nvalid_meas) = sqrt( Q_fit(1:nvalid_meas)*Q_fit(1:nvalid_meas) + &
                                          U_fit(1:nvalid_meas)*U_fit(1:nvalid_meas) ) / I_fit(1:nvalid_meas)         
        end select
        call print_meas_fit ( iu_main_output, meas_type, meas_type_char, nvalid_meas, &
                              sza, vis(1:nvalid_meas), fis(1:nvalid_meas), arg(1:nvalid_meas), &
                              temp_meas(1:nvalid_meas), temp_fit(1:nvalid_meas) )
      case(meas_type_U)
        temp_meas(1:nvalid_meas) = segment_meas%pixels(ipix)%meas(iw)%U(1:nvalid_meas)
        temp_fit(1:nvalid_meas) = segment_fit%pixels(ipix)%meas(iw)%U(1:nvalid_meas)
        select case(RIN%iPOBS)
        case(1)
          meas_type_char = 'U'
        case(2)
          meas_type_char = 'U/I'
          I_meas(1:nvalid_meas) = segment_meas%pixels(ipix)%meas(iw)%I(1:nvalid_meas)
          I_fit(1:nvalid_meas) = segment_fit%pixels(ipix)%meas(iw)%I(1:nvalid_meas)
          temp_meas(1:nvalid_meas) = temp_meas(1:nvalid_meas) / I_meas(1:nvalid_meas)
          temp_fit(1:nvalid_meas) = temp_fit(1:nvalid_meas) / I_meas(1:nvalid_meas)
        end select
        call print_meas_fit ( iu_main_output, meas_type, meas_type_char, nvalid_meas, &
                              sza, vis(1:nvalid_meas), fis(1:nvalid_meas), arg(1:nvalid_meas), &
                              temp_meas(1:nvalid_meas), temp_fit(1:nvalid_meas) )
      case(meas_type_P)
        temp_meas(1:nvalid_meas) = segment_meas%pixels(ipix)%meas(iw)%P(1:nvalid_meas)
        temp_fit(1:nvalid_meas) = segment_fit%pixels(ipix)%meas(iw)%P(1:nvalid_meas)
        select case(RIN%iPOBS)
        case(3)
          meas_type_char = 'P'
        case(5)
! iPOBS=5 In case when number of angeles or angles for I differe from ones for P
! linear polarization has been divided by I at the time of setting 
! modeled measurements in subroutine set_pixel_Stokes_vec_fit
          meas_type_char = 'P/I'
        case(4)
          meas_type_char = 'P/I'
          I_meas(1:nvalid_meas) = segment_meas%pixels(ipix)%meas(iw)%I(1:nvalid_meas)
          I_fit(1:nvalid_meas) = segment_fit%pixels(ipix)%meas(iw)%I(1:nvalid_meas)
          temp_meas(1:nvalid_meas) = temp_meas(1:nvalid_meas) / I_meas(1:nvalid_meas)
          temp_fit(1:nvalid_meas) = temp_fit(1:nvalid_meas) / I_fit(1:nvalid_meas)
        end select
        call print_meas_fit ( iu_main_output, meas_type, meas_type_char, nvalid_meas, &
                              sza, vis(1:nvalid_meas), fis(1:nvalid_meas), arg(1:nvalid_meas), &
                              temp_meas(1:nvalid_meas), temp_fit(1:nvalid_meas) )
      end select
enddo loop_meas_type
enddo loop_wl
enddo loop_pixel

      end subroutine print_fitting
! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine print_fitting_FS ( iu_main_output, RIN, segment_meas, &
                                    segment_vec_meas, segment_vec_fit )
	  
      use mod_par_inv, only : KW, KNBVM, NBVM
      use mod_sdata_meas_type
      use mod_sdata_derived_type
      use mod_inversion_utils
      use mod_retr_settings_derived_type
      use mod_sdata, only : get_pixel_wl, get_vert_prof_h
            										      
      implicit none  

!---------------------------------------------------------------------------
! IN :
      integer, intent(in) :: iu_main_output
      type(retr_input_settings), intent(in) :: RIN
      type(segment_data), intent(inout) :: segment_meas
      type(pixel_vector),dimension(KIMAGE), intent(in) :: segment_vec_meas
      type(pixel_vector),dimension(KIMAGE), intent(in) :: segment_vec_fit


!---------------------------------------------------------------------------
! LOCAL :
      integer ::  nw, npixels
      real, dimension(KW) :: wave
      integer :: ipix, iw, ip, iv, IWB, JJS
      real, dimension(KNBVM) :: temp_meas, temp_fit
      character(len=8) :: meas_type_char
      real :: sza
      real, dimension(NBVM) :: vis,fis                     
      real, dimension(KNBVM) :: arg
      integer :: nmeas_type, meas_type, nvalid_meas 
      logical :: status_funct 
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
! iPOBS = 1    I, Q, U
!         2    I, Q/I, U/I
!         3    I, P    or sqrt(Q*Q+U*U)
!         4    I, P/I  or sqrt(Q*Q+U*U)/I
!         5    I, P/I  P/I provided as measurements in input
!
! iIOBS = 1    I
!       = 2    I/sum(I(:))
!---------------------------------------------------------------------------
      npixels = segment_meas%npixels
      write(iu_main_output,*)
      write(iu_main_output,'(a)') '*** FITTING from meas vector ***'
loop_pixel: do ipix=1,npixels
      call get_pixel_wl ( segment_meas%pixels(ipix), nw, wave )      
loop_wl: do iw=1,nw
      write(iu_main_output,'(a)') '--------------------------------------------------------------------------'
      write(iu_main_output,'(a,i4,a,i3,f10.3,a)') 'pixel # ',ipix,'   wavelength # ',iw,wave(iw),' (um)'
      write(iu_main_output,'(a)') '--------------------------------------------------------------------------'
      nmeas_type = segment_meas%pixels(ipix)%meas(iw)%NIP
      IWB=SUM(segment_vec_meas(ipix)%nFS(1:iw))-segment_vec_meas(ipix)%nFS(iw)
      ! IWD - number of elements in FS before current (iw) wavelength
      JJS=IWB
loop_meas_type: do ip=1,nmeas_type
      meas_type_char(:) = ' '
      meas_type   = segment_meas%pixels(ipix)%meas(iw)%meas_type(ip)
      if(RIN%iPOBS .ge. 3 .and. meas_type .eq. meas_type_U) &
      exit loop_meas_type
      nvalid_meas = segment_meas%pixels(ipix)%meas(iw)%NBVM(ip)
      do iv=1,nvalid_meas
        JJS = JJS+1
        temp_meas(iv) = exp(segment_vec_meas(ipix)%FS(JJS))
        temp_fit (iv) = exp(segment_vec_fit (ipix)%FS(JJS))
        if( (meas_type .ge. meas_type_Q   .and. meas_type .le. meas_type_P) .or. &
             meas_type .eq. meas_type_p12 .or.        &
            (meas_type .ge. meas_type_p33 .and. meas_type .le. meas_type_p44) &
          ) then
          temp_meas(iv) = temp_meas(iv) - RIN%SHIFT
          temp_fit (iv) = temp_fit (iv) - RIN%SHIFT
        endif
      enddo  ! iv
      if(meas_type .ge. meas_type_LS .and. meas_type .le. meas_type_RL) then
        call get_vert_prof_h ( iw, ip, segment_meas%pixels(ipix), &
                                    segment_meas%pixels(ipix)%meas(iw)%NBVM(ip), &
                                    arg(:) )
      else
        do iv=1,nvalid_meas
          sza = segment_meas%pixels(ipix)%meas(iw)%sza
          vis(iv) = segment_meas%pixels(ipix)%meas(iw)%thetav(iv,ip)
          fis(iv) = segment_meas%pixels(ipix)%meas(iw)%phi(iv,ip)
          status_funct = geom2scat_angl ( sza,  vis(iv), fis(iv), arg(iv) )
        enddo
      endif
      select case(meas_type)
! Phase matrix
      case(meas_type_tod)
        if(meas_type .eq. meas_type_tod) meas_type_char = 'tod'
      case(meas_type_aod)
        if(meas_type .eq. meas_type_aod) meas_type_char = 'aod'
      case(meas_type_p11)
        meas_type_char = 'p11'
      case(meas_type_p12)
        if(RIN%iPOBS .eq. 1) then
          meas_type_char = 'p12'
        elseif(RIN%iPOBS .eq. 2) then
          meas_type_char = '-p12/p11'
        elseif(RIN%iPOBS .eq. 5) then
          meas_type_char = '-p12/p11'
        endif
      case(meas_type_p22)
        if(RIN%iPOBS .eq. 1) then
          meas_type_char = 'p22'
        elseif(RIN%iPOBS .eq. 2) then
          meas_type_char = 'p22/p11'
        endif
      case(meas_type_p33)
        if(RIN%iPOBS .eq. 1) then
          meas_type_char = 'p33'
        elseif(RIN%iPOBS .eq. 2) then
          meas_type_char = 'p33/p11'
        endif
      case(meas_type_p34)
        if(RIN%iPOBS .eq. 1) then
          meas_type_char = 'p34'
        elseif(RIN%iPOBS .eq. 2) then
          meas_type_char = 'p34/p11'
        endif
      case(meas_type_p44)
        if(RIN%iPOBS .eq. 1) then
          meas_type_char = 'p44'
        elseif(RIN%iPOBS .eq. 2) then
          meas_type_char = 'p44/p11'
        endif
! Lidar
      case(meas_type_LS)  
        meas_type_char = 'LS'
      case(meas_type_DP)
        meas_type_char = 'DP'
      case(meas_type_RL)
        meas_type_char = 'RL'
! Stokes vector
      case(meas_type_I)
        if(RIN%iIOBS .eq. 1) then
          meas_type_char = 'I'
        elseif(RIN%iIOBS .eq. 2) then
          meas_type_char = 'I/sum'
        endif
      case(meas_type_Q)
        select case(RIN%iPOBS)
        case(1)
          meas_type_char = 'Q'
        case(2)
          meas_type_char = 'Q/I'
        case(3)
          meas_type_char = 'P'
        case(4)
          meas_type_char = 'P/I'
        end select
      case(meas_type_U)
        select case(RIN%iPOBS)
        case(1)
          meas_type_char = 'U'
        case(2)
          meas_type_char = 'U/I'
        end select
      case(meas_type_P)
        select case(RIN%iPOBS)
        case(3)
          meas_type_char = 'P'
        case(5)
! iPOBS=5 In case when number of angeles or angles for I differe from ones for P
! linear polarization has been divided by I at the time of setting 
! modeled measurements in subroutine set_pixel_Stokes_vec_fit
          meas_type_char = 'P/I'
        case(4)
          meas_type_char = 'P/I'
        end select
      end select
      call print_meas_fit ( iu_main_output, meas_type, meas_type_char, nvalid_meas, &
                            sza, vis(1:nvalid_meas), fis(1:nvalid_meas), arg(1:nvalid_meas), &
                            temp_meas(1:nvalid_meas), temp_fit(1:nvalid_meas) )
enddo loop_meas_type
enddo loop_wl
enddo loop_pixel

      end subroutine print_fitting_FS

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine print_meas_fit ( iu_main_output, meas_type, meas_type_char, & 
                                  nvalid_meas, sza, vis, fis, arg, temp_meas, temp_fit )	        										      
      use mod_sdata_meas_type
      use mod_sdata_derived_type
      
      implicit none  
!---------------------------------------------------------------------------
! IN :
      integer, intent(in) :: iu_main_output, nvalid_meas, meas_type
!      character(len=8), intent(in) :: meas_type_char
      character(*), intent(in) :: meas_type_char
      real, intent(in) :: sza
      real, dimension(nvalid_meas), intent(in) :: vis, fis, arg
      real, dimension(nvalid_meas), intent(in) :: temp_meas, temp_fit
!---------------------------------------------------------------------------
! LOCAL :
      integer ::  iv
!---------------------------------------------------------------------------
      if(meas_type .eq. meas_type_tod .or. meas_type .eq. meas_type_aod) then
        write(iu_main_output,'(2(7x,a))') &
        'meas_'//trim(meas_type_char),' fit_'//trim(meas_type_char)
        do iv=1,nvalid_meas
          write(iu_main_output,'(2e15.5)') temp_meas(iv),temp_fit(iv)
        enddo ! IV                        
      elseif(meas_type .ge. meas_type_LS .and. meas_type .le. meas_type_RL) then
        write(iu_main_output,'(3x,a,4x,a,7x,a,8x,a)') '#','hight[m]', &
        'meas_'//trim(meas_type_char),' fit_'//trim(meas_type_char)
        do iv=1,nvalid_meas
          write(iu_main_output,'(i4,f11.2,2e15.5)') iv,arg(iv),temp_meas(iv),temp_fit(iv)
        enddo ! IV                  
      elseif(meas_type .ge. meas_type_p11 .and. meas_type .le. meas_type_p44) then
        write(iu_main_output,'(3x,a,2x,a,2a15)') '#','sca_ang', &
        'meas_'//trim(meas_type_char),' fit_'//trim(meas_type_char)
        do iv=1,nvalid_meas
          write(iu_main_output,'(i4,f9.2,2e15.5)') iv,arg(iv),temp_meas(iv),temp_fit(iv)
        enddo ! IV                  
      elseif(meas_type .ge. meas_type_I .and. meas_type .le. meas_type_P) then
        write(iu_main_output,'(3x,a,3(6x,a),2x,a,2a15)') '#','sza','vis','fis','sca_ang', &
        'meas_'//trim(meas_type_char),' fit_'//trim(meas_type_char)
        do iv=1,nvalid_meas
          write(iu_main_output,'(i4,4f9.2,2e15.5)') &
          iv,sza,vis(iv),fis(iv),arg(iv),temp_meas(iv),temp_fit(iv)
        enddo ! IV                        
      endif
      
      end subroutine print_meas_fit

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine print_error_estimates ( iu_main_output, RIN, segment_meas, GOUT )

      use mod_retr_settings_derived_type
      use mod_retr_general_output_derived_type
      use mod_time_utils
      use mod_sdata_derived_type
      use mod_stop_report

      implicit none
!	------------------------------------------------------------------------------------------------------
! IN :
      integer, intent(in) :: iu_main_output
      type(retr_input_settings), intent(in) :: RIN
      type(segment_data), intent(in) :: segment_meas
      type(output_segment_general), intent(in) :: GOUT
!	------------------------------------------------------------------------------------------------------
! LOCAL :
      integer  ::  ipix,IW,IOF,ISD
      integer  ::  NW,npixels
      integer  ::  alloc_stat
      real,dimension(:,:,:),allocatable     :: TEMP
      character(len=255)                    :: opt_funct
      character(len=12),dimension(KIMAGE)   :: pdate,ptime 
!	------------------------------------------------------------------------------------------------------
! Print error estimates of retrieved parameter logarithms
      if( RIN%products%errest%par ) then
        call print_PAR_ERR_estimates( iu_main_output, RIN, segment_meas, GOUT%errest%par )
        if ( stop_report%status ) return
      endif

! Print error estimates of retrieved optical characteristic logarithms
      if( RIN%products%errest%aerosol%opt .or. RIN%products%errest%aerosol%lidar) then
        !write(iu_main_output,'(a,i2)') 'INVSING =',RIN%INVSING
        call print_OPT_ERR_estimates ( iu_main_output, RIN, segment_meas, GOUT )
        if ( stop_report%status ) return
      endif

      return
      end subroutine print_error_estimates

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine print_OPT_ERR_estimates ( iu_main_output, RIN, segment_meas, GOUT )
      use mod_retr_settings_derived_type
      use mod_retr_general_output_derived_type
      use mod_time_utils
      use mod_sdata_derived_type
      use mod_stop_report

      implicit none
!	------------------------------------------------------------------------------------------------------
! IN :
      integer, intent(in) :: iu_main_output
      type(retr_input_settings), intent(in) :: RIN
      type(segment_data), intent(in) :: segment_meas
      type(output_segment_general), intent(in) :: GOUT

! NOF  - number of optical functions for error estimates
! NSD1 - number of SD modes + 1(total)
!	------------------------------------------------------------------------------------------------------
! LOCAL :
      integer :: ipix,IW,IOF,ISD
      integer :: NW,npixels
      integer :: alloc_stat
      real,dimension(:,:,:), allocatable :: TEMP
      character(len=255) :: opt_funct
      character(len=12),dimension(KIMAGE) :: pdate, ptime
      integer :: NSD1, NOF_beg, NOF_end
! NOF  - number of optical functions for error estimates
! NSD1 - number of SD modes + 1(total)
!	------------------------------------------------------------------------------------------------------
      npixels = segment_meas%npixels
      NW = RIN%nw
      if( GOUT%products%errest%aerosol%opt .and. GOUT%products%errest%aerosol%lidar ) then
         NOF_beg = 1
         NOF_end = 3 ! ext, ssa, lr
      elseif( .not. GOUT%products%errest%aerosol%opt .and. GOUT%products%errest%aerosol%lidar ) then
         NOF_beg = 3
         NOF_end = 3 ! lr
      elseif( .not. GOUT%products%errest%aerosol%lidar .and. GOUT%products%errest%aerosol%opt ) then
         NOF_beg = 1
         NOF_end = 2 ! ext, ssa
      endif !
      select case(RIN%NSD)
      case(1)
         NSD1=RIN%NSD
      case(2)
         NSD1=RIN%NSD+1  ! fine, coarse, total    
      end select

      !write(*,*) npixels,NW,NOF,NSD1,'  npixels,NW,NOF,NSD1'

      call yyyymmdd_hhmmss ( segment_meas,pdate,ptime )
      
      allocate(TEMP(NW,NSD1,npixels),stat=alloc_stat)
      if (alloc_stat /= 0) stop 'error while trying to allocate TEMP in write_OPT_ERR'

! *** RANDOM
      write(iu_main_output,'(a)')    '--------------------------------------------------------------------------------------'
      if(RIN%KL .eq. 1) then
         write(iu_main_output,'(a)') 'Standard deviations of retrieved optical characteristic logarithms (~relative errors) :'
      else
         write(iu_main_output,'(a)') 'Standard deviations of retrieved optical characteristic  (absolute errors) :'      
      endif
      write(iu_main_output,'(a)')    '--------------------------------------------------------------------------------------'
      write(iu_main_output,'(a5,10x,20000(3x,A10,x))')  &
                  "Date:",(pdate(ipix),ipix=1,npixels)
      write(iu_main_output,'(a5,12x,20000(3x,A8,3x))')  & 
                  "Time:",(ptime(ipix),ipix=1,npixels)

      do IOF=NOF_beg,NOF_end
      opt_funct = ''
      select case(IOF)
      case(index_ext)
         opt_funct = 'Wavelength (um), Aerosol Optical Depth (Random)'
         do ipix=1,npixels
         if(NSD1 .eq. 3) then
            TEMP(1:NW,NSD1,ipix) = GOUT%errest%aerosol%opt%pixel(ipix)%ERR_extt(1:NW)
         endif
         do ISD=1,RIN%NSD
            TEMP(1:NW,ISD,ipix) = GOUT%errest%aerosol%opt%pixel(ipix)%ERR_ext(1:NW,ISD)
         enddo ! ISD
         enddo ! ipix
      case(index_ssa)
         opt_funct = 'Wavelength (um), Single Scattering Albedo (Random)'      
         do ipix=1,npixels
         if(NSD1 .eq. 3) then
            TEMP(1:NW,NSD1,ipix) = GOUT%errest%aerosol%opt%pixel(ipix)%ERR_ssat(1:NW)
         endif
         do ISD=1,RIN%NSD
            TEMP(1:NW,ISD,ipix) = GOUT%errest%aerosol%opt%pixel(ipix)%ERR_ssa(1:NW,ISD)
         enddo ! ISD
         enddo ! ipix 
      case(index_lr)
         opt_funct =  'Wavelength (um), Lidar Ratio (Random)'      
         do ipix=1,npixels
         if(NSD1 .eq. 3) then
            TEMP(1:NW,NSD1,ipix) = GOUT%errest%aerosol%lidar%pixel(ipix)%ERR_lrt(1:NW)
         endif
         do ISD=1,RIN%NSD
            TEMP(1:NW,ISD,ipix) = GOUT%errest%aerosol%lidar%pixel(ipix)%ERR_lr(1:NW,ISD)
         enddo ! ISD
         enddo ! ipix       
      case default
         write(tmp_message,'(a,i0,a)') 'IOF = ',IOF,'value is not valid'
         G_ERROR(trim(tmp_message))
      end select

      if(NSD1 .eq. 3) then
        write(iu_main_output,'(2a)') trim(opt_funct),',   Total '
        do IW=1,NW
          write(iu_main_output,'(f14.4,1000e14.5)')  &
                  RIN%wave(IW),(TEMP(IW,NSD1,ipix),ipix=1,npixels)
        enddo ! IW
      endif
      do ISD=1,RIN%NSD
        write(iu_main_output,'(2a,i3)') trim(opt_funct),',   fraction - ',ISD
        do IW=1,NW
          write(iu_main_output,'(f14.4,1000e14.5)')  &
                  RIN%wave(IW),(TEMP(IW,ISD,ipix),ipix=1,npixels)

        enddo ! IW
      enddo ! ISD      

      enddo ! IOF

! *** BIAS

      write(iu_main_output,'(a)') & 
      '-----------------------------------------------------------------------------------------------'
      if(RIN%KL .eq. 1) then
         write(iu_main_output,'(a)') &
         'BIAS - Standard deviations of systematic errors of retrieved optical characteristic logarithms :'
      else
         write(iu_main_output,'(a)') & 
         'BIAS - Standard deviations of systematic errors of retrieved optical characteristics :'
      endif
      write(iu_main_output,'(a)') &
      '-----------------------------------------------------------------------------------------------'
      write(iu_main_output,'(a5,10x,20000(3x,A10,x))')  &
                  "Date:",(pdate(ipix),ipix=1,npixels)
      write(iu_main_output,'(a5,12x,20000(3x,A8,3x))')  & 
                  "Time:",(ptime(ipix),ipix=1,npixels)

      do IOF=NOF_beg,NOF_end
      opt_funct = ''
      select case(IOF)
      case(index_ext)
         opt_funct = 'Wavelength (um), Aerosol Optical Depth (Bias)'
         do ipix=1,npixels
         if(NSD1 .eq. 3) then
            TEMP(1:NW,NSD1,ipix) = GOUT%errest%aerosol%opt%pixel(ipix)%BIAS_extt(1:NW)
         endif
         do ISD=1,RIN%NSD
            TEMP(1:NW,ISD,ipix) = GOUT%errest%aerosol%opt%pixel(ipix)%BIAS_ext(1:NW,ISD)
         enddo ! ISD
         enddo ! ipix 
      case(index_ssa)
         opt_funct = 'Wavelength (um), Single Scattering Albedo (Bias)'            
         do ipix=1,npixels
         if(NSD1 .eq. 3) then
            TEMP(1:NW,NSD1,ipix) = GOUT%errest%aerosol%opt%pixel(ipix)%BIAS_ssat(1:NW)
         endif
         do ISD=1,RIN%NSD
            TEMP(1:NW,ISD,ipix) = GOUT%errest%aerosol%opt%pixel(ipix)%BIAS_ssa(1:NW,ISD)
         enddo ! ISD
         enddo ! ipix 
      case(index_lr)
         opt_funct =  'Wavelength (um), Lidar Ratio (Bias)'      
         do ipix=1,npixels
         if(NSD1 .eq. 3) then
            TEMP(1:NW,NSD1,ipix) = GOUT%errest%aerosol%lidar%pixel(ipix)%BIAS_lrt(1:NW)
         endif
         do ISD=1,RIN%NSD
            TEMP(1:NW,ISD,ipix) = GOUT%errest%aerosol%lidar%pixel(ipix)%BIAS_lr(1:NW,ISD)
         enddo ! ISD
         enddo ! ipix       
      case default
         write(tmp_message,'(a,i0,a)') 'IOF = ',IOF,'value is not valid'
         G_ERROR(trim(tmp_message))
      end select

      if(NSD1 .eq. 3) then
        write(iu_main_output,'(2a)') trim(opt_funct),',   Total '
        do IW=1,NW
          write(iu_main_output,'(f14.4,1000e14.5)')  &
                  RIN%wave(IW),(TEMP(IW,NSD1,ipix),ipix=1,npixels)
        enddo ! IW
      endif
      do ISD=1,RIN%NSD
        write(iu_main_output,'(2a,i3)') trim(opt_funct),',   fraction - ',ISD
        do IW=1,NW
          write(iu_main_output,'(f14.4,1000e14.5)')  &
                  RIN%wave(IW),(TEMP(IW,ISD,ipix),ipix=1,npixels)
        enddo ! IW
      enddo ! ISD      

      enddo ! IOF

      deallocate(TEMP,stat=alloc_stat)
      if (alloc_stat /= 0) then
        write(tmp_message,'(a)') 'error while trying to deallocate TEMP'
        G_ERROR(trim(tmp_message))
      endif

      return
      end subroutine print_OPT_ERR_estimates

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      subroutine print_PAR_ERR_estimates( iu_main_output, RIN, segment_meas, GOUT_errest_par )

      use mod_retr_settings_derived_type
      use mod_retr_general_output_derived_type
      use mod_time_utils
      
      implicit none
!	------------------------------------------------------------------------------------------------------
! IN :
      integer, intent(in) :: iu_main_output
      type(retr_input_settings), intent(in) :: RIN
      type(segment_data), intent(in) :: segment_meas
      type(output_segment_err_estim_par), intent(in) :: GOUT_errest_par
!	------------------------------------------------------------------------------------------------------
! LOCAL :
      integer :: npixels, INVSING, KNSINGF
      integer :: ipix, i
      character(len=12), dimension(KIMAGE) :: pdate, ptime
!	------------------------------------------------------------------------------------------------------
      npixels = segment_meas%npixels
      INVSING = RIN%IPFP%INVSING
      KNSINGF = RIN%KNSINGF

      pdate(:) = ' '
      ptime(:) = ' '
      do ipix=1,npixels
         call convert_time_to_string(segment_meas%pixels(ipix)%t, "%F", pdate(ipix))
         call convert_time_to_string(segment_meas%pixels(ipix)%t, "%X", ptime(ipix))
      enddo 

      if(RIN%IPRI_verbose) then
      
         write(iu_main_output,'(a,i2)') ' INVSING =',INVSING
         write(iu_main_output,'(a)')    '-------------------------------------------------------------------------'
         if(RIN%KL .eq. 1) then
            write(iu_main_output,'(a)') 'Standard deviations of retrieved parameter logarithms (~relative errors) :'
         else
            write(iu_main_output,'(a)') 'Standard deviations of retrieved parameters  (absolute errors) :'      
         endif
         write(iu_main_output,'(a)')    '-------------------------------------------------------------------------'
         write(iu_main_output,'(a5,10x,20000(3x,A10,x))')  & 
                     "Date:",(pdate(ipix),ipix=1,npixels)
         write(iu_main_output,'(a5,12x,20000(3x,A8,3x))')  & 
                     "Time:",(ptime(ipix),ipix=1,npixels)

         do i=1,KNSINGF
         write(iu_main_output,'(i14,1000e14.5)')  &
                     i,(GOUT_errest_par%pixel(ipix)%ERRP(i),ipix=1,npixels)
         enddo ! i

         write(iu_main_output,'(a)')    '----------------------------------------------------------------------------------'      
         if(RIN%KL .eq. 1) then
            write(iu_main_output,'(a)') 'BIAS - Standard deviation of systematic errors of retrieved parameter logarithms :'
         else
            write(iu_main_output,'(a)') 'BIAS - Standard deviation of systematic errors of retrieved parameters :'
         endif
         write(iu_main_output,'(a)')    '----------------------------------------------------------------------------------'
         write(iu_main_output,'(a5,10x,20000(3x,A10,x))')  & 
                     "Date:",(pdate(ipix),ipix=1,npixels)
         write(iu_main_output,'(a5,12x,20000(3x,A8,3x))')  & 
                     "Time:",(ptime(ipix),ipix=1,npixels)

         do i=1,KNSINGF
         write(iu_main_output,'(i14,1000e14.5)')  &
                     i,(GOUT_errest_par%pixel(ipix)%BIASP(i),ipix=1,npixels)
         enddo ! i
      
      endif ! RIN%IPRI_verbose

      return
      end subroutine print_PAR_ERR_estimates

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! Print phase matrix
      subroutine print_phmx ( iu_main_output, RIN, segment_meas, GOUT )

      use mod_retr_settings_derived_type
      use mod_retr_general_output_derived_type
      use mod_time_utils
      use mod_sdata_derived_type
      	  
      implicit none
! -----------------------------------------------------------------------------------
      integer, intent(in) :: iu_main_output
      type(retr_input_settings),    intent(in) :: RIN
      type(output_segment_general), intent(in) :: GOUT
      type(segment_data),           intent(in) :: segment_meas
! -----------------------------------------------------------------------------------
      integer :: II, ipix, IW, ISD, IDIM1, NDIM3, npixels, par_type
      character(LEN=12), dimension(KIMAGE) :: pdate, ptime
      real,dimension(KIDIM3,KIDIM2,KIMAGE) :: sph
! -----------------------------------------------------------------------------------
      npixels = segment_meas%npixels

      do IDIM1=1,RIN%NDIM%n1
        par_type = RIN%NDIM%par_type(IDIM1)        
        if(par_type .gt. par_type_SHD_beg  .and. par_type .lt. par_type_SHD_end) then
          call unpack_retr_param_to_print ( RIN,GOUT%retrieval%par,IDIM1,npixels,sph )
          if(RIN%NDIM%n2(IDIM1) .eq. 1 .and. RIN%NDIM%n3(1,IDIM1) .eq. 2) & 
          sph(2,:,:) = sph(1,:,:)     	  
        exit
        endif ! par_type .gt. par_type_SHD_beg  .and. 
      enddo ! IDIM1

      write(iu_main_output,*)
      write(iu_main_output,*)	"  Phase Matrix"
      write(iu_main_output,*)

      call yyyymmdd_hhmmss ( segment_meas,pdate,ptime )

      select case (RIN%DLSF%keyEL)
      case(1)
        do ipix=1,segment_meas%npixels
          write(iu_main_output,'(a)') 	"--------------------------------------------------------------------------"
          write(iu_main_output,'(a,i0,a,a10,a,a8)') 'ipix=',ipix,'  yymmdd = ',pdate(ipix),'  hhmmss = ',ptime(ipix)
          write(iu_main_output,'(a)') 	"--------------------------------------------------------------------------"
          do IW=1,RIN%nw
            do ISD=1,RIN%NSD
              NDIM3 = RIN%NDIM%n3(ISD,IDIM1)
              write(iu_main_output,10) &
              'wl=',RIN%wave(IW),'  isd=',ISD, &
              '  sca=',GOUT%aerosol%opt%pixel(ipix)%wl(IW)%ssa(ISD)* &
                       GOUT%aerosol%opt%pixel(ipix)%wl(IW)%ext(ISD), &
              '  ssa=',GOUT%aerosol%opt%pixel(ipix)%wl(IW)%ssa(ISD), &
              '  n=',GOUT%aerosol%rind%pixel(ipix)%wl(IW)%mreal(ISD), &
              '  k=',GOUT%aerosol%rind%pixel(ipix)%wl(IW)%mimag(ISD), &
              '  fsphere=',sph(1:NDIM3,ISD,ipix)
              write(iu_main_output,'(2a)')	"         # 	    angle        ", &
              "ph11/sca"
              do II = 1,GOUT%aerosol%phmx%nangle
                write(iu_main_output,11) II,GOUT%aerosol%phmx%angle(II), &
                                  GOUT%aerosol%phmx%pixel(ipix)%wl(IW)%ph11(II,ISD)
              enddo ! II
            enddo ! ISD
          enddo ! IW
        enddo ! ipix
      case(2)
        do ipix=1,segment_meas%npixels
          write(iu_main_output,'(a)') 	"--------------------------------------------------------------------------"
          write(iu_main_output,'(a,i0,a,a10,a,a8)') 'ipix=',ipix,'  yymmdd = ',pdate(ipix),'  hhmmss = ',ptime(ipix)
          write(iu_main_output,'(a)') 	"--------------------------------------------------------------------------"
          do IW=1,RIN%nw
            do ISD=1,RIN%NSD
              NDIM3 = RIN%NDIM%n3(ISD,IDIM1)
              write(iu_main_output,10) &
              'wl=',RIN%wave(IW),'  isd=',ISD, &
              '  sca=',GOUT%aerosol%opt%pixel(ipix)%wl(IW)%ssa(ISD)* &
                       GOUT%aerosol%opt%pixel(ipix)%wl(IW)%ext(ISD), &
              '  ssa=',GOUT%aerosol%opt%pixel(ipix)%wl(IW)%ssa(ISD), &
              '  n=',GOUT%aerosol%rind%pixel(ipix)%wl(IW)%mreal(ISD), &
              '  k=',GOUT%aerosol%rind%pixel(ipix)%wl(IW)%mimag(ISD), &
              '  fsphere=',sph(1:NDIM3,ISD,ipix)
              write(iu_main_output,'(2a)')	"         # 	    angle        ", &
              "ph11/sca        ph12/sca"
              do II = 1,GOUT%aerosol%phmx%nangle
                write(iu_main_output,11) II,GOUT%aerosol%phmx%angle(II), &
                                  GOUT%aerosol%phmx%pixel(ipix)%wl(IW)%ph11(II,ISD), &
                                  GOUT%aerosol%phmx%pixel(ipix)%wl(IW)%ph12(II,ISD)
              enddo ! II
            enddo ! ISD
          enddo ! IW
        enddo ! ipix
      case(3)
        do ipix=1,segment_meas%npixels
          write(iu_main_output,'(a)') 	"--------------------------------------------------------------------------"
          write(iu_main_output,'(a,i0,a,a10,a,a8)') 'ipix=',ipix,'  yymmdd = ',pdate(ipix),'  hhmmss = ',ptime(ipix)
          write(iu_main_output,'(a)') 	"--------------------------------------------------------------------------"
          do IW=1,RIN%nw
            do ISD=1,RIN%NSD
              NDIM3 = RIN%NDIM%n3(ISD,IDIM1)
              write(iu_main_output,10) &
              'wl=',RIN%wave(IW),'  isd=',ISD, &
              '  sca=',GOUT%aerosol%opt%pixel(ipix)%wl(IW)%ssa(ISD)* &
                       GOUT%aerosol%opt%pixel(ipix)%wl(IW)%ext(ISD), &
              '  ssa=',GOUT%aerosol%opt%pixel(ipix)%wl(IW)%ssa(ISD), &
              '  n=',GOUT%aerosol%rind%pixel(ipix)%wl(IW)%mreal(ISD), &
              '  k=',GOUT%aerosol%rind%pixel(ipix)%wl(IW)%mimag(ISD), &
              '  fsphere=',sph(1:NDIM3,ISD,ipix)
              write(iu_main_output,'(2a)')	"         # 	    angle        ", &
              "ph11/sca        ph12/sca        ph22/sca"
              do II = 1,GOUT%aerosol%phmx%nangle
                write(iu_main_output,11) II,GOUT%aerosol%phmx%angle(II),                    &
                                  GOUT%aerosol%phmx%pixel(ipix)%wl(IW)%ph11(II,ISD), &
                                  GOUT%aerosol%phmx%pixel(ipix)%wl(IW)%ph12(II,ISD), &
                                  GOUT%aerosol%phmx%pixel(ipix)%wl(IW)%ph22(II,ISD)

              enddo ! II
            enddo ! ISD
          enddo ! IW
        enddo ! ipix
      case(4)
        do ipix=1,segment_meas%npixels
          write(iu_main_output,'(a)') 	"--------------------------------------------------------------------------"
          write(iu_main_output,'(a,i0,a,a10,a,a8)') 'ipix=',ipix,'  yymmdd = ',pdate(ipix),'  hhmmss = ',ptime(ipix)
          write(iu_main_output,'(a)') 	"--------------------------------------------------------------------------"
          do IW=1,RIN%nw
            do ISD=1,RIN%NSD
              NDIM3 = RIN%NDIM%n3(ISD,IDIM1)
              write(iu_main_output,10) &
              'wl=',RIN%wave(IW),'  isd=',ISD, &
              '  sca=',GOUT%aerosol%opt%pixel(ipix)%wl(IW)%ssa(ISD)* &
                       GOUT%aerosol%opt%pixel(ipix)%wl(IW)%ext(ISD), &
              '  ssa=',GOUT%aerosol%opt%pixel(ipix)%wl(IW)%ssa(ISD), &
              '  n=',GOUT%aerosol%rind%pixel(ipix)%wl(IW)%mreal(ISD), &
              '  k=',GOUT%aerosol%rind%pixel(ipix)%wl(IW)%mimag(ISD), &
              '  fsphere=',sph(1:NDIM3,ISD,ipix)
              write(iu_main_output,'(2a)')	"         # 	    angle        ", &
              "ph11/sca        ph12/sca        ph22/sca        ph33/sca"
              do II = 1,GOUT%aerosol%phmx%nangle
                write(iu_main_output,11) II,GOUT%aerosol%phmx%angle(II),                    &
                                  GOUT%aerosol%phmx%pixel(ipix)%wl(IW)%ph11(II,ISD), &
                                  GOUT%aerosol%phmx%pixel(ipix)%wl(IW)%ph12(II,ISD), &
                                  GOUT%aerosol%phmx%pixel(ipix)%wl(IW)%ph22(II,ISD), &
                                  GOUT%aerosol%phmx%pixel(ipix)%wl(IW)%ph33(II,ISD)
              enddo ! II
            enddo ! ISD
          enddo ! IW
        enddo ! ipix
      case(5)
        do ipix=1,segment_meas%npixels
          write(iu_main_output,'(a)') 	"--------------------------------------------------------------------------"
          write(iu_main_output,'(a,i0,a,a10,a,a8)') 'ipix=',ipix,'  yymmdd = ',pdate(ipix),'  hhmmss = ',ptime(ipix)
          write(iu_main_output,'(a)') 	"--------------------------------------------------------------------------"
          do IW=1,RIN%nw
            do ISD=1,RIN%NSD
              NDIM3 = RIN%NDIM%n3(ISD,IDIM1)
              write(iu_main_output,10) &
              'wl=',RIN%wave(IW),'  isd=',ISD, &
              '  sca=',GOUT%aerosol%opt%pixel(ipix)%wl(IW)%ssa(ISD)* &
                       GOUT%aerosol%opt%pixel(ipix)%wl(IW)%ext(ISD), &
              '  ssa=',GOUT%aerosol%opt%pixel(ipix)%wl(IW)%ssa(ISD), &
              '  n=',GOUT%aerosol%rind%pixel(ipix)%wl(IW)%mreal(ISD), &
              '  k=',GOUT%aerosol%rind%pixel(ipix)%wl(IW)%mimag(ISD), &
              '  fsphere=',sph(1:NDIM3,ISD,ipix)
              write(iu_main_output,'(2a)')	"         # 	    angle        ", &
              "ph11/sca        ph12/sca        ph22/sca        ph33/sca        ph34/sca"
              do II = 1,GOUT%aerosol%phmx%nangle
                write(iu_main_output,11) II,GOUT%aerosol%phmx%angle(II),                    &
                                  GOUT%aerosol%phmx%pixel(ipix)%wl(IW)%ph11(II,ISD), &
                                  GOUT%aerosol%phmx%pixel(ipix)%wl(IW)%ph12(II,ISD), &
                                  GOUT%aerosol%phmx%pixel(ipix)%wl(IW)%ph22(II,ISD), &
                                  GOUT%aerosol%phmx%pixel(ipix)%wl(IW)%ph33(II,ISD), &
                                  GOUT%aerosol%phmx%pixel(ipix)%wl(IW)%ph34(II,ISD)
              enddo ! II
            enddo ! ISD
          enddo ! IW
        enddo ! ipix
      case(6)
        do ipix=1,segment_meas%npixels
          write(iu_main_output,'(a)') 	"--------------------------------------------------------------------------"
          write(iu_main_output,'(a,i0,a,a10,a,a8)') 'ipix=',ipix,'  yymmdd = ',pdate(ipix),'  hhmmss = ',ptime(ipix)
          write(iu_main_output,'(a)') 	"--------------------------------------------------------------------------"
          do IW=1,RIN%nw
            do ISD=1,RIN%NSD
              NDIM3 = RIN%NDIM%n3(ISD,IDIM1)
              write(iu_main_output,10) &
              'wl=',RIN%wave(IW),'  isd=',ISD, &
              '  sca=',GOUT%aerosol%opt%pixel(ipix)%wl(IW)%ssa(ISD)* &
                       GOUT%aerosol%opt%pixel(ipix)%wl(IW)%ext(ISD), &
              '  ssa=',GOUT%aerosol%opt%pixel(ipix)%wl(IW)%ssa(ISD), &
              '  n=',GOUT%aerosol%rind%pixel(ipix)%wl(IW)%mreal(ISD), &
              '  k=',GOUT%aerosol%rind%pixel(ipix)%wl(IW)%mimag(ISD), &
              '  fsphere=',sph(1:NDIM3,ISD,ipix)
              write(iu_main_output,'(2a)')	"         # 	   angle        ", &
              "ph11/sca        ph12/sca        ph22/sca        ph33/sca        ph34/sca        ph44/sca"
              do II = 1,GOUT%aerosol%phmx%nangle
                write(iu_main_output,11) II,GOUT%aerosol%phmx%angle(II),                    &
                                  GOUT%aerosol%phmx%pixel(ipix)%wl(IW)%ph11(II,ISD), &
                                  GOUT%aerosol%phmx%pixel(ipix)%wl(IW)%ph12(II,ISD), &
                                  GOUT%aerosol%phmx%pixel(ipix)%wl(IW)%ph22(II,ISD), &
                                  GOUT%aerosol%phmx%pixel(ipix)%wl(IW)%ph33(II,ISD), &
                                  GOUT%aerosol%phmx%pixel(ipix)%wl(IW)%ph34(II,ISD), &
                                  GOUT%aerosol%phmx%pixel(ipix)%wl(IW)%ph44(II,ISD)
              enddo ! II
            enddo ! ISD
          enddo ! IW
        enddo ! ipix
      end select

      write(iu_main_output,'(a)') 	"--------------------------------------------------------------------------"
10    format(a,f9.6,a,i0,4(a,e12.5),a,25e12.5)
11    format(6x,I4,3x,f12.3,2x,6(e14.5,2x))

      return
      end subroutine print_phmx

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! Print broadband flux

      subroutine print_forcing_bbflux (iu,RIN,segment_meas,GOUT)

      use mod_retr_general_output_derived_type
      use mod_print_array
      use mod_sdata_derived_type
            
      implicit none  
!----------------------------------------------------------------------------------------
! IN :
      integer,                     intent(in)  ::  iu
      type(retr_input_settings),   intent(in)  ::  RIN
      type(segment_data),          intent(in)  ::  segment_meas
      type(output_segment_general),intent(in)  ::  GOUT
!----------------------------------------------------------------------------------------
! LOCAL :
      integer                              :: npixels, IHLV, NHLV
!----------------------------------------------------------------------------------------
      npixels  = segment_meas%npixels
      NHLV     = GOUT%forcing%bbflux%pixel(1)%NHLV

      write(iu,*) ''
      write(iu,'(a)') 'Upward Broadband Flux without Aerosol (w/m^2):'
      do IHLV=1,NHLV
        write(iu,'(f9.3,a5)', advance='no') GOUT%forcing%bbflux%pixel(1)%HLV(IHLV), ' km  '
        call rprint_array_real (iu,GOUT%forcing%bbflux%pixel(:)%BBUFX0(IHLV), '(e14.5)', 1, npixels)
      enddo ! IHLV

      write(iu,*) ''
      write(iu,'(a)') 'Downward Broadband Flux without Aerosol (w/m^2):'
      do IHLV=1,NHLV
        write(iu,'(f9.3,a5)', advance='no') GOUT%forcing%bbflux%pixel(1)%HLV(IHLV), ' km  '
        call rprint_array_real (iu,GOUT%forcing%bbflux%pixel(:)%BBDFX0(IHLV), '(e14.5)', 1, npixels)
      enddo ! IHLV

      write(iu,*) ''
      write(iu,'(a)') 'Upward Broadband Flux with Aerosol (w/m^2):'
      do IHLV=1,NHLV
        write(iu,'(f9.3,a5)', advance='no') GOUT%forcing%bbflux%pixel(1)%HLV(IHLV), ' km  '
        call rprint_array_real (iu,GOUT%forcing%bbflux%pixel(:)%BBUFXA(IHLV), '(e14.5)', 1, npixels)
      enddo ! IHLV

      write(iu,*) ''
      write(iu,'(a)') 'Downward Broadband Flux with Aerosol (w/m^2):'
      do IHLV=1,NHLV
        write(iu,'(f9.3,a5)', advance='no') GOUT%forcing%bbflux%pixel(1)%HLV(IHLV), ' km  '
        call rprint_array_real (iu,GOUT%forcing%bbflux%pixel(:)%BBDFXA(IHLV), '(e14.5)', 1, npixels)
      enddo ! IHLV

      return
      end subroutine print_forcing_bbflux

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! Print net aerosol radiative forcing

      subroutine print_forcing_forcing (iu,RIN,segment_meas,GOUT)

      use mod_retr_general_output_derived_type
      use mod_print_array
      use mod_sdata_derived_type
            
      implicit none  
!----------------------------------------------------------------------------------------
! IN :
      integer,                     intent(in)  ::  iu
      type(retr_input_settings),   intent(in)  ::  RIN
      type(segment_data),          intent(in)  ::  segment_meas
      type(output_segment_general),intent(in)  ::  GOUT
!----------------------------------------------------------------------------------------
! LOCAL :
      integer                              :: npixels, IHLV, NHLV
!----------------------------------------------------------------------------------------
      npixels  = segment_meas%npixels
      NHLV     = GOUT%forcing%bbflux%pixel(1)%NHLV

      write(iu,*) ''
      write(iu,'(a)') 'Net Aerosol Forcing (w/m^2):'
      do IHLV=1,NHLV
        write(iu,'(f9.3,a5)', advance='no') GOUT%forcing%forcing%pixel(1)%HLV(IHLV), ' km  '
        call rprint_array_real (iu,GOUT%forcing%forcing%pixel(:)%NETFORC(IHLV), '(e14.5)', 1, npixels)
      enddo ! IHLV

      write(iu,*) ''
      write(iu,'(a)') 'Aerosol Forcing Efficiency (w/m^2):'
      do IHLV=1,NHLV
        write(iu,'(f9.3,a5)', advance='no') GOUT%forcing%forcing%pixel(1)%HLV(IHLV), ' km  '
        call rprint_array_real (iu,GOUT%forcing%forcing%pixel(:)%FORCEFF(IHLV), '(e14.5)', 1, npixels)
      enddo ! IHLV

      return
      end subroutine print_forcing_forcing

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

