! **
! **  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.  
! **  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
! **

#include "../constants_set/mod_globals.inc"
module mod_sdata

      use mod_time_utils
      use mod_index_cloud
      use mod_sdata_meas_type
      use mod_sdata_derived_type
      use mod_stop_report
      !use mod_par_OS,  only : NBVM
      !use mod_par_inv, only : KWM,KW,KVERTM,KNBVM,KIP,KSURF, &
      !                        KITIME,KIX,KIY,KMESS,KKNOISE 

      implicit none

contains

!	sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

   subroutine print_array_int(array, element_format, from_element, to_element)
     integer, intent(in) :: array(:)
     character(*), intent(in) :: element_format
     integer, intent(in), optional :: from_element
     integer, intent(in), optional :: to_element
     integer :: i
     integer :: from_element_
     integer :: to_element_
     integer :: nelements

     if (present(from_element)) then
        from_element_ = from_element
     else
        from_element_ = 1
     end if

     if (present(to_element)) then
        to_element_ = to_element
     else
        to_element_ = size(array)
     end if

     nelements = to_element_ - from_element_ + 1
     write(*, '("[",I0,"/",I0,"]")', advance='no') nelements, size(array)

     
     do i = from_element_, to_element_
        write(*,'(" ")', advance='no')
        write(*,fmt=element_format, advance='no') array(i)
        write(*,'(" ")', advance='no')
     end do
     write(*,*)
     
   end subroutine print_array_int

   subroutine print_array_real(array, element_format, from_element, to_element)
     real(kind=C_FLOAT), intent(in) :: array(:)
     character(*), intent(in) :: element_format
     integer, intent(in), optional :: from_element
     integer, intent(in), optional :: to_element
     integer :: i
     integer :: from_element_
     integer :: to_element_
     integer :: nelements

     if (present(from_element)) then
        from_element_ = from_element
     else
        from_element_ = 1
     end if

     if (present(to_element)) then
        to_element_ = to_element
     else
        to_element_ = size(array)
     end if

     nelements = to_element_ - from_element_ + 1
     write(*, '("[",I0,"/",I0,"]")', advance='no') nelements, size(array)
     do i = from_element_, to_element_
        write(*,'(" ")', advance='no')
        write(*,fmt=element_format, advance='no') array(i)
        write(*,'(" ")', advance='no')
     end do
     write(*,*)
     
   end subroutine print_array_real

   subroutine print_data_wl(label, dwl)
     character(*), intent(in) :: label
     type(data_wl), intent(in) :: dwl
     integer :: ip
     integer :: ivm

     write(*,'(A,"%NIP: ",I0)') trim(label), dwl%NIP

     write(*,'(A,"%NBVM(ip=1:",I0,"): ")', advance='no') trim(label), dwl%NIP
     call print_array_int(dwl%NBVM, "(I0)", 1, dwl%NIP)

     
     write(*,'(A,"%meas_type(ip=1:",I0,"): ")', advance='no') trim(label), dwl%NIP
     call print_array_int(dwl%meas_type,"(I0)", 1, dwl%NIP)

     write(*,'(A,"%wl: ",F9.3)') trim(label), dwl%wl
     write(*,'(A,"%ind_wl: ",I0)') trim(label), dwl%ind_wl
     write(*,'(A,"%sza: ",F9.3)') trim(label), dwl%sza

     do ip = 1, dwl%NIP
        if(dwl%meas_type(ip) .lt. meas_type_lid_beg .or. dwl%meas_type(ip) .gt. meas_type_lid_end ) then
          write(*,'(A,"%thetav(ivm=1:",I0,",ip=",I0,"): ")', advance='no') trim(label), dwl%NBVM(ip), ip
          call print_array_real(dwl%thetav(:,ip),"(E12.3)", 1, dwl%NBVM(ip))
        !else
          !write(*,'(A,"%HVP(ivm=1:",I0,",ip=",I0,"): ")', advance='no') trim(label), dwl%NBVM(ip), ip
          !call print_array_real(dwl%thetav(:,ip),"(E12.3)", 1, dwl%NBVM(ip))
        end if        
     end do     

     do ip = 1, dwl%NIP
        if(dwl%meas_type(ip) .lt. meas_type_lid_beg .or. dwl%meas_type(ip) .gt. meas_type_lid_end ) then
          write(*,'(A,"%phi(ivm=1:",I0,",ip=",I0,"): ")', advance='no') trim(label), dwl%NBVM(ip), ip
          call print_array_real(dwl%phi(:,ip),"(E12.3)", 1, dwl%NBVM(ip))
        endif
     end do

     write(*,'(A,"%Nsurf: ",I0)') trim(label), dwl%Nsurf
     write(*,'(A,"%groundpar(isurf=1:",I0,"): ")', advance='no') trim(label), dwl%Nsurf
     call print_array_real(dwl%groundpar,"(E12)", 1, dwl%Nsurf)
     
     write(*,'(A,"%gaspar: ",F9.3)') trim(label), dwl%gaspar

     do ip = 1, dwl%NIP
        select case (dwl%meas_type(ip))
        case(MEAS_TYPE_TOD)
           write(*,'(A,"%tod: ")', advance='no') trim(label)
           call print_array_real(dwl%tau, "(E12.3)", 1, dwl%NBVM(ip))
        case(MEAS_TYPE_AOD)
           write(*,'(A,"%aod: ")', advance='no') trim(label)
           call print_array_real(dwl%tau, "(E12.3)", 1, dwl%NBVM(ip))
        case(MEAS_TYPE_P11)
           write(*,'(A,"%p11: ")', advance='no') trim(label)
           call print_array_real(dwl%p11, "(E12.3)", 1, dwl%NBVM(ip))
        case(MEAS_TYPE_P12)
           write(*,'(A,"%p12: ")', advance='no') trim(label)
           call print_array_real(dwl%p12, "(E12.3)", 1, dwl%NBVM(ip))
        case(MEAS_TYPE_P22)
           write(*,'(A,"%p22: ")', advance='no') trim(label)
           call print_array_real(dwl%p22, "(E12.3)", 1, dwl%NBVM(ip))
        case(MEAS_TYPE_P33)
           write(*,'(A,"%p33: ")', advance='no') trim(label)
           call print_array_real(dwl%p33, "(E12.3)", 1, dwl%NBVM(ip))
        case(MEAS_TYPE_P34)
           write(*,'(A,"%p34: ")', advance='no') trim(label)
           call print_array_real(dwl%p34, "(E12.3)", 1, dwl%NBVM(ip))
        case(MEAS_TYPE_P44)
           write(*,'(A,"%p44: ")', advance='no') trim(label)
           call print_array_real(dwl%p44, "(E12.3)", 1, dwl%NBVM(ip))
        case(MEAS_TYPE_LS)
           write(*,'(A,"%LS: ")', advance='no') trim(label)
           call print_array_real(dwl%LS, "(E12.3)", 1, dwl%NBVM(ip))
        case(MEAS_TYPE_DP)
           write(*,'(A,"%DP: ")', advance='no') trim(label)
           call print_array_real(dwl%DP, "(E12.3)", 1, dwl%NBVM(ip))
        case(MEAS_TYPE_RL)
           write(*,'(A,"%RL: ")', advance='no') trim(label)
           call print_array_real(dwl%RL, "(E12.3)", 1, dwl%NBVM(ip))
        case(MEAS_TYPE_I)
           write(*,'(A,"%I: ")', advance='no') trim(label)
           call print_array_real(dwl%I, "(E12.3)", 1, dwl%NBVM(ip))
        case(MEAS_TYPE_Q)
           write(*,'(A,"%Q: ")', advance='no') trim(label)
           call print_array_real(dwl%Q, "(E12.3)", 1, dwl%NBVM(ip))
        case(MEAS_TYPE_U)
           write(*,'(A,"%U: ")', advance='no') trim(label)
           call print_array_real(dwl%U, "(E12.3)", 1, dwl%NBVM(ip))
        case(MEAS_TYPE_P)
           write(*,'(A,"%P: ")', advance='no') trim(label)
           call print_array_real(dwl%P, "(E12.3)", 1, dwl%NBVM(ip))   
        case default
           write(tmp_message,'(a,i4)') &
           "Fatal error. Unexpected value for meas_type (I0). Please check the code.", &
           dwl%meas_type(ip)
           G_ERROR(trim(tmp_message))
        end select
     end do ! ip
     
     write(*,'(A,"%IFCOV(ip=1:",I0,"): ")', advance='no') trim(label), dwl%NIP
     call print_array_int(dwl%IFCOV, "(I0)", 1, dwl%NIP)

     write(*,'(A,"%IFMP(ip=1:",I0,"): ")', advance='no') trim(label), dwl%NIP
     call print_array_int(dwl%IFMP, "(I0)", 1, dwl%NIP)

     do ip = 1, dwl%NIP
        write(*,'(A,"%CMTRX(ivm=1:",I0,",ip=",I0,"): ")', advance='no') trim(label), &
             dwl%NBVM(ip)*dwl%IFCOV(ip), ip
        call print_array_real(dwl%CMTRX(:,ip),"(E12.3)", 1, dwl%NBVM(ip)*dwl%IFCOV(ip))
     end do

     do ip = 1, dwl%NIP
        write(*,'(A,"%MPROF(ivm=1:",I0,",ip=",I0,"): ")', advance='no') trim(label), &
             dwl%NBVM(ip)*dwl%IFMP(ip), ip
        call print_array_real(dwl%MPROF(:,ip),"(E12.3)", 1, dwl%NBVM(ip)*dwl%IFMP(ip))
     end do

   end subroutine print_data_wl

   subroutine print_pixel(label, one_pixel)
     character(*), intent(in) :: label
     type(pixel), intent(in) :: one_pixel
     integer :: iwl
     character(256) :: dwl_label
     character(20) :: readable_time
     
     call convert_time_to_string(one_pixel%t, '%FT%H:%M:%SZ', readable_time)

     write(*,'(A,"%HOBS: ",F11.3)') trim(label), one_pixel%HOBS
     write(*,'(A,"%nwl: ",I0)') trim(label), one_pixel%nwl
     write(*,'(A,"%cloudy: ",I0)') trim(label), one_pixel%cloudy
     write(*,'(A,"%x: ",F9.3)') trim(label), one_pixel%x
     write(*,'(A,"%y: ",F9.3)') trim(label), one_pixel%y
     write(*,'(A,"%t: ",I0," (",A,")")') trim(label), one_pixel%t, trim(readable_time)
     write(*,'(A,"%ix: ",I0)') trim(label), one_pixel%ix
     write(*,'(A,"%iy: ",I0)') trim(label), one_pixel%iy
     write(*,'(A,"%it: ",I0)') trim(label), one_pixel%it
     write(*,'(A,"%MASL: ",F9.3)') trim(label), one_pixel%MASL
     write(*,'(A,"%land_percent: ",F9.3)') trim(label), one_pixel%land_percent
     write(*,'(A,"%irow: ",I0)') trim(label), one_pixel%irow
     write(*,'(A,"%icol: ",I0)') trim(label), one_pixel%icol
     write(*,'(A,"%IFGAS: ",I0)') trim(label), one_pixel%IFGAS
     write(*,'(A,"%meas(iwl=1:",I0,"): ")', advance='no') trim(label), one_pixel%nwl
     write(*,*)
     do iwl = 1, one_pixel%nwl
        write(dwl_label, '(A,"%meas(",I0,")")') trim(label), iwl
        call print_data_wl(dwl_label, one_pixel%meas(iwl))
     end do
     write(*,'(A,"%HVP(KVERTM=1:",I0,"): ")', advance='no') trim(label), KVERTM
     call print_array_real(one_pixel%HVP, "(F0.1)")
     
   end subroutine print_pixel

   subroutine print_segment(label, segment)
     character(*), intent(in) :: label
     type(segment_data), intent(in) :: segment
     integer :: ipixel
     character(256) :: pixel_label

     write(*,'(A," npixels: ",I0," NX: ",I0," NY: ",I0," NT: ",I0," [KIMAGE: ",I0,"]")') &
          label, segment%npixels, segment%NX, segment%NY, segment%NT, KIMAGE
     do ipixel = 1, segment%npixels
        write(pixel_label, '(A,"%pixels(",I0,")")') trim(label), ipixel
        call print_pixel(pixel_label, segment%pixels(ipixel))
     end do
     
   end subroutine print_segment

!	sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

    subroutine read_input_sdata ( iu_main_output, &
                                  sdata_file,     &                              
                                  segment		      &
                                )
		implicit none
!  ------------------------------------------------------------------------------------------------------
		integer,intent(in)		                 ::	iu_main_output
    character(*),intent(in)                :: sdata_file
		type(segment_data),intent(out)         ::	segment
!	------------------------------------------------------------------------------------------------------
!	------------------------------------------------------------------------------------------------------
!        character(LEN=100)      ::  msg
		integer		      ::	npixels
		integer		      ::	iu_file
		integer					::	ii
		integer					::	IT,IPIX
		integer					::	IW,IV,IP,IBRF,k
		integer					::	io_status
		integer					::	max_npixels 
!	------------------------------------------------------------------------------------------------------
		integer					     ::	NPIXELT
		character(LEN=20)		 ::	iso8601_string      
		character(LEN=16)		 ::	data_filename
		integer					     ::	Nsurf,IFGAS
 		integer					     :: NX,NY,NT  !,NMT,IMT(14)
		real                 ::	HOBS   ! observation Height
		integer(KIND_TIME)   ::	TIME_sec

!	------------------------------------------------------------------------------------------------------
		integer              ::	IX
		integer              ::	IY
		integer              ::	cloudy
		integer              ::	irow
		integer              ::	icol
		integer              ::	nwave_lengths
		real                 ::	y_coord_latitude
		real                 ::	x_coord_longitude
		real                 ::	MASL     ! Metres Above Sea Level 
		real                 ::	land_percent

    integer,dimension(KWM)    ::  num_meas_type
		real,   dimension(KWM)    ::  wave_length
		real,   dimension(KWM)    ::  sza ! solar zenith angle

		real,dimension(KNBVM,KIP,KWM)  ::  thetav
		real,dimension(KNBVM,KIP,KWM)  ::  phi
		real,dimension(KSURF,KWM)      ::  groundpar !??? ask
		real,dimension(KWM)            ::  gaspar    !??? ask
!	------------------------------------------------------------------------------------------------------
    integer                            ::  meas
    integer,dimension(KIP,KWM)         ::  IFCOV,IFMP      
		real,dimension(KNBVM,KIP,KWM)      ::  FPSW,CMTRX
		real,dimension(KVERTM,KIP,KWM)     ::  MPROF
    integer,dimension(KIP,KWM)         ::  meas_type
		integer,dimension(KIP,KWM)         ::  nvalid_meas

    real                               ::  tau,I,Q,U,Q2,U2,P,LS,DP,RL,  &
                                           p11, p12, p22, p33, p34, p44
    logical		                         ::	 status_funct
!	------------------------------------------------------------------------------------------------------
      integer    ::  ind_cloudy
		ind_cloudy = 0
!	------------------------------------------------------------------------------------------------------
      max_npixels = size(segment%pixels(:))
      npixels = 0
      ii = 0
!	------------------------------------------------------------------------------------------------------
      open (newunit=iu_file, FILE=TRIM(sdata_file), STATUS='old')

      READ (iu_file,*,IOSTAT= io_status) 
      READ (iu_file,*,IOSTAT= io_status) NX,NY,NT

      segment%NT = NT
      segment%NX = NX
      segment%NY = NY
      
      !if(NT .gt. KITIME .or. NX .gt. KIX .or. NY .gt. KIY .or. NMT .gt. KIP) then
         !write(*,*) 'KITIME,KIX,KIY,KIP : ',KITIME,KIX,KIY,KIP
         !write(*,*) 'NT,    NX, NY, NMT : ',NT,NX,NY,NMT
         !write(*,*) 'Check parameters in mod_par_inv'
         !stop 'stop in read_input_sdata'
      !endif
      if(NT .gt. KITIME .or. NX .gt. KIX .or. NY .gt. KIY) then
         write(*,*) 'KITIME,KIX,KIY : ',KITIME,KIX,KIY
         write(*,*) 'NT,    NX, NY  : ',NT,NX,NY
         write(*,*) 'Check parameters in mod_par_inv'
         stop 'stop in read_input_sdata'
      endif

      DO IT=1,NT 
      READ(iu_file,*,IOSTAT= io_status) NPIXELT,iso8601_string,  &
                                        HOBS,Nsurf,IFGAS
                                        

      TIME_sec = convert_string_to_time(iso8601_string, status_funct)

      if (.not. status_funct) then
        write(0,*) 'iso8601_string=',iso8601_string
        stop 'in read_input_sdata: error on convert_string_to_time'
      end if

      if(IFGAS .gt. 1) then
         write(*,*) 'IFGAS=',IFGAS,' > 1 value is not valid '
         stop 'stop in read_input_sdata'
      endif ! IFGAS .gt. 1
      
!write(*,*) 'read pixel before loop IPIX, IT=',IT
      DO IPIX=1, NPIXELT
      ii = ii + 1 
      if(ii .GT. max_npixels) then
         write(*,*) 'Too many pixels to read'
         stop 'stop in read_input_sdata'
      endif
      READ(iu_file,*,IOSTAT= io_status)   &
            IX,         &
            IY,			    &
            cloudy,	  	&
            irow,       &
            icol,       &
            x_coord_longitude,   &
            y_coord_latitude,		 &
            MASL,                &
            land_percent,        &
            
            nwave_lengths,       &
            (wave_length(IW),          &               
               IW=1,nwave_lengths),    &

            (num_meas_type(IW),         &
               IW=1,nwave_lengths),     & 
                          
            ((meas_type(IP,IW),         &
               IP=1,num_meas_type(IW)), &                     
               IW=1,nwave_lengths),     & 

            ((nvalid_meas(IP,IW),       &
               IP=1,num_meas_type(IW)), &                     
               IW=1,nwave_lengths),     &
               
            (sza(IW),                      &
               IW=1,nwave_lengths),	       &

            (((thetav(IV,IP,IW),           &
               IV=1,nvalid_meas(IP,IW)), 	 &
               IP=1,num_meas_type(IW)),	   &                                       
               IW=1,nwave_lengths),	       &

            (((phi(IV,IP,IW),              &
               IV=1,nvalid_meas(IP,IW)), 	 &
               IP=1,num_meas_type(IW)),	   &                                       
               IW=1,nwave_lengths),	       &
                  
            (((FPSW(IV,IP,IW),		         &
               IV=1,nvalid_meas(IP,IW)),   &
               IP=1,num_meas_type(IW)),	   &                                       
               IW=1,nwave_lengths),        &

            ((groundpar(IBRF,IW),       &
               IBRF=1,Nsurf),           &
               IW=1,nwave_lengths),		  &
               
            ((gaspar(IW),               &
               k=1,IFGAS),	            &
               IW=1,nwave_lengths),     &
               
            ((IFCOV(IP,IW),             &
               IP=1,num_meas_type(IW)), &                     
               IW=1,nwave_lengths),     & 

            (((CMTRX(IV,IP,IW),                       & 
               IV=1,IFCOV(IP,IW)*nvalid_meas(IP,IW)), &
               IP=1,num_meas_type(IW)),  &
               IW=1,nwave_lengths),		   &
                  

            ((IFMP(IP,IW),              &
               IP=1,num_meas_type(IW)), &                     
               IW=1,nwave_lengths),     & 
               
            (((MPROF(IV,IP,IW),                & 
               IV=1,IFMP(IP,IW)*nvalid_meas(IP,IW)),  &
               IP=1,num_meas_type(IW)),        &
               IW=1,nwave_lengths)            
                                                           
      if( io_status .NE. 0 ) then
         write(*,*) 'Can not read data file'
         stop 'stop in read_input_sdata'
      endif

!write(*,*) "IFCOV:"
!do IW=1,nwave_lengths
!write(*,'(10i5)')  (IFCOV(IP,IW),IP=1,num_meas_type(IW))
!enddo

!write(*,*) "meas:"
!do IW=1,nwave_lengths
!do IP=1,num_meas_type(IW)
!write(*,'(10e16.5)')  (FPSW(IV,IP,IW),IV=1,nvalid_meas(IP,IW))
!enddo ! IP
!enddo ! IW 
      
!if(ii .eq. 2) then
!do IW=1,nwave_lengths
!do IP=1,num_meas_type(IW)
   !write(*,'(/,2(a,i4))') 'CMTRX: iw=',iw,'  ip=',ip
   !write(*,'(10e16.5)') (CMTRX(IV,IP,IW),IV=1,IFCOV(IP,IW)*nvalid_meas(IP,IW))
!enddo 
!enddo 
!stop
!endif

! CHECK MEAS DATA

      !write(*,'(a,i4)') 'pixel# =',ii !cloudy,'  - 1: cloudy'
      do iw=1,nwave_lengths
      do ip=1,num_meas_type(iw)
      meas = meas_type(ip,iw)
      !write(*,*) iw,ip,meas,'  - iw,ip,meas'
      !write(*,*) 'FPSW(1:NBVM)  ',(FPSW(IV,IP,IW),IV=1,nvalid_meas(IP,IW)) 	 
      !write(*,*) 'MNOISEI(IP,IW) =',MNOISEI(IP,IW) 	 

      do iv=1,nvalid_meas(ip,iw)
      select case(meas)
      case(meas_type_tod)
! tod
         tau = FPSW(IV,IP,IW)
         IF(tau .LE. 0.0) THEN
         WRITE(iu_main_output,*) 'PROBLEM WITH DATA !!!'
         WRITE(iu_main_output,*) IT,IX,IY,IV,IW,' IT,IX,IY,IV,IW'
         WRITE(iu_main_output,*) tau,' - TOD'
         cloudy=0 
		   ENDIF ! tau .LE. 0.0        
      case(meas_type_aod)
! aod
         tau = FPSW(IV,IP,IW)
         IF(tau .LE. 0.0) THEN
         WRITE(iu_main_output,*) 'PROBLEM WITH DATA !!!'
         WRITE(iu_main_output,*) IT,IX,IY,IV,IW,' IT,IX,IY,IV,IW'
         WRITE(iu_main_output,*) tau,' - AOD'
         cloudy=0 
		   ENDIF ! tau .LE. 0.0
      case(meas_type_p11)
! p11         
         p11 = FPSW(IV,IP,IW)
         IF(p11 .LE. 1e-5) THEN
         WRITE(iu_main_output,*) 'PROBLEM WITH DATA !!!'
         WRITE(iu_main_output,*) IT,IX,IY,IV,IW,' IT,IX,IY,IV,IW'
         WRITE(iu_main_output,*) p11,' - p11'
         cloudy=0 
		   ENDIF ! tau .LE. 0.0        
      case(meas_type_p12)
! p12         
         p12 = FPSW(IV,IP,IW)
         IF(p12/p11 .LE. -1.0) THEN
         WRITE(iu_main_output,*) 'PROBLEM WITH DATA !!!'
         WRITE(iu_main_output,*) IT,IX,IY,IV,IW,' IT,IX,IY,IV,IW'
         WRITE(iu_main_output,*) p11,p12,p12/p11,' - p11,p12,p12/p11'
         cloudy=0 
		   ENDIF ! tau .LE. 0.0    
      case(meas_type_p22)
! p22         
         p22 = FPSW(IV,IP,IW)
         IF(p22 .LE. 1e-5) THEN
         WRITE(iu_main_output,*) 'PROBLEM WITH DATA !!!'
         WRITE(iu_main_output,*) IT,IX,IY,IV,IW,' IT,IX,IY,IV,IW'
         WRITE(iu_main_output,*) p22,' - p22'
         cloudy=0 
		   ENDIF ! tau .LE. 0.0     
      case(meas_type_LS)
! LS
         LS = FPSW(IV,IP,IW)
         IF(LS .LE. 0.0) THEN
         WRITE(iu_main_output,*) 'PROBLEM WITH DATA !!!'
         WRITE(iu_main_output,*) IT,IX,IY,IV,IW,' IT,IX,IY,IV,IW'
         WRITE(iu_main_output,*) LS,' - Lidar Signal'
         cloudy=0 
		   ENDIF ! LR .LE. 0.0        

      case(meas_type_DP)
! DP
         DP = FPSW(IV,IP,IW)
         IF(DP .LE. 0.0) THEN
         WRITE(iu_main_output,*) 'PROBLEM WITH DATA !!!'
         WRITE(iu_main_output,*) IT,IX,IY,IV,IW,' IT,IX,IY,IV,IW'
         WRITE(iu_main_output,*) LS,' - Depolarization'
         cloudy=0 
		   ENDIF ! DP .LE. 0.0        

      case(meas_type_RL)
! RL
         RL = FPSW(IV,IP,IW)
         IF(RL .LE. 0.0) THEN
         WRITE(iu_main_output,*) 'PROBLEM WITH DATA !!!'
         WRITE(iu_main_output,*) IT,IX,IY,IV,IW,' IT,IX,IY,IV,IW'
         WRITE(iu_main_output,*) RL,' - Roman Lidar signal'
         cloudy=0 
		   ENDIF ! DP .LE. 0.0        

      case(meas_type_I)
! I
         I = FPSW(IV,IP,IW)
         IF(I .LE. 1e-5) THEN
         WRITE(iu_main_output,*) 'PROBLEM WITH DATA !!!'
         WRITE(iu_main_output,*) IT,IX,IY,IV,IW,' IT,IX,IY,IV,IW'
         WRITE(iu_main_output,*) I,' - RADIANCE'
         cloudy=0 
		   ENDIF ! I .LE. 1e-5        

      case(meas_type_Q)
! Q
         I = FPSW(IV,IP-1,IW)
         Q = FPSW(IV,IP,  IW)
         U = FPSW(IV,IP+1,IW)         
         IF(abs(Q/I) .GE. 1.0) THEN
         WRITE(iu_main_output,*) 'PROBLEM WITH DATA !!!'
         WRITE(iu_main_output,*) IT,IX,IY,IV,IW,' IT,IX,IY,IV,IW'
         WRITE(iu_main_output,*) I,Q,U,'  - I,Q,U (Q)'
         cloudy=0 
		   ENDIF ! Q .eq. 0.0 .or.         

      case(meas_type_U)
! U
         I = FPSW(IV,IP-2,IW)
         Q = FPSW(IV,IP-1,IW)
         U = FPSW(IV,IP  ,IW)         
         IF(abs(U/I) .GE. 1.0) THEN
         WRITE(iu_main_output,*) 'PROBLEM WITH DATA !!!'
         WRITE(iu_main_output,*) IT,IX,IY,IV,IW,' IT,IX,IY,IV,IW'
         WRITE(iu_main_output,*) I,Q,U,'  - I,Q,U (U)'
         cloudy=0 
		   ENDIF ! U .eq. 0.0 .or.         
         Q2 = Q*Q
         U2 = U*U
         P  = SQRT(Q2+U2)
         IF(abs(P/I) .GE. 1.0) THEN
         WRITE(iu_main_output,*) 'PROBLEM WITH DATA !!!'
         WRITE(iu_main_output,*) IT,IX,IY,IV,IW,' IT,IX,IY,IV,IW'
         WRITE(iu_main_output,*) P/I,'  - Degree of Linear Polarization'
         WRITE(iu_main_output,*) I,Q,U,'  - I,Q,U'
         cloudy=0 
		   ENDIF ! P/I .LE. -1.0        

      case(meas_type_P)
! P
         I = FPSW(IV,IP-1,IW)
         P = FPSW(IV,IP,  IW)
         IF(abs(P/I) .GE. 1.0) THEN
         WRITE(iu_main_output,*) 'PROBLEM WITH DATA !!!'
         WRITE(iu_main_output,*) IT,IX,IY,IV,IW,' IT,IX,IY,IV,IW'
         WRITE(iu_main_output,*) P/I,'  - Degree of Linear Polarization'
         WRITE(iu_main_output,*) P,I,'  - P,I'
         write(*,*) 'DATA STILL BE USED!!!!'
         !cloudy=0 
		   ENDIF ! P/I .LE. -1.0        
         
      case default
         write(*,*) 'iw=',iw,'  ip=',ip,' meas_type=',meas_type,' - Unknown value of meas_type'
         stop 'stop in read_input_sdata'         
      end select

      enddo ! iv      
      enddo ! ip
      enddo ! iw
      
      !write(*,*) cloudy,'  - 2: cloudy'      
      if(cloudy.eq.1) then
		  npixels = npixels+1
      if(size(segment%pixels(:)) .lt. npixels) then
         write(*,*) 'size(segment)=',size(segment%pixels(:)),'  npixels=',npixels
         write(*,*) 'Too many pixels to be set'
         stop 'stop in read_input_sdata'
      endif
!write(*,*) 'npixels=',npixels,'  before set_pixel' 
!write(*,*) 'num_meas_type:',num_meas_type(1:nwave_lengths)
!write(*,*) 'nvalid_meas 1:',nvalid_meas(1,1:nwave_lengths) 
!write(*,*) 'nvalid_meas 2:',nvalid_meas(2,1:nwave_lengths) 

      call set_pixel (                          &
                        HOBS,                   &
                        IX,							        &
                        IY,							        &
                        IT,							        &
                        cloudy,						      &
                        irow,						        &
                        icol,						        &
                        x_coord_longitude,		  &
                        y_coord_latitude,			  &
                        TIME_sec,					      &
                        MASL,						        &
                        land_percent,				    &
                        nwave_lengths,				  &
                        nvalid_meas, 		        &
                        wave_length, 			      &

                        sza,                 	  &
                        thetav, 					      &
                        phi, 						        &
                        Nsurf, 						      &
                        groundpar, 					    &
                        IFGAS, 						      &
                        gaspar, 					      &
                        FPSW,                   &
                        num_meas_type,          &
                        meas_type,IFCOV,IFMP,   &
                        CMTRX,MPROF,			      &
                        segment%pixels(npixels)	&					
                     )

                ind_cloudy = ind_cloudy+cloudy

				
      endif ! cloudy

      ENDDO  ! IPIX

		ENDDO  ! IT 
      
      segment%npixels = npixels
      
      CLOSE(iu_file)
!stop
				
!	------------------------------------------------------------------------------------------------------

      if(ind_cloudy.eq.0) then
         write(*,*) 'All pixels are cloudy.'
         stop 'stop in read_input_sdata'
      endif
	
!	------------------------------------------------------------------------------------------------------
      end subroutine read_input_sdata

!	sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine init_pixel ( pixel_cont )
	
!  ------------------------------------------------------------------------------------------------------
      implicit none
!  ------------------------------------------------------------------------------------------------------
      type(pixel),intent(inout)  ::  pixel_cont
!  ------------------------------------------------------------------------------------------------------
      integer                    ::  iw
!  ------------------------------------------------------------------------------------------------------

      do iw=1,KWM
        pixel_cont%meas(iw)%meas_type(:) = 0
        pixel_cont%meas(iw)%wl	         = 0.0
        pixel_cont%meas(iw)%sza          = 0.0

        pixel_cont%meas(iw)%thetav(:,:)  = 0.0
        pixel_cont%meas(iw)%phi(:,:)	   = 0.0
        pixel_cont%meas(iw)%Nsurf	       = 0
        pixel_cont%meas(iw)%groundpar(:) = 0.0
        pixel_cont%meas(iw)%gaspar       = 0.0
        pixel_cont%meas(iw)%NBVM(:)      = 0
        pixel_cont%meas(iw)%NIP          = 0
        pixel_cont%meas(iw)%tau(:)	     = 0.0
        pixel_cont%meas(iw)%p11(:)	     = 0.0
        pixel_cont%meas(iw)%p12(:)	     = 0.0
        pixel_cont%meas(iw)%p22(:)	     = 0.0
        pixel_cont%meas(iw)%p33(:)	     = 0.0
        pixel_cont%meas(iw)%p34(:)	     = 0.0
        pixel_cont%meas(iw)%p44(:)	     = 0.0

        pixel_cont%meas(iw)%I(:)	      = 0.0
        pixel_cont%meas(iw)%Q(:)	      = 0.0
        pixel_cont%meas(iw)%U(:)	      = 0.0
        pixel_cont%meas(iw)%P(:)	      = 0.0
        pixel_cont%meas(iw)%LS(:)	      = 0.0
        pixel_cont%meas(iw)%DP(:)	      = 0.0

        pixel_cont%meas(iw)%CMTRX(:,:)	= 1.0 ! 0.0
        pixel_cont%meas(iw)%MPROF(:,:)	= 0.0

        pixel_cont%meas(iw)%IFCOV(:)    = 0
        pixel_cont%meas(iw)%IFMP(:)     = 0
      enddo ! iw
      
      pixel_cont%IFGAS	      = 0	
      pixel_cont%HOBS         = 0.0
      pixel_cont%nwl          = 0
      pixel_cont%cloudy       = 0
      pixel_cont%x            = 0.0
      pixel_cont%y            = 0.0
      pixel_cont%t            = 0
      pixel_cont%ix           = 0	
      pixel_cont%iy           = 0
      pixel_cont%it           = 0				

      pixel_cont%MASL	        = 0.0
      pixel_cont%land_percent = 0.0
      pixel_cont%irow         = 0
      pixel_cont%icol         = 0
      pixel_cont%HVP(:)       = 0.0

!	------------------------------------------------------------------------------------------------------

      end subroutine init_pixel

!	sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine init_segment_vec ( npixels, segment_vec  )
	
!  ------------------------------------------------------------------------------------------------------
      implicit none
!  ------------------------------------------------------------------------------------------------------
      integer,                             intent(in)     ::  npixels
      type(pixel_vector),dimension(KIMAGE),intent(inout)  ::  segment_vec
!  ------------------------------------------------------------------------------------------------------
      integer                    ::  ipix
!  ------------------------------------------------------------------------------------------------------
      do ipix=1,npixels
         segment_vec(ipix)%nFS(:)  = 0
         segment_vec(ipix)%FS(:)   = 0.0
         segment_vec(ipix)%KMIMAGE = 0
      enddo ! iw

      return
      end subroutine init_segment_vec

!	sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss


      subroutine set_pixel (							          &
                                    HOBS,           &         
                                    IX,						  &
                                    IY,						  &
                                    IT,						  &
                                    cloudy,					&
                                    irow,					  &
                                    icol,						&
                                    x_coord_longitude,	&
                                    y_coord_latitude,   &
                                    time_sec,			      &
                                    MASL,					      &
                                    land_percent,			  &
                                    nwave_lengths,			&
                                    nvalid_meas, 	      &
                                    wave_length, 		    &
                                    sza,                &
                                    thetav, 					  &
                                    phi, 						    &
                                    Nsurf, 					    &
                                    groundpar, 				  &
                                    IFGAS, 					    &
                                    gaspar, 					  &
                                    FPSW,						    &
                                    num_meas_type,      &
                                    meas_type,IFCOV,IFMP,  &
                                    CMTRX,MPROF,           &
                                    pixel_cont				     &					
                                   )
	
      implicit none
!  ------------------------------------------------------------------------------------------------------
		type(pixel),intent(out)	::	pixel_cont
!  ------------------------------------------------------------------------------------------------------
		real,intent(in)			    ::	HOBS
		integer,intent(in)			::	IX
		integer,intent(in)			::	IY
		integer,intent(in)			::	IT
		integer,intent(in)			::	cloudy
		integer,intent(in)			::	irow
		integer,intent(in)			::	icol
		integer,intent(in)			::	nwave_lengths
		real,intent(in)			    ::	y_coord_latitude
		real,intent(in)			    ::	x_coord_longitude
		integer(KIND_TIME),intent(in)	::	time_sec
		real,intent(in)			    ::	MASL
		real,intent(in)			    ::	land_percent

		integer,dimension(:,:),intent(in) ::  nvalid_meas
		integer,dimension(:,:),intent(in) ::  meas_type
		integer,dimension(:),intent(in)	  ::  num_meas_type
		integer,dimension(:,:),intent(in) ::  IFCOV,IFMP
		real,dimension(:),intent(in)	    ::  wave_length
		real,dimension(:),intent(in)	    ::  sza
		real,dimension(:,:,:),intent(in)  ::  thetav
		real,dimension(:,:,:),intent(in)  ::  phi

		integer,intent(in)			          ::  Nsurf
		integer,intent(in)			          ::  IFGAS

		real,dimension(:,:),intent(in)	     ::  groundpar
		real,dimension(:),intent(in)	       ::  gaspar
		real,dimension(:,:,:),intent(inout)  ::  FPSW
		real,dimension(:,:,:),intent(in)     ::  CMTRX,MPROF
!	------------------------------------------------------------------------------------------------------
		integer					             ::  iw  !,ii
!	------------------------------------------------------------------------------------------------------

		call init_pixel( pixel_cont)
			
    pixel_cont%HOBS = HOBS
      
		pixel_cont%icol = icol
		pixel_cont%irow = irow

		pixel_cont%y = y_coord_latitude
		pixel_cont%x = x_coord_longitude
		pixel_cont%t = time_sec

		pixel_cont%it = it
		pixel_cont%iy = iy
		pixel_cont%ix = ix
		pixel_cont%cloudy = cloudy

		pixel_cont%land_percent = land_percent
		pixel_cont%MASL = MASL

		pixel_cont%nwl = nwave_lengths

		do iw = 1,nwave_lengths
      call set_data_wl (iw,                            &
                        wave_length(iw),		           &
                        sza(iw),	                     &
                        thetav(:,:,iw),		             &
                        phi(:,:,iw),		               &
                        Nsurf,		                     &
                        groundpar(:,iw),		           &
                        IFGAS,	                       &
                        gaspar(iw),		                 &
                        nvalid_meas(:,iw),		         &
                        num_meas_type(iw),             &
                        meas_type(:,iw),               & 
                        FPSW(:,:,iw),		               &
                        IFCOV(:,iw),IFMP(:,iw),        &
                        CMTRX(:,:,iw),MPROF(:,:,iw),   &
                        pixel_cont                     &
                       )

      enddo ! iw

   end subroutine set_pixel
	
!	sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

   subroutine set_data_wl (iw,                 &
                           wl_value,		       &
                           sza,	               &
                           thetav,		         &
                           phi,		             &
                           Nsurf,		           &
                           groundpar,		       &
                           IFGAS,	             &
                           gaspar,		         &
                           nvalid_meas,		     &
                           num_meas_type,	     &
                           meas_type,          &
                           FPSW,               &
                           IFCOV,IFMP,         &
                           CMTRX,MPROF,		     &
                           pixel_cont          &
                          )
	
		implicit none
!  ------------------------------------------------------------------------------------------------------
		type(pixel),intent(out)	::	pixel_cont
!  ------------------------------------------------------------------------------------------------------
		integer,intent(in)					::	IFGAS
		integer,intent(in)					::	iw
		real,intent(in)						::	wl_value
		real,intent(in)		            ::	sza
		real,dimension(:,:),intent(in)	::	thetav
		real,dimension(:,:),intent(in)	::	phi
		integer,intent(in)					::	Nsurf ! ***Number of surface components ???
		real,dimension(:),intent(in)		::	groundpar
		real,             intent(in)		::	gaspar
		integer,dimension(:),intent(in)	::	nvalid_meas
		integer,intent(in)					::	num_meas_type
      integer,dimension(:),intent(in)  :: meas_type
		real,dimension(:,:),intent(inout)::	FPSW
      integer,dimension(:),intent(in)  :: IFCOV,IFMP      
		real,dimension(:,:),intent(in)   :: CMTRX,MPROF

!	------------------------------------------------------------------------------------------------------
      integer      ::  ii,IP
      integer      ::  ind ! ! ind=1 meas => pixel;  ind=2 pixel => meas
!	------------------------------------------------------------------------------------------------------

      !if(size(thetav).LT.nvalid_meas) then
        !write(*,*) 'Size(thetav)=',Size(thetav),' .LT.  NBVM=',NBVM 
        !stop 'stop in set_data_wl' 
	   !endif
				
      pixel_cont%meas(iw)%wl   = wl_value
      pixel_cont%meas(iw)%NIP  = num_meas_type		  
      pixel_cont%meas(iw)%sza  = sza
		
      do IP=1,num_meas_type
        pixel_cont%meas(iw)%meas_type(IP) = meas_type(IP)        
        pixel_cont%meas(iw)%NBVM(IP)      = nvalid_meas(IP)		
        pixel_cont%meas(iw)%IFCOV(IP)     = IFCOV(IP)
        pixel_cont%meas(iw)%IFMP(IP)      = IFMP(IP)      		
        do ii = 1,nvalid_meas(IP)
          if(meas_type(IP) .lt. meas_type_lid_beg .or. meas_type(IP) .gt. meas_type_lid_end) then
            !write(*,*) '1: IP=',IP,'  meas_type=', meas_type(IP)
            pixel_cont%meas(iw)%thetav(ii,IP) = thetav(ii,IP)
            pixel_cont%meas(iw)%phi(ii,IP)    = phi(ii,IP)
          else
            !write(*,*) '2: IP=',IP,'  meas_type=', meas_type(IP)
! For Lidar measurements thetav(ii,IP) contains values of altidudes in meters
! In current code version a set of heights should be the same for different pixels if 
! pixel contains lidar measurements.
! The order of heights and variables depending on heights have to be in order HVP(1) > HVP(nvalid_meas)
            pixel_cont%HVP(ii) = thetav(ii,IP) 
            if(thetav(ii,IP) .lt. thetav(nvalid_meas(IP),IP)) then
              write(*,*) 'IP=',IP,'  HVP(1)=',thetav(ii,IP),' .lt. HVP(nvalid_meas)=', &
                                               thetav(nvalid_meas(IP),IP)
              write(*,*) 'Heights should be decreasing order array'
              stop 'stop in set_data_wl'        
            endif      
            pixel_cont%meas(iw)%MPROF(ii,IP) = MPROF(ii,IP)
          endif ! meas_type(IP) .lt. meas_type_lid_beg .or.          
          pixel_cont%meas(iw)%CMTRX(ii,IP) = CMTRX(ii,IP)
        enddo ! ii
      enddo ! IP
      
      ind = 1 ! meas => pixel
      do IP=1,num_meas_type         
      call set_pixel_meas (                             & ! IN
                           iw,                          &
                           nvalid_meas(IP),             &
                           meas_type(IP),               &
                           ind,                         &
                           FPSW(1:nvalid_meas(IP),IP),  & ! INOUT
                           pixel_cont                   &
                          )
      enddo ! IP
      
      pixel_cont%meas(iw)%Nsurf = Nsurf		  
      do ii = 1,Nsurf
		pixel_cont%meas(iw)%groundpar(ii) = groundpar(ii)
      enddo ! ii
		  
      pixel_cont%IFGAS = IFGAS		  
      if(IFGAS .eq. 1) &
		pixel_cont%meas(iw)%gaspar = gaspar

   end subroutine set_data_wl

!	sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

   subroutine set_pixel_meas (            & ! IN
                              iw,         &
                              NBVM,       &
                              meas_type,  &
                              ind,        &
                              meas,	      & ! INOUT
                              pixel_cont  &
                             )
	
      implicit none
!  ---------------------------------------------------------------------------------------
		integer,intent(in)					::	iw
		integer,intent(in)					::	NBVM
		integer,intent(in)					::	meas_type
		integer,intent(in)					::	ind ! ind=1 meas  => pixel (set);  
                                        ! ind=2 pixel => meas (get)      
!  ---------------------------------------------------------------------------------------
		type(pixel),intent(inout)	            :: pixel_cont
		real,dimension(NBVM),intent(inout)	  :: meas
!	----------------------------------------------------------------------------------------
      integer      ::  iv
!	----------------------------------------------------------------------------------------
	
	  !if(size(meas).LT.NBVM) then
        !write(*,*) 'Size(meas)=',Size(meas),' .LT.  NBVM=',NBVM 
        !stop 'stop in set_pixel_meas' 
	  !endif

      select case(meas_type)
      case(meas_type_aod)
         if(ind .eq. 1) then         
         pixel_cont%meas(iw)%tau(1:NBVM) = meas(1:NBVM)      
         else
         meas(1:NBVM) = pixel_cont%meas(iw)%tau(1:NBVM)      
         endif ! ind
      case(meas_type_tod)
         if(ind .eq. 1) then         
         pixel_cont%meas(iw)%tau(1:NBVM) = meas(1:NBVM)      
         else
         meas(1:NBVM) = pixel_cont%meas(iw)%tau(1:NBVM)      
         endif ! ind
      case(meas_type_p11)
         if(ind .eq. 1) then         
         pixel_cont%meas(iw)%p11(1:NBVM) = meas(1:NBVM)      
         else
         meas(1:NBVM) = pixel_cont%meas(iw)%p11(1:NBVM)      
         endif ! ind
      case(meas_type_p12)
         if(ind .eq. 1) then         
         pixel_cont%meas(iw)%p12(1:NBVM) = meas(1:NBVM)      
         else
         meas(1:NBVM) = pixel_cont%meas(iw)%p12(1:NBVM)      
         endif ! ind
      case(meas_type_p22)
         if(ind .eq. 1) then         
         pixel_cont%meas(iw)%p22(1:NBVM) = meas(1:NBVM)      
         else
         meas(1:NBVM) = pixel_cont%meas(iw)%p22(1:NBVM)      
         endif ! ind
      case(meas_type_p33)
         if(ind .eq. 1) then         
         pixel_cont%meas(iw)%p33(1:NBVM) = meas(1:NBVM)      
         else
         meas(1:NBVM) = pixel_cont%meas(iw)%p33(1:NBVM)      
         endif ! ind
      case(meas_type_p34)
         if(ind .eq. 1) then         
         pixel_cont%meas(iw)%p34(1:NBVM) = meas(1:NBVM)      
         else
         meas(1:NBVM) = pixel_cont%meas(iw)%p34(1:NBVM)      
         endif ! ind
      case(meas_type_p44)
         if(ind .eq. 1) then         
         pixel_cont%meas(iw)%p44(1:NBVM) = meas(1:NBVM)      
         else
         meas(1:NBVM) = pixel_cont%meas(iw)%p44(1:NBVM)      
         endif ! ind
      case(meas_type_LS)
         if(ind .eq. 1) then 
         pixel_cont%meas(iw)%LS(1:NBVM)  = meas(1:NBVM)      
         else
         meas(1:NBVM) = pixel_cont%meas(iw)%LS(1:NBVM)      
         endif ! ind
      case(meas_type_DP)
         if(ind .eq. 1) then 
         pixel_cont%meas(iw)%DP(1:NBVM)  = meas(1:NBVM)      
         else
         meas(1:NBVM) = pixel_cont%meas(iw)%DP(1:NBVM)      
         endif ! ind
      case(meas_type_RL)
         if(ind .eq. 1) then 
         pixel_cont%meas(iw)%RL(1:NBVM)  = meas(1:NBVM)      
         else
         meas(1:NBVM) = pixel_cont%meas(iw)%RL(1:NBVM)      
         endif ! ind
      case(meas_type_I)
         if(ind .eq. 1) then 
         pixel_cont%meas(iw)%I(1:NBVM)   = meas(1:NBVM)      
         else
         meas(1:NBVM) = pixel_cont%meas(iw)%I(1:NBVM)         
         endif ! ind
      case(meas_type_Q)
         if(ind .eq. 1) then 
         pixel_cont%meas(iw)%Q(1:NBVM)   = meas(1:NBVM)      
         else
         meas(1:NBVM) = pixel_cont%meas(iw)%Q(1:NBVM)      
         endif ! ind
      case(meas_type_U)
         if(ind .eq. 1) then 
         pixel_cont%meas(iw)%U(1:NBVM)   = meas(1:NBVM)      
         else
         meas(1:NBVM) = pixel_cont%meas(iw)%U(1:NBVM)      
         endif ! ind
      case(meas_type_P)
         if(ind .eq. 1) then 
         pixel_cont%meas(iw)%P(1:NBVM)   = meas(1:NBVM)      
         else
         meas(1:NBVM) = pixel_cont%meas(iw)%P(1:NBVM)      
         endif ! ind
      case default
         write(tmp_message,'(2(a,i0,3x),a)') &
         'iw = ',iw,'meas_type = ',meas_type,'- Unknown value of meas_type'
         G_ERROR(trim(tmp_message))
      end select
            
   end subroutine set_pixel_meas

!	sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! Assign general wave length indices to pixel wave lengths
      subroutine set_segment_pixel_wl_index ( iu_main_output, RIN, segment )

      use mod_retr_settings_derived_type

      implicit none
!	---------------------------------------------------------------------------------------
      integer,                   intent(in)     :: iu_main_output
      type(retr_input_settings), intent(in)     :: RIN
      type(segment_data),        intent(inout)	:: segment
!  --------------------------------------------------------------------------------------
!  --------------------------------------------------------------------------------------
      integer :: nwl
      !real,   dimension(KWM)           :: wl   
      !integer,dimension(KWM)           :: ind_wl   
      real,   dimension(KW) :: wl
      integer,dimension(KW) :: ind_wl
      integer :: iw,ipix
      real :: tiny = 1e-3
      logical :: status
!  --------------------------------------------------------------------------------------

      do ipix=1,segment%npixels

      call get_pixel_wl (  segment%pixels(ipix),  & ! IN
                           nwl,wl                 & ! out
                        )

      ! search if pixel wavelengths present in wavelength set for
      ! retrieved parameters RIN%WAVE(1:RIN%NW) and define their indices 
      ! ind_wl(1:nwl) in RIN%WAVE(1:RIN%NW) array
      call R_subset_index ( tiny,                       &
                            RIN%NW, RIN%WAVE(1:RIN%NW), &
                            nwl, wl(1:nwl),             &
                            ind_wl(1:nwl), status       & ! OUT
                          )
      if(.not. status) then
        write(tmp_message,'(a,i0,a,a,20f7.4,2a,20f7.4)') &
        'Can not fill in indices for pixel wavelengths, ipix = ',ipix, &
        NEW_LINE('A'), &
        'WL: ',RIN%WAVE(1:RIN%NW), &
        NEW_LINE('A'), &
        'wl: ',wl(1:nwl)
        G_ERROR(trim(tmp_message))
      endif

      do iw=1,nwl
         segment%pixels(ipix)%meas(iw)%ind_wl = ind_wl(iw)
      enddo ! iw                              

      enddo ! ipix
                                              
      end subroutine set_segment_pixel_wl_index

!	sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

	subroutine get_pixel_wl (                  &
                              pixel_cont,    & ! IN
                              nwl,           & ! OUT
                              wl_val,        & 
                              ind_wl         &   
                           )
      implicit none
!	---------------------------------------------------------------------------------------
      type(pixel),intent(in)	         :: pixel_cont
!  --------------------------------------------------------------------------------------
      integer,             intent(out)            :: nwl
      real,dimension(:),   intent(out), optional  :: wl_val   
      integer,dimension(:),intent(out), optional  :: ind_wl   

!	----------------------------------------------------------------------------------------
      integer             ::  iw
!	  character(LEN=100)  ::  msg
!	----------------------------------------------------------------------------------------

      nwl = pixel_cont%nwl
      !if(size(wl_val).LT.nwl) then
         !write(*,*) 'Size(wl_val) .LT. NWL'
         !stop 'stop in get_pixel_wl'
      !endif

      do iw=1,nwl						
         if(present(wl_val)) wl_val(iw) = pixel_cont%meas(iw)%wl
         if(present(ind_wl)) ind_wl(iw) = pixel_cont%meas(iw)%ind_wl
      enddo ! iw
      
	end subroutine get_pixel_wl

!	sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

	subroutine get_pixel_geom (  iw,ip,        & ! IN
	                             pixel_cont,   &
	                             NBVM,         & ! OUT
                               sza,          & 
	                             thetav,       & 
                               phi           &  
                             )
      implicit none
!	---------------------------------------------------------------------------------------
      type(pixel),intent(in)   :: pixel_cont
      integer,intent(in)               :: iw,ip
!  --------------------------------------------------------------------------------------
      integer,intent(out)              :: NBVM
      real,   intent(out)              :: sza
      real,dimension(:),intent(out)    :: thetav   
      real,dimension(:),intent(out)    :: phi   
!	---------------------------------------------------------------------------------------
      integer             ::  iv
!	    character(LEN=100)  ::  msg
!	---------------------------------------------------------------------------------------
      thetav(:) = 0.
      phi(:)    = 0.
	  
      NBVM = pixel_cont%meas(iw)%NBVM(ip)	  
      !if(size(thetav).LT.NBVM) then
         !write(*,*) 'Size(thetav) .LT. NBVML'
         !stop 'stop in get_pixel_geom'
      !endif
      !if(size(phi).LT.NBVM) then
         !write(*,*) 'Size(phi) .LT. NBVML'
         !stop 'stop in get_pixel_geom'
      !endif
	  
      sza = pixel_cont%meas(iw)%sza
      do iv=1,NBVM
      thetav(iv) = pixel_cont%meas(iw)%thetav(iv,ip)
      phi(iv)    = pixel_cont%meas(iw)%phi(iv,ip)
      enddo ! iv

	end subroutine get_pixel_geom

!	sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! TL
      subroutine get_pixel_txy (                   &
	                              pixel_cont,        &
	                              irow,icol,x,y,t    &  ! X - lon, Y - lat
							                 )
      implicit none
!	---------------------------------------------------------------------------------------
      type(pixel),intent(in)	:: pixel_cont
!  --------------------------------------------------------------------------------------
      integer,   intent(out)    :: irow,icol
      real,      intent(out)    :: x,y  ! X - lon, Y - lat 
      integer(KIND_TIME), intent(out)    :: t   !universal date and time stamp in sec
!	---------------------------------------------------------------------------------------
!	---------------------------------------------------------------------------------------

      irow  = pixel_cont%irow						
      icol  = pixel_cont%icol
      x     = pixel_cont%x
      y     = pixel_cont%y
      t     = pixel_cont%t

      end subroutine get_pixel_txy

!	sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! TL
!      subroutine radiance_correction ( iu_main_output,RIN,segment )

!      use mod_retr_settings_derived_type

!      implicit none
!	---------------------------------------------------------------------------------------
!      integer,                       intent(in)    ::  iu_main_output
!      type(segment_data),            intent(inout) ::  segment
!      type(retr_input_settings),     intent(in)    ::  RIN
  
!  --------------------------------------------------------------------------------------  
!      integer  ::  ipix, ip, iw
!	---------------------------------------------------------------------------------------

!      ip=1
!      iw = RIN%I_corr_iwl
!      do ipix=1,segment%npixels
!        segment%pixels(ipix)%meas(iw)%I(1:segment%pixels(ipix)%meas(iw)%NBVM(ip)) =  &
!        segment%pixels(ipix)%meas(iw)%I(1:segment%pixels(ipix)%meas(iw)%NBVM(ip))    &
!                                                                * RIN%coeff_corr
!      enddo ! ipix
      
!      if(.not. RIN%IPRI_silent_mode)  write(iu_main_output,*) 'Warning !!!'
!      if(.not. RIN%IPRI_silent_mode)  write(iu_main_output,'(a,i3,a,f5.2)')   &
!         'PARASOL radiance measurements at #iw=',iw,  &
!         ' have been corrected with coeff =',RIN%coeff_corr
!      return
!      end subroutine radiance_correction

!	ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

! Fill in measurement vector for segment

      subroutine set_segment_meas_vector_FS ( RIN, INVSING, segment_meas, segment_vec_meas )
	  
      use mod_retr_settings_derived_type

      implicit none
!	----------------------------------------------------------------------------
! IN :	
      type(retr_input_settings), intent(in) ::  RIN
      integer,                   intent(in) ::  INVSING
      type(segment_data),        intent(in) ::  segment_meas
!	----------------------------------------------------------------------------
! INOUT :	
!	----------------------------------------------------------------------------
! OUT :
      type(pixel_vector),dimension(KIMAGE),intent(out) ::  segment_vec_meas
!	-----------------------------------------------------------------------------
! LOCAL :
      integer                              ::  IWb,IWe,ipix
      integer                              ::  npixels,KMIMAGE
!	-----------------------------------------------------------------------------
      npixels = segment_meas%npixels
      
      call init_segment_vec ( npixels, segment_vec_meas )

      IWB = 1
      do ipix=1,npixels
        IWE = segment_meas%pixels(ipix)%nwl
! Fill in pixel measurement vector
        call set_pixel_meas_vector_FS ( RIN,IWB,IWE,ipix,           & ! IN
                                        segment_meas%pixels(ipix),  & 
                                        segment_vec_meas(ipix)      & ! INOUT
                                      )
        if ( stop_report%status ) return
! Calculate the number of measurements in each pixel
        segment_vec_meas(ipix)%KMIMAGE = SUM(segment_vec_meas(ipix)%nFS(IWB:IWE))           
      enddo ! ipix

      return
      end subroutine set_segment_meas_vector_FS

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! Add random noise to segment measurement vector
! INOISE  - the number of different noise sources              
! SGMS(I) - std of noise in i -th source                      
! INN(I)  - EQ.1.THEN error is absolute with                  
!         - EQ.0 THEN error assumed relative
! DNN(I)  - variation of the noise of the I-th source

      subroutine add_rnoise_segment ( RIN,deep_random_switch,  & ! IN
                                      segment_meas,            &  
                                      segment_vec_meas,        & ! INOUT
                                      MNOISEI                  & ! IN                                      
                                    ) 
	  
      use mod_retr_settings_derived_type
      
      implicit none
!	----------------------------------------------------------------------------
! IN :	
      type(retr_input_settings),        intent(in) ::  RIN
      logical,                          intent(in) ::  deep_random_switch
      type(segment_data),               intent(in) ::  segment_meas
      integer,dimension(KIP,KWM,KIMAGE),intent(in) ::  MNOISEI   
!	----------------------------------------------------------------------------
! INOUT :	
      type(pixel_vector),dimension(KIMAGE), intent(inout) ::  segment_vec_meas
!	----------------------------------------------------------------------------
! LOCAL :
      integer                              ::  ipix,JJS
      integer                              ::  npixels
!	------------------------------------------------------------------------------------------
      npixels = segment_meas%npixels

! Disturb measurements with random noise
      do ipix=1,npixels
        JJS = SUM(segment_vec_meas(ipix)%nFS(1:segment_meas%pixels(ipix)%nwl))      
        if(RIN%KL .eq. 1) &
        segment_vec_meas(ipix)%FS(1:JJS) = EXP(segment_vec_meas(ipix)%FS(1:JJS))
        call add_rnoise_pixel ( deep_random_switch,          & ! IN
                                RIN, MNOISEI(:,:,ipix),      &
                                segment_meas%pixels(ipix),   & 
                                segment_vec_meas(ipix)       & ! INOUT
                              )
        if(RIN%KL .eq. 1) &
        segment_vec_meas(ipix)%FS(1:JJS) = LOG(segment_vec_meas(ipix)%FS(1:JJS))
      enddo ! ipix

      return
      end subroutine add_rnoise_segment

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss 

! ***************************************************
! **       Calculating the number of 
! **       the measurements KM,
! **       preparing vector of measurements FP 
! ***************************************************

      subroutine set_pixel_meas_vector_FS ( RIN, IWb, IWe, ipix, & ! IN
                                            pixel_cont,          & ! IN
                                            pixel_vec            & ! INOUT
                                          )
	  
      use mod_retr_settings_derived_type
      use mod_inversion_utils, only  : check_nan

      implicit none
!	----------------------------------------------------------------------------
! IN :	
      type(retr_input_settings),  intent(in) :: RIN
      integer,                    intent(in) :: IWb,IWe,ipix
      type(pixel),                intent(in) :: pixel_cont
!	----------------------------------------------------------------------------
! INOUT :	
      type(pixel_vector),         intent(inout) :: pixel_vec
!	----------------------------------------------------------------------------
! OUT :		
!	-----------------------------------------------------------------------------
! LOCAL :
      integer             ::  NWL, NIP, meas_type
      integer             ::  JJS, IW, IV, IP, JJSb
      real                ::  meas, I, Q2, U2 
      integer             ::  nFS
      logical             ::  LERR
!	-----------------------------------------------------------------------------
! iPOBS = 
!           1    I, Q, U  or P11,P12,P22,P33,P34,P44
!           2    I, Q/I, U/I  or p11,-P12/P11,P22/p11,P33/P11,P34/P11,P44/P11
!           3    I, P    or sqrt(Q*Q+U*U)
!           4    I, P/I  or sqrt(Q*Q+U*U)/I
!           5    I, P/I meas
! iIOBS =
!           1    I
!           2    I/sum(I(:))

! Types of measurements :
! 11 - tod(wl)       - total optical depth (aer+mol+gas)
! 12 - aod(wl)       - aerosol optical depth
! 21 - p11(angle,wl)  - phase matrix element p11
! 22 - p12(angle,wl)  - phase matrix element p12
! 23 - p22(angle,wl)  - phase matrix element p22
! 24 - p33(angle,wl)  - phase matrix element p33
! 25 - p34(angle,wl)  - phase matrix element p34
! 26 - p44(angle,wl)  - phase matrix element p44

! 31  - LS(height,wl)  - lidar signal
! 32  - DP(height,wl)  - depolarization ratio
! 33  - RL(height,wl)  - Raman lidar signal

! 41 - I(angle,wl)    - Stokes parameter I
! 42 - Q(angle,wl)    - Stokes parameter Q
! 43 - U(angle,wl)    - Stokes parameter U
! 44 - P(angle,wl)    - polarization sqrt(Q*Q+U*U) or P/I(iPOBS=5)
!	-----------------------------------------------------------------------------      

      JJS  = SUM(pixel_vec%nFS(1:IWb))-pixel_vec%nFS(IWb)
      JJSb = JJS

      NWL = pixel_cont%NWL
      LOOP_WL : do iw=IWb,IWe
      nFS = 0
      NIP = pixel_cont%meas(iw)%NIP
      LOOP_meas_type : do IP=1,NIP
         meas_type = pixel_cont%meas(iw)%meas_type(IP)
         if(RIN%iPOBS .ge. 3 .and. meas_type .eq. meas_type_U) &
         exit LOOP_meas_type
         do IV=1,pixel_cont%meas(iw)%NBVM(IP)
         JJS=JJS+1              
         select case(meas_type)
         case(meas_type_aod)
            pixel_vec%FS(JJS) = pixel_cont%meas(iw)%tau(IV)
            nFS = nFS+1
         case(meas_type_tod)
            pixel_vec%FS(JJS) = pixel_cont%meas(iw)%tau(IV)
            nFS = nFS+1
         case(meas_type_p11)
            pixel_vec%FS(JJS) = pixel_cont%meas(iw)%p11(IV)  
            nFS = nFS+1
         case(meas_type_p12) 
            select case(RIN%iPOBS)
            case(1,5)
              pixel_vec%FS(JJS) = pixel_cont%meas(iw)%p12(IV) + RIN%SHIFT
            case(2)
              pixel_vec%FS(JJS) = -pixel_cont%meas(iw)%p12(IV)/pixel_cont%meas(iw)%p11(IV) + RIN%SHIFT
            end select
            nFS = nFS+1
         case(meas_type_p22)
            select case(RIN%iPOBS)
            case(1)
              pixel_vec%FS(JJS) = pixel_cont%meas(iw)%p22(IV)
            case(2)
              pixel_vec%FS(JJS) = pixel_cont%meas(iw)%p22(IV)/pixel_cont%meas(iw)%p11(IV)
            end select
            nFS = nFS+1
         case(meas_type_p33)
            select case(RIN%iPOBS)
            case(1)
              pixel_vec%FS(JJS) = pixel_cont%meas(iw)%p33(IV) + RIN%SHIFT
            case(2)
              pixel_vec%FS(JJS) = pixel_cont%meas(iw)%p33(IV)/pixel_cont%meas(iw)%p11(IV) + RIN%SHIFT
            end select
            nFS = nFS+1
         case(meas_type_p34)
            select case(RIN%iPOBS)
            case(1)
              pixel_vec%FS(JJS) = pixel_cont%meas(iw)%p33(IV) + RIN%SHIFT
            case(2)
              pixel_vec%FS(JJS) = pixel_cont%meas(iw)%p33(IV)/pixel_cont%meas(iw)%p11(IV) + RIN%SHIFT
            end select
            nFS = nFS+1
         case(meas_type_p44)
            select case(RIN%iPOBS)
            case(1)
              pixel_vec%FS(JJS) = pixel_cont%meas(iw)%p44(IV) + RIN%SHIFT
            case(2)
              pixel_vec%FS(JJS) = pixel_cont%meas(iw)%p44(IV)/pixel_cont%meas(iw)%p11(IV) + RIN%SHIFT
            end select
            nFS = nFS+1
         case(meas_type_LS)
            pixel_vec%FS(JJS) = pixel_cont%meas(iw)%LS(IV)
            nFS = nFS+1
         case(meas_type_DP) 
            pixel_vec%FS(JJS) = pixel_cont%meas(iw)%DP(IV)
            nFS = nFS+1
         case(meas_type_RL) 
            pixel_vec%FS(JJS) = pixel_cont%meas(iw)%RL(IV)
            nFS = nFS+1
         case(meas_type_I) 
            select case(RIN%iIOBS)
            case(1)
               pixel_vec%FS(JJS) = pixel_cont%meas(iw)%I(IV)
            case(2)
               pixel_vec%FS(JJS) = pixel_cont%meas(iw)%I(IV)/ &
               sum(pixel_cont%meas(iw)%I(1:pixel_cont%meas(iw)%NBVM(IP)))
            end select
            nFS = nFS+1
         case(meas_type_Q)
            select case(RIN%iPOBS)
            case(1) 
               pixel_vec%FS(JJS) = pixel_cont%meas(iw)%Q(IV) + RIN%SHIFT
            case(2) 
               I = pixel_cont%meas(iw)%I(IV)
               pixel_vec%FS(JJS) = pixel_cont%meas(iw)%Q(IV)/I + RIN%SHIFT
            case(3) 
               meas = pixel_cont%meas(iw)%Q(IV)
               Q2 = meas*meas
               meas = pixel_cont%meas(iw)%U(IV)
               U2 = meas*meas
               pixel_vec%FS(JJS) = sqrt(Q2+U2)
            case(4)
               I = pixel_cont%meas(iw)%I(IV)
               meas = pixel_cont%meas(iw)%Q(IV)
               Q2 = meas*meas
               meas = pixel_cont%meas(iw)%U(IV)
               U2 = meas*meas
               pixel_vec%FS(JJS) = sqrt(Q2+U2)/I
            end select
            nFS = nFS+1
         case(meas_type_U)
            select case(RIN%iPOBS)
            case(1) 
               pixel_vec%FS(JJS) = pixel_cont%meas(iw)%U(IV) + RIN%SHIFT
            case(2) 
               I = pixel_cont%meas(iw)%I(IV)
               pixel_vec%FS(JJS) = pixel_cont%meas(iw)%U(IV)/I + RIN%SHIFT
            end select
            nFS = nFS+1
         case(meas_type_P)
            select case(RIN%iPOBS)
            case(3,5)
!tl 25-02-2015
! iPOBS=5 - in case when number of angles or angle values for I differe from ones for P
! linear polarization has been divided by I at the time of setting
! modeled measurements in subroutine set_pixel_Stokes_vec_fit
               pixel_vec%FS(JJS) = pixel_cont%meas(iw)%P(IV)
            case(4) 
               I = pixel_cont%meas(iw)%I(IV)
               pixel_vec%FS(JJS) = pixel_cont%meas(iw)%P(IV)/I
            case default
               write(tmp_message,'(a,i0,a)') 'RIN%iPOBS = ',RIN%iPOBS,'  - value is not valid'
               G_ERROR(trim(tmp_message))
            end select
            nFS = nFS+1
         case default
            write(tmp_message,'(a,i0,a)') 'meas_type = ',meas_type,'  - unknown value of meas_type'
            G_ERROR(trim(tmp_message))
         end select ! meas_type
         if(RIN%KL .eq. 1) then
            if(pixel_vec%FS(JJS) .le. 0.0) then
              if( (meas_type .ge. meas_type_Q   .and. meas_type .le. meas_type_U .and. RIN%iPOBS .le. 2) .or.  &
                  (meas_type .eq. meas_type_p12) .or.                                   &
                  (meas_type .ge. meas_type_p33 .and. meas_type .le. meas_type_p44)     &
                ) then
                write(tmp_message,'(a,2(2x,a,es11.4),a,4(a,i0))') &
                '!!! Measurement vector problem  FS <=0 !!!','  FS=F+shift =',pixel_vec%FS(JJS), &
                'shift for applying logarithm to negative measurements =',RIN%SHIFT, &
                NEW_LINE('A'), &
                'pixel # ',ipix,'  wl # ',iw,'  meas type # ',ip,'  meas # ',iv
                G_ERROR(trim(tmp_message))
              else
                write(tmp_message,'(a,2x,a,es11.4,a,4(a,i0))') '!!! Measurement vector problem FS<=0 !!!', &
                'FS =',pixel_vec%FS(JJS), &
                NEW_LINE('A'), &
                'pixel # ',ipix,'  wl # ',iw,'  meas type # ',ip,'  meas # ',iv
                G_ERROR(trim(tmp_message))
              endif
            endif
         endif ! KL .eq. 1
         enddo ! IV
      enddo LOOP_meas_type
      pixel_vec%nFS(iw) = nFS
      enddo LOOP_WL
      if(RIN%KL .eq. 1) then 
        !if(any(pixel_vec%FS(JJSb+1:JJS) .le. 0.0)) then
          !write(*,*) 'FS:'
          !write(*,'(10e13.5)') pixel_vec%FS(JJSb+1:JJS)
          !write(*,*) '!!! Meas vector FS<=0 problem !!! ipix=',ipix,'  in set_pixel_meas_vector_FS'
          !stop 'stop in set_pixel_meas_vector_FS'
        !endif
        iv = SUM(pixel_vec%nFS(IWb:IWe))
        LERR = check_nan(iv,pixel_vec%FS(JJSb+1:JJS))
        if(.not. LERR) then
          if(RIN%IPRI_verbose .eqv. .true.) then
            write(*,*) 'FS:'
            write(*,'(10e13.5)') pixel_vec%FS(JJSb+1:JJS)
          endif
          write(tmp_message,'(a,i0,2x,a)') '!!! Meas vector ln(FS)=NaN problem !!! ipix=',ipix, &
          'in set_pixel_meas_vector_FS'
          G_ERROR(trim(tmp_message))
        endif
      pixel_vec%FS(JJSb+1:JJS) = LOG(pixel_vec%FS(JJSb+1:JJS))
      endif ! RIN%KL .eq. 1

      return
      end subroutine set_pixel_meas_vector_FS

!	sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! Add random noise to measurements (simulated measurements)
      subroutine add_rnoise_pixel ( deep_random_switch, & ! IN
                                    RIN, MNOISEI,       &
                                    pixel_cont,         & 
                                    pixel_vec           & ! INOUT
                                  )
      use mod_retr_settings_derived_type

      implicit none
! -----------------------------------------------------------------
! IN :
      logical,                   intent(in) :: deep_random_switch	   
      type(retr_input_settings), intent(in) :: RIN
      type(pixel),               intent(in) :: pixel_cont
      integer,dimension(KIP,KWM),intent(in) :: MNOISEI   

! -----------------------------------------------------------------
! INOUT :	  	  
      type(pixel_vector), intent(inout) :: pixel_vec
! -----------------------------------------------------------------
! LOCAL :
      integer      :: iw,IP,iv,iMN,NWL,NIP,meas_type,JJS
! ----------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------
! Modeling RANDOM NOISE                      
! INOISE  - the number of different noise sources           
! SGMS(I) - std of noise in i -th source                    

      JJS = 0
      NWL = pixel_cont%nwl
LOOP_WL : do iw=1,NWL
      NIP = pixel_cont%meas(iw)%NIP
LOOP_meas_type : do IP=1,NIP
         meas_type = pixel_cont%meas(iw)%meas_type(IP)
         iMN = MNOISEI(IP,iw)  
         if(iMN .lt. 1 .or. iMN .gt. RIN%NOISE%INOISE) then
            write(tmp_message,'(2(a,i0),2a)') 'iMN = ',iMN,'  INOISE = ',RIN%NOISE%INOISE, &
            NEW_LINE('A'), &
            'Noise index is not in valid range 1<iMN<INOISE'
            G_ERROR(trim(tmp_message))
         endif !
         if(RIN%iPOBS .ge. 3 .and. meas_type .eq. meas_type_U) & ! fit sqrt(Q*Q+U*U), sqrt(Q*Q+U*U)P/I
         exit LOOP_meas_type

         if(RIN%NOISE%SGMS(iMN) .gt. 0.0) then
           do IV=1,pixel_cont%meas(iw)%NBVM(IP)
           JJS = JJS+1
             if( (meas_type .ge. meas_type_p33 .and. meas_type .le. meas_type_p44) .or. &
                  meas_type .eq. meas_type_p12 ) then
               pixel_vec%FS(JJS) = pixel_vec%FS(JJS) - RIN%SHIFT
               call rnoise (deep_random_switch,RIN%NOISE%SGMS(iMN),RIN%NOISE%INN(iMN),pixel_vec%FS(JJS))
               pixel_vec%FS(JJS) = pixel_vec%FS(JJS) + RIN%SHIFT
             elseif( meas_type .eq. meas_type_Q .and. RIN%iPOBS .lt. 3) then
               pixel_vec%FS(JJS) = pixel_vec%FS(JJS) - RIN%SHIFT
               call rnoise (deep_random_switch,RIN%NOISE%SGMS(iMN),RIN%NOISE%INN(iMN),pixel_vec%FS(JJS))
               pixel_vec%FS(JJS) = pixel_vec%FS(JJS) + RIN%SHIFT
             else
               call rnoise (deep_random_switch,RIN%NOISE%SGMS(iMN),RIN%NOISE%INN(iMN),pixel_vec%FS(JJS))
             endif
           enddo ! IV
         else
           do IV=1,pixel_cont%meas(iw)%NBVM(IP)
           JJS = JJS+1
           enddo ! IV
         endif ! RIN%NOISE%SGMS(iMN) .gt. 0.0
enddo LOOP_meas_type
enddo LOOP_WL

      return
      end subroutine add_rnoise_pixel

!	sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine rnoise ( deep_random_switch, & ! IN
                          SGMS,INN,           &
                          FS                  & ! INOUT
                        )
      implicit none
! -----------------------------------------------------------------
! IN :
      logical, intent(in) :: deep_random_switch
      real,    intent(in) :: SGMS
      integer, intent(in) :: INN	   
! -----------------------------------------------------------------
! INOUT :
      real,    intent(inout) :: FS
! -----------------------------------------------------------------
! LOCAL :
      real              :: EMG=0.
      real,dimension(1) :: RDM
      real*8            :: RDM0(1)=0.4
      real              :: FS0
! ----------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------
! Modeling   RANDOM NOISE                       
! SGMS - std of noise in i -th source                 
! INN  - EQ.1.THEN error is absolute with             

         if (deep_random_switch) then
!  Random number generator to provide more deep randomization	
            call RDMU2(1,RDM0(1))
         endif ! end of if (deep_random_switch)
            
         call RDMG(RDM0(1), 1, EMG, SGMS, RDM(1))
         FS0 = FS
         select case(INN)
         case(0)
            FS = FS * (1.0+RDM(1))
         case(1)
            FS = FS + RDM(1)
         end select
         !write(*,'(2(a,i3),4(a,e12.5))') 'in rnoise: INN=',INN, &
         !               '  RDM0=',RDM0(1),'  RDM=',RDM(1),'  FS0=',FS0,'  FS=',FS
         
      return
      end subroutine rnoise

!	sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! cmatrix
!	sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

   subroutine write_simul_pixels (  iu_main_output, &
                                    sdata_sim_file, &
                                    npixels,        &
                                    pixel_cont,     &
                                    index_clouds    &
                                 )

      use mod_par_inv,  only : KIMAGE

      implicit none
! -----------------------------------------------------------------------------------
! IN : 

      integer,                       intent(in)    ::  iu_main_output,npixels
      character(*),                  intent(in)    ::  sdata_sim_file
      type(pixel),dimension(KIMAGE), intent(inout) ::  pixel_cont
      type(ind_clouds),              intent(in)	   ::  index_clouds			
	  
! -----------------------------------------------------------------------------------
! LOCAL :
      integer             :: iu_SDATA
      integer             :: IT,IX,IY,IW,IV,Nsurf,IFGAS,npixels1,k,IBRF,IP
      character           :: file_name*20,NWLC*20,CFMT*150
      real                :: HOBS,tau,I,Q,U,Q2,U2,P,LS,DP,RL, &
                             p11, p12, p22, p33, p34, p44
      integer             :: NBVM,NIP
      integer             :: NWL,IIMAGE,NT,NX,NY
      integer             :: meas_type,ind
      integer,dimension(KIP,KWM)         ::  IFMP     
      real,dimension(KNBVM,KIP,KWM)   ::  FPSW
      real,dimension(KNBVM,KIP,KWM)   ::  thetav
      real,dimension(KNBVM,KIP,KWM)   ::  phi
      character(len=20)               ::  iso8601_string
      logical                         ::  status_funct
!	------------------------------------------------------------------------------------------------------
      open (newunit=iu_SDATA, FILE=trim(sdata_sim_file), STATUS='replace', &
                                                   FORM='formatted',RECL=2000000)      
!	------------------------------------------------------------------------------------------------------

      NT = index_clouds%NT
      NX = index_clouds%NX
      NY = index_clouds%NY
      npixels1 = NT*NX*NY
	   
      if(npixels1 .ne. npixels) THEN
         write(tmp_message,'(2(a,i0))') &
         'npixels1 = ',npixels1,' .ne. npixels1 = ',npixels
         G_ERROR(trim(tmp_message))
      endif
	  
      Nsurf    = pixel_cont(1)%meas(1)%Nsurf
      IFGAS     = pixel_cont(1)%IFGAS
	  
#ifdef WARN_DRY
#warning "__SDATA_VERSION__ duplicated"
#endif                 
      WRITE(iu_SDATA,*) 'SDATA version 2.0'
      write (iu_SDATA,'(i3,2i4,a)') NX,NY,NT,'  : NX NY NT'

!C*****************************************************
!C***  The parameters define TIME and SPATIAL 
!C***   distribution of PIXELS considered in the retrieval    
!C*****************************************************
!C*** NT - number of different images for the same area
!C*** NX - number of the pixel raws in each image
!C*** NY - number of the pixels in each pixel raw
!C*** X(KITIME,KX,KY),y(KITIME,KX,KY) - X and Y coordinates of each pixel
!C*** ICLOUD(KITIME,KX,KY) - cloud index for each pixel
!C***        =0  - cloud free
!C***        =1  - cloudy 
!C*****************************************************      

	   DO IT=1,NT     
      IIMAGE = index_clouds%INIMAGE(1,1,IT)
      
      call convert_time_to_string(pixel_cont(IIMAGE)%T, '%FT%H:%M:%SZ', iso8601_string)
      !write(*,'("time1: ", A)') iso8601_string
      
      write(NWLC,*) pixel_cont(IIMAGE)%nwl

      WRITE(iu_SDATA,'(/,i3,a25,f9.2,2i4,3x,a)')   &
	    NX*NY,trim(iso8601_string),    &
      pixel_cont(IIMAGE)%HOBS,pixel_cont(IIMAGE)%meas(1)%NSURF,pixel_cont(IIMAGE)%IFGAS, &
	   '    : NPIXELS TIMESTAMP  HEIGHT_OBS(m)  NSURF  IFGAS'    

!       write(*,*) 'NWLC : ',trim(adjustl(NWLC))
                    
    DO IY=1,NY
    DO IX=1,NX
      IIMAGE=index_clouds%INIMAGE(IX,IY,IT)

	    IF(pixel_cont(IIMAGE)%cloudy .EQ. 0) THEN
         write(tmp_message,'(2(a,i0),a)') 'ipixel = ',IIMAGE, &
         '  cloudy = ',pixel_cont(IIMAGE)%cloudy,'  - Pixel is cloudy'
         G_ERROR(trim(tmp_message))
	    ENDIF
      
      ind = 2 ! pixel => meas

      do IW=1,pixel_cont(IIMAGE)%nwl
      do IP=1,pixel_cont(IIMAGE)%meas(IW)%NIP 
      NBVM =  pixel_cont(IIMAGE)%meas(IW)%NBVM(IP)
      meas_type = pixel_cont(IIMAGE)%meas(IW)%meas_type(IP)      
      call set_pixel_meas (                      & ! IN
                           iw,                   &
                           NBVM,                 &
                           meas_type,            &
                           ind,                  &
                           FPSW(1:NBVM,IP,IW),   & ! INOUT
                           pixel_cont(IIMAGE)    &
                          )
                          
      do IV=1,pixel_cont(IIMAGE)%meas(IW)%NBVM(IP)
      if(meas_type .ge. meas_type_lid_beg .and. meas_type .le. meas_type_lid_end) then
      thetav(IV,IP,IW) = pixel_cont(IIMAGE)%HVP(IV)
      phi(IV,IP,IW)    = 0.0      
      else
      thetav(IV,IP,IW) = pixel_cont(IIMAGE)%meas(IW)%thetav(IV,IP)
      phi(IV,IP,IW)    = pixel_cont(IIMAGE)%meas(IW)%phi(IV,IP)
      endif ! meas_type .ge. meas_type_lid_beg .and.
      enddo ! IV


      enddo ! IP
      enddo ! IW


      !CFMT= '(i3,2i4,2i6,3f10.3,f8.2,i4,'//trim(adjustl(NWLC))// &
            !'i4,2('//trim(adjustl(NWLC))//'f10.3),'//trim(adjustl(NWLC))// &
	        !'i4,5000e14.6)'	  
     ! WRITE(iu_SDATA,TRIM(CFMT)) 

      WRITE(iu_SDATA,*)                      &
         IX,                                 &
         IY,                                 &
         pixel_cont(IIMAGE)%cloudy,          &
         pixel_cont(IIMAGE)%irow,            &
         pixel_cont(IIMAGE)%icol,            &
         pixel_cont(IIMAGE)%X,               &
         pixel_cont(IIMAGE)%Y,               & ! X - lon, Y - lat
         pixel_cont(IIMAGE)%MASL,            &
         pixel_cont(IIMAGE)%land_percent,    &
         
         pixel_cont(IIMAGE)%nwl,             &                                 
         (pixel_cont(IIMAGE)%meas(IW)%wl,    &
            IW=1,pixel_cont(IIMAGE)%nwl),    &
            
         (pixel_cont(IIMAGE)%meas(IW)%NIP,         &
            IW=1,pixel_cont(IIMAGE)%nwl),          &

         ((pixel_cont(IIMAGE)%meas(IW)%meas_type(IP),  &
            IP=1,pixel_cont(IIMAGE)%meas(IW)%NIP),     &                     
            IW=1,pixel_cont(IIMAGE)%nwl),              & 
            
			   ((pixel_cont(IIMAGE)%meas(IW)%NBVM(IP),   &
            IP=1,pixel_cont(IIMAGE)%meas(IW)%NIP), &                     
            IW=1,pixel_cont(IIMAGE)%nwl),          &                           
                                    
                                                                                    
			    (pixel_cont(IIMAGE)%meas(IW)%sza, 		   &
            IW=1,pixel_cont(IIMAGE)%nwl),		       &


         (((thetav(IV,IP,IW),                           &
            IV=1,pixel_cont(IIMAGE)%meas(IW)%NBVM(IP)), &
            IP=1,pixel_cont(IIMAGE)%meas(IW)%NIP),	    &
            IW=1,pixel_cont(IIMAGE)%nwl),	              &
            
         (((phi(IV,IP,IW),                              &
            IV=1,pixel_cont(IIMAGE)%meas(IW)%NBVM(IP)), &
            IP=1,pixel_cont(IIMAGE)%meas(IW)%NIP),	    &
            IW=1,pixel_cont(IIMAGE)%nwl),	              &

         (((FPSW(IV,IP,IW),                             &
            IV=1,pixel_cont(IIMAGE)%meas(IW)%NBVM(IP)), &
            IP=1,pixel_cont(IIMAGE)%meas(IW)%NIP),      &
            IW=1,pixel_cont(IIMAGE)%nwl),               &
            
         ((pixel_cont(IIMAGE)%meas(IW)%groundpar(IBRF), &
            IBRF=1,pixel_cont(IIMAGE)%meas(IW)%Nsurf),  &
            IW=1,pixel_cont(IIMAGE)%nwl),		            &
            
         ((pixel_cont(IIMAGE)%meas(IW)%gaspar,          & 
            k=1,pixel_cont(IIMAGE)%IFGAS),	            &
            IW=1,pixel_cont(IIMAGE)%nwl),               &
            
         ((pixel_cont(IIMAGE)%meas(IW)%IFCOV(IP),       & 
            IP=1,pixel_cont(IIMAGE)%meas(IW)%NIP),      &
            IW=1,pixel_cont(IIMAGE)%nwl),		            &

         (((pixel_cont(IIMAGE)%meas(IW)%CMTRX(IV,IP),   & 
            IV=1,pixel_cont(IIMAGE)%meas(IW)%IFCOV(IP)* &
                 pixel_cont(IIMAGE)%meas(IW)%NBVM(IP)), &
            IP=1,pixel_cont(IIMAGE)%meas(IW)%NIP),      &
            IW=1,pixel_cont(IIMAGE)%nwl),		            &
            
         ((pixel_cont(IIMAGE)%meas(IW)%IFMP(IP),        & 
            IP=1,pixel_cont(IIMAGE)%meas(IW)%NIP),      &
            IW=1,pixel_cont(IIMAGE)%nwl),		            &

         (((pixel_cont(IIMAGE)%meas(IW)%MPROF(IV,IP),   & 
            IV=1,pixel_cont(IIMAGE)%meas(IW)%IFMP(IP)*  &
                 pixel_cont(IIMAGE)%meas(IW)%NBVM(IP)), &
            IP=1,pixel_cont(IIMAGE)%meas(IW)%NIP),      &
            IW=1,pixel_cont(IIMAGE)%nwl)

! CHECK MEAS DATA

      do iw=1,pixel_cont(IIMAGE)%nwl
      do ip=1,pixel_cont(IIMAGE)%meas(IW)%NIP
      meas_type = pixel_cont(IIMAGE)%meas(IW)%meas_type(IP)
      do iv=1,pixel_cont(IIMAGE)%meas(IW)%NBVM(IP)
      select case(meas_type)
      case(meas_type_tod)
! tau         
         tau = FPSW(IV,IP,IW)
         IF(tau .LE. 0.0) THEN
         WRITE(tmp_message,'(a,a,5i4,a,a,es11.4,a)') 'PROBLEM WITH DATA !!!', &
         NEW_LINE('A'), &
         IT,IX,IY,IV,IW,' IT,IX,IY,IV,IW', &
         NEW_LINE('A'), &
         tau,' - TOD'
         G_ERROR(trim(tmp_message))
         ENDIF ! tau .LE. 0.0

      case(meas_type_aod)
! tau         
         tau = FPSW(IV,IP,IW)
         IF(tau .LE. 0.0) THEN
         WRITE(tmp_message,'(a,5i4,a,es11.4,a)') 'PROBLEM WITH DATA !!!', &
         NEW_LINE('A'), &
         IT,IX,IY,IV,IW,' IT,IX,IY,IV,IW', &
         NEW_LINE('A'), &
         tau,' - AOD'
         G_ERROR(trim(tmp_message))
         ENDIF ! tau .LE. 0.0

      case(meas_type_p11)
! p11         
         p11 = FPSW(IV,IP,IW)
         IF(p11 .LE. 0.0) THEN
         WRITE(tmp_message,'(a,a,5i4,a,a,es11.4,a)') 'PROBLEM WITH DATA !!!', &
         NEW_LINE('A'), &
         IT,IX,IY,IV,IW,' IT,IX,IY,IV,IW', &
         NEW_LINE('A'), &
         p11,' - p11'
         G_ERROR(trim(tmp_message))
         ENDIF ! p11 .LE. 0.0

      case(meas_type_p12)
! p12         
         p11 = FPSW(IV,IP-1,IW)
         p12 = FPSW(IV,IP,IW)
         IF(abs(p12/p11) .GE. 1.0) THEN
         WRITE(tmp_message,'(a,a,5i4,a,a,2e12.4,a)') 'PROBLEM WITH DATA !!!', &
         NEW_LINE('A'), &
         IT,IX,IY,IV,IW,' IT,IX,IY,IV,IW', &
         NEW_LINE('A'), &
         p12,p11,' - p12,p11'
         G_ERROR(trim(tmp_message))
         ENDIF ! p12 .LE. 0.0

      case(meas_type_p22)
! p22         
         p22 = FPSW(IV,IP,IW)
         IF(p22 .LE. 0.0) THEN
         WRITE(tmp_message,'(a,a,5i4,a,a,e11.4,a)') 'PROBLEM WITH DATA !!!', &
         NEW_LINE('A'), &
         IT,IX,IY,IV,IW,' IT,IX,IY,IV,IW', &
         NEW_LINE('A'), &
         p22,' - p22'
         G_ERROR(trim(tmp_message))
         ENDIF ! p22 .LE. 0.0

      case(meas_type_LS)
! LS
         LS = FPSW(IV,IP,IW)
         IF(LS .LE. 0.0) THEN
         WRITE(tmp_message,'(a,a,5i4,a,a,e11.4,a)') 'PROBLEM WITH DATA !!!', &
         NEW_LINE('A'), &
         IT,IX,IY,IV,IW,' IT,IX,IY,IV,IW', &
         NEW_LINE('A'), &
         LS,' - Lidar Signal'
         G_ERROR(trim(tmp_message))
         ENDIF ! LR .LE. 0.0

      case(meas_type_DP)
! DP
         DP = FPSW(IV,IP,IW)
         IF(DP .LE. 0.0) THEN
         WRITE(tmp_message,'(a,a,5i4,a,a,e11.4,a)') 'PROBLEM WITH DATA !!!', &
         NEW_LINE('A'), &
         IT,IX,IY,IV,IW,' IT,IX,IY,IV,IW', &
         NEW_LINE('A'), &
         LS,' - Depolarization'
         G_ERROR(trim(tmp_message))
         ENDIF ! DP .LE. 0.0

      case(meas_type_RL)
! RL
         RL = FPSW(IV,IP,IW)
         IF(RL .LE. 0.0) THEN
         WRITE(tmp_message,'(a,a,5i4,a,a,e11.4,a)') 'PROBLEM WITH DATA !!!', &
         NEW_LINE('A'), &
         IT,IX,IY,IV,IW,' IT,IX,IY,IV,IW', &
         NEW_LINE('A'), &
         RL,' - Roman Lidar signal'
         G_ERROR(trim(tmp_message))
         ENDIF ! RL .LE. 0.0

      case(meas_type_I)
! I
         I = FPSW(IV,IP,IW)
         IF(I .LE. 1e-5) THEN
         WRITE(tmp_message,'(a,a,5i4,a,a,e11.4,a)') 'PROBLEM WITH DATA !!!', &
         NEW_LINE('A'), &
         IT,IX,IY,IV,IW,' IT,IX,IY,IV,IW', &
         NEW_LINE('A'), &
         I,' - RADIANCE'
         G_ERROR(trim(tmp_message))
         ENDIF ! I .LE. 1e-5

      case(meas_type_Q)
! Q
         I = FPSW(IV,IP-1,IW)
         Q = FPSW(IV,IP,  IW)
         U = FPSW(IV,IP+1,IW)         
         IF(Q .eq. 0.0 .or. abs(Q/I) .GE. 1.0) THEN
         WRITE(tmp_message,'(a,a,5i4,a,a,3e13.4,a)') 'PROBLEM WITH DATA !!!', &
         NEW_LINE('A'), &
         IT,IX,IY,IV,IW,' IT,IX,IY,IV,IW', &
         NEW_LINE('A'), &
         I,Q,U,'  - I,Q,U'
         G_ERROR(trim(tmp_message))
         ENDIF ! Q .eq. 0.0 .or.

      case(meas_type_U)
! U
         I = FPSW(IV,IP-2,IW)
         Q = FPSW(IV,IP-1,IW)
         U = FPSW(IV,IP  ,IW)         
         IF(U .eq. 0.0 .or. abs(U/I) .GE. 1.0) THEN
         WRITE(tmp_message,'(a,a,5i4,a,a,3e13.4,a)') 'PROBLEM WITH DATA !!!', &
         NEW_LINE('A'), &
         IT,IX,IY,IV,IW,' IT,IX,IY,IV,IW', &
         NEW_LINE('A'), &
         I,Q,U,'  - I,Q,U'
         G_ERROR(trim(tmp_message))
         ENDIF ! U .eq. 0.0 .or.
         Q2 = Q*Q
         U2 = U*U
         P  = SQRT(Q2+U2)
         IF(abs(P/I) .GE. 1.0) THEN
         WRITE(tmp_message,'(a,a,5i4,a,a,es11.4,a,a,3es12.4,a)') &
         'PROBLEM WITH DATA !!!', &
         NEW_LINE('A'), &
         IT,IX,IY,IV,IW,' IT,IX,IY,IV,IW', &
         NEW_LINE('A'), &
         P/I,'  - Degree of Linear Polarization', &
         NEW_LINE('A'), &
         I,Q,U,'  - I,Q,U'
         G_ERROR(trim(tmp_message))
         ENDIF ! P/I .LE. -1.0

      case(meas_type_P)
! P
         I = FPSW(IV,IP-1,IW)
         P = FPSW(IV,IP,  IW)
         IF(abs(P/I) .GE. 1.0) THEN
         WRITE(tmp_message,'(a,a,5i4,a,a,es11.4,a,a,es11.4,a)') &
         'PROBLEM WITH DATA !!!', &
         NEW_LINE('A'), &
         IT,IX,IY,IV,IW,' IT,IX,IY,IV,IW', &
         NEW_LINE('A'), &
         P/I,'  - Degree of Linear Polarization', &
         NEW_LINE('A'), &
         P,I,'  - P,I'
         G_ERROR(trim(tmp_message))
         ENDIF ! P/I .LE. -1.0

      case default
         write(tmp_message,'(3(a,i0),a)') 'iw = ',iw,' ip = ',ip, &
         ' meas_type = ',meas_type,' - unknown value of meas_type'
         G_ERROR(trim(tmp_message))
      end select

      enddo ! iv      
      enddo ! ip
      enddo ! iw

      ENDDO ! IY
      ENDDO ! IX
      ENDDO ! IT	  	

      CLOSE (iu_SDATA)

   return
   end subroutine write_simul_pixels 

!	sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

	function write_sdata_pixels (							      &
                                 sdata_sim_file,	&
                                 npixels,         &
                                 segment					&
                               )
		use m_inssor
!	**************************************************************************************
		character(*),            intent(in)		  ::	sdata_sim_file
		integer,                 intent(in)		  ::	npixels
		type(segment_data),      intent(inout)	::	segment
		logical		         ::	write_sdata_pixels
		integer					   ::	ii,jj,kk	
		integer					   ::	io_status
		integer,dimension(npixels)	                ::  ix_pool,iy_pool,it_pool 
		integer,dimension(npixels,npixels,npixels)  :: 	index_3d ! *** (IT,IX,IY)

		character(LEN=25)		::	iso8601_string
    integer(KIND_TIME)   :: time_sec
		integer					::	Nsurf,IFGAS
 		integer					::  NX,NY,NT
		logical         ::  status_funct 
!		real,dimension(max_num_wave_length_v2,max_num_valid_directions_v2)		::	FPSW
!	------------------------------------------------------------------------------------------------------
      integer  ::  id_sdata_like_file
		integer,dimension(:),allocatable			::	num_el_time_gr		! *** number of elements in each timegroup
!	------------------------------------------------------------------------------------------------------
		integer			::	it_group_begin,it_group_end
		integer			::	ix_group_begin,ix_group_end
		integer			::	iy_group_begin,iy_group_end

		integer			::	it,it_shift,icount_time_gr
!	------------------------------------------------------------------------------------------------------
		character(LEN=255),parameter 	::	timegroup_header =	&
						"   : NPIXELS  TIMESTAMP  HEIGHT_OBS(m)  NSURF  IFGAS"
		character(LEN=64) 	::	format_timegroup_hdr =	&
						'(/,i3,3x,a20,f11.2,2i4,3x,a,i5)'

!			trim(format_timegroup_hdr)
!	------------------------------------------------------------------------------------------------------
		integer 		:: dbg_msg_level
!	------------------------------------------------------------------------------------------------------
		dbg_msg_level 	= 12
!	------------------------------------------------------------------------------------------------------
		ix_pool	= 0
		iy_pool	= 0
		it_pool	= 0
!	------------------------------------------------------------------------------------------------------
		Nsurf = segment%pixels(1)%meas(1)%Nsurf
		IFGAS  = segment%pixels(1)%IFGAS

    NT = segment%NT
    if(.NOT. allocated(num_el_time_gr)) then
			allocate(num_el_time_gr(1:NT))
		endif
    NX = segment%NX
    NY = segment%NY
!	------------------------------------------------------------------------------------------------------
      open (newunit=id_sdata_like_file, FILE=trim(sdata_sim_file), STATUS='replace', & 
                                                      FORM='formatted',RECL=2000000)      
!	------------------------------------------------------------------------------------------------------
!tl If time group consists of 2 pixels (ix=2 iy=1; ix=1 iy=2) sorting indexes does not
!tl work properly. The problem is not fixed (goto 333 - 333 continue).

		ii = 0
    num_el_time_gr(:) = 0
    icount_time_gr    = 1  
      
		do  ii = 1,segment%npixels
			it_pool(ii) = segment%pixels(ii)%it
			ix_pool(ii) = segment%pixels(ii)%ix
			iy_pool(ii) = segment%pixels(ii)%iy
			if (segment%pixels(ii)%it .EQ. 0 ) then
				write(*,*) "Error. it(",ii,") = ",segment%pixels(ii)%it
			endif
			if (segment%pixels(ii)%ix .EQ. 0 ) then
				write(*,*) "Error. ix(",ii,") = ",segment%pixels(ii)%ix
			endif
			if (segment%pixels(ii)%iy .EQ. 0 ) then
				write(*,*) "Error. iy(",ii,") = ",segment%pixels(ii)%iy
			endif
			! *** forming it,ix,iy matrix with numbers of corresponding
			! 	linear pixels set. Then we will use it,ix,iy from corresponding 
			!	i?_pool arrays to address pixels in any desired it,ix,iy order	
			index_3d (                        &
                segment%pixels(ii)%it,  &
                segment%pixels(ii)%ix,  &
                segment%pixels(ii)%iy   &
               ) = ii 
!write(*,*) ii,segment%pixels(ii)%it,segment%pixels(ii)%ix,segment%pixels(ii)%iy,index_3d &
!											  (					     &
!												 segment%pixels(ii)%it,	 &
!												 segment%pixels(ii)%ix,	 &
!												 segment%pixels(ii)%iy	 &
!											  ),					   &
!'  ii,it,ix,iy,index_3d'
!write(*,*) ii,it_pool(ii),ix_pool(ii),iy_pool(ii),'  ii,it_pool(ii),ix_pool(ii),iy_pool(ii)'
!write(*,*) ii,segment%pixels(ii)%it,segment%pixels(ii)%ix,segment%pixels(ii)%iy,index_3d &
!											  (					     &
!												 it_pool(ii),	 &
!												 ix_pool(ii),	 &
!												 iy_pool(ii)	 &
!											  ),					   &
!'  ii,it(ii),ix(ii),iy(ii),index_3d'

      if(icount_time_gr .ne. segment%pixels(ii)%it) icount_time_gr = icount_time_gr + 1
      num_el_time_gr(icount_time_gr) = num_el_time_gr(icount_time_gr) + 1
!write(*,*) ii,icount_time_gr,num_el_time_gr(icount_time_gr),'  ii,icount_time_gr,num_el_time_gr(icount_time_gr)'
		enddo  !  ii
!write(*,*) 
!stop		

goto 333
		it_group_begin = 1
		it = 0
		do  ii = 1,segment%npixels - 1
			if(it_pool(ii+1) .gt. it_pool(ii) .or. ii+1 .eq. segment%npixels) then ! *** looking for new timegroup or for the last element of array
				if(ii+1 .eq. segment%npixels .and. it_pool(ii+1) .eq. it_pool(ii)) then
					it_group_end = ii + 1
				else
					it_group_end = ii
				endif
				if(allocated(num_el_time_gr)) then
					it = it + 1
					num_el_time_gr(it) = it_group_end - it_group_begin +1 
					if (dbg_msg_level .ge. 12) then
						write (*,*) "num_el_time_gr(",it,") =", num_el_time_gr(it) 
						write (*,*) "it_group_end = ",it_group_end
						write (*,*) "it_group_begin = ",it_group_begin
					endif
				endif
				call inssor(iy_pool(it_group_begin:it_group_end)) ! *** sorting part of iy indexes' array
				iy_group_begin = it_group_begin
				do jj = it_group_begin,it_group_end -1
					if(iy_pool(jj+1) .gt. iy_pool(jj) .or. jj+1 .eq. it_group_end) then ! *** looking for new timegroup or for the last element of array
						if(jj+1 .eq. it_group_end ) then
							iy_group_end = jj + 1
						else
							iy_group_end = jj
						endif
						call inssor(iy_pool(iy_group_begin:iy_group_end)) ! *** sorting part of iy indexes' array
						ix_group_begin = iy_group_begin
						do kk = iy_group_begin,iy_group_end -1
							if(iy_pool(kk+1) .gt. iy_pool(kk) .or. kk+1 .eq. iy_group_end) then ! *** looking for new timegroup or for the last element of array
								ix_group_end = kk
								if(kk+1 .eq. iy_group_end ) then
									ix_group_end = kk + 1
								else
									ix_group_end = kk
								endif
								call inssor(ix_pool(ix_group_begin:ix_group_end)) ! *** sorting part of ix indexes' array
								ix_group_begin = ix_group_end + 1
							endif
							iy_group_begin = iy_group_end + 1
						enddo ! end of do  kk = iy_group_begin,iy_group_end -1
						it_group_begin = it_group_end + 1
					endif ! end of if	(	iy_pool(jj+1) .GT.iy_pool(jj) .OR.jj+1 .EQ. it_group_end ) then ! *** looking for new timegroup or for the last element of array
				enddo ! end of do  jj = it_group_begin,it_group_end -1
!				if(it_pool(ii+1) .gt. it_pool(ii) .and.	ii+1 .eq. segment%npixels) then ! *** in this case we have ONLY ONE PIXEL in timgroup i.e. no need to sort anything over ix and iy. We need just to output as it is.
				if(it_pool(ii+1) .gt. it_pool(ii)) then ! *** in this case we have ONLY ONE PIXEL in timgroup i.e. no need to sort anything over ix and iy. We need just to output as it is.
					it_group_begin = it_group_end + 1
					if(ii+1 .eq. segment%npixels) then
						if(allocated(num_el_time_gr)) then
							it = it + 1
							num_el_time_gr(it) = 1 
							if(dbg_msg_level .ge. 12) then
								write (*,*) "num_el_time_gr(",it,") =", num_el_time_gr(it) 
								write (*,*) "Near the bottom"
							endif
						endif
					endif
				endif  !  if	(	it_pool(ii+1) .GT.it_pool(ii) .AND. ii+1 .EQ. segment%npixels) then ! *** looking for new timegroup or for the last element of array
			endif ! end of if	(it_pool(ii+1) .GT.it_pool(ii).OR.ii+1 .EQ. segment%npixels ) then ! *** looking for new timegroup or for the last element of array
		enddo ! end of do  ii = 1,segment%npixels - 1

		if  (segment%npixels .EQ. 1) then
			num_el_time_gr(1) = 1 ! in case there is only one pixel in the set
		endif
333 continue

!do ii=1,segment%npixels
!write(*,*) ii,it_pool(ii),ix_pool(ii),iy_pool(ii),'  ii,it_pool(ii),ix_pool(ii),iy_pool(ii)'
!enddo 

#ifdef WARN_DRY
#warning "__SDATA_VERSION__ duplicated"
#endif           

		write (id_sdata_like_file,'(a)',IOSTAT= io_status) 'SDATA version 2.0'
		write (id_sdata_like_file,'(i3,2i4,a)',IOSTAT= io_status) NX,NY,NT,'  : NX NY NT'
		it = 1
!write(*,*) it,it_pool(it),ix_pool(it),iy_pool(it),index_3d &
!											  (					     &
!												 it_pool(it),	 &
!												 ix_pool(it),	 &
!												 iy_pool(it)	 &
!											  ),					   &
!'  1: ii,it_pool(ii),ix_pool(ii),iy_pool(ii),index_3d'

			time_sec = segment%pixels            &	
                         (                 &
                           index_3d        &
                            (              &
                              it_pool(1),  &
                              ix_pool(1),  &
                              iy_pool(1)   &
                            )              &
                         )%t				


!		write(id_sdata_like_file,* ,IOSTAT= io_status) 	&
		!character(LEN=64) 	::	format_timegroup_hdr =	&
		!				'(/,i3,i8,3x,a,3x,a,f9.2,3i4,a)'
    call convert_time_to_string(time_sec, '%FT%H:%M:%SZ', iso8601_string)
    !write(*,'("time2: ", A)') iso8601_string
		write(id_sdata_like_file,trim(format_timegroup_hdr) ,IOSTAT= io_status) &
            num_el_time_gr(it),                 &
            trim(adjustl(iso8601_string)),		  &
            segment%pixels(1)%HOBS,             & 
            segment%pixels(1)%meas(1)%Nsurf,    &
            segment%pixels(1)%IFGAS,            & 
            trim(adjustl(timegroup_header)),    &
            1
               				
!				segment%pixels(1)%time_string,segment%pixels(1)%data_filename
		it_shift = 0
		do  ii = 1,segment%npixels

			if(num_el_time_gr(it) + it_shift + 1 .EQ. ii ) then

!	------------------------------------------------------------------------------------------------------
! *** unpacking date/time and checking them if they are present
!write(*,*) ii,it_pool(ii),ix_pool(ii),iy_pool(ii),index_3d &
!											  (					     &
!												 it_pool(ii),	 &
!												 ix_pool(ii),	 &
!												 iy_pool(ii)	 &
!											  ),					   &
!'  2: ii,it_pool(ii),ix_pool(ii),iy_pool(ii),index_3d'
            time_sec = segment%pixels	 &	
								 (							       &
									 index_3d 				   &
											  (					     &
												 it_pool(ii),	 &
												 ix_pool(ii),	 &
												 iy_pool(ii)	 &
											  )					     &
								 )%t		
								                   
! *** END unpacking date/time and checking them if they are present
!	------------------------------------------------------------------------------------------------------
				it_shift = it_shift + num_el_time_gr(it)
				it = it + 1
				write(id_sdata_like_file,*,IOSTAT= io_status) 
!				write(id_sdata_like_file,*,IOSTAT= io_status) 		&
            
            call convert_time_to_string(time_sec, '%FT%H:%M:%SZ', iso8601_string)
            !write(*,'("time3: ",A)') iso8601_string
				write(id_sdata_like_file,trim(format_timegroup_hdr) ,IOSTAT= io_status) &
				num_el_time_gr(it),                    &
            iso8601_string,                    &
            segment%pixels(ii)%HOBS,           & 
            segment%pixels(ii)%meas(1)%Nsurf,  &
            segment%pixels(ii)%IFGAS,          & 
				trim(adjustl(timegroup_header)),       &
				it_pool(ii)
!						segment%pixels(ii)%time_string,segment%pixels(ii)%data_filename
			endif ! end of if(num_el_time_gr(it) + it_shift + 1 .EQ. ii )
		
			status_funct =                &
				write_one_pixel_sdata       &
						(                       &
						id_sdata_like_file,     &
						segment%pixels          &
							(                     &	
                index_3d            &
                  (                 &
                    it_pool(ii),    &
                    ix_pool(ii),    &
                    iy_pool(ii)     &
                  )                 &
							)                     &
						)

		enddo ! ii

!		READ (id_sdata_like_file,*,IOSTAT= io_status) NX,NY,NT,HOBS,IPLANE,Nsurf,IFGAS,IP_ANGLE
      close(id_sdata_like_file)			

!	------------------------------------------------------------------------------------------------------
		if(allocated(num_el_time_gr)) then
			deallocate(num_el_time_gr)
		endif
!	------------------------------------------------------------------------------------------------------
		write_sdata_pixels 	= .true. 

	end function write_sdata_pixels

!ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

	function write_one_pixel_sdata (							        &
                                    id_sim_sdata_file,	&
                                    one_pixel				    &
                                 )
!	**************************************************************************************
!	------------------------------------------------------------------------------------------------------
		integer,intent(in)		    ::  id_sim_sdata_file
		type(pixel),intent(inout) ::  one_pixel
!	------------------------------------------------------------------------------------------------------
		logical		                ::  write_one_pixel_sdata
!	------------------------------------------------------------------------------------------------------
		integer					          ::  IW,IV,IP,IBRF,k,ind, &
                                  nvalid_meas,meas_type
!	------------------------------------------------------------------------------------------------------
      real,dimension(KNBVM,KIP,KWM)  ::  thetav
		real,dimension(KNBVM,KIP,KWM)    ::  phi
      real,dimension(KNBVM,KIP,KWM)  ::  FPSW
		character(LEN=255)        	     ::  CFMT
		character(LEN=16)                ::  nwl_string
!	------------------------------------------------------------------------------------------------------
!	------------------------------------------------------------------------------------------------------
                    
      ind = 2 ! pixel => meas
      do IW=1,one_pixel%nwl
      do IP=1,one_pixel%meas(IW)%NIP 
      meas_type   = one_pixel%meas(IW)%meas_type(IP)
      nvalid_meas = one_pixel%meas(IW)%NBVM(IP)       
      call set_pixel_meas (                                  & ! IN
                           iw,                               &
                           nvalid_meas,                      &
                           one_pixel%meas(IW)%meas_type(IP), &
                           ind,                              &
                           FPSW(1:nvalid_meas,IP,IW),        & ! INOUT
                           one_pixel                         &
                          )
      if(meas_type .gt. meas_type_lid_beg .and. meas_type .lt. meas_type_lid_end) then
         thetav(1:nvalid_meas,IP,IW) = one_pixel%HVP(1:nvalid_meas)
         phi   (1:nvalid_meas,IP,IW) = 0.0
      else
         thetav(1:nvalid_meas,IP,IW) = one_pixel%meas(IW)%thetav(1:nvalid_meas,IP)
         phi   (1:nvalid_meas,IP,IW) = one_pixel%meas(IW)%phi   (1:nvalid_meas,IP)         
      endif ! meas_type

      enddo ! IP
      enddo ! IW

		write(nwl_string,*) one_pixel%nwl
		!CFMT= '(i3,2i4,2i6,3f10.3,f8.2,i4,'//trim(adjustl(nwl_string))// &
				!'i4,2('//trim(adjustl(nwl_string))//'f10.3),'//trim(adjustl(nwl_string))// &
				!'i4,5000e14.6,2x)'
      !write(id_sim_sdata_file,TRIM(CFMT)) 	&

      write(id_sim_sdata_file,*) &
			one_pixel%ix,			   &
			one_pixel%iy,			   &
			one_pixel%cloudy,		   &
			one_pixel%irow,		   &
			one_pixel%icol,		   &
			one_pixel%x,		      &
			one_pixel%y,		      &
			one_pixel%MASL,       &
			one_pixel%land_percent,	&
         
			one_pixel%nwl,		          &
			(one_pixel%meas(IW)%wl, 	 &
            IW=1,one_pixel%nwl),		 &
            
         (one_pixel%meas(IW)%NIP,         &
            IW=1,one_pixel%nwl),          & 
                                   
! integer
         ((one_pixel%meas(IW)%meas_type(IP),  &
            IP=1,one_pixel%meas(IW)%NIP),     &                     
            IW=1,one_pixel%nwl),              & 

			((one_pixel%meas(IW)%NBVM(IP),   &
            IP=1,one_pixel%meas(IW)%NIP), &                     
            IW=1,one_pixel%nwl),          &                           
! real
			(one_pixel%meas(IW)%sza, 		   &
            IW=1,one_pixel%nwl),		      &
            
         (((thetav(IV,IP,IW),                  &
            IV=1,one_pixel%meas(IW)%NBVM(IP)), &
            IP=1,one_pixel%meas(IW)%NIP),	     &                                       
            IW=1,one_pixel%nwl),	              &
            
         (((phi(IV,IP,IW),                     &
            IV=1,one_pixel%meas(IW)%NBVM(IP)), &
            IP=1,one_pixel%meas(IW)%NIP),	     &                                       
            IW=1,one_pixel%nwl),	              &
                  
         (((FPSW(IV,IP,IW),		              &
            IV=1,one_pixel%meas(IW)%NBVM(IP)), &
            IP=1,one_pixel%meas(IW)%NIP),	     &                                       
            IW=1,one_pixel%nwl),	              &
                  
         ((one_pixel%meas(IW)%groundpar(IBRF), &
            IBRF=1,one_pixel%meas(IW)%Nsurf),  &
            IW=1,one_pixel%nwl),		           &
            
         ((one_pixel%meas(IW)%gaspar,          & 
            k=1,one_pixel%IFGAS),	              &
            IW=1,one_pixel%nwl),               &
            
         ((one_pixel%meas(IW)%IFCOV(IP),       & 
            IP=1,one_pixel%meas(IW)%NIP),      &
            IW=1,one_pixel%nwl),		           &
                        
         (((one_pixel%meas(IW)%CMTRX(IV,IP),    & 
            IV=1,one_pixel%meas(IW)%IFCOV(IP)*  &
                 one_pixel%meas(IW)%NBVM(IP)),  &
            IP=1,one_pixel%meas(IW)%NIP),       &
            IW=1,one_pixel%nwl),		            &
            
         ((one_pixel%meas(IW)%IFMP(IP),         & 
            IP=1,one_pixel%meas(IW)%NIP),       &
            IW=1,one_pixel%nwl),		            &
                        
         (((one_pixel%meas(IW)%MPROF(IV,IP),    & 
            IV=1,one_pixel%meas(IW)%IFMP(IP)*   &
                 one_pixel%meas(IW)%NBVM(IP)),  &
            IP=1,one_pixel%meas(IW)%NIP),       &
            IW=1,one_pixel%nwl)                
              
!	------------------------------------------------------------------------------------------------------
      write_one_pixel_sdata = .true. 
!	------------------------------------------------------------------------------------------------------
      end function write_one_pixel_sdata

!	sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine get_vert_prof_h (  iw,ip,    	 &
                                    pixel_cont,  & ! IN
                                    NBVM,HVP		 & ! OUT
                                 )

      implicit none
!	------------------------------------------------------------------------------------------------------
      integer,           intent(in)  :: iw,ip
      type(pixel),       intent(in)  :: pixel_cont
!  ------------------------------------------------------------------------------------------------------
      integer,           intent(out) :: NBVM
      real,dimension(:), intent(out) :: HVP   
!	------------------------------------------------------------------------------------------------------
      integer  ::  iv
!	------------------------------------------------------------------------------------------------------
      HVP(:) = 0.0
	  
      NBVM = pixel_cont%meas(iw)%NBVM(ip)	  
      if(size(HVP) .LT. NBVM) then
         write(*,*) 'Size(HVP)=',Size(HVP),' .LT. NBVM=',NBVM
         stop 'stop in get_vert_prof_h'
      endif
	  
      do iv=1,NBVM
      HVP(iv) = pixel_cont%HVP(iv)
      enddo ! iv

	end subroutine get_vert_prof_h

! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine get_HVP_lidar ( segment,NHVP,HVP )

      implicit none
!	------------------------------------------------------------------------------------------------------
      type(segment_data), intent(in)      ::  segment
      integer,            intent(inout)   ::  NHVP
      real,dimension(:),  intent(inout)   ::  HVP       
!	------------------------------------------------------------------------------------------------------
      integer  ::  meas_type,ipix,iw,ip
!	------------------------------------------------------------------------------------------------------
      NHVP = 1
      HVP(:) = 0.0
                    
      do ipix=1,segment%npixels
      do iw=1,segment%pixels(ipix)%nwl
      do ip=1,segment%pixels(ipix)%meas(iw)%NIP
         meas_type = segment%pixels(ipix)%meas(iw)%meas_type(ip)
         if(meas_type .ge. meas_type_lid_beg .and. meas_type .le. meas_type_lid_end) then
            NHVP = segment%pixels(ipix)%meas(iw)%NBVM(ip) 
            HVP(1:NHVP) = segment%pixels(ipix)%HVP(1:NHVP)
            !write(*,*) 'In get_HVP: IW,IP, NBVM(IP)',IW,IP,NHVP
            return 
         endif ! meas_type .ge. meas_type_lid_beg .and.
      enddo ! ip
      enddo ! iw
      enddo ! ipix
      
      return
      end subroutine get_HVP_lidar

! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

! Assign noise number to every measurement types in each pixel

      subroutine assign_noise_index(RIN,segment_meas,MNOISEI)
      
      use mod_par_inv,    only : KW,KIMAGE
      use mod_retr_settings_derived_type
      
      implicit none
!	-------------------------------------------------------------------------
      type(retr_input_settings),        intent(in)     ::  RIN
      type(segment_data),               intent(in)     ::  segment_meas
      integer,dimension(KIP,KWM,KIMAGE),intent(inout)  ::  MNOISEI
!	-------------------------------------------------------------------------            
      integer                :: nwl
      real,   dimension(KW)  :: wl
      integer,dimension(KW)  :: ind_wl
      integer                :: ipix,ip,iw,ip1,iw1,iin
!	-------------------------------------------------------------------------
!	-------------------------------------------------------------------------
      do ipix=1,segment_meas%npixels
      call get_pixel_wl ( segment_meas%pixels(ipix),  &
                          nwl,wl,ind_wl               &
                        )
if(ipix .eq. -5) then
write(*,*) 'ind_wl: ',ind_wl(1:nwl)
stop
endif
      do iw=1,segment_meas%pixels(ipix)%nwl
      do ip=1,segment_meas%pixels(ipix)%meas(iw)%NIP     
         do iin=1,RIN%NOISE%INOISE
         do ip1=1,RIN%NOISE%NMT(iin)
         do iw1=1,RIN%NOISE%NWLP(ip1,iin)               
if(iin .eq. -2) then
write(*,'(2a,6i5,a)') 'ip1, iw1, ind_wl(iw), IWLP(iw1,ip1,iin),', &
' meas_type(ip), MT(ip1,iin) - ',ip1,iw1,ind_wl(iw),RIN%NOISE%IWLP(iw1,ip1,iin), &
segment_meas%pixels(ipix)%meas(iw)%meas_type(ip),RIN%NOISE%MT(ip1,iin), &
'  -  in assign_noise_index'
endif
            if ( ind_wl(iw) .eq. RIN%NOISE%IWLP(iw1,ip1,iin)  .and.  &
                 segment_meas%pixels(ipix)%meas(iw)%meas_type(ip) .eq. RIN%NOISE%MT(ip1,iin) &
               ) then
               MNOISEI(ip,iw,ipix)=iin

!               write(*,*) 'ipix=',ipix,'  ind_wl(iw)=',ind_wl(iw),'  ip=',ip, &
!               '  iin=',iin,'  ip1=',ip1,'  iw1=',iw1,  &
!               '  MNOISE=',MNOISEI(ip,iw,ipix)
            endif
         enddo ! iw1
         enddo ! ip1
         enddo ! iin
      enddo ! ip
      enddo ! iw
      if(RIN%IPRI_verbose) then
         if(ipix .eq. 1) write(*,'(a)') 'in assign_noise_index:'
         write(*,'(4x,a,i0,3x,a,100i5)') 'pixel # ',ipix,'MNOISE:  ',((MNOISEI(ip,iw,ipix),  &
         ip=1,segment_meas%pixels(ipix)%meas(iw)%NIP),iw=1,segment_meas%pixels(ipix)%nwl)
      endif !  RIN%IPRI_verbose
      enddo ! ipix

      return
      end subroutine assign_noise_index

! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss 
! Fill out an array of wave lengths for retrieved parameters (important if pixeles have different neither numbers 
!                                                        or values of wave lengths)

      subroutine set_RIN_wave_length_array (  iu_main_output,  & ! IN
                                              segment_meas,    &
                                              RIN              & ! INOUT
                                           )
      
      use mod_retr_settings_derived_type
      use m_inssor
      use mod_par_inv, only: KW

      implicit none
!	-------------------------------------------------------------------------
      integer,                      intent(in)     :: iu_main_output
      type(segment_data),           intent(in)     :: segment_meas
      type(retr_input_settings),    intent(inout)  :: RIN
!	-------------------------------------------------------------------------
      integer :: ipix,iw,iw1,nwl_temp,nwl,IDIM1,IDIM2
      real    :: tiny
      logical :: add_wl
!	-------------------------------------------------------------------------
      tiny  = 1e-5
      
! Calculate RIN%NW - number of wavelengths from wavelength dependent characteristics ( refractive index or surface )
      do IDIM1=1,RIN%NDIM%n1 
        if(RIN%NDIM%par_type(IDIM1) .eq. par_type_RERI_spect) then  
          RIN%NW = 0
          do IDIM2=1,RIN%NDIM%n2(IDIM1)
            RIN%NW = max( RIN%NW,RIN%NDIM%n3(IDIM2,IDIM1) )
          enddo
        endif
        if(RIN%NDIM%par_type(IDIM1) .gt. par_type_SURF1_land_beg .and. & 
           RIN%NDIM%par_type(IDIM1) .lt. par_type_SURF_water_end) then
          RIN%NW = 0
          do IDIM2=1,RIN%NDIM%n2(IDIM1)
            RIN%NW = max( RIN%NW,RIN%NDIM%n3(IDIM2,IDIM1) )
          enddo       
        endif ! RIN%NDIM%par_type(IDIM1) .eq. par_type_RERI_spect
      enddo ! IDIM1

! A set of wavelengths from measurement data
      nwl = 0
      RIN%WAVE(:) = 0.0      
      nwl = segment_meas%pixels(1)%nwl
      if(nwl .gt. KW) then
        write(tmp_message,'(2(a,i0),a)') &
        'Number of wavelengths in segment data nwl = ',nwl, &
        '  is bigger than constant KW = ',KW,'  in module mod_par_inv'
        G_ERROR(trim(tmp_message))
      endif
      RIN%WAVE(1:nwl) = segment_meas%pixels(1)%meas(1:nwl)%wl
      nwl_temp = nwl
      
loop_pix :  do ipix=2,segment_meas%npixels

      do iw=1,segment_meas%pixels(ipix)%nwl
      add_wl = .true.
         do iw1=1,nwl_temp
            if (abs(segment_meas%pixels(ipix)%meas(iw)%wl-RIN%WAVE(iw1))/RIN%WAVE(iw1) .le. tiny) then 	
             add_wl = .false.	
             exit
            endif ! segment_meas%pixels(ipix)%meas(iw)%wl .eq. wl(iw1)
         enddo ! iw1
         if(add_wl) then
         nwl = nwl+1
            if(nwl .gt. KW) then
              write(tmp_message,'(2(a,i0),a)') &
              'Number of wavelengths in segment data nwl = ',nwl, &
              '  is bigger than constant KW = ',KW,'  in module mod_par_inv'
              G_ERROR(trim(tmp_message))
            endif
            if(nwl .gt. RIN%NW) then
              write(tmp_message,'(2(a,i0),a)') &
              'Number of wavelengths in segment data nwl = ',nwl, &
              '  is bigger than number of wavelengths NW = ',RIN%NW,'  in settings'
              G_ERROR(trim(tmp_message))
            endif ! nwl .gt. RIN%NW
         RIN%WAVE(nwl) = segment_meas%pixels(ipix)%meas(iw)%wl                  
         nwl_temp = nwl         
         endif ! add_wl
      enddo ! iw      
      
      enddo loop_pix

      if(nwl .ne. RIN%NW) then
          write(tmp_message,'(2(a,i0),a)') &
          'Number of wavelengths in segment data nwl = ',nwl, &
          ' must be equal to number of wavelengths NW = ',RIN%NW,'  in settings'
          G_ERROR(trim(tmp_message))
      endif ! nwl .ne. RIN%NW

! Sort wavelengths into increasing order (Insertion sort)

      call inssor(RIN%WAVE(1:nwl))
      
      !write(*,*) 'wave lengths:'
      !write(*,*) (wl(iw),iw=1,nwl)
      !write(*,*) 

      !write(*,*) 'sub set_RIN_wave_length_array: wave lengths from SDATA   : ',RIN%WAVE(1:RIN%NW)
      if(RIN%IPRI_additional_info) then
         write(iu_main_output,*) 'NW,(RIN%WAVE(IW),IW=1,NW) in set_RIN_wave_length_array'
         write(iu_main_output,'(i4,20f8.4)') RIN%NW,(RIN%WAVE(IW),IW=1,RIN%NW) 
      endif

      return
      end subroutine set_RIN_wave_length_array

! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss 

    subroutine set_index_clouds(RIN, segment, index_clouds)

        use mod_retr_settings_derived_type

        implicit none
!	-------------------------------------------------------------------------        
		  type(retr_input_settings),  intent(in)  :: RIN
        type(segment_data),       intent(in)  :: segment
        type(ind_clouds),         intent(out) :: index_clouds
!	-------------------------------------------------------------------------        
        integer :: i,ipix,IX,IY,IT
!	-------------------------------------------------------------------------
        
        ! Initialize
         index_clouds%NX=segment%NX
         index_clouds%NY=segment%NY
         index_clouds%NT=segment%NT
         index_clouds%ITIMAGE(:)   = 0
         index_clouds%IXIMAGE(:)   = 0
         index_clouds%IYIMAGE(:)   = 0
         index_clouds%INIMAGE(:,:,:) = 0
         index_clouds%T(:,:,:)       = 0.0
         index_clouds%X(:,:,:)       = 0.0
         index_clouds%Y(:,:,:)       = 0.0
         index_clouds%ICLOUD(:,:,:)  = 0

        
        ! Set values
        DO ipix=1, segment%npixels
            IX=segment%pixels(ipix)%IX
            IY=segment%pixels(ipix)%IY
            IT=segment%pixels(ipix)%IT
            index_clouds%INIMAGE(IX,IY,IT)=ipix
            !index_clouds%ICLOUD(IX,IY,IT)=segment%pixels(ipix)%cloudy
            index_clouds%ICLOUD(IX,IY,IT)=1
            index_clouds%ITIMAGE(ipix)=IT
            index_clouds%IYIMAGE(ipix)=IY
            index_clouds%IXIMAGE(ipix)=IX
            index_clouds%X(IX,IY,IT)=segment%pixels(ipix)%X
            index_clouds%Y(IX,IY,IT)=segment%pixels(ipix)%Y
            index_clouds%T(IX,IY,IT)=segment%pixels(ipix)%T
        ENDDO        

      if(RIN%IPRI_verbose) then
         write(*,'(a)') 'in set_index_clouds:'
         write(*,'(4x,3(a,i4))') 'NT=',index_clouds%NT,'  NX=',index_clouds%NX,'  NY=',index_clouds%NY
         i=0
         do IT=1,index_clouds%NT
         do IY=1,index_clouds%NY
         do IX=1,index_clouds%NX
          if(index_clouds%INIMAGE(IX,IY,IT) .gt. 0) then 	  
            i=i+1
            !write(*,'(2(a,i5),a,i16,2(a,f10.3),4(a,i4))') 'i=',i,'  IIMAGE=',  &
            write(*,'(4x,a,i0,3x,a,i0,2(3x,a,f9.3),3x,a,f13.1,7(3x,a,i0))') 'i = ',i,'pixel # ',  &
            index_clouds%INIMAGE(IX,IY,IT),        &
            'X = ',index_clouds%X(IX,IY,IT),       &
            'Y = ',index_clouds%Y(IX,IY,IT),       &
            'T = ',real(index_clouds%T(IX,IY,IT)), &
            'ix = ',index_clouds%IXIMAGE(index_clouds%INIMAGE(IX,IY,IT)),    &
            'iy = ',index_clouds%IYIMAGE(index_clouds%INIMAGE(IX,IY,IT)),    &
            'it = ',index_clouds%ITIMAGE(index_clouds%INIMAGE(IX,IY,IT)),    &
            'icloud = ',index_clouds%ICLOUD(IX,IY,IT),                       &
            'out_x = ',segment%pixels(index_clouds%INIMAGE(IX,IY,IT))%out_x, &
            'out_y = ',segment%pixels(index_clouds%INIMAGE(IX,IY,IT))%out_y, &
            'out_t = ',segment%pixels(index_clouds%INIMAGE(IX,IY,IT))%out_t
          endif
         enddo
         enddo
         enddo
      endif ! IPRI_mode .eq. 1

    return
    end subroutine set_index_clouds

!	sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

end module mod_sdata

!	sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! Assign general wave length indices to pixel wave lengths.
! For real Parasol data radiance measurements at wavelength #iw can be corrected.
      subroutine set_segment_meas ( iu_main_output,    & ! IN
                                    RIN,               & 
                                    segment_meas       & ! INOUT
                                  )
      use mod_sdata
      use mod_retr_settings_derived_type

      implicit none
!	---------------------------------------------------------------------------------------
      integer,                   intent(in)     :: iu_main_output
      type(retr_input_settings), intent(in)     :: RIN
      type(segment_data),        intent(inout)	:: segment_meas
! -----------------------------------------------------------------------------------------
! Assign general wave length indices to pixel wave lengths.
      call set_segment_pixel_wl_index (  iu_main_output, & ! IN
                                         RIN,               & 
                                         segment_meas       & ! INOUT
                                      )
!! For real Parasol data radiance measurements at wavelength #iw can be corrected.
!      if(RIN%I_corr)  &
!      call radiance_correction (  iu_main_output,  & ! IN
!                                  RIN,             &
!                                  segment_meas     & ! INOUT
!                               )

      return
      end subroutine set_segment_meas

!	sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

              
!> This function has to be called just before the inversion subroutine in order to
!> prepare segment and settings. Settings will be filled with some segment information 
!> like for example, wavelengths
subroutine prepare_segment_settings(iu_main_output, & ! IN
                                    segment_meas,   & ! INOUT
                                    RIN             ) ! INOUT
      use mod_sdata
      use mod_retr_settings_derived_type
      
      implicit none
! -----------------------------------------------------------------------------------------
      integer,                    intent(in)    ::  iu_main_output
      type(segment_data),         intent(inout) ::  segment_meas
      type(retr_input_settings),  intent(inout) ::  RIN
              
! Set RIN fields (radii, surface model flags for RT_OSH, number of retrieved parameters ).
! Set RIN with wave length array for retrieved parameters (combined array of wavelengths in case  
!     different pixels contain niether diff numbers or diff values of wavelengths).
! Set RIN with NDVI wavelength indices.
      call set_input_settings ( iu_main_output, & ! IN
                                segment_meas,   &
                                RIN             & ! INOUT
                              )
      if ( stop_report%status ) return
! Set segment_meas. Assign general wave length indices to pixel wave lengths.
! Set segment_meas. For real Parasol data, radiance measurements at wavelength 
!                   RIN%I_corr_iwl can be corrected.
      call set_segment_meas ( iu_main_output, & ! IN
                              RIN,            & 
                              segment_meas    & ! INOUT
                            )              
              
end subroutine prepare_segment_settings

