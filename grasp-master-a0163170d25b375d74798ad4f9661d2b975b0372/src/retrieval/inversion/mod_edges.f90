! **
! **  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.  
! **  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
! **

!> @file mod_edges.f90

        !> File contains module with edge structure definition and edge related routines.
        !>
        !> @author Tatsiana Lapionak
        !> @date 20 APR 2015 
        !>
        !> @brief Module of structure (derived type) for edges data.  

MODULE mod_edges

      use iso_c_binding
      use mod_time_utils
      use mod_par_inv, only : KIX,KIY,KITIME,KIEDGE,KPARS             
      	  
      implicit none
! -----------------------------------------------------------------------
! Structure contains information for edges 
!
!> X groups (1 - before, 2 - after)
    type, bind(C) :: edges_group_X
        ! nx = input_settings%edges%nx
        ! ny = segment_meas%NY
        ! nt = number of times with edge data; nt <= segment%NT
        integer(kind=C_INT)   ::  nx !< nx - size of X edges; provided in retrieval settings
        integer(kind=C_INT)   ::  ny !< ny - Y size of main measurement segment
        integer(kind=C_INT)   ::  nt !< nt - size of T edges with data; nt <= segment%NT   
        integer(kind=C_INT)   ::  icloud(KIEDGE,KIY,KITIME) !< icloud(KIEDGE,KIY,KITIME) - cloud present indicator      
        real(kind=C_FLOAT)    ::  x(KIEDGE,KIY,KITIME)      !< x(KIEDGE,KIY,KITIME) - pixel longitude (X)
        real(kind=C_FLOAT)    ::  y(KIEDGE,KIY,KITIME)      !< y(KIEDGE,KIY,KITIME) - pixel latitude  (Y)
        integer(KIND_TIME)    ::  t(KITIME)                 !< t(KITIME)            - pixel time      (T)
        integer(kind=C_INT)   ::  out_x(KIEDGE,KIY,KITIME)  !< x(KIEDGE,KIY,KITIME) - absolute position of pixel inside the tile in x dimension
        integer(kind=C_INT)   ::  out_y(KIEDGE,KIY,KITIME)  !< y(KIEDGE,KIY,KITIME) - absolute position of pixel inside the tile in y dimension
        integer(kind=C_INT)   ::  out_t(KITIME)             !< t(KITIME)            - absolute position of pixel inside the tile in t dimension
        integer(kind=C_INT)   ::  it(KITIME)                !< it(KITIME) - array contains time indices for time with data
        real(kind=C_FLOAT)    ::  AP(KPARS,KIEDGE,KIY,KITIME) !<  AP(KPARS,KIEDGE,KIY,KITIME) - parameters in edge data
    end type edges_group_X

!> Y groups (1 - before, 2 - after)
    type, bind(C) :: edges_group_Y
        ! nx = segment_meas%NX
        ! ny = input_settings%edges%ny
        ! nt = number of times with edge data; nt <= segment%NT
        integer(kind=C_INT)   ::  nx !< nx - X size of main measurement segment
        integer(kind=C_INT)   ::  ny !< ny - size of Y edges; provided in retrieval settings
        integer(kind=C_INT)   ::  nt !< nt - size of T edges with data; nt <= segment%NT     
        integer(kind=C_INT)   ::  icloud(KIX,KIEDGE,KITIME) !< icloud(KIEDGE,KIY,KITIME) - cloud present indicator
        real(kind=C_FLOAT)    ::  x(KIX,KIEDGE,KITIME)      !< x(KIX,KIEDGE,KITIME) - pixel longitude (X)
        real(kind=C_FLOAT)    ::  y(KIX,KIEDGE,KITIME)      !< y(KIX,KIEDGE,KITIME) - pixel latitude  (Y)
        integer(KIND_TIME)    ::  t(KITIME)                 !< t(KITIME)            - pixel time      (T)
        integer(kind=C_INT)   ::  out_x(KIX,KIEDGE,KITIME)  !< x(KIX,KIEDGE,KITIME) - absolute position of pixel inside the tile in x dimension
        integer(kind=C_INT)   ::  out_y(KIX,KIEDGE,KITIME)  !< y(KIX,KIEDGE,KITIME) - absolute position of pixel inside the tile in y dimension
        integer(kind=C_INT)   ::  out_t(KITIME)             !< t(KITIME)            - absolute position of pixel inside the tile in t dimension
        integer(kind=C_INT)   ::  it(KITIME)                !< it(KITIME) - array contains time indices for time with data
        real(kind=C_FLOAT)    ::  AP(KPARS,KIX,KIEDGE,KITIME) !< AP(KPARS,KIX,KIEDGE,KITIME) - parameters in edge data
    end type edges_group_Y
                
!> T groups (1 - before, 2 - after)
    type, bind(C) :: edges_group_T
        ! nx = segment_meas%NX
        ! ny = segment_meas%NY
        ! nt = input_settings%edges%nt
        integer(kind=C_INT)   ::  nx !< nx - X size of main measurement segment
        integer(kind=C_INT)   ::  ny !< ny - Y size of main measurement segment
        integer(kind=C_INT)   ::  nt !< nt - size of T edges with data; nt <= input_settings%edges%nt
        integer(kind=C_INT)   ::  icloud(KIX,KIY,KIEDGE) !< icloud(KIEDGE,KIY,KITIME) - cloud present indicator
        real(kind=C_FLOAT)    ::  x(KIX,KIY,KIEDGE)      !< x(KIX,KIY,KIEDGE) - pixel longitude (X)
        real(kind=C_FLOAT)    ::  y(KIX,KIY,KIEDGE)      !< y(KIX,KIY,KIEDGE) - pixel latitude  (Y)
        integer(KIND_TIME)    ::  t(KIEDGE)              !< t(KIEDGE)         - pixel time      (T)
        integer(kind=C_INT)   ::  out_x(KIX,KIY,KIEDGE)  !< x(KIX,KIY,KIEDGE) - absolute position of pixel inside the tile in x dimension
        integer(kind=C_INT)   ::  out_y(KIX,KIY,KIEDGE)  !< y(KIX,KIY,KIEDGE) - absolute position of pixel inside the tile in y dimension
        integer(kind=C_INT)   ::  out_t(KIEDGE)          !< t(KIEDGE)         - absolute position of pixel inside the tile in t dimension
        integer(kind=C_INT)   ::  it(KIEDGE)             !< it(KIEDGE) - array contains time indices for time with data
        real(kind=C_FLOAT)    ::  AP(KPARS,KIX,KIY,KIEDGE) !< AP(KPARS,KIX,KIY,KIEDGE) - parameters in edge data
    end type edges_group_T

!> segment edge structure 
    type, bind(C) :: segment_edges
        integer(kind=C_INT)   ::  N_I_EDGE   !< N_I_EDGE - number of edge groups with data (maximum 6: X before/after, Y before/after, T before/after)
        integer(kind=C_INT)   ::  I_EDGE(6)  !< I_EDGE(6) - present group index array I_EDGE(1:N_I_EDGE)
        type(edges_group_X)   ::  group_X(2) !< group_X(2) - array of present group indices (group(1:N_I_EDGE))
        type(edges_group_Y)   ::  group_Y(2) !< group_Y(2) - array of present group indices (group(1:N_I_EDGE))
        type(edges_group_T)   ::  group_T(2) !< group_T(2) - array of present group indices (group(1:N_I_EDGE))      
    end type segment_edges
! -----------------------------------------------------------------------

END MODULE mod_edges

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

        !> Routine initializes matrices edges structure
        !>
        !> @author Tatsiana Lapionak
        !> @date 20 APR 2015
        !> 
        !> @param[inout]    edges - structure of edges data

      subroutine edges_initialization (edges) 

      use mod_edges

      implicit none
! -----------------------------------------------------------------------------------------	  
      type(segment_edges), intent(inout) :: edges
      integer   :: i	  
! -----------------------------------------------------------------------------------------
      edges%N_I_EDGE  = 0
      edges%I_EDGE(:) = 0  
      do i=1,2
        edges%group_X(i)%nx = 0
        edges%group_X(i)%ny = 0
        edges%group_X(i)%nt = 0
        edges%group_X(i)%it(:) = 0
        edges%group_X(i)%icloud(:,:,:) = 0 
        edges%group_X(i)%x(:,:,:) = 0.0
        edges%group_X(i)%y(:,:,:) = 0.0
        edges%group_X(i)%t(:)     = 0.0
        edges%group_X(i)%out_x(:,:,:) = -1 ! -1 is invalid value (not set). Valid index are in range 0..tile_with-1
        edges%group_X(i)%out_y(:,:,:) = -1
        edges%group_X(i)%out_t(:)     = -1
        edges%group_X(i)%AP(:,:,:,:) = -999.0
      enddo ! i
      do i=1,2
        edges%group_Y(i)%nx = 0
        edges%group_Y(i)%ny = 0
        edges%group_Y(i)%nt = 0
        edges%group_Y(i)%it(:) = 0
        edges%group_Y(i)%icloud(:,:,:) = 0 
        edges%group_Y(i)%x(:,:,:) = 0.0
        edges%group_Y(i)%y(:,:,:) = 0.0
        edges%group_Y(i)%t(:)     = 0.0
        edges%group_Y(i)%out_x(:,:,:) = -1
        edges%group_Y(i)%out_y(:,:,:) = -1
        edges%group_Y(i)%out_t(:)     = -1        
        edges%group_Y(i)%AP(:,:,:,:) = -999.0
      enddo ! i
      do i=1,2
        edges%group_T(i)%nx = 0
        edges%group_T(i)%ny = 0
        edges%group_T(i)%nt = 0
        edges%group_T(i)%it(:) = 0
        edges%group_T(i)%icloud(:,:,:) = 0 
        edges%group_T(i)%x(:,:,:) = 0.0
        edges%group_T(i)%y(:,:,:) = 0.0
        edges%group_T(i)%t(:)     = 0.0
        edges%group_T(i)%out_x(:,:,:) = -1
        edges%group_T(i)%out_y(:,:,:) = -1
        edges%group_T(i)%out_t(:)     = -1        
        edges%group_T(i)%AP(:,:,:,:) = -999.0
      enddo ! i

      end subroutine edges_initialization

      ! C interface for above function
    subroutine grasp_input_edges_initialization(edges) bind(C)   
        use mod_edges
        type(segment_edges),               intent(inout)  :: edges   
         
        call edges_initialization(edges)
    end subroutine grasp_input_edges_initialization
! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

        !> Routine reads matrices edges data from "input_edge.dat" file.
        !>
        !> @author Tatsiana Lapionak
        !> @date 20 APR 2015 
        !>
        !> @param[in]     iu_main_output - unit number of main output file of retrieval
        !> @param[in]     RIN   - structure of settings for retrieval
        !> @param[inout]  edges - structure of edges data

      subroutine read_edges(iu_main_output,RIN,edges)

      use mod_edges
      use mod_retr_settings_derived_type
      !use iso_c_binding
      use mod_par_inv, only : KIX,KIY,KITIME,KIEDGE,KPARS
                  
      implicit none
!	------------------------------------------------------------------------------------------
      integer,                    intent(in)    ::  iu_main_output
      type(retr_input_settings),  intent(in)    ::  RIN
      type(segment_edges),        intent(inout) ::  edges
!	------------------------------------------------------------------------------------------
      integer ::  i,II_E,IT_E,IX,IY,IT,ipix,npar
      integer ::  iEDGE_input,NPIXT,N_I_EDGE,I_EDGE(6),NX_edge,NY_edge,NT_edge
      character(LEN=20)  ::	iso8601_string      
		  integer(KIND_TIME) ::	TIME_edge 
      logical            :: status_funct 

      integer,dimension(:,:,:),allocatable :: ICLOUD_edge
      real,dimension(:,:,:),   allocatable :: X_edge,Y_edge
      real,dimension(:,:,:,:),   allocatable :: AP_edge
      integer,parameter        :: KIX_E = max(KIEDGE,KIX)
      integer,parameter        :: KIY_E = max(KIEDGE,KIY)
      integer,parameter        :: KIT_E = max(KIEDGE,KITIME)
      integer                  :: alloc_stat
!	------------------------------------------------------------------------------------------      

      if(.not. allocated(ICLOUD_edge)) then
        allocate(ICLOUD_edge(KIX_E,KIY_E,KIT_E),stat=alloc_stat)       
        if (alloc_stat /= 0) then
        stop 'error while trying to allocate ICLOUD_edge(:,:,:) in read_edges'
        endif              
      endif
      if(.not. allocated(X_edge)) then
        allocate(X_edge(KIX_E,KIY_E,KIT_E),stat=alloc_stat)       
        if (alloc_stat /= 0) then
        stop 'error while trying to allocate X_edge(:,:,:) in read_edges'
        endif              
      endif
      if(.not. allocated(Y_edge)) then
        allocate(Y_edge(KIX_E,KIY_E,KIT_E),stat=alloc_stat)       
        if (alloc_stat /= 0) then
        stop 'error while trying to allocate Y_edge(:,:,:) in read_edges'
        endif              
      endif
      if(.not. allocated(AP_edge)) then
        allocate(AP_edge(KPARS,KIX_E,KIY_E,KIT_E),stat=alloc_stat)       
        if (alloc_stat /= 0) then
        stop 'error while trying to allocate AP_edge(:,:,:) in read_edges'
        endif              
      endif

      open (newunit=iEDGE_input, FILE=trim(RIN%DLSF%internal_file_path)//"input_edge.dat",status='old')
!* Read parameters for edge groups:
!*      NX_before NX_after NY_before NY_after NT_before NT_after
        read(iEDGE_input,*) N_I_EDGE, npar        
        write(iu_main_output,*) N_I_EDGE, npar,'  - N_I_EDGE, npar'
        if(N_I_EDGE .eq. 0) return
        
LOOP_I_EDGE: do II_E=1,N_I_EDGE ! MAX of N_I_EDGE = 6
!*  X - longitude, Y - latitude, T - time
!*  I_EDGE = 1 - X_before 
!*         = 2 - X_after 
!*         = 3 - Y_before
!*         = 4 - Y_after
!*         = 5 - T_before
!*         = 6 - T_after
          read(iEDGE_input,*) I_EDGE(II_E),NX_edge,NY_edge,NT_edge
          write(iu_main_output,*)
          write(iu_main_output,*) II_E,I_EDGE(II_E),NX_edge,NY_edge,NT_edge,  &
                '  - II_E,I_EDGE(II_E),NX_edge,NY_edge,NT_edge; in read_edges'
          do IT_E=1,NT_edge
            read(iEDGE_input,*) NPIXT,IT,iso8601_string
            write(iu_main_output,*) NPIXT,IT,iso8601_string,'  - NPIXT,IT,iso8601_string'
            TIME_edge = convert_string_to_time(iso8601_string, status_funct)
            if (.not. status_funct) then
              write(iu_main_output,*) 'iso8601_string=',iso8601_string
              stop 'read_edges: error on convert_string_to_time'
            end if

            do ipix = 1,NPIXT
              read(iEDGE_input,*) IX,IY,ICLOUD_edge(IX,IY,IT),X_edge(IX,IY,IT),Y_edge(IX,IY,IT),  &
                                                       AP_edge(1:npar,IX,IY,IT)
              write(iu_main_output,*) IX,IY,ICLOUD_edge(IX,IY,IT),X_edge(IX,IY,IT),Y_edge(IX,IY,IT), &
                    '  - IX,IY,ICLOUD_edge(IX,IY,IT),X_edge(IX,IY,IT),Y_edge(IX,IY,IT)'
              select case(II_E)
              case(1,2) ! X groups
                edges%group_X(II_E)%nx = NX_edge
                edges%group_X(II_E)%ny = NY_edge
                edges%group_X(II_E)%nt = NT_edge
                edges%group_X(II_E)%t(IT) = TIME_edge
                edges%group_X(II_E)%icloud(IX,IY,IT) = ICLOUD_edge(IX,IY,IT)           
                edges%group_X(II_E)%x(IX,IY,IT) = X_edge(IX,IY,IT)
                edges%group_X(II_E)%y(IX,IY,IT) = Y_edge(IX,IY,IT)
                edges%group_X(II_E)%AP(1:npar,IX,IY,IT) = AP_edge(1:npar,IX,IY,IT)
              case(3,4) ! Y groups
                edges%group_Y(II_E-2)%nx = NX_edge
                edges%group_Y(II_E-2)%ny = NY_edge
                edges%group_Y(II_E-2)%nt = NT_edge
                edges%group_Y(II_E-2)%t(IT) = TIME_edge
                edges%group_Y(II_E-2)%icloud(IX,IY,IT) = ICLOUD_edge(IX,IY,IT)           
                edges%group_Y(II_E-2)%x(IX,IY,IT) = X_edge(IX,IY,IT)
                edges%group_Y(II_E-2)%y(IX,IY,IT) = Y_edge(IX,IY,IT)
                edges%group_Y(II_E-2)%AP(1:npar,IX,IY,IT) = AP_edge(1:npar,IX,IY,IT)            
              case(5,6) ! T groups
                edges%group_T(II_E-4)%nx = NX_edge
                edges%group_T(II_E-4)%ny = NY_edge
                edges%group_T(II_E-4)%nt = NT_edge
                edges%group_T(II_E-4)%t(IT) = TIME_edge
                edges%group_T(II_E-4)%icloud(IX,IY,IT) = ICLOUD_edge(IX,IY,IT)           
                edges%group_T(II_E-4)%x(IX,IY,IT) = X_edge(IX,IY,IT)
                edges%group_T(II_E-4)%y(IX,IY,IT) = Y_edge(IX,IY,IT)
                edges%group_T(II_E-4)%AP(1:npar,IX,IY,IT) = AP_edge(1:npar,IX,IY,IT)
              end select
            enddo ! ipix                        
          enddo ! IT_E
enddo LOOP_I_EDGE
      close(iEDGE_input) 

      if(allocated(ICLOUD_edge)) then
        deallocate(ICLOUD_edge,stat=alloc_stat)       
        if (alloc_stat /= 0) then
        stop 'error while trying to deallocate ICLOUD_edge(:,:,:) in read_edges'
        endif              
      endif
      if(allocated(X_edge)) then
        deallocate(X_edge,stat=alloc_stat)       
        if (alloc_stat /= 0) then
        stop 'error while trying to deallocate X_edge(:,:,:) in read_edges'
        endif              
      endif
      if(.not. allocated(Y_edge)) then
        deallocate(Y_edge,stat=alloc_stat)       
        if (alloc_stat /= 0) then
        stop 'error while trying to deallocate Y_edge(:,:,:) in read_edges'
        endif              
      endif
      if(allocated(AP_edge)) then
        deallocate(AP_edge,stat=alloc_stat)       
        if (alloc_stat /= 0) then
        stop 'error while trying to deallocate AP_edge(:,:,:,:) in read_edges'
        endif              
      endif
         
      end subroutine read_edges
          
! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss          

        !> Function checks if data for edges are present.
        !> If present, function sets the data into edge structure and returns edges_present = .true.
        !>
        !> @author Tatsiana Lapionak
        !> @date 20 APR 2015
        !> 
        !> @param[in]     iu_main_output - unit number of main output file of retrieval
        !> @param[in]     nx_segm - max number of pixels with different X coordinates for retrieved segment
        !> @param[in]     ny_segm - max number of pixels with different Y coordinates for retrieved segment
        !> @param[in]     nt_segm - max number of pixels with different T coordinates for retrieved segment
        !> @param[in]     RIN   - structure of settings for retrieval
        !> @param[inout]  edges - structure of edges data

      function edges_present(iu_main_output,nx_segm,ny_segm,nt_segm,RIN,edges)

      use mod_edges
      use mod_retr_settings_derived_type
      
      implicit none
!	------------------------------------------------------------------------------------------
      integer,                    intent(in)    ::  iu_main_output,nx_segm,ny_segm,nt_segm
      type(retr_input_settings),  intent(in)    ::  RIN
      type(segment_edges),        intent(inout) ::  edges
      logical ::  edges_present
!	------------------------------------------------------------------------------------------
      integer ::  i,II_E,IX,IY,IT,ipix,it_counter 
      integer ::  NPIXT,N_I_EDGE,I_EDGE(6) 
      integer ::  nx_edge_max,ny_edge_max,nt_edge_max
      integer ::  nx_temp,ny_temp,nt_temp
      integer ::  icloud_X(KIEDGE,KIY,KITIME),  &
                  icloud_Y(KIX,KIEDGE,KITIME),  &
                  icloud_T(KIX,KIY,KIEDGE)
      logical ::  lstop = .false.     
      logical ::  ltest = .false.
!	------------------------------------------------------------------------------------------
      edges_present = .false.

      if( edges%group_X(1)%nx .eq. 0 .and. edges%group_X(2)%nx .eq. 0 .and. &
          edges%group_Y(1)%ny .eq. 0 .and. edges%group_Y(2)%ny .eq. 0 .and. &
          edges%group_T(1)%nt .eq. 0 .and. edges%group_T(2)%nt .eq. 0 ) then
! Retrieval without edges was asked by settings
        return
      endif

      if(RIN%IPRI_verbose) then
          do II_E=1,6
            select case(II_E) 
            case(1,2)
              i = II_E
              if(edges%group_X(i)%nx .ne. RIN%edges%nx) then
                write(iu_main_output,'(4(a,i3))') 'edge group #',II_E,'  i=',i,  &
                '  edges%group_X(i)%nx=',edges%group_X(i)%nx,' .ne. RIN%edges%nx=',RIN%edges%nx
                lstop = .true.
              elseif(edges%group_X(i)%ny .ne. ny_segm) then
                write(iu_main_output,'(4(a,i3))') 'edge group #',II_E,'  i=',i,  &
                '  edges%group_X(i)%ny=',edges%group_X(i)%ny,  ' .ne. ny_segm=',ny_segm
                lstop = .true.
              elseif(edges%group_X(i)%nt .gt. nt_segm) then  
                write(iu_main_output,'(4(a,i3))') 'edge group #',II_E,'  i=',i,  &
                '  edges%group_X(i)%nt=',edges%group_X(i)%nt,' .gt. nt_segm=',nt_segm
                lstop = .true.
              endif
              write(iu_main_output,'(2(a,i3),3i3,a)') 'edge group #',II_E,'  i=',i, &
              edges%group_X(i)%nx, edges%group_X(i)%ny, edges%group_X(i)%nt,'  nx, ny, nt - edges'
              write(iu_main_output,'(2(a,i3),3i3,a,/)') 'edge group #',II_E,'  i=',i, &
              RIN%edges%nx, ny_segm, nt_segm,'  RIN%edges%nx, ny_segm, nt_segm'
            case(3,4)
              i = II_E-2
              if(edges%group_Y(i)%nx .ne. nx_segm)          then
                write(iu_main_output,'(4(a,i3))') 'edge group #',II_E,'  i=',i,  &
                '  edges%group_Y(i)%nx=',edges%group_Y(i)%nx,' .ne. nx_segm=',nx_segm
                lstop = .true.
              elseif(edges%group_Y(i)%ny .ne. RIN%edges%ny) then
                write(iu_main_output,'(4(a,i3))') 'edge group #',II_E,'  i=',i,  &
                '  edges%group_Y(i)%ny=',edges%group_Y(i)%ny,' .ne. RIN%edges%ny=',RIN%edges%ny
                lstop = .true.
              elseif(edges%group_Y(i)%nt .gt. nt_segm)      then  
                write(iu_main_output,'(4(a,i3))') 'edge group #',II_E,'  i=',i,  &
                '  edges%group_Y(i)%nt=',edges%group_Y(i)%nt,' .gt. nt_segm=',nt_segm
                lstop = .true.
              endif
              write(iu_main_output,'(2(a,i3),3i3,a)') 'edge group #',II_E,'  i=',i, &
              edges%group_Y(i)%nx, edges%group_Y(i)%ny, edges%group_Y(i)%nt,'  nx, ny, nt - edges'
              write(iu_main_output,'(2(a,i3),3i3,a,/)') 'edge group #',II_E,'  i=',i, &
              nx_segm, RIN%edges%ny, nt_segm,'  nx_segm, RIN%edges%ny, nt_segm'
            case(5,6)
              i = II_E-4
              if(edges%group_T(i)%nx .ne. nx_segm)          then
                write(iu_main_output,'(4(a,i3))') 'edge group #',II_E,'  i=',i,  &
                '  edges%group_T(i)%nx=',edges%group_T(i)%nx,' .ne. nx_segm=',nx_segm
                lstop = .true.
              elseif(edges%group_T(i)%ny .ne. ny_segm)      then
                write(iu_main_output,'(4(a,i3))') 'edge group #',II_E,'  i=',i,  &
                '  edges%group_T(i)%ny=',edges%group_T(i)%ny,' .ne. ny_segm=',ny_segm
                lstop = .true.
              elseif(edges%group_T(i)%nt .gt. RIN%edges%nt) then  
                write(iu_main_output,'(4(a,i3))') 'edge group #',II_E,'  i=',i,  &
                '  edges%group_T(i)%nt=',edges%group_T(i)%nt,' .gt. RIN%edges%nt=',RIN%edges%nt
                lstop = .true.
              endif
              write(iu_main_output,'(2(a,i3),3i3,a)') 'edge group #',II_E,'  i=',i, &
              edges%group_T(i)%nx, edges%group_T(i)%ny, edges%group_T(i)%nt,'  nx, ny, nt - edges'
              write(iu_main_output,'(2(a,i3),3i3,a,/)') 'edge group #',II_E,'  i=',i, &
              nx_segm, ny_segm, RIN%edges%nt,'  nx_segm, ny_segm, RIN%edges%nt'
            end select
          enddo ! II_E
          if(lstop) then
            write(iu_main_output,'(a)') 'STOP 1 in edges_present; execution has to be interrupted'    
            stop 
          endif
      endif ! .not. RIN%IPRI_silent_mode
      
      if(ltest) then
          ! Print edge retrieved parameters
          write(iu_main_output,'(a)') 'AP_X_before:'
          do it=1,2
            write(iu_main_output,'(10e16.5)') edges%group_X(1)%AP(1:RIN%KNSING,1,1,it)
          enddo
          write(iu_main_output,'(a)') 'AP_X_after :'
          do it=1,2
            write(iu_main_output,'(10e16.5)') edges%group_X(2)%AP(1:RIN%KNSING,1,1,it)
          enddo
          write(iu_main_output,'(a)') 'AP_Y_before:'
          do it=1,2
            write(iu_main_output,'(10e16.5)') edges%group_Y(1)%AP(1:RIN%KNSING,1,1,it)
          enddo
          write(iu_main_output,'(a)') 'AP_Y_after :'
          do it=1,2
            write(iu_main_output,'(10e16.5)') edges%group_Y(2)%AP(1:RIN%KNSING,1,1,it)
          enddo                  
          write(iu_main_output,'(a)') 'AP_T_before:'
          do it=1,2
            write(iu_main_output,'(10e16.5)') edges%group_T(1)%AP(1:RIN%KNSING,1,1,it)
          enddo
          write(iu_main_output,'(a)') 'AP_T_after :'
          do it=1,2
            write(iu_main_output,'(10e16.5)') edges%group_T(2)%AP(1:RIN%KNSING,1,1,it)
          enddo
          ! Print size of edges provided by framework    
          do II_E=1,6
            select case(II_E) 
            case(1,2)
              i = II_E
              write(*,'(a,2i3,a,3i3,a)') 'II_E,edge_group_X: ',II_E,i,  &
              '    nx,ny,nt: ',edges%group_X(i)%nx,edges%group_X(i)%ny,edges%group_X(i)%nt,'  - provided by framework'
            case(3,4)
              i = II_E-2
              write(*,'(a,2i3,a,3i3,a)') 'II_E,edge_group_Y: ',II_E,i,  &
              '    nx,ny,nt: ',edges%group_Y(i)%nx,edges%group_Y(i)%ny,edges%group_Y(i)%nt,'  - provided by framework'
            case(5,6)
              i = II_E-4
              write(*,'(a,2i3,a,3i3,a)') 'II_E,edge_group_T: ',II_E,i,  &
              '    nx,ny,nt: ',edges%group_T(i)%nx,edges%group_T(i)%ny,edges%group_T(i)%nt,'  - provided by framework'
            end select
          enddo ! II_E
          ! Print edge pixels provided by framework 
          do II_E=1,6
            select case(II_E) 
            case(1,2)
              i = II_E
              do IT=1,edges%group_X(i)%nt
              do IY=1,edges%group_X(i)%ny
              do IX=1,edges%group_X(i)%nx
                write(*,'(a,5i3,a,3i5,a)') 'edge_group_X,IX,IY,IT,icloud: ',i,IX,IY,IT,edges%group_X(i)%icloud(IX,IY,IT),     &
                                                 '     out_x,out_y,out_t: ',edges%group_X(i)%out_x(IX,IY,IT),                 & 
                                                                            edges%group_X(i)%out_y(IX,IY,IT),                 & 
                                                                            edges%group_X(i)%out_t(IT),'  - in edges_present' 
              enddo ! IX
              enddo ! IY
              enddo ! IT  
            case(3,4)
              i = II_E-2
              do IT=1,edges%group_Y(i)%nt
              do IY=1,edges%group_Y(i)%ny
              do IX=1,edges%group_Y(i)%nx
                write(*,'(a,5i3,a,3i5,a)') 'edge_group_Y,IX,IY,IT,icloud: ',i,IX,IY,IT,edges%group_Y(i)%icloud(IX,IY,IT),     &
                                                 '     out_x,out_y,out_t: ',edges%group_Y(i)%out_x(IX,IY,IT),                 & 
                                                                            edges%group_Y(i)%out_y(IX,IY,IT),                 & 
                                                                            edges%group_Y(i)%out_t(IT),'  - in edges_present' 
              enddo ! IX
              enddo ! IY
              enddo ! IT  
            case(5,6)
              i = II_E-4
              do IT=1,edges%group_T(i)%nt
              do IY=1,edges%group_T(i)%ny
              do IX=1,edges%group_T(i)%nx
                write(*,'(a,5i3,a,3i5,a)') 'edge_group_T,IX,IY,IT,icloud: ',i,IX,IY,IT,edges%group_T(i)%icloud(IX,IY,IT),     &
                                                 '     out_x,out_y,out_t: ',edges%group_T(i)%out_x(IX,IY,IT),                 & 
                                                                            edges%group_T(i)%out_y(IX,IY,IT),                 & 
                                                                            edges%group_T(i)%out_t(IT),'  - in edges_present' 
              enddo ! IX
              enddo ! IY
              enddo ! IT 
            end select
          enddo ! II_E
      endif ! ltest

      nx_edge_max = RIN%edges%nx
      ny_edge_max = RIN%edges%ny        
      nt_edge_max = RIN%edges%nt 
      if(RIN%IPRI_verbose) write(iu_main_output,'(a,3i3,a)')   & 
      'NX_edge,NY_edge,NT_edge: ',RIN%edges%nx,RIN%edges%ny,RIN%edges%nt,'  - sizes for edges provided by settings '
      N_I_EDGE    = 0
      I_EDGE(1:6) = 0 
 
LOOP_group: do II_E=1,6

! Counting groups and fill in edges structure with number of edges groups and array of group indices 
      select case(II_E)
      case(1,2)
! X groups (1 - before; 2 - after)
        i = II_E
        nt_temp = 0
        edges%group_X(i)%it(:) = 0
        do IT=1,nt_segm
          it_counter = 0
          do IY=1,ny_segm
            do IX=1,nx_edge_max
              if(edges%group_X(i)%icloud(IX,IY,IT) .eq. 1) then
                it_counter = it_counter+1
                !write(*,'(a,2i3,a,3i3,a,i3)') 'II_E,edge_group_X: ',II_E,i,'  ix,iy,it: ',ix,iy,it,  &
                !'  icloud(framework)=',edges%group_X(i)%icloud(IX,IY,IT)
                !write(iu_main_output,'(a)') 'AP_X:'
                !write(*,'(10e16.4)') edges%group_X(i)%AP(1:RIN%KNSING,IX,IY,IT)
              endif
            enddo ! IX
          enddo ! IY
          if(it_counter .gt. 0) then
            nt_temp = nt_temp+1
            edges%group_X(i)%it(nt_temp) = IT
          endif
        enddo ! IT
        if(nt_temp .gt. 0) then
          N_I_EDGE  = N_I_EDGE+1
          I_EDGE(N_I_EDGE) = II_E 
        endif
        if(RIN%IPRI_verbose) write(iu_main_output,'(2(a,2i3),a)') 'II_E,edge_group_X: ',II_E,i,  &
        '    edges%group_X(i)%nt,nt_temp: ',edges%group_X(i)%nt,nt_temp,'  - test time in edges_present'
      case(3,4)
! Y groups (3 - before; 4 - after)
        i = II_E-2
        nt_temp = 0
        edges%group_Y(i)%it(:) = 0
        do IT=1,nt_segm
          it_counter = 0
          do IY=1,ny_edge_max
            do IX=1,nx_segm
              if(edges%group_Y(i)%icloud(IX,IY,IT) .eq. 1) then
                it_counter = it_counter+1
                !write(*,'(a,2i3,a,3i3,a,i3)') 'II_E,edge_group_Y: ',II_E,i,'  ix,iy,it: ',ix,iy,it,  &
                !'  icloud(framework)=',edges%group_Y(i)%icloud(IX,IY,IT)
                !write(iu_main_output,'(a)') 'AP_Y:'
                !write(*,'(10e16.4)') edges%group_Y(i)%AP(1:RIN%KNSING,IX,IY,IT)
              endif
            enddo ! IX
          enddo ! IY
          if(it_counter .gt. 0) then
            nt_temp = nt_temp+1
            edges%group_Y(i)%it(nt_temp) = IT
          endif
        enddo ! IT
        if(nt_temp .gt. 0) then
          N_I_EDGE  = N_I_EDGE+1
          I_EDGE(N_I_EDGE) = II_E 
        endif
        if(RIN%IPRI_verbose) write(iu_main_output,'(2(a,2i3),a)') 'II_E,edge_group_Y: ',II_E,i,  &
        '    edges%group_Y(i)%nt,nt_temp: ',edges%group_Y(i)%nt,nt_temp,'  - test time in edges_present'
      case(5,6) ! T edges (5 - before; 6 - after)
! T groups (5 - before; 6 - after)
        i = II_E-4
        nt_temp = 0
        edges%group_T(i)%it(:) = 0
        do IT=1,nt_edge_max
          it_counter = 0
          do IY=1,ny_segm
            do IX=1,nx_segm
              if(edges%group_T(i)%icloud(IX,IY,IT) .eq. 1) then
                it_counter = it_counter+1
                !write(*,'(a,2i3,a,3i3,a,i3)') 'II_E,edge_group_T: ',II_E,i,'  ix,iy,it: ',ix,iy,it,  &
                !'  icloud(framework)=',edges%group_T(i)%icloud(IX,IY,IT)
                !write(iu_main_output,'(a)') 'AP_T:'
                !write(*,'(10e16.4)') edges%group_T(i)%AP(1:RIN%KNSING,IX,IY,IT)
              endif
            enddo ! IX
          enddo ! IY
          if(it_counter .gt. 0) then
            nt_temp = nt_temp+1
            edges%group_T(i)%it(nt_temp) = IT
          endif
        enddo ! IT
        if(nt_temp .gt. 0) then
          N_I_EDGE  = N_I_EDGE+1
          I_EDGE(N_I_EDGE) = II_E 
        endif
        if(RIN%IPRI_verbose) write(iu_main_output,'(2(a,2i3),a)') 'II_E,edge_group_T: ',II_E,i,  &
        '    edges%group_T(i)%nt,nt_temp: ',edges%group_T(i)%nt,nt_temp,'  - test time in edges_present' 
      case default
        if(RIN%IPRI_verbose) &
        write(iu_main_output,*) 'group# =',II_E,' is not exist (1 <= group# <= 6)'
      end select

enddo LOOP_group

      if(N_I_EDGE .gt. 0) then
        edges_present = .true.
        edges%N_I_EDGE = N_I_EDGE
        edges%I_EDGE(1:N_I_EDGE) = I_EDGE(1:N_I_EDGE) 
      endif ! N_I_EDGE .ge. 1
      if(RIN%IPRI_verbose) then
        write(iu_main_output,'(a,i4)') 'in edges_groups_present: N_I_EDGE=',N_I_EDGE                
        if(N_I_EDGE .gt. 0)  &
        write(iu_main_output,'(a,6i4)') 'in edges_groups_present: I_EDGE: ',I_EDGE(1:N_I_EDGE)        
        write(iu_main_output,'(a,l2)') 'in edges_groups_present: edges_present=',edges_present
      endif

!      if(lstop) then
      if(lstop .and. N_I_EDGE .gt. 0) then
        write(iu_main_output,'(a)') 'STOP 2 in edges_present; execution has to be interrupted'    
        stop 
      endif

!stop 'stop test in edges_present'                             
      end function edges_present
          
! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss 