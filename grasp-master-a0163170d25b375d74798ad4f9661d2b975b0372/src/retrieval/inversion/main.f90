! **
! **  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.  
! **  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
! **

program Main
  use iso_fortran_env, only : error_unit
  use getoptions ! thanks to Dominik Epple
  use mod_globals, only : GBL_FILE_PATH_LEN, GBL_APPNAME, globals
  use mod_grasp
  use mod_stop_report
  
  implicit none
  
  character(*), parameter :: version = '0.1.0'

  character :: opt_key

  integer, parameter :: NARGS_MAX = 255
  character(GBL_FILE_PATH_LEN) :: arguments(NARGS_MAX)
  integer :: nargs
  
  character(GBL_FILE_PATH_LEN) :: sdata_file
  character(GBL_FILE_PATH_LEN) :: main_output_file = '-'            ! by default, main results will be sent to the standard output (screen)
  character(GBL_FILE_PATH_LEN) :: plotting_output_file = ''         ! by default, the results for plotting won't be written
  character(GBL_FILE_PATH_LEN) :: inversion_input_config_file = '-' ! by default, the configuration for inversion will be read from the standard input
                                                                    ! (like it was done before)
  character(GBL_FILE_PATH_LEN) :: sdata_sim_file = ''
  logical :: debug_mode = .false.
  logical :: force_mode = .false.
  
  nargs = 0

  do
     opt_key = getopt('c:dfo:p:s:')
     if (opt_key .eq. '>') then
        exit
     else if (opt_key .eq. 'd') then
        debug_mode = .true.
     else if (opt_key .eq. 'f') then
        force_mode = .true.
     else if (opt_key .eq. 'o') then
        main_output_file = optarg
     else if (opt_key .eq. 'p') then
        plotting_output_file = optarg
     else if (opt_key .eq. 's') then
	    sdata_sim_file = optarg
     else if (opt_key .eq. 'c') then
        inversion_input_config_file = optarg
     else if (opt_key .eq. '.') then
        if (nargs == NARGS_MAX) then
           write(error_unit,*) GBL_APPNAME, ': too many arguments, max allowed is ', NARGS_MAX
           write(error_unit,*)
           call usage()
        end if
        nargs = nargs + 1
        arguments(nargs) = optarg
     else
        write(error_unit,*) GBL_APPNAME, ': unknown option: ', trim(optarg)
        stop
     end if

  end do
  
  if (nargs /= 1) then
     call usage()
  end if
  
  globals%debug_mode = debug_mode

  if (debug_mode) then
    write(*,'(A,": debug mode is enabled")') GBL_APPNAME
  
    if (trim(inversion_input_config_file) .eq. '-') then
      write(*,'(A,": the inversion input configuration settings &
      &will be read from the standard input (redirection with <)")') GBL_APPNAME
    else
      write(*,'(A,": the inversion input configuration settings &
      &will be read from a file: ", A)') GBL_APPNAME, trim(inversion_input_config_file)
    end if
  
    if (trim(main_output_file) .eq. '-') then
      write(*,'(A,": main results will be sent to the standard output")') GBL_APPNAME
    else
      write(*,'(A,": main results will be sent to a file: ", A)') GBL_APPNAME, trim(main_output_file)
    end if
  
    if (trim(plotting_output_file) .ne. '') then
      write(*,'(A,": plotting results will be sent to a file: ", A)') GBL_APPNAME, trim(plotting_output_file)
    end if
	
  end if
    
  sdata_file = arguments(1)
  
  call check_files(inversion_input_config_file, sdata_file, main_output_file,  &
                   plotting_output_file,sdata_sim_file                         &
				  )

  call grasp_main(inversion_input_config_file, sdata_file, main_output_file,   &
                  plotting_output_file,sdata_sim_file                          &
				 )  
  if ( stop_report%status ) then
    call print_stop_report()
  endif
  stop

contains
 
  subroutine check_files(inversion_input_config_file, sdata_file, main_output_file,   &
                         plotting_output_file,sdata_sim_file)
    character(*), intent(in) :: inversion_input_config_file
    character(*), intent(in) :: sdata_file
    character(*), intent(in) :: main_output_file
    character(*), intent(in) :: plotting_output_file
    character(*), intent(in) :: sdata_sim_file
    
    if (.not. file_exists(sdata_file)) then
      write(error_unit, '(A, ": ",A," not found")') GBL_APPNAME, trim(sdata_file)
      stop
    end if
	
    if (trim(inversion_input_config_file) .ne. '-') then
      if (.not. file_exists(inversion_input_config_file)) then
        write(error_unit, '(A, ": ",A," not found")') GBL_APPNAME, trim(inversion_input_config_file)
        stop
      end if
    end if
    
    if (trim(main_output_file) .ne. '-') then
      if (file_exists(trim(main_output_file))) then
        write(error_unit, '(A, ": ",A," already exists")') GBL_APPNAME, trim(main_output_file)
        if (.not. force_mode) then
          stop
        else
          write(error_unit, '(A, ": ",A," will be overwritten (force mode enabled)")') GBL_APPNAME, trim(main_output_file)
        end if
      end if
    end if
	
    if (trim(plotting_output_file) .ne. '') then
      if (file_exists(trim(plotting_output_file))) then
        write(error_unit, '(A,": ",A," already exists")') GBL_APPNAME, trim(plotting_output_file)
        if (.not. force_mode) then
          stop
        else
          write(error_unit, '(A,": ",A," will be overwritten (force mode enabled)")') GBL_APPNAME, trim(plotting_output_file)
        end if
      end if
    end if
    if (trim(sdata_sim_file) .ne. '') then
      if (file_exists(trim(sdata_sim_file))) then
        write(error_unit, '(A,": ",A," already exists")') GBL_APPNAME, trim(sdata_sim_file)
        if (.not. force_mode) then
          stop
        else
          write(error_unit, '(A,": ",A," will be overwritten (force mode enabled)")') GBL_APPNAME, trim(sdata_sim_file)
        end if
      end if
    end if

  end subroutine check_files

  subroutine usage()
    write(error_unit, '("usage: ", A, " sdata_file [< inversion_input_config_file]")') GBL_APPNAME
    write(error_unit, '(A)') ! put an empty line
    write(error_unit, '("OPTIONS:")')
    write(error_unit, '("  -o main_output_file             sets the output file where to &
                        &send the inversion results (standard output by default)")')
    write(error_unit, '("  -p plotting_output_file         sets the output file where to &
                        &send the inversion results for plotting (none by default)")')
    write(error_unit, '("  -c inversion_input_config_file  sets the inversion input configuration &
                        &file (will be read from the standard input by default)")')
    write(error_unit, '("  -s sdata_sim_file               sets the output file where to &
                        &send simulated data (none by default)")')
    write(error_unit, '("  -d                              sets the debug mode")')
    write(error_unit, '("  -f                              forces the application to overwrite &
                        &existing output files with the same name (will abort by default)")')
    write(error_unit, '(A)') ! put an empty line

    stop
  end subroutine usage
  
  function file_exists(pathname)
        character(*), intent(in) :: pathname
        logical :: file_exists
        
        inquire(file=trim(pathname), exist=file_exists)
  end function

! unfortunately, the intrinsic routine inquire's behaviour concerning directories is known to be not portable.
! Till a propoer solution is found, one has to write specific codes for each compiler to be supported.
!#ifdef __GFORTRAN__
!  function directory_exists(pathname)
!    character(*), intent(in) :: pathname
!    logical :: directory_exists

!    inquire(file=trim(pathname), exist=directory_exists)
!  end function
!#elif defined (__INTEL_COMPILER)
!  function directory_exists(pathname)
!    character(*), intent(in) :: pathname
!    logical :: directory_exists
    
!    inquire(directory=trim(pathname), exist=directory_exists)
!  end function
!#else
!#error "The code has been tested so far with ifort or gfortran. For other compilers' support, please contact the GRASP team."
!#endif


end program Main
