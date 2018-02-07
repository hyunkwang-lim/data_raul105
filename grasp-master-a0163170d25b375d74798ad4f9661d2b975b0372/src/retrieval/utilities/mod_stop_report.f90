! **
! **  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.  
! **  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
! **

        !> @file mod_stop_report.f90
        !> Module contains derived types and routines related to report grasp execution stops
        !>

#include "../constants_set/mod_globals.inc"
module mod_stop_report

  use iso_fortran_env
  use mod_sdata_derived_type
  use mod_globals, only : GBL_FILE_PATH_LEN

  implicit none

  integer, parameter :: KLINEpar = 5
  integer, parameter :: ERROR_MESSAGE_LEN = 1024

  character (len=ERROR_MESSAGE_LEN) :: tmp_message

  type :: stop_report_type
    logical :: status
    character (len=ERROR_MESSAGE_LEN) :: message
    character (len=GBL_FILE_PATH_LEN) :: file
    integer :: line
  end type

  type :: crash_report_type
    type(segment_data), pointer :: segment
  end type

  type(crash_report_type) :: crash_report
  type(stop_report_type) :: stop_report

  contains

! **********************************************************************
        !> @brief Routine assigns Pointer to segment Target
        !> 
        !> @author Tatsiana Lapionak
        !> @date 10 APR 2017
        !>
        !> @param[in]  segment - inverted segment data

    subroutine set_segment_for_stop_report(segment)
      type(segment_data), target, intent(in) :: segment

      crash_report%segment => segment

    end subroutine set_segment_for_stop_report

! **********************************************************************
        !> @brief Routine prints first pixel information of segment
        !>
        !> @author Tatsiana Lapionak
        !> @date 15 DEC 2015
        !>

    subroutine print_crash_report()
        type(segment_data), pointer :: segment

        segment => crash_report%segment

        write(error_unit,'(a,i0,2(2x,a,f9.4),2(2x,a,i0),/)') &
          'id_inversion = ', segment%id_inversion,  &
          'coordinates for pixel # 1 : lon = ',segment%pixels(1)%x, &
          'lat = ',segment%pixels(1)%y,     &
          'irow = ',segment%pixels(1)%irow, &
          'icol = ',segment%pixels(1)%icol
    end subroutine print_crash_report

! **********************************************************************
        !> @brief Routine initializes stop report structure
        !> 
        !> @author Tatsiana Lapionak
        !> @date 15 DEC 2015
        !>

    subroutine initialize_stop_report()
      tmp_message = ''
      stop_report%status = .false.
      stop_report%message = ''
      stop_report%file = ''
      stop_report%line = 0

    end subroutine initialize_stop_report

! **********************************************************************
    subroutine print_stop_report()

        write(error_unit,'(a)') 'Execution has to be terminated.'
        call print_crash_report()

        !write(error_unit,'(2a,3x,a,i0)') 'Problem in file = ',trim(stop_report%file), &
        !'line = ',stop_report%line
        !write(error_unit,'(a)') trim(stop_report%message)

    end subroutine print_stop_report

! **********************************************************************
    logical function error_present()
      error_present = .false.
      if(stop_report%status .eqv. .true.) then
        error_present = .true.
      endif

    end function error_present

! **********************************************************************

end module mod_stop_report
