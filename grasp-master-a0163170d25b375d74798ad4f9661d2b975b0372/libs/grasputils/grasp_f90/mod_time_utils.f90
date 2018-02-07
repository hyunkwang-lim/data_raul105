!@Author: Fabrice Ducos <fabrice.ducos@univ-lille1.fr>


module mod_time_utils
  ! Author: Fabrice Ducos <fabrice.ducos@univ-lille1.fr>
  ! Date  : 29. november 2013
  
  use iso_c_binding
  use mod_c_utils
  implicit none
  private
  
#ifndef SIZEOF_TIME_T
#error "The size of time_t is not known on your system, I can't build a reliable version of the library."
#endif

  integer, parameter :: KIND_TIME = c_int64_t

#if (SIZEOF_TIME_T == 32)
  integer, parameter :: KIND_TIME_T = c_int32_t
#elif (SIZEOF_TIME_T == 64)
  integer, parameter :: KIND_TIME_T = c_int64_t
#else
#error "Your platform uses a definition of time_t that is not supported. Please contact the maintainer <fabrice.ducos@univ-lille1.fr>"
#endif

  integer(kind=KIND_TIME), parameter :: TIME_CONV_ERROR = -1_KIND_TIME

  public :: KIND_TIME
  public :: TIME_CONV_ERROR
  public :: convert_string_to_time
  public :: convert_time_to_string

  interface
     integer(c_int) function convert_string_to_time_c(cstr_time, ctime, time_format) bind(c, name='convert_string_to_time')
       use iso_c_binding
       import KIND_TIME_T
       character(kind=c_char, len=1), intent(in) :: cstr_time
       integer(KIND_TIME_T), intent(out) :: ctime
       integer(c_int), intent(in), value :: time_format
     end function convert_string_to_time_c

     type(c_ptr) function time_to_string_c(ctime, time_format) bind(c, name='time_to_string')
       use iso_c_binding
       import KIND_TIME_T
       integer(KIND_TIME_T), intent(in), value :: ctime
       character(kind=c_char, len=1), intent(in) :: time_format
     end function time_to_string_c

     type(c_ptr) function time_to_string_r_c(ctime, time_format, cstr_size, cstr_time) bind(c, name='time_to_string_r')
       use iso_c_binding
       import KIND_TIME_T
       integer(KIND_TIME_T), intent(in), value :: ctime
       character(kind=c_char, len=1), intent(in) :: time_format
       integer(kind=c_size_t), intent(in), value :: cstr_size
       character(kind=c_char,len=1), intent(in) :: cstr_time
     end function time_to_string_r_c

  end interface

  contains

    ! version of convert_time_to_string based on a non thread-safe version
    ! of time_to_string_c (deprecated and not published anymore)
    subroutine convert_time_to_string_static(time, time_format, str_time, success)
      integer(KIND_TIME) :: time
      character(*), intent(in)  :: time_format
      character(*), intent(out) :: str_time

      logical, intent(out), optional :: success
      type(c_ptr) :: cstr_time
      character(len(str_time)), pointer :: str_time_ptr
      integer(KIND_TIME_T) :: ctime
      
      if (time < 0) then
         if (present(success)) success = .false.
         str_time = ''
         return
      end if
      ctime = int(time, kind=KIND_TIME_T)
      
      cstr_time = time_to_string_c(ctime, time_format)
      call c_f_pointer(cstr_time, str_time_ptr)
      
      str_time = str_time_ptr

      if (len_trim(str_time) > 0) then
         if (present(success)) success = .true.
      else
         if (present(success)) success = .false.
      end if

    end subroutine convert_time_to_string_static

    subroutine convert_time_to_string(time, time_format, str_time, success)
      integer(KIND_TIME) :: time
      character(*), intent(in)  :: time_format
      character(*), intent(out) :: str_time

      logical, intent(out), optional :: success
      type(c_ptr) :: cstr_time_ptr
      integer(kind=c_size_t) :: cstr_size
      character(len(str_time)), pointer :: str_time_ptr
      character(kind=c_char, len=64) :: cstr_time
      integer(KIND_TIME_T) :: ctime
 
      if (time < 0) then
         if (present(success)) success = .false.
         str_time = ''
         return
      end if
      ctime = int(time, kind=KIND_TIME_T)

      cstr_size = len(str_time) + 1
      call fstring2cstring(str_time, cstr_size, cstr_time)

      cstr_time_ptr = time_to_string_r_c(ctime, time_format, cstr_size, cstr_time)
      call c_f_pointer(cstr_time_ptr, str_time_ptr)
      
      str_time = str_time_ptr

      if (len_trim(str_time) > 0) then
         if (present(success)) success = .true.
      else
         if (present(success)) success = .false.
      end if

    end subroutine convert_time_to_string



    function convert_string_to_time(str_time, success) result(time)
      character(*), intent(in) :: str_time
      logical, intent(out), optional :: success
      integer(KIND_TIME) :: time

      integer :: retval_c
      integer, parameter :: TIMEFMT_ISO8601 = 0
      character(kind=c_char, len=64) :: cstr_time
      integer(kind=c_size_t) :: cstr_size
      integer(KIND_TIME_T) :: ctime
      
      if (len(str_time) == 0) then
#ifdef DEBUG
        write(error_unit, '("convert_string_to_time(""",A,"""): &
	  &empty string was given")') trim(str_time)
#endif
        time = TIME_CONV_ERROR
        if (present(success)) success = .false.
        return
      end if

      cstr_size = len(cstr_time)
      call fstring2cstring(str_time, cstr_size, cstr_time, .true.)

      retval_c = convert_string_to_time_c(cstr_time, ctime, TIMEFMT_ISO8601)
      if (retval_c /= 0) then
#ifdef DEBUG
        write(error_unit, '("convert_string_to_time(""",A,"""): &
	  &failed in convert_string_to_time_c")') str_time
#endif
        time = TIME_CONV_ERROR
        if (present(success)) success = .false.
        return
      end if

      if (ctime < 0) then
         time = TIME_CONV_ERROR
         success = .false.
      end if

      time = ctime
      if (present(success)) success = .true.
    end function convert_string_to_time


end module mod_time_utils

#ifdef UNIT_TEST

program test
  use iso_c_binding
  use mod_time_utils
  implicit none
  
  character(128) :: a_long_string
  integer, parameter :: sizeof_time_t_ = SIZEOF_TIME_T
  logical :: retval
  character(*), parameter :: time_format = '%FT%H:%M:%SZ' ! strftime specification for ISO8601 (ref: strftime man page)
  
  a_long_string = repeat("0123456789ABCDEF", 8)

  retval = .true.
  retval = retval .and. test_string_to_time('1970-01-01T00:00:00Z', 0_KIND_TIME)
  retval = retval .and. test_string_to_time('1970-01-01T00:00:01Z', 1_KIND_TIME)
  retval = retval .and. test_string_to_time('2000-01-01T00:00:00Z', 946684800_KIND_TIME)
  retval = retval .and. test_string_to_time('1969-12-31T23:59:59Z', TIME_CONV_ERROR, expected_to_fail = .true.)
  retval = retval .and. test_string_to_time('XXX', TIME_CONV_ERROR, expected_to_fail = .true.)
  retval = retval .and. test_string_to_time(a_long_string, TIME_CONV_ERROR, expected_to_fail = .true.)
  retval = retval .and. test_string_to_time(' ', TIME_CONV_ERROR, expected_to_fail = .true.)
  retval = retval .and. test_string_to_time('', TIME_CONV_ERROR, expected_to_fail = .true.)
  
  retval = retval .and. test_time_to_string(0_KIND_TIME, time_format)
  retval = retval .and. test_time_to_string(1_KIND_TIME, time_format)
  retval = retval .and. test_time_to_string(946684800_KIND_TIME, time_format)

  if (retval .eqv. .true.) then
     stop 0
  else
     write(*,'("Something unexpected occured. sizeof(time_t) was assumed to be ",I0," bits. It appears to be wrong. &
          &Please run the get_sizeof_time_t command to find the size of time_t in bits.")') sizeof_time_t_
     stop 1
  end if

contains

  function test_time_to_string(time, time_format, expected_to_fail) result(retval)
    integer(KIND_TIME), intent(in) :: time
    character(*), intent(in) :: time_format
    logical, intent(in), optional :: expected_to_fail

    character(20) :: str_time
    logical :: success
    logical :: retval
    logical :: expected_to_fail_

    if (present(expected_to_fail)) then
       expected_to_fail_ = expected_to_fail
    else
       expected_to_fail_ = .false.
    end if

    call convert_time_to_string(time, time_format, str_time, success)
    
    if (.not. success) then
       if (expected_to_fail_ .eqv. .true.) then
          write(*,'("convert_time_to_string(",I0,",",""",A,""",") failed, but this was expected")') &
               time, time_format
          retval = .true.
       else
          write(*,'("convert_time_to_string(",I0,",",""",A,""",") failed")') &
               time, time_format
          retval = .false.
       end if
       return
    end if

    write(*,*) str_time
    retval = .true.
    ! TODO: validate the output strings with an assert (validated by visual control for the moment)

  end function test_time_to_string

  function test_string_to_time(date, control_value, expected_to_fail) result(retval)
    character(*), intent(in) :: date
    integer(KIND_TIME), intent(in) :: control_value
    logical, intent(in), optional :: expected_to_fail
    
    integer(KIND_TIME) :: time
    logical :: success
    logical :: retval
    logical :: expected_to_fail_

    if (present(expected_to_fail)) then
       expected_to_fail_ = expected_to_fail
    else
       expected_to_fail_ = .false.
    end if

    time = convert_string_to_time(date, success)
    
    if (.not. success) then
       if (expected_to_fail_ .eqv. .true.) then
          write(*,'("convert_string_to_time(""",A,""") failed, but this was expected")') date
          retval = .true.
       else
          write(*,'("convert_string_to_time(""",A,""") failed")') date
          retval = .false.
       end if
       return
    end if

    retval = assert_equals(trim(date), time, control_value)
    
  end function test_string_to_time

  ! not a substitute for a real unit-test library, but allows some
  ! quick tests without depending on any external library
  function assert_equals(expression_to_evaluate, evaluated_expression, &
    control_value, expected_to_fail) result(retval)
    character(*), intent(in) :: expression_to_evaluate
    integer(KIND_TIME), intent(in) :: evaluated_expression
    integer(KIND_TIME), intent(in) :: control_value
    logical, intent(in), optional  :: expected_to_fail
    logical :: retval

    logical :: expected_to_fail_
    
    if (present(expected_to_fail)) then
       expected_to_fail_ = expected_to_fail
    else
       expected_to_fail_ = .false.
    end if

    if (evaluated_expression == control_value) then
       if (.not. expected_to_fail_) then
          write(*,'(A,": ",I0," ... test passed")') expression_to_evaluate, evaluated_expression
          retval = .true.
       else
          write(*,'(A,": ",I0," ... test passed, but this was not expected !")') &
               expression_to_evaluate, evaluated_expression
          retval = .false.
       end if
    else
       if (.not. expected_to_fail_) then
          write(*,'(A,": ",I0," ... test failed (",I0," was expected)")') &
               expression_to_evaluate, evaluated_expression, control_value
          retval = .false.
       else
          write(*,'(A,": ",I0," ... test failed, but this was expected")') &
               expression_to_evaluate, evaluated_expression, control_value
          retval = .true.
       end if
    end if
  end function assert_equals

end program test

#endif

