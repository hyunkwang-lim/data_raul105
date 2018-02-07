! **
! **  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.  
! **  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
! **


module mod_c_utils

  use iso_fortran_env
  use iso_c_binding
  implicit none

  ! interface to some C standard routines (shall be extended later)
  interface
     function strlen(cstring) bind(c)
       use iso_c_binding
       character(kind=C_CHAR),  intent(in) :: cstring(*)
       integer(kind=C_SIZE_T) :: strlen
     end function strlen
  end interface

  contains
    ! Copyright (C) 2012, 2013 Fabrice Ducos, fabrice.ducos@univ-lille1.fr
    ! F. Ducos, 30 August 2013: fstring is no longer allocatable in cstring2fstring    
    ! a routine to convert cstrings (char *) into fortran strings
    subroutine cstring2fstring(cstring, fstring)
      character(kind=C_CHAR),  intent(in) :: cstring(*)
      character(kind=C_CHAR,len=*), intent(out) :: fstring

      integer(kind=C_SIZE_T) :: length
      integer(kind=C_SIZE_T) :: i

      length = strlen(cstring)
      
      do i = 1, length
         fstring(i:i) = cstring(i)
      end do

      do i = length + 1, len(fstring)
         fstring(i:i) = ' '
      end do
      
    end subroutine cstring2fstring
    
    subroutine fstring2cstring(fstring,cstring_size, cstring, keep_trailing_spaces)
      character(kind=C_CHAR, len=*), intent(in) :: fstring
      integer(kind=C_SIZE_T), intent(in), value :: cstring_size ! includes the terminating NULL char
      character(kind=C_CHAR), intent(out) :: cstring(*)
      logical, intent(in), optional :: keep_trailing_spaces

      integer(kind=C_SIZE_T) :: i
      integer(kind=C_SIZE_T) :: fstring_len
      integer(kind=C_SIZE_T) :: min_len
      logical :: keep_trailing_spaces_
      
      if (.not. present(keep_trailing_spaces)) then
         keep_trailing_spaces_ = .false.
      else
         keep_trailing_spaces_ = keep_trailing_spaces
      end if

      if (keep_trailing_spaces_ .eqv. .false.) then
         fstring_len = len_trim(fstring)
      else
         fstring_len = len(fstring)
      end if

      min_len=min(fstring_len, cstring_size - 1)

      do i = 1, min_len
         cstring(i) = fstring(i:i)
      end do
      
      cstring(min_len + 1) = C_NULL_CHAR
      
    end subroutine fstring2cstring    


end module mod_c_utils
