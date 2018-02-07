! **
! **  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.  
! **  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
! **

#include "constants_set/mod_globals.inc" 
module mod_globals

  implicit none
  
  integer, parameter :: GBL_FILE_PATH_LEN = _GBL_FILE_PATH_LEN
  character(*), parameter :: GBL_APPNAME = 'grasp_retrieval'
  
  ! for data shared at the application scope
  ! (e.g. user options) ;
  ! should not be abused to create unnecessary
  ! dependencies. Mainly intended for debugging
  ! purposes.
  
  ! only the main program should have
  ! writing access to this structure
  
  type globals_type
    logical      :: debug_mode = .false.
	! other flags or user settings should be put here later
  end type globals_type
  
  type(globals_type), save :: globals

end module mod_globals
