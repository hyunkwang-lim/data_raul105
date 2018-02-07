! **
! **  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.  
! **  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
! **

module mo_grasp_output_segment_print_classic
    use iso_c_binding
    use mod_retr_settings_derived_type
    use mod_sdata
    use mod_retr_general_output_derived_type
    use mod_print_array
    use mod_globals, only  : GBL_FILE_PATH_LEN
    use mod_c_utils

    !private
    implicit none   

    contains

    subroutine grasp_output_segment_print_classic_f90_screen(RIN,sdata, ROUT) bind(C)                              
        type(retr_input_settings), intent(in)       ::  RIN
        type(segment_data),        intent(in)       ::  sdata
        type(output_segment_general), intent(in)  ::  ROUT

        call print_output_results(6,RIN,sdata, ROUT)
    end subroutine grasp_output_segment_print_classic_f90_screen

    subroutine grasp_output_segment_print_classic_f90_file(filename,RIN,sdata, ROUT) bind(C)  
        character(kind=C_CHAR)  ::  filename(GBL_FILE_PATH_LEN)
        type(retr_input_settings), intent(in)       ::  RIN
        type(segment_data),        intent(in)       ::  sdata
        type(output_segment_general), intent(in)  ::  ROUT
        
        integer :: out_unit = 20
        character (len=GBL_FILE_PATH_LEN) :: fortranfilename
    
        call cstring2fstring(filename,fortranfilename)
        
        open (unit=out_unit,file=trim(fortranfilename),action="write",status="replace")
        
        call print_output_results(out_unit,RIN,sdata, ROUT)
        
        close (out_unit)
    end subroutine grasp_output_segment_print_classic_f90_file
 

end module mo_grasp_output_segment_print_classic

