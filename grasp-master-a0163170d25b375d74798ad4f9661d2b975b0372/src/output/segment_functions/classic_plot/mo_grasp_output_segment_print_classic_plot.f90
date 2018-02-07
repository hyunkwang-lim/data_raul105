! **
! **  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.  
! **  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
! **

module grasp_output_segment_print_classic_plot
    use iso_c_binding
!    use mo_retr_input
    use mod_sdata
    use mod_retr_general_output_derived_type
    use mod_print_array
    !use mod_grasp

    !private
    implicit none   

    contains

    subroutine grasp_output_segment_print_classic_plot_f90(RIN,sdata, ROUT) bind(C)                              
        type(retr_input_settings), intent(in)       ::  RIN
        type(segment_data),        intent(in)       ::  sdata
        type(output_segment_general), intent(in)  ::  ROUT

        call print_classic_plot(RIN, sdata, ROUT)
    end subroutine grasp_output_segment_print_classic_plot_f90


end module grasp_output_segment_print_classic_plot

