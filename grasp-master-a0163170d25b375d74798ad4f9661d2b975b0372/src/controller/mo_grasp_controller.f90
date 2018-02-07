! **
! **  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.  
! **  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
! **

module mo_grasp_controller
    use iso_fortran_env ! input_unit, output_unit, error_unit
    use iso_c_binding
    use mod_grasp
    use mod_retr_settings_derived_type
    !use mod_retr_general_output
    use mod_sdata
    use mod_alloc_arrays
    use mod_par_inv
    use mod_grasp
    use inversion_subsystem
    use mod_edges
    use mod_stop_report
    !private
    implicit none   

    type(kernels_triangle_bin)      ::  KERNELS1
    type(kernels_lognormal_bin)     ::  KERNELS2   
    real,dimension(:,:,:),allocatable ::  US

    contains
        

    subroutine grasp_prepare_segment_settings(RIN, segment_meas) bind(C)
      type(retr_input_settings), intent(inout) :: RIN
      type(segment_data), intent(inout) :: segment_meas
      
      RIN%main_output_file = "-"
      RIN%plotting_output_file = ""

      ! TEMPORAL SOLUTION: Retrieval code need two paths to kernels but control 
      ! unit reads just one. Temporally here it is duplicated the path read in settings
      ! Also we add a '/' symbol at the end. Now it is optional to put it in settings file
      RIN%DLSF%distname_O=trim(RIN%DLSF%distname_O)//'/'
      RIN%DLSF%distname_N=RIN%DLSF%distname_O
      
      call prepare_segment_settings(6, segment_meas, RIN )
    end subroutine grasp_prepare_segment_settings
    
    subroutine grasp_init_inversion(RIN) bind(C)
        type(retr_input_settings),         intent(inout)  :: RIN
        LOGICAL :: constant_errors = .FALSE.

        ! Validate constants
        if (KNBVM .ne. max(KVERTM,NBVM)) then
            write(*,*) 'Constant bad defined: KNBVM must be the maximum between KVERT and NBVM'
            constant_errors = .TRUE.
        endif
        if (KDIF .ne. max( KPARS, KITIME, KIX, KIY)) then
            write(*,*) 'Constant bad defined: KDIF must be the maximum between KVERT, KPARS, KITIME, KIX and KIY'
            constant_errors = .TRUE.
        endif
                
        if (constant_errors) then
            stop
        endif
                
        ! Init inversion                  
        call alloc_arrays(RIN,KERNELS1,KERNELS2,US)
        ! Initialize inversion subsystem
        call inversion_init(RIN)
    end subroutine grasp_init_inversion
    
    
    function grasp_input_inversion(RIN,sdata, iguess, edges, ROUT) bind(C)                              
        type(retr_input_settings),         intent(inout)  :: RIN
        type(segment_data),                intent(inout)  :: sdata
        real(C_FLOAT),                     intent(in)     :: iguess(KPARS,KIMAGE)
        type(segment_edges),               intent(in)     :: edges   
        type(output_segment_general),      intent(out)    :: ROUT
        integer                                           :: grasp_input_inversion
        
        grasp_input_inversion=0
        ! here is how to display a segment from the FORTRAN side (for debugging)
        call initialize_stop_report()
        call set_segment_for_stop_report(sdata)
        call grasp_prepare_segment_settings(RIN, sdata)
        if (error_present()) then
            call print_stop_report()
            grasp_input_inversion=-2
            return
        endif
        call inversion(RIN, 6, sdata, iguess, edges, ROUT, KERNELS1, KERNELS2, US)
        if (error_present()) then
            call print_stop_report()
            grasp_input_inversion=-1
            return
        endif
        return 
    end function grasp_input_inversion
    
    
    subroutine grasp_finalize_inversion(RIN) bind(C)
        type(retr_input_settings),     intent(in)  :: RIN
        ! Finalize inversion subsystem
        call inversion_finalize()
        ! dealloc kernels and Jacobian matrix
        call dealloc_arrays(RIN,KERNELS1,KERNELS2,US)
    end subroutine grasp_finalize_inversion
    

end module mo_grasp_controller


    

