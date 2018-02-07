# Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
# Licensed under the GRASP Open Source License V1.0 (see LICENSE file)

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp")

add_library(output_segment_function_classic_plot
	mo_grasp_output_segment_print_classic_plot.f90
	grasp_output_segment_classic_plot.c
)
target_link_libraries(output_segment_function_classic_plot grasp_retrieval_lib grasp_settings)
