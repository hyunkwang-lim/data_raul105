# Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
# Licensed under the GRASP Open Source License V1.0 (see LICENSE file)

add_library(output_segment_function_none
	grasp_output_segment_none.c
	grasp_output_segment_none.h
)
target_link_libraries(output_segment_function_none grasp_settings)

