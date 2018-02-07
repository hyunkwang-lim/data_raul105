# Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
# Licensed under the GRASP Open Source License V1.0 (see LICENSE file)

add_library(output_tile_function_ascii
	grasp_output_tile_ascii.c
)
target_link_libraries(output_tile_function_ascii grasp_settings grasp_output)
