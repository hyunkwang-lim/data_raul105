# Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
# Licensed under the GRASP Open Source License V1.0 (see LICENSE file)

add_library(output_current_function_none grasp_output_current_none.c)
target_link_libraries(output_current_function_none grasp_settings)

