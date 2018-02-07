# Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
# Licensed under the GRASP Open Source License V1.0 (see LICENSE file)

add_library(input_transformer_none grasp_input_transformer_none.c)
target_link_libraries(input_transformer_none grasp_settings)
