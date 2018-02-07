# Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
# Licensed under the GRASP Open Source License V1.0 (see LICENSE file)


add_library(grasp_input_util
  grasp_input_preloader.c
  grasp_driver_settings.c
  grasp_box.c
  grasp_reindexer3D.c
)
target_link_libraries(grasp_input_util
  grasp_global
  grasp_settings
)