# Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
# Licensed under the GRASP Open Source License V1.0 (see LICENSE file)

include_directories(${INPUT_SRC} ${INPUT_BIN} ${OUTPUT_SRC} ${OUTPUT_BIN})

add_library(grasp_settings
	grasp_settings.c
	grasp_settings_assign_data.c
	grasp_settings_validators.c
	grasp_settings_data_types.c
)


target_link_libraries(grasp_settings yamlsettings)
