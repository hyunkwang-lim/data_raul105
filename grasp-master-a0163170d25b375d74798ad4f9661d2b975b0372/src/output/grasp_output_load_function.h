/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

#include "grasp_output.h"
#include "yamlsettings/yamlsettings_dictionary.h"
#include "../input/grasp_input.h"

grasp_output_segment_function_t grasp_output_get_segment_output_funtion(char segment_output_function_name[_GBL_FILE_PATH_LEN]);

grasp_output_tile_function_t grasp_output_get_tile_output_funtion(char tile_output_function_name[_GBL_FILE_PATH_LEN]);

grasp_output_current_function_t grasp_output_get_current_output_funtion(char current_output_function_name[_GBL_FILE_PATH_LEN]);

grasp_settings_parameter_array **grasp_output_segment_functions_get_settings_parameters(grasp_settings *settings);

grasp_settings_parameter_array **grasp_output_tile_functions_get_settings_parameters(grasp_settings *settings);

grasp_settings_parameter_array **grasp_output_current_functions_get_settings_parameters(grasp_settings *settings);

int grasp_output_segment_nfunctions();

int grasp_output_tile_nfunctions();

int grasp_output_current_nfunctions();

const char *grasp_output_segment_function_get_name(int n);

const char *grasp_output_tile_function_get_name(int n);

const char *grasp_output_current_function_get_name(int n);
