/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

#include "../input/grasp_input.h"
#include "yamlsettings/yamlsettings_dictionary.h"


// For drivers
grasp_input_driver_t grasp_input_driver_get_function(char *name);

grasp_settings_parameter_array **grasp_input_driver_get_settings_parameters(grasp_settings *settings);

int grasp_input_ndrivers();

const char *grasp_input_driver_get_name(int n);

// For transformers
grasp_input_transformer_t grasp_input_transformer_get_function(char *name);

grasp_settings_parameter_array **grasp_input_transformer_get_settings_parameters(grasp_settings *settings);

int grasp_input_ntransformers();

const char *grasp_input_transformer_get_name(int n);
