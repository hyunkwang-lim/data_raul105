/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

/**
 * @file   grasp_settings.h
 * @author David Fuertes (david.fuertes@univ-lille1.fr)
 * @date   August, 2013
 * @brief  Main interface binded for input settings 
 *
 * Here is defined option parameters and main functions for retrieved it
 * from a configuration file
 * 
 */
#ifndef _GRASP_SETTINGS_H
#define _GRASP_SETTINGS_H


#ifdef	__cplusplus
extern "C" {
#endif
    
#include <stdbool.h>
#include "yamlsettings/yamlsettings_dictionary.h"
#include "mod_par_OS.inc"
#include "mod_par_inv.inc"
#include "../retrieval/constants_set/mod_globals.inc"
#include "../input/grasp_input_settings.h"
#include "../controller/grasp_controller_settings.h"
#include "../output/grasp_output_settings.h"    
#include "../settings/grasp_products.h"
#include "../global/grasp_parameters.h"
    
#ifndef RESOURCES_PREFIX
#define RESOURCES_PREFIX "/usr/local/grasp-resources/"
#endif
    
#include "grasp_settings_t.h"
#include "../output/grasp_output_stream_t.h"

typedef struct grasp_settings_parameter_array_{
    int nparameters;
    yamlsettings_parameter *parameters;
} grasp_settings_parameter_array;

// Return a dictionary assig
yamlsettings_dictionary_t *grasp_settings_dictionary_get(grasp_settings *settings);


// Function that read  Retrieval INput structure (RIN) from a Yaml file. Arguments:
// grasp_settings **settings pointer to settings that it will be allocated and returned
// inputFile filename of input YML configuration
// Return -1 if file can't be open or is a malformed YAML
// Return -2 if Reached maximum of imported files. This can occur with a recursive importation
// Return  0 if everything is ok
// Return a number greater than 0 thats means the number of errors when data is validated
// nparameter is the number of extra parameters
// parameters is an array of parameters in string format (example: retrieval.number_wavelengths=3)
int grasp_settings_read(grasp_output_stream *controller_stream, yamlsettings_dictionary_t **dictionary, int nparameters, char const *parameters[]);

// Return an array of parameter allocated
grasp_settings_parameter_array *grasp_settings_parameter_array_allocate(int nelements);

// Deallocate current x parameter array
void grasp_settings_parameter_array_destroy(grasp_settings_parameter_array *x);

// This function will print help information (parameter name + description) in a defined stream
void grasp_settings_help(grasp_output_stream *controller_stream, const char *filter);

// Free memory from settings
void grasp_settings_destroy(grasp_settings *settings);

// Function which is called after process settings. This function will call other functions in order to set properly all settings structure
void grasp_settings_postprocess_function(yamlsettings_dictionary_t *dictionary);


#ifdef	__cplusplus
}
#endif

#endif /*_GRASP_SETTINGS_H */
