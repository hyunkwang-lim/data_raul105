/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

/**
 * @file   grasp_settings_validators.h
 * @author David Fuertes (david.fuertes@univ-lille1.fr)
 * @date   August, 2013
 * @brief  Validator functions used in dictionary
 *
 * This functions read yamls files in a file and return a GTree
 */

#ifndef GRASP_SETTINGS_VALIDATORS_H
#define	GRASP_SETTINGS_VALIDATORS_H

#ifdef	__cplusplus
extern "C" {
#endif


#include "yamlsettings/yamlsettings_dictionary.h"


// The parameter is a directory and it must be exist (stored like fortran string)
int grasp_settings_validator_directory_fortran(int param_index,yamlsettings_dictionary_t *dictionary, char arguments[YAMLSETTINGS_VALIDATOR_MAX_ARGUMENTS][YAMLSETTINGS_VALIDATOR_MAX_ARGUMENT_SIZE]);

// This function call to grasp_output_stream_filename_validation to know if the stream is well formed
int grasp_settings_validator_stream(int param_index,yamlsettings_dictionary_t *dictionary, char arguments[YAMLSETTINGS_VALIDATOR_MAX_ARGUMENTS][YAMLSETTINGS_VALIDATOR_MAX_ARGUMENT_SIZE]);

// This validator checks if the input driver is compiled within current framework compilation
int grasp_settings_validator_input_driver(int param_index,yamlsettings_dictionary_t *dictionary, char arguments[YAMLSETTINGS_VALIDATOR_MAX_ARGUMENTS][YAMLSETTINGS_VALIDATOR_MAX_ARGUMENT_SIZE]);

// This validator checks if the input transformer name is available in the data transformers compiled with the system.
int grasp_settings_validator_input_transformer(int param_index,yamlsettings_dictionary_t *dictionary, char arguments[YAMLSETTINGS_VALIDATOR_MAX_ARGUMENTS][YAMLSETTINGS_VALIDATOR_MAX_ARGUMENT_SIZE]);

// This validator checks if output segment function is compiled within current framework compilation
int grasp_settings_validator_output_segment_function(int param_index,yamlsettings_dictionary_t *dictionary, char arguments[YAMLSETTINGS_VALIDATOR_MAX_ARGUMENTS][YAMLSETTINGS_VALIDATOR_MAX_ARGUMENT_SIZE]);

// This validator  checks that bins have to be defined only in one way, or bin bins or by min and max
int graspsettings_validator_bins (int param_index,yamlsettings_dictionary_t *dictionary, char arguments[YAMLSETTINGS_VALIDATOR_MAX_ARGUMENTS][YAMLSETTINGS_VALIDATOR_MAX_ARGUMENT_SIZE]);

// This validator check if output tile function is compiled within current framework compilation
int grasp_settings_validator_output_tile_function(int param_index,yamlsettings_dictionary_t *dictionary, char arguments[YAMLSETTINGS_VALIDATOR_MAX_ARGUMENTS][YAMLSETTINGS_VALIDATOR_MAX_ARGUMENT_SIZE]);

// This validator checks if output current function is compiled within current framework compilation
int grasp_settings_validator_output_current_function(int param_index,yamlsettings_dictionary_t *dictionary, char arguments[YAMLSETTINGS_VALIDATOR_MAX_ARGUMENTS][YAMLSETTINGS_VALIDATOR_MAX_ARGUMENT_SIZE]);

// This validator checks that all not retrieved characteristics are at the end of NDIM arrays (iguess).
int graspsettings_validator_characteristic_retrieved(int param_index,yamlsettings_dictionary_t *dictionary, char arguments[YAMLSETTINGS_VALIDATOR_MAX_ARGUMENTS][YAMLSETTINGS_VALIDATOR_MAX_ARGUMENT_SIZE]);

// This validator checks mandatory parameters depending on characteristic types defined
int graspsettings_validator_characteristic_type(int param_index,yamlsettings_dictionary_t *dictionary, char arguments[YAMLSETTINGS_VALIDATOR_MAX_ARGUMENTS][YAMLSETTINGS_VALIDATOR_MAX_ARGUMENT_SIZE]);

// This validator checks if the use of models is valid for current configuration
int graspsettings_validator_models(int param_index,yamlsettings_dictionary_t *dictionary, char arguments[YAMLSETTINGS_VALIDATOR_MAX_ARGUMENTS][YAMLSETTINGS_VALIDATOR_MAX_ARGUMENT_SIZE]);

// This validator checks if two single (not an array) parameters readed are divisibles (mod=0)
int graspsettings_validator_divisible(int param_index,yamlsettings_dictionary_t *dictionary, char arguments[YAMLSETTINGS_VALIDATOR_MAX_ARGUMENTS][YAMLSETTINGS_VALIDATOR_MAX_ARGUMENT_SIZE]);

// It check if a wevalength index is defined between 1 and the bigest of wavelength indexes
int graspsettings_validator_indexes_of_wavelengths(int param_index,yamlsettings_dictionary_t *dictionary, char arguments[YAMLSETTINGS_VALIDATOR_MAX_ARGUMENTS][YAMLSETTINGS_VALIDATOR_MAX_ARGUMENT_SIZE]);

// The parameter has the samennumber of elements that other element or 0 (don't defined). Argument 1=name of second element 
int graspsettings_validator_same_nelements_or_zero(int param_index,yamlsettings_dictionary_t *dictionary, char arguments[YAMLSETTINGS_VALIDATOR_MAX_ARGUMENTS][YAMLSETTINGS_VALIDATOR_MAX_ARGUMENT_SIZE]);

// To validate if kernels_folder is a folder which exists
int graspsettings_validator_kernelpath(int param_index,yamlsettings_dictionary_t *dictionary, char arguments[YAMLSETTINGS_VALIDATOR_MAX_ARGUMENTS][YAMLSETTINGS_VALIDATOR_MAX_ARGUMENT_SIZE]);

#ifdef	__cplusplus
}
#endif

#endif	/* GRASP_SETTINGS_VALIDATORS_H */
