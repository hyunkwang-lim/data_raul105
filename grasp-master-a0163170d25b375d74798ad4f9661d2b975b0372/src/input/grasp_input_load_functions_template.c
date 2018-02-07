/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <grasp/utils.h>
#include "grasp_input_load_functions.h"
#include "yamlsettings/yamlsettings_dictionary.h"
// Driver list

@GRASP_INPUT_LOAD_DRIVER_INCLUDE@

// Find a driver and return it.
grasp_input_driver_t grasp_input_driver_get_function(char *name){
    grasp_input_driver_t none;

    if(0){
	@GRASP_INPUT_LOAD_DRIVER_FUNCTION@
    }else{
        printf("ERROR: Unknown driver\n");
        exit(-1);
    }        

    return none;
}

// Load settings
grasp_settings_parameter_array **grasp_input_driver_get_settings_parameters(grasp_settings *settings){

    grasp_settings_parameter_array **result;
    result = (grasp_settings_parameter_array **) trackmem_malloc(sizeof (grasp_settings_parameter_array *)*@GRASP_INPUT_LOAD_DRIVER_NDRIVERS@);
	@GRASP_INPUT_LOAD_DRIVER_PARAMETER@

    return result;
}

int grasp_input_ndrivers(){
    return @GRASP_INPUT_LOAD_DRIVER_NDRIVERS@;
}

const char *grasp_input_driver_get_name(int n){
    assert(n>=0 && n<grasp_input_ndrivers());
    
    switch(n){
@GRASP_INPUT_LOAD_DRIVER_NAMES@            default: return NULL;
    }
}

// Transformer list

@GRASP_INPUT_LOAD_TRANSFORMER_INCLUDE@


// Find a transformer and return it.
grasp_input_transformer_t grasp_input_transformer_get_function(char *name){
    grasp_input_transformer_t none;

    if(0){
      @GRASP_INPUT_LOAD_TRANSFORMER_FUNCTION@
    }else{
        printf("ERROR: Unknown transformer %s\n", name);
        exit(-1);
    }        

    return none;
} 

// Load transformer settings
grasp_settings_parameter_array **grasp_input_transformer_get_settings_parameters(grasp_settings *settings){

    grasp_settings_parameter_array **result;
    result = (grasp_settings_parameter_array **) trackmem_malloc(sizeof (grasp_settings_parameter_array *)*@GRASP_INPUT_LOAD_TRANSFORMER_NTRANSFORMERS@);
    @GRASP_INPUT_LOAD_TRANSFORMER_PARAMETER@

    return result;
}

int grasp_input_ntransformers(){
    return @GRASP_INPUT_LOAD_TRANSFORMER_NTRANSFORMERS@;
}

const char *grasp_input_transformer_get_name(int n){
    assert(n>=0 && n<grasp_input_ntransformers());
    
    switch(n){
@GRASP_INPUT_LOAD_TRANSFORMER_NAMES@            default: return NULL;
    }
}



