/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <grasp/utils.h>
#include "grasp_output_load_function.h"
#include "yamlsettings/yamlsettings_dictionary.h"
// Segment function list
@GRASP_OUTPUT_LOAD_FUNCTION_SEGMENT_INCLUDE@

// Tile function list
@GRASP_OUTPUT_LOAD_FUNCTION_TILE_INCLUDE@

// Current function list
@GRASP_OUTPUT_LOAD_FUNCTION_CURRENT_INCLUDE@


// Find a function and return it.
grasp_output_segment_function_t grasp_output_get_segment_output_funtion(char *name){
    grasp_output_segment_function_t none;

    if(0){
	@GRASP_OUTPUT_LOAD_FUNCTION_SEGMENT@
    }else{
        printf("ERROR: Unknown output segment function: %s\n", name);
        exit(-1);
    }        

    none.init = NULL;
    none.function = NULL;
    none.close = NULL;

    return none;
}

grasp_output_tile_function_t grasp_output_get_tile_output_funtion(char *name){
    grasp_output_tile_function_t none;

    if(0){
	@GRASP_OUTPUT_LOAD_FUNCTION_TILE@
    }else{
        printf("ERROR: Unknown output tile function: %s\n", name);
        exit(-1);
    }        


    none.init = NULL;
    none.function = NULL;
    none.close = NULL;

    return none;
}

grasp_output_current_function_t grasp_output_get_current_output_funtion(char *name){
    grasp_output_current_function_t none;

    if(0){
	@GRASP_OUTPUT_LOAD_FUNCTION_CURRENT@
    }else{
        printf("ERROR: Unknown output current function: %s\n", name);
        exit(-1);
    }        

    none.init = NULL;
    none.function = NULL;
    none.close = NULL;

    return none;
}



// Load settings
grasp_settings_parameter_array **grasp_output_segment_functions_get_settings_parameters(grasp_settings *settings){

    grasp_settings_parameter_array **result;
    result = (grasp_settings_parameter_array **) trackmem_malloc(sizeof (grasp_settings_parameter_array *)*@GRASP_OUTPUT_LOAD_FUNCTION_SEGMENT_SIZE@);
    @GRASP_OUTPUT_LOAD_FUNCTION_SEGMENT_PARAMETERS@

    return result;
}

grasp_settings_parameter_array **grasp_output_tile_functions_get_settings_parameters(grasp_settings *settings){

    grasp_settings_parameter_array **result;
    result = (grasp_settings_parameter_array **) trackmem_malloc(sizeof (grasp_settings_parameter_array *)*@GRASP_OUTPUT_LOAD_FUNCTION_TILE_SIZE@);
    @GRASP_OUTPUT_LOAD_FUNCTION_TILE_PARAMETERS@

    return result;
}

grasp_settings_parameter_array **grasp_output_current_functions_get_settings_parameters(grasp_settings *settings){

    grasp_settings_parameter_array **result;
    result = (grasp_settings_parameter_array **) trackmem_malloc(sizeof (grasp_settings_parameter_array *)*@GRASP_OUTPUT_LOAD_FUNCTION_CURRENT_SIZE@);
    @GRASP_OUTPUT_LOAD_FUNCTION_CURRENT_PARAMETERS@

    return result;
}


int grasp_output_segment_nfunctions(){
    return @GRASP_OUTPUT_LOAD_FUNCTION_SEGMENT_SIZE@;
}

int grasp_output_tile_nfunctions(){
    return @GRASP_OUTPUT_LOAD_FUNCTION_TILE_SIZE@;
}

int grasp_output_current_nfunctions(){
    return @GRASP_OUTPUT_LOAD_FUNCTION_CURRENT_SIZE@;
}

const char *grasp_output_segment_function_get_name(int n){
    assert(n>=0 && n<grasp_output_segment_nfunctions());
    
    switch(n){
@GRASP_OUTPUT_LOAD_FUNCTION_SEGMENT_NAMES@            default: return NULL;
    }
}

const char *grasp_output_tile_function_get_name(int n){
    assert(n>=0 && n<grasp_output_tile_nfunctions());
    
    switch(n){
@GRASP_OUTPUT_LOAD_FUNCTION_TILE_NAMES@            default: return NULL;
    }
}

const char *grasp_output_current_function_get_name(int n){
    assert(n>=0 && n<grasp_output_current_nfunctions());
    
    switch(n){
@GRASP_OUTPUT_LOAD_FUNCTION_CURRENT_NAMES@            default: return NULL;
    }
}