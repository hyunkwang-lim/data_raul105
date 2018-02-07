/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

#include "grasp_input_transformer_none.h"
#include "../../../input/grasp_input.h"



grasp_input_transformer_t grasp_input_transformer_none(){
    grasp_input_transformer_t x;
    
    x.init=grasp_input_transformer_none_init;
    x.function=grasp_input_transformer_none_function;
    x.close=grasp_input_transformer_none_close;
    
    return x;
}

int grasp_input_transformer_none_init(grasp_settings *settings, grasp_tile_description_t *input_information){
    return 0;
}

int grasp_input_transformer_none_function(grasp_settings *settings,grasp_segment_t *segment){
    return 0;
}

int grasp_input_transformer_none_close(void){
    return 0;
}

grasp_settings_parameter_array *grasp_input_transformer_settings_none(grasp_settings *settings){
    return NULL;
}