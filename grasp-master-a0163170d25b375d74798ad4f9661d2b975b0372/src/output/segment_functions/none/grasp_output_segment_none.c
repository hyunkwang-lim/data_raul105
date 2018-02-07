/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

#include "grasp_output_segment_none.h"


grasp_output_segment_function_t grasp_output_segment_function_none(){
    grasp_output_segment_function_t x;
    
    x.init=grasp_output_segment_function_none_init;
    x.function=grasp_output_segment_function_none_process;
    x.close=grasp_output_segment_function_none_close;
    
    return x;
}

int grasp_output_segment_function_none_init(grasp_settings *settings, grasp_tile_description_t *input_information){
    return 0;
}

int grasp_output_segment_function_none_close(void){
    return 0;
}

int grasp_output_segment_function_none_process(grasp_output_stream *stream, grasp_settings *settings, grasp_segment_t *segment, output_segment_general *output, grasp_tile_description_t *tile_description,int icol,int irow,int itime){
    return 0;
}

grasp_settings_parameter_array *grasp_output_segment_function_settings_none(grasp_settings *settings){
    return NULL;
}

