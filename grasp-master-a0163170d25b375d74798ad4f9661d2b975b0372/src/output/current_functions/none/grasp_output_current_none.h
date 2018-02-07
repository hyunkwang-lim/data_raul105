/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

/* 
 * File:   grasp_output_current_print_none.h
 * Author: fuertes
 *
 * Created on September 7, 2014, 12:24 PM
 */

#ifndef GRASP_OUTPUT_CURRENT_PRINT_NONE_H
#define	GRASP_OUTPUT_CURRENT_PRINT_NONE_H

#ifdef	__cplusplus
extern "C" {
#endif

#include "../../grasp_output_stream_t.h"
#include "../../../settings/grasp_settings.h"
#include "../../../input/grasp_input.h"
#include "../../grasp_output.h"    

                               
grasp_output_current_function_t grasp_output_current_function_none();

int grasp_output_current_function_none_init(grasp_settings *settings, grasp_tile_description_t *input_information);

int grasp_output_current_function_none_close(void);    
    
int grasp_output_current_function_none_process(grasp_output_stream *stream, grasp_settings *settings, grasp_segment_t *current, output_segment_general *output, grasp_results_t *results, grasp_tile_description_t *tile_description,int icol,int irow,int itime);

grasp_settings_parameter_array *grasp_output_current_function_settings_none(grasp_settings *settings);


#ifdef	__cplusplus
}
#endif

#endif	/* GRASP_OUTPUT_CURRENT_PRINT_NONE_H */
