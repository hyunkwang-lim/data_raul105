/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

/* 
 * File:   grasp_output_tile_print_ascii.h
 * Author: fuertes
 *
 * Created on September 8, 2014, 12:41 PM
 */

#ifndef GRASP_OUTPUT_TILE_PRINT_ASCII_H
#define	GRASP_OUTPUT_TILE_PRINT_ASCII_H

#ifdef	__cplusplus
extern "C" {
#endif

#include "../../grasp_output_stream.h"
#include "../../../settings/grasp_settings.h"
#include "../../../input/grasp_input.h"
    
    
grasp_output_tile_function_t grasp_output_tile_function_ascii();

int grasp_output_tile_function_ascii_init(grasp_settings *settings, grasp_tile_description_t *input_information);

int grasp_output_tile_function_ascii_close(void);

int grasp_output_tile_function_ascii_process(grasp_output_stream *stream, grasp_settings *settings, grasp_tile_description_t *tile_description, grasp_results_t *results); 

grasp_settings_parameter_array *grasp_output_tile_function_settings_ascii(grasp_settings *settings);


#ifdef	__cplusplus
}
#endif

#endif	/* GRASP_OUTPUT_TILE_PRINT_ASCII_H */

