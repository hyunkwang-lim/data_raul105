/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

/* 
 * File:   grasp_input_functions.h
 * Author: fuertes
 *
 * Created on 3 de octubre de 2014, 14:35
 */

#ifndef GRASP_INPUT_FUNCTIONS_H
#define	GRASP_INPUT_FUNCTIONS_H

#ifdef	__cplusplus
extern "C" {
#endif

#include "../settings/grasp_settings.h"
#include "../input/grasp_input.h"
    
#include "grasp_input_segment.h"    
#include "../settings/grasp_settings_t.h"
#include "grasp_input_tile_description.h"
    
// Predefined pointer to functions

// This function retrieve a segment from a tile and return number of pixels retrieved. If number of pixels is lower than 0 means that there is an error retrieving pixel.
// Segement is an input/output argument. You must allocate the memory before
typedef int (*grasp_input_driver_function_t)(grasp_settings *settings, grasp_segment_t *segment, int col, int row, int itime);

// Driver type its a generic function that initialize a driver and return iterator function
typedef int (*grasp_input_driver_init_t)(grasp_settings *settings, grasp_tile_description_t *input_information);     

typedef int (*grasp_input_driver_close_t)(void);


typedef struct grasp_input_driver_t_{
    grasp_input_driver_init_t  init;
    grasp_input_driver_function_t  get_segment;
    grasp_input_driver_close_t close;
}grasp_input_driver_t;    
    
// This type define a generic function to print output information
typedef int (*grasp_input_transformer_function_t)(grasp_settings *settings,grasp_segment_t *segment);   
typedef int (*grasp_input_transformer_init_t)(grasp_settings *settings, grasp_tile_description_t *input_information);     
typedef int (*grasp_input_transformer_close_t)(void);

typedef struct grasp_input_transformer_t_{
    grasp_input_transformer_init_t  init;
    grasp_input_transformer_function_t  function;
    grasp_input_transformer_close_t close;
}grasp_input_transformer_t;


#ifdef	__cplusplus
}
#endif

#endif	/* GRASP_INPUT_FUNCTIONS_H */

