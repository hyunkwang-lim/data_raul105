/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

/* 
 * File:   grasp_output_functions.h
 * Author: fuertes
 *
 * Created on 3 de octubre de 2014, 14:36
 */

#ifndef GRASP_OUTPUT_FUNCTIONS_H
#define	GRASP_OUTPUT_FUNCTIONS_H

#ifdef	__cplusplus
extern "C" {
#endif

#include "grasp_output.h"    
#include "grasp_output_segment_result.h"
#include "../controller/grasp_controller_functions.h"
#include "grasp_output_tile_result.h"
   
// Generic function for initialize segment output function
typedef int (*grasp_output_segment_function_init_t)(grasp_settings *settings, grasp_tile_description_t *input_information);   
// This type define a generic function to print output information
typedef int (*grasp_output_segment_function_process_t)(grasp_output_stream *stream, grasp_settings *settings,grasp_segment_t *segment,output_segment_general *output,grasp_tile_description_t *tile_description,int icol,int irow,int itime);     
// Generic function for close segment output function
typedef int (*grasp_output_segment_function_close_t)(void);

typedef struct grasp_output_segment_function_t_{
    grasp_output_segment_function_init_t  init;
    grasp_output_segment_function_process_t  function;
    grasp_output_segment_function_close_t close;
}grasp_output_segment_function_t;


// Generic function for initialize tile output function
typedef int (*grasp_output_tile_function_init_t)(grasp_settings *settings, grasp_tile_description_t *input_information);   
// This type define a generic function to print output information
typedef int (*grasp_output_tile_function_process_t)(grasp_output_stream *stream, grasp_settings *settings, grasp_tile_description_t *tile_description, grasp_results_t *results);   
// Generic function for close tile output function
typedef int (*grasp_output_tile_function_close_t)(void);

typedef struct grasp_output_tile_function_t_{
    grasp_output_tile_function_init_t  init;
    grasp_output_tile_function_process_t  function;
    grasp_output_tile_function_close_t close;
}grasp_output_tile_function_t;    


// Generic function for initialize current output function
typedef int (*grasp_output_current_function_init_t)(grasp_settings *settings, grasp_tile_description_t *input_information);   
// This type define a generic function to print output information
typedef int (*grasp_output_current_function_process_t)(grasp_output_stream *stream, grasp_settings *settings,grasp_segment_t *current,output_segment_general *output, grasp_results_t *results, grasp_tile_description_t *tile_description,int icol,int irow,int itime);     
// Generic function for close current output function
typedef int (*grasp_output_current_function_close_t)(void);

typedef struct grasp_output_current_function_t_{
    grasp_output_current_function_init_t  init;
    grasp_output_current_function_process_t  function;
    grasp_output_current_function_close_t close;
}grasp_output_current_function_t;

#ifdef	__cplusplus
}
#endif

#endif	/* GRASP_OUTPUT_FUNCTIONS_H */

