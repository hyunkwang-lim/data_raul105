/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

/* 
 * File:   grasp_controller_settings.h
 * Author: fuertes
 *
 * Created on 4 de diciembre de 2013, 11:27
 */

#ifndef GRASP_OUTPUT_SETTINGS_H
#define	GRASP_OUTPUT_SETTINGS_H

#ifdef	__cplusplus
extern "C" {
#endif
    
#include "../retrieval/constants_set/mod_globals.inc"  
#include "grasp_output_functions_settings.h"    
    
#define GRASP_MAX_OUTPUT_FUNC 5
    
typedef struct output_settings_ {
    int  nsegment_output_function;
    char segment_output_function[GRASP_MAX_OUTPUT_FUNC][_GBL_FILE_PATH_LEN];
    int  ntile_output_function;
    char tile_output_function   [GRASP_MAX_OUTPUT_FUNC][_GBL_FILE_PATH_LEN];
    int  nsegment_stream;
    char segment_stream       [GRASP_MAX_OUTPUT_FUNC][_GBL_FILE_PATH_LEN];
    int  ntile_stream;
    char tile_stream          [GRASP_MAX_OUTPUT_FUNC][_GBL_FILE_PATH_LEN];
    int  ncurrent_output_function;
    char current_output_function[GRASP_MAX_OUTPUT_FUNC][_GBL_FILE_PATH_LEN];    
    int  ncurrent_stream;
    char current_stream       [GRASP_MAX_OUTPUT_FUNC][_GBL_FILE_PATH_LEN];    
    grasp_output_segment_functions_settings segment_function;
    grasp_output_tile_functions_settings tile_function;
    grasp_output_current_functions_settings current_function;
} output_settings_t;


#ifdef	__cplusplus
}
#endif

#endif	/* GRASP_OUTPUT_SETTINGS_H */

