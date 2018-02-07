/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

/**
 * @file grasp_controller_functions.h
 * @author David Fuertes
 * @date 3 Oct 2013
 * @brief Structure to load extension functions
 *
 * The behavior of GRAPS control unit can be extended through external functions
 * which will be called in different stages of the process. Following file defines
 * the functions which will be used.
 */

#ifndef GRASP_CONTROLLER_FUNCTIONS_H
#define	GRASP_CONTROLLER_FUNCTIONS_H

#ifdef	__cplusplus
extern "C" {
#endif

#include "../input/grasp_input_functions.h"
#include "../output/grasp_output_functions.h"

/**
 * This structure contains pointers to extension functions that will be used in the process
 */   
typedef struct grasp_processing_functions_t_{
    grasp_input_driver_t driver; /**< @brief Driver which is used */   
    
    int ntransformers; /**< @brief Number of transformer function will be used */   
    grasp_input_transformer_t *transformers;   /**< @brief  Transformer extension functions */   
    
    int nsegment_output_functions; /**< @brief Number of extension function to print the output from a single segment */   
    grasp_output_segment_function_t *segment_output_functions; /**< @brief  Output segment function pointers */   
    
    int ntile_output_functions; /**< @brief  Number of extension function to print the output from a entire tile */   
    grasp_output_tile_function_t *tile_output_functions;  /**< @brief Output tile function pointers */   
    
    int ncurrent_output_functions; /**< @brief Number of extension function to print the output from the part of the tile already processed (partial tile) */   
    grasp_output_current_function_t *current_output_functions; /**< @brief Output current function pointers */   
}grasp_processing_functions_t;


#ifdef	__cplusplus
}
#endif

#endif	/* GRASP_CONTROLLER_FUNCTIONS_H */

