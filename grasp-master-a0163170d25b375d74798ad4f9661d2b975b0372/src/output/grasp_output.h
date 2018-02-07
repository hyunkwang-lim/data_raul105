/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

/* 
 * File:   grasp_output.h
 * Author: david
 *
 * Created on 5 de octubre de 2013, 14:30
 */

#ifndef GRASP_OUTPUT_H
#define	GRASP_OUTPUT_H

#include "mod_par_DLS.inc"
#include "../settings/grasp_settings.h"
#include "../input/grasp_input.h"
#include "mod_par_OS.inc"
#include "mod_par_inv.inc"
#include "mod_par_DLS_bin.inc"
#include "../retrieval/constants_set/mod_globals.inc"
#include "grasp_output_stream_t.h"

#ifdef	__cplusplus
extern "C" {
#endif
     
#include <stdint.h>
    
#include "grasp_output_segment_result.h" 
#include "grasp_output_tile_result.h"

    
    // This funtion initilize results arrays
    int grasp_output_initialize_results(grasp_settings *settings, const grasp_tile_description_t *tile_description, grasp_results_t *results);
    
    // This function assign the output from a inversion to a tile result
    void grasp_output_process_output(grasp_settings *settings,grasp_segment_t *segment,output_segment_general *output,const grasp_tile_description_t *tile_description, grasp_results_t *results,int icol,int irow,int itime);
    
    // This function deallocate input pixel results...
    void grasp_output_destroy_result(const grasp_tile_description_t *tile_description, grasp_results_t *results);     
    
#ifdef	__cplusplus
}
#endif

#endif	/* GRASP_OUTPUT_H */

