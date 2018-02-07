/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

/**
 * @file mo_grasp_controller.h
 * @author David Fuertes
 * @date 25 Oct 2013
 * @brief C interface to access some Fortran functions from retrieval algorithm
 *
 * GRASP algorithm is developed in Fortran. Some functions from scientific core 
 * are used by the control unit. Controller is the responsible to organize the
 * workflow for processing a big tile (many segments). In this workflow the controller
 * will need to call inversion subroutine (implemented in Fortran) and other 
 * functions in order to prepare inputs and outputs. This file contains a C interface
 * for some functions implemented in scientific core but needed to be called from
 * the controller
 */

#ifndef MO_GRASP_CONTROLLER_H
#define	MO_GRASP_CONTROLLER_H

#ifdef	__cplusplus
extern "C" {
#endif
    
#include "mod_par_inv.inc" 
#include "../settings/grasp_settings_t.h"
#include "../input/grasp_input_segment.h"
#include "../output/grasp_output_segment_result.h"
#include "../retrieval/mod_par_inv.inc"
    
/**
 * Before calling inversion scientific core needs to do some changes and set
 * some input parameters based on input settings and input measures
 * @param RIN Input settings read which will be checked and prepared to process a specific segment
 * @param segment_meas Input measures will be checked and prepared to process with specific behavior (settings)
 */
extern void grasp_prepare_segment_settings(retr_input *RIN, sensor_data_t *segment_meas);

/**
 * Before calling inversion algorithm it needs to initialize some global variables.
 * It take some time but once it is initialized you can process many segments so it is only
 * necessary to call at the beginning of first inversion.
 * @param RIN Input retrieval settings
 */
extern void grasp_init_inversion(retr_input *RIN);

/**
 * C interface for inversion Fortran function.
 * @param RIN Input retrieval settings
 * @param sdata Input measures
 * @param iguess Initial guess for each pixel. If it is -999 it will replace the value read from settings otherwise this value will be used for a specific pixel
 * @param edges Edges information of current segment
 * @param ROUT Output structure filled with results from inversion
 */
extern int grasp_input_inversion(const retr_input *RIN, const sensor_data_t *sdata, const float iguess[_KIMAGE][_KPARS],  const segment_edges *edges, output_segment_general *ROUT);

/**
 * Init funtion allocate some global structures which need to be deallocated.
 * It has to be called after inverting all segments
 * @param RIN Input retrieval settings.
 */
extern void grasp_finalize_inversion(retr_input *RIN);

#ifdef	__cplusplus
}
#endif

#endif	/* MO_GRASP_CONTROLLER_H */

