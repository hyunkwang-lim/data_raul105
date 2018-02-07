/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

#ifndef GRASP_INPUT_FUNCTIONS_SETTINGS_H
#define	GRASP_INPUT_FUNCTIONS_SETTINGS_H

#ifdef	__cplusplus
extern "C" {
#endif

@GRASP_INPUT_DRIVER_SETTINGS_DRIVER_INCLUDE@

    typedef struct grasp_input_driver_settings_t_{
	@GRASP_INPUT_DRIVER_SETTINGS_STRUCT@
    }grasp_input_driver_settings_t;


// For transformers

@GRASP_INPUT_TRANSFORMER_SETTINGS_TRANSFORMER_INCLUDE@

    typedef struct grasp_input_transformer_settings_t_{
      @GRASP_INPUT_TRANSFORMER_SETTINGS_STRUCT@
    }grasp_input_transformer_settings_t;



#ifdef	__cplusplus
}
#endif

#endif	/* GRASP_INPUT_FUNCTIONS_SETTINGS_H */
