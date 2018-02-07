/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

#ifndef GRASP_OUTPUT_FUNCTIONS_SETTINGS_H
#define	GRASP_OUTPUT_FUNCTIONS_SETTINGS_H

#ifdef	__cplusplus
extern "C" {
#endif

@GRASP_OUTPUT_FUNCTIONS_SETTINGS_SEGMENT_INCLUDES@
@GRASP_OUTPUT_FUNCTIONS_SETTINGS_TILE_INCLUDES@
@GRASP_OUTPUT_FUNCTIONS_SETTINGS_CURRENT_INCLUDES@

	typedef struct grasp_output_segment_functions_settings_{@GRASP_OUTPUT_FUNCTIONS_SETTINGS_SEGMENT_STRUCT@
	}grasp_output_segment_functions_settings;

	typedef struct grasp_output_tile_functions_settings_{@GRASP_OUTPUT_FUNCTIONS_SETTINGS_TILE_STRUCT@
	}grasp_output_tile_functions_settings;
	
  typedef struct grasp_output_current_functions_settings_{@GRASP_OUTPUT_FUNCTIONS_SETTINGS_CURRENT_STRUCT@
	}grasp_output_current_functions_settings;

#ifdef	__cplusplus
}
#endif

#endif	/* GRASP_OUTPUT_FUNCTIONS_SETTINGS_H */
