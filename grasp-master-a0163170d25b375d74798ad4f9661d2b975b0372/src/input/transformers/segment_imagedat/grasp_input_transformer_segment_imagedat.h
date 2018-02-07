/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

#ifndef GRASP_INPUT_TRANSFORMER_SEGMENT_IMAGEDAT_H
#define	GRASP_INPUT_TRANSFORMER_SEGMENT_IMAGEDAT_H

#ifdef	__cplusplus
extern C {
#endif

#include "../../../settings/grasp_settings.h"
#include "../../../input/grasp_input.h"
#include "../../../input/grasp_input_functions.h"   

grasp_input_transformer_t grasp_input_transformer_segment_imagedat();

int grasp_input_transformer_segment_imagedat_init(grasp_settings *settings, grasp_tile_description_t *input_information);

int grasp_input_transformer_segment_imagedat_close(void);

int grasp_input_transformer_segment_imagedat_function(grasp_settings *settings,grasp_segment_t *segment);

grasp_settings_parameter_array *grasp_input_transformer_settings_segment_imagedat(grasp_settings *settings);

#ifdef	__cplusplus
}
#endif

#endif	/* GRASP_INPUT_TRANSFORMER_SEGMENT_IMAGEDAT_H */
