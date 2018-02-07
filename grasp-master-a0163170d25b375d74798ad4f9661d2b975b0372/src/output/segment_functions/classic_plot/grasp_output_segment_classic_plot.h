/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

/* 
 * File:   grasp_output_segment_print_plot_classic.h
 * Author: fuertes
 *
 * Created on September 7, 2014, 12:23 PM
 */

#ifndef GRASP_OUTPUT_SEGMENT_PRINT_PLOT_CLASSIC_H
#define	GRASP_OUTPUT_SEGMENT_PRINT_PLOT_CLASSIC_H

#ifdef	__cplusplus
extern "C" {
#endif

#include "../../grasp_output_stream.h"

grasp_output_segment_function_t grasp_output_segment_function_classic_plot();

int grasp_output_segment_function_classic_plot_init(grasp_settings *settings, grasp_tile_description_t *input_information);

int grasp_output_segment_function_classic_plot_close(void);  
    
int grasp_output_segment_function_classic_plot_process(grasp_output_stream *stream, grasp_settings *settings, grasp_segment_t *segment, output_segment_general *output, grasp_tile_description_t *tile_description,int icol,int irow,int itime);

grasp_settings_parameter_array *grasp_output_segment_function_settings_classic_plot(grasp_settings *settings);

extern void grasp_output_segment_print_classic_plot_f90(retr_input *RIN, sensor_data_t *sdata,output_segment_general *ROUT);


#ifdef	__cplusplus
}
#endif

#endif	/* GRASP_OUTPUT_SEGMENT_PRINT_PLOT_CLASSIC_H */

