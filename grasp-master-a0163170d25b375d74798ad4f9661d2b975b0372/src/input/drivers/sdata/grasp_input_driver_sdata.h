/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */


grasp_input_driver_t grasp_input_driver_sdata();

int grasp_input_driver_sdata_init(grasp_settings *settings, grasp_tile_description_t *input_information);

int grasp_input_driver_sdata_get_segment(grasp_settings *settings, grasp_segment_t *segment, int col, int row, int itime);

int grasp_input_driver_sdata_close(void);

grasp_settings_parameter_array *grasp_input_driver_settings_sdata(grasp_settings *settings);