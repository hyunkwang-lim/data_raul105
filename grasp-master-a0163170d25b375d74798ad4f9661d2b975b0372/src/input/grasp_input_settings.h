/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

/* 
 * File:   grasp_input_segment.h
 * Author: fuertes
 *
 * Created on May 30, 2014, 11:36 PM
 */

#ifndef GRASP_INPUT_SETTINGS_H
#define	GRASP_INPUT_SETTINGS_H

#ifdef	__cplusplus
extern "C" {
#endif

#include "../retrieval/constants_set/mod_globals.inc" 
#include "grasp_input_functions_settings.h"
    
#define GRASP_INPUT_MAX_FILES 3000
#define GRASP_INPUT_MAX_TRANSFORMERS 5    
    
typedef struct {
	/// the longitude of the pixel, i.e. around the globe
	float lon;
	/// the latitude of the pixel, i.e. towards the poles
	float lat;
} coord_t;

typedef struct {
    // Row coordinate
    int row;
    // Column coordinate
    int col;
} coord_grid_t;


typedef struct input_settings_t_{
	/// the the width and heigth of a pixel in geo-coords
	coord_t resolution;
        // the center of tile
        coord_t coordinates;
        // corner in row cols
        coord_grid_t coordinates_grid;
	// the offset of the segmentation grid taking 0,0 like reference. Output will be refered to 0,0.
        // If you want that the output will be refered to set 0,0 in this parameter
	coord_grid_t grid_offset;         
        // Reference for coordinates: center or corner
        char coordinates_refence[25]; 
        // type of coordinates: latlon or rowcol
        char coordinates_type[25];
	/// the width of the covered area
	int area_width;
	/// the height of the covered area
	int area_height;
        // Initial date for data processing
        char time_from[255];
        // Final date for data processing
        char time_to[255];
	/// the width of a segment
	int segment_nx;
	/// the height of a segment
	int segment_ny;
	// the number of different times that a segment can have
	int segment_nt;
        // Number of files to process
        int nfiles;
        // Name of file if it is necessary
        char files[GRASP_INPUT_MAX_FILES][_GBL_FILE_PATH_LEN];
        // Name of driver for retrieve data
        char driver_name[_GBL_FILE_PATH_LEN];
        // If this option is true raw segment information will be printed in screen
        char print_raw_segment[_GBL_FILE_PATH_LEN];;
        // If this option is true segment information after clean NaN values will be printed in screen
        char print_clean_segment[_GBL_FILE_PATH_LEN];;
        // If it is true the name of used files will be dump in screen
        char print_used_files[_GBL_FILE_PATH_LEN];;
        // Stream where dump sdata information
        char sdata_stream[_GBL_FILE_PATH_LEN];
        // Stream where dump imagedat (initial guess) information
        char imagedat_stream[_GBL_FILE_PATH_LEN];
        // Number of pixels loaded in block for X
        int preload_nsegmentx;
        // Number of pixels loaded in block for Y
        int preload_nsegmenty;
        // Number of pixels loaded in block for T
        int preload_nsegmentt;
        // Settings from drivers
        grasp_input_driver_settings_t driver;
        // Settings from transformers
        grasp_input_transformer_settings_t transformer;
        // Number of transformers to be used
        int ntransformers;
        // Name of transformers to be used
        char transformers_name[GRASP_INPUT_MAX_TRANSFORMERS][_GBL_FILE_PATH_LEN];
} input_settings_t;

#ifdef	__cplusplus
}
#endif

#endif	/* GRASP_INPUT_SETTINGS_H */

