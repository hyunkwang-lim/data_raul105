/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

/* 
 * File:   grasp_input_tile_description.h
 * Author: fuertes
 *
 * Created on 3 de octubre de 2014, 14:53
 */

#ifndef GRASP_INPUT_TILE_DESCRIPTION_H
#define	GRASP_INPUT_TILE_DESCRIPTION_H

#ifdef	__cplusplus
extern "C" {
#endif

    
typedef struct grasp_tile_dimensions_t_{
    int segment_nrows;  // Number of rows in SEGMENT
    int segment_ncols;  // Number of cols in SEGMENT
    int segment_ntimes; // Number of times in SEGMENT
    int tile_nt; // Number of total different t's in tile
    int tile_nx; // Number of total different x's in tile
    int tile_ny; // Number of total different y's in tile        
    int npixel_estimated;
}grasp_tile_dimensions_t;

typedef struct grasp_tile_description_t_{
    // Information set by framework
    int ninput_files; // Number of input files
    char **input_files; // Name of each input file
    
    grasp_tile_dimensions_t dimensions;
    
    // This information should be provided by the drivers and it should be available at the end of the tile process (Except for polder that use this information in output segment.)
    int nused_files; // Number of used files
    char **used_files; // Name of used files  
}grasp_tile_description_t;


#ifdef	__cplusplus
}
#endif

#endif	/* GRASP_INPUT_TILE_DESCRIPTION_H */

