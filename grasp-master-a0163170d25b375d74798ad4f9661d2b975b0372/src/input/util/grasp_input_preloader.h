/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

/* 
 * File:   grasp_input_preloader.h
 * Author: fuertes
 *
 * Created on 16 de julio de 2014, 9:37
 */

#ifndef GRASP_INPUT_PRELOADER_H
#define	GRASP_INPUT_PRELOADER_H

#ifdef	__cplusplus
extern "C" {
#endif

#include <stdbool.h>
#include "grasp_input.h"
#include "../settings/grasp_settings.h"
#include "../input/grasp_input.h"
    
typedef struct grasp_input_preloader_t_{
    int nsubtilecols; // number of subtiles in cols
    int nsubtilerows; // number of subtiles in rows
    int nsubtiletimes; // number of subtiles in times
    int nsegmentspercol; // number of segments that it has each col
    int nsegmentsperrow; // number of segments that it has each row
    int nsegmentspertime; // number of segments that it has each times  
    int nsegmentcols; // Number of total segment in X range
    int nsegmentrows; // Number of total segment in Y range
    int nsegmenttimes; // Number of total segment in T range
    bool *subtiles; // matrix of subtiles where false means without load (t,x,y)
    bool *segment_read; // matrix of segments where false means that segment in position (t,x,y) is not read yet
}grasp_input_preloader_t;

// Initilize current grasp_input_preloader_t *preload with settings an tile
void grasp_input_preloader_init(grasp_input_preloader_t *preload, grasp_settings *settings,grasp_tile_description_t *tile_description);

// This function will return true if the segment in the position (segment_t_position, segment_x_position, segment_y_position) is already allocated
bool grasp_input_preload_is_allocated(grasp_input_preloader_t *preload, int segment_t_position, int segment_x_position, int segment_y_position);

// This function set a subtile like loaded (ready in memory)
void grasp_input_preload_allocate(grasp_input_preloader_t *preload, int segment_t_position, int segment_x_position, int segment_y_position);

// This funtion set a subtile like deallocated
void grasp_input_preload_deallocate(grasp_input_preloader_t *preload, int segment_t_position, int segment_x_position, int segment_y_position);

// This function mark a segment like read
void grasp_input_preload_set_read(grasp_input_preloader_t *preload, int segment_t_position, int segment_x_position, int segment_y_position);

// Return true if all segments in the subtile where the segment is placed is read.
bool grasp_input_preload_is_read(grasp_input_preloader_t *preload, int segment_t_position, int segment_x_position, int segment_y_position);

// Deallocate preloader
void grasp_input_preloader_destroy(grasp_input_preloader_t *preload);    



/* Example of use in a bridge:
 * 
records_type *records;
grasp_input_preloader_t preload;

  
int grasp_input_aeronet_get_segment(grasp_settings *settings, grasp_segment_t *segment, int col, int row, int itime){    
    if(grasp_input_preload_is_allocate(&preload, itime, row, col)==FALSE){ // If the segment is not allocate, we'll allocate a block of segments
        // Allocate records ...
 
        //  Set subtile like allocate
        grasp_input_preload_allocate(&preload, itime, row, col); 
    }
  
   // Set the segment as read
   grasp_input_preload_set_read(&preload, itime, row, col);
 
   // Get and set segment information ...   

   if(grasp_input_preload_is_read(&preload, itime, row, col)==TRUE){ // If the block os segments is completely read we'll deallocate it
      // Deallocate the subtile...
   
      // Marking it like deallocated 
      grasp_input_preload_deallocate(&preload, itime, row, col); 
   }
  
   return ipixel;
}

int grasp_input_aeronet_close_driver(void){
    trackmem_free(records);
    grasp_input_preloader_destroy (&preload);
    return 0; 
}

int grasp_input_aeronet_driver(grasp_settings *settings, input_tile_t *input_tile){
    // Some validations...    
 
    // Set tile information...

    // Init preload 
    grasp_input_preloader_init(&preload, settings, tile);
       
    // Return iterator
    return 0;    
}

*/
  
#ifdef	__cplusplus
}
#endif

#endif	/* GRASP_INPUT_PRELOADER_H */

