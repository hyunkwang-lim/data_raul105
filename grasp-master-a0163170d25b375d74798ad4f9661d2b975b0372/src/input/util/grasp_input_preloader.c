/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

#include <assert.h>
#include "grasp_input_preloader.h"
#include <grasp/utils.h>
#include <math.h>


void grasp_input_preloader_init(grasp_input_preloader_t *preload, grasp_settings *settings, grasp_tile_description_t *tile_description){ 
    int i,ncols, nrows, ntimes;
    
    // Calculating sizes
    // x
    if(settings->input.preload_nsegmentx<1){ // If in settings the number of segment for be preloader is negative, we are considerer that we are not going to preload, everithing will be load in first call.
        ncols=1;                             // So, this means that we are going to have only 1 subtile (the same size than the tile)
    }else{
        ncols=ceil((float)tile_description->dimensions.segment_ncols/(float)settings->input.preload_nsegmentx); 
    }
    // y
    if(settings->input.preload_nsegmenty<1){
        nrows=1;
    }else{
        nrows=ceil((float)tile_description->dimensions.segment_nrows/(float)settings->input.preload_nsegmenty);
    }
    // t
    if(settings->input.preload_nsegmenty<1){
        ntimes=1;
    }else{
        ntimes=ceil((float)tile_description->dimensions.segment_ntimes/(float)settings->input.preload_nsegmentt); 
    }
    
    
    // Allocating
    preload->subtiles = (bool *) trackmem_malloc(sizeof (bool)*ncols*nrows*ntimes);
    assert(preload->subtiles!=NULL);
    preload->segment_read = (bool *) trackmem_malloc(sizeof (bool )*tile_description->dimensions.segment_ncols*tile_description->dimensions.segment_nrows*tile_description->dimensions.segment_ntimes); 
    assert(preload->subtiles!=NULL);
    
    // Initializating values
    preload->nsubtilecols=ncols;
    preload->nsubtilerows=nrows;
    preload->nsubtiletimes=ntimes;
    
    // x
    if(settings->input.preload_nsegmentx<1){    // If the preloader is not going to be used, all segments will be loaded in the tile
       preload->nsegmentspercol=tile_description->dimensions.segment_ncols;
    }else{
       preload->nsegmentspercol=settings->input.preload_nsegmentx; 
    }
    // y
    if(settings->input.preload_nsegmenty<1){
        preload->nsegmentsperrow=tile_description->dimensions.segment_nrows;
    }else{
        preload->nsegmentsperrow=settings->input.preload_nsegmenty;
    }
    // t
    if(settings->input.preload_nsegmentt<1){
        preload->nsegmentspertime=tile_description->dimensions.segment_ntimes;
    }else{
        preload->nsegmentspertime=settings->input.preload_nsegmentt;
    }
    
    preload->nsegmentcols=tile_description->dimensions.segment_ncols;
    preload->nsegmentrows=tile_description->dimensions.segment_nrows;
    preload->nsegmenttimes=tile_description->dimensions.segment_ntimes;
    
    for (i = 0; i < ncols*nrows*ntimes; i++) {
        preload->subtiles[i]=FALSE;
    }
    for (i = 0; i < tile_description->dimensions.segment_ncols*tile_description->dimensions.segment_nrows*tile_description->dimensions.segment_ntimes; i++) {
        preload->segment_read[i]=FALSE;
    }    
}

int grasp_preload_subitle_x_position(grasp_input_preloader_t *preload, int segment_x_position){
    return  floor((float)segment_x_position/(float)preload->nsubtilecols);
}

int grasp_preload_subitle_y_position(grasp_input_preloader_t *preload, int segment_y_position){
    return  floor((float)segment_y_position/(float)preload->nsubtilerows);
}

int grasp_preload_subitle_t_position(grasp_input_preloader_t *preload, int segment_t_position){
    return  floor((float)segment_t_position/(float)preload->nsubtiletimes);
}

bool grasp_input_preload_is_allocated(grasp_input_preloader_t *preload, int segment_t_position, int segment_x_position, int segment_y_position){
    int x,y,t;
    assert(segment_x_position<preload->nsegmentcols);
    assert(segment_y_position<preload->nsegmentrows);
    assert(segment_t_position<preload->nsegmenttimes);
    
    x=grasp_preload_subitle_x_position(preload, segment_x_position);
    y=grasp_preload_subitle_y_position(preload, segment_y_position);
    t=grasp_preload_subitle_t_position(preload, segment_t_position);
    
    return preload->subtiles[index3D(t, x, y, preload->nsubtilecols, preload->nsubtilerows)];
}

void grasp_input_preload_allocate(grasp_input_preloader_t *preload, int segment_t_position, int segment_x_position, int segment_y_position){
    int x,y,t;
    assert(segment_x_position<preload->nsegmentcols);
    assert(segment_y_position<preload->nsegmentrows);
    assert(segment_t_position<preload->nsegmenttimes);
    
    x=grasp_preload_subitle_x_position(preload, segment_x_position);
    y=grasp_preload_subitle_y_position(preload, segment_y_position);
    t=grasp_preload_subitle_t_position(preload, segment_t_position);
    
    preload->subtiles[index3D(t, x, y, preload->nsubtilecols, preload->nsubtilerows)]=TRUE;
}

void grasp_input_preload_deallocate(grasp_input_preloader_t *preload, int segment_t_position, int segment_x_position, int segment_y_position){
    int x,y,t;
    assert(segment_x_position<preload->nsegmentcols);
    assert(segment_y_position<preload->nsegmentrows);
    assert(segment_t_position<preload->nsegmenttimes);
    
    x=grasp_preload_subitle_x_position(preload, segment_x_position);
    y=grasp_preload_subitle_y_position(preload, segment_y_position);
    t=grasp_preload_subitle_t_position(preload, segment_t_position);
    
    preload->subtiles[index3D(t, x, y, preload->nsubtilecols, preload->nsubtilerows)]=FALSE; 
}

void grasp_input_preload_set_read(grasp_input_preloader_t *preload, int segment_t_position, int segment_x_position, int segment_y_position){
    assert(segment_x_position<preload->nsegmentcols);
    assert(segment_y_position<preload->nsegmentrows);
    assert(segment_t_position<preload->nsegmenttimes);
    
    preload->segment_read[index3D(segment_t_position, segment_x_position, segment_y_position, preload->nsegmentcols, preload->nsegmentrows)]=TRUE;
}

bool grasp_input_preload_is_read(grasp_input_preloader_t *preload, int segment_t_position, int segment_x_position, int segment_y_position){
    int x, y, t, minx, maxx, miny, maxy, mint, maxt;
    assert(segment_x_position<preload->nsegmentcols);
    assert(segment_y_position<preload->nsegmentrows);
    assert(segment_t_position<preload->nsegmenttimes);
    
    minx=grasp_preload_subitle_x_position(preload, segment_x_position)*preload->nsegmentspercol;
    miny=grasp_preload_subitle_y_position(preload, segment_y_position)*preload->nsegmentsperrow;
    mint=grasp_preload_subitle_t_position(preload, segment_t_position)*preload->nsegmentspertime;
    
    if(grasp_preload_subitle_x_position(preload, segment_x_position)-1==preload->nsubtilecols){ // If it's the last subtile
        maxx=preload->nsegmentcols;
    }else{ // If it's an intermedian subtile
        maxx=((grasp_preload_subitle_x_position(preload, segment_x_position)+1)*preload->nsegmentcols)-1;
    }
    
    if(grasp_preload_subitle_y_position(preload, segment_y_position)-1==preload->nsubtilerows){ // If it's the last subtile
        maxy=preload->nsegmentrows;
    }else{ // If it's an intermedian subtile
        maxy=((grasp_preload_subitle_y_position(preload, segment_y_position)+1)*preload->nsegmentrows)-1;
    }    
    
    if(grasp_preload_subitle_t_position(preload, segment_t_position)-1==preload->nsubtiletimes){ // If it's the last subtile
        maxt=preload->nsegmenttimes;
    }else{ // If it's an intermedian subtile
        maxt=((grasp_preload_subitle_t_position(preload, segment_t_position)+1)*preload->nsegmenttimes)-1;
    }    
    
    for (t = mint; t < maxt; t++) {
        for (x = minx; x < maxx; x++) {
            for (y = miny; y < maxy; y++) {
                if(preload->segment_read[index3D(t, x, y, preload->nsegmentcols, preload->nsegmentrows)]==FALSE){
                    return FALSE;
                }
            }
        }
    }
    
    return TRUE;
    
}

void grasp_input_preloader_destroy(grasp_input_preloader_t *preload){
    trackmem_free(preload->subtiles);
    trackmem_free(preload->segment_read);    
}