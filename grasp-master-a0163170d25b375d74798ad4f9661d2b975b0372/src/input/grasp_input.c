/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "grasp_input.h"
#include "grasp_input_load_functions.h"
#include <assert.h>
#include "mod_par_inv.inc"
#include "../global/grasp_retrieval_meas_type.h"
#include "drivers/sdata/sdata.h"
#include "drivers/sdata/sdata-impl.h" /* SDATA_RECORD */
#include "../output/grasp_output.h"
#include "../output/grasp_output_load_function.h"
#include "../output/grasp_output_stream.h"

int grasp_input_initialize_tile_description(grasp_settings *settings, grasp_tile_description_t *tile_description){   
    tile_description->nused_files=0;
    tile_description->used_files=NULL;    
    
    tile_description->dimensions.segment_ncols=0;
    tile_description->dimensions.segment_nrows=0;
    tile_description->dimensions.segment_ntimes=0;
    tile_description->dimensions.tile_nx=0;
    tile_description->dimensions.tile_ny=0;
    tile_description->dimensions.tile_nt=0;
   
    tile_description->ninput_files=grasp_input_expand_files(settings,&(tile_description->input_files));  
 
    if(tile_description->ninput_files==0){
        fprintf(stderr, "ERROR: Drivers need at least one valid data file. Please, check input.file parameter in settings file\n");
        exit(-1);        
    }   
    
    return 0;
}

int grasp_input_position_of_inversion(const grasp_tile_dimensions_t *tile_dimensions, int id_inversion, int *icol, int *irow, int *itime){
    int indexes[3];
    int maximums[3];
    
    if(id_inversion<0){
        return -1;
    }
    
    maximums[0]=tile_dimensions->segment_ntimes;
    maximums[1]=tile_dimensions->segment_nrows;
    maximums[2]=tile_dimensions->segment_ncols;
             
    if(id_inversion > maximums[0]*maximums[1]*maximums[2]){ // If the index is out bound it return -1
        return -1;
    }
    
    indexesofND(id_inversion, indexes, 3,  maximums);
    
    *itime=indexes[0];
    *irow=indexes[1];
    *icol=indexes[2];
    
    return 0;
}

int grasp_input_extract_segment(grasp_settings *settings, grasp_input_driver_t *driver, int ntransformers,grasp_input_transformer_t *transformers, grasp_segment_t *segment, grasp_results_t *results, grasp_tile_dimensions_t *tile_dimensions, int id_inversion){
    int i,npixels;
    grasp_output_stream grasp_stream;
    FILE *f;
    int icol, irow,itime;
        
    if(settings->input.driver.sdata.debug){
        sdata_set_debug_stream(stderr);
    }
    // Set to 0 all segment
    memset(&(segment->sdata),0,sizeof(sensor_data_t));
    grasp_parameters_initialize(segment->iguess);
    
    grasp_input_position_of_inversion(tile_dimensions,id_inversion, &icol, &irow,&itime);
    
    segment->sdata.id_inversion=id_inversion;    
    
    // Retrieve segment with the specific driver
    npixels = driver->get_segment(settings, segment, icol, irow, itime);    

    grasp_input_edges_find(settings, segment, tile_dimensions ,results, icol, irow, itime);       
                
#ifdef DEBUG2
    fprintf(stderr, "%s:%d: %s: ", __FILE__, __LINE__, __func__);
    grasp_input_print_segment(stderr, "after get_segment():\n", segment);
#endif
    
    // Printing raw segment if it is requested
    grasp_output_stream_initialize(settings->input.print_raw_segment, &grasp_stream);
    f = grasp_output_stream_open(&grasp_stream, settings, segment, NULL, tile_dimensions, icol, irow, itime);
    if (grasp_output_stream_writable(&grasp_stream)) {
        fprintf(f, "Raw segment data:\n");
        grasp_input_print_segment(f, "segment", segment);
    }
    grasp_output_stream_close(&grasp_stream);
    grasp_output_stream_destroy(&grasp_stream);
                  
    // Calling segment transformers
    for (i = 0; i < ntransformers; i++) {                 
        transformers[i].function(settings,segment);
    }    
    
    // Here we'll apply some data filters like "only land" or similar...
    
    // Here we'll call new clean segment function
        
    grasp_input_clean_segment(segment, -999);
    
    // Printing raw segment if it is requested
    grasp_output_stream_initialize(settings->input.print_clean_segment, &grasp_stream);
    f = grasp_output_stream_open(&grasp_stream, settings, segment, NULL, tile_dimensions, icol, irow, itime);
    if (grasp_output_stream_writable(&grasp_stream)) {
        fprintf(f, "Clean segment data:\n");
        grasp_input_print_segment(f, "segment", segment);
    }
    grasp_output_stream_close(&grasp_stream);
    grasp_output_stream_destroy(&grasp_stream);

    return npixels;
}

// This function return true if a time (unixtimestamp) is present in the segement
bool segment_has_time_index(grasp_segment_t *segment, int it){
    int i;
    
    for (i = 0; i < segment->sdata.npixels; i++) {
        if(segment->sdata.pixel[i].out_t==it){
            return true;
        }
    }
    
    return false;
}

void grasp_input_edges_loop_and_copy_data_x(grasp_settings *settings, grasp_segment_t *segment, grasp_tile_dimensions_t *tile_dimensions, grasp_results_t *results, int minit, int maxit,int miniy, int maxiy, int minix, int maxix, int before_after){   
    int it,ix,iy; //index of loops
    int ipix; // index of pixel in result map
    int ipar, npar; // number or parameters and index of parameters in loops
    int segment_it;//, iit;
    bool timeblock_has_data;
    assert(before_after>=0 && before_after<=1);
        
    // total parameter per pixel
    npar=settings->retrieval.KNSING;

    assert(
       minit >=0 && minit <=tile_dimensions->tile_nt && maxit >=0 && maxit <=tile_dimensions->tile_nt && 
       miniy >=0 && miniy <=tile_dimensions->tile_ny && maxiy >=0 && maxiy <=tile_dimensions->tile_ny && 
       minix >=0 && minix <=tile_dimensions->tile_nx && maxix >=0 && maxix <=tile_dimensions->tile_nx 
    );
    segment_it=0;
    for (it = minit; it < maxit; it++) {
        timeblock_has_data=false;
        for (iy = miniy; iy < maxiy; iy++) {
            for (ix = minix; ix < maxix; ix++) {
                // Getting pixel position in result map
                ipix=index3D(it,ix,iy, tile_dimensions->tile_nx,tile_dimensions->tile_ny);
                assert(ipix>=0 && ipix< tile_dimensions->tile_nt*tile_dimensions->tile_nx*tile_dimensions->tile_ny);
                if(results->tile_result_map[ipix]!=NULL && segment_has_time_index(segment, it)){
                    timeblock_has_data=true;
                    // Copying the information
                    assert(segment_it>=0 && segment_it<_KITIME);
                    assert(iy-miniy>=0 && iy-miniy<_KIY);
                    assert(ix-minix>=0 && ix-minix<_KIEDGE);
                    segment->edges.group_x[before_after].t[segment_it]=grasp_output_tile_pixel_information_time(results,it,ix,iy);
                    segment->edges.group_x[before_after].x [segment_it][iy-miniy][ix-minix]=grasp_output_tile_pixel_information_longitude(results,it,ix,iy);
                    segment->edges.group_x[before_after].y [segment_it][iy-miniy][ix-minix]=grasp_output_tile_pixel_information_latitude(results,it,ix,iy);
                    segment->edges.group_x[before_after].icloud [segment_it][iy-miniy][ix-minix]=1;
                    segment->edges.group_x[before_after].out_t[segment_it]=grasp_output_tile_pixel_information_out_t(results,it,ix,iy);
                    segment->edges.group_x[before_after].out_x [segment_it][iy-miniy][ix-minix]=grasp_output_tile_pixel_information_out_x(results,it,ix,iy);
                    segment->edges.group_x[before_after].out_y [segment_it][iy-miniy][ix-minix]=grasp_output_tile_pixel_information_out_y(results,it,ix,iy);
                    for (ipar = 0; ipar < npar; ipar++) { 
                        assert(ipar>=0 && ipar<_KPARS);
                        segment->edges.group_x[before_after].AP[segment_it][iy-miniy][ix-minix][ipar]=grasp_output_tile_retrieval_par_parameters(results,it,ix,iy,ipar);
                    }          
                }
            }
        }
        if(timeblock_has_data){
            segment_it++;
        }
    }
    segment->edges.group_x[before_after].nt=segment_it; 
}

void grasp_input_edges_loop_and_copy_data_y(grasp_settings *settings, grasp_segment_t *segment, grasp_tile_dimensions_t *tile_dimensions, grasp_results_t *results, int minit, int maxit,int miniy, int maxiy, int minix, int maxix, int before_after){   
    int it,ix,iy; //index of loops
    int ipix; // index of pixel in result map
    int ipar, npar; // number or parameters and index of parameters in loops
    int segment_it;//, iit; 
    bool timeblock_has_data;
    assert(before_after>=0 && before_after<=1);
    
    // total parameter per pixel
    npar=settings->retrieval.KNSING;
  
    assert(
       minit >=0 && minit <=tile_dimensions->tile_nt && maxit >=0 && maxit <=tile_dimensions->tile_nt && 
       miniy >=0 && miniy <=tile_dimensions->tile_ny && maxiy >=0 && maxiy <=tile_dimensions->tile_ny && 
       minix >=0 && minix <=tile_dimensions->tile_nx && maxix >=0 && maxix <=tile_dimensions->tile_nx 
    );
    segment_it=0;
    for (it = minit; it < maxit; it++) {
        timeblock_has_data=false;
        for (iy = miniy; iy < maxiy; iy++) {
            for (ix = minix; ix < maxix; ix++) {
                // Getting pixel position in result map
                ipix=index3D(it,ix,iy, tile_dimensions->tile_nx,tile_dimensions->tile_ny);  
                assert(ipix>=0 && ipix< tile_dimensions->tile_nt*tile_dimensions->tile_nx*tile_dimensions->tile_ny);
                if(results->tile_result_map[ipix]!=NULL && segment_has_time_index(segment, it)){
                    timeblock_has_data=true;
                    // Copying the information
                    assert(segment_it>=0 && segment_it<_KITIME);
                    assert(iy-miniy>=0 && iy-miniy<_KIEDGE);
                    assert(ix-minix>=0 && ix-minix<_KIX);    
                    segment->edges.group_y[before_after].t[segment_it]=grasp_output_tile_pixel_information_time(results,it,ix,iy);
                    segment->edges.group_y[before_after].x [segment_it][iy-miniy][ix-minix]=grasp_output_tile_pixel_information_longitude(results,it,ix,iy);
                    segment->edges.group_y[before_after].y [segment_it][iy-miniy][ix-minix]=grasp_output_tile_pixel_information_latitude(results,it,ix,iy);
                    segment->edges.group_y[before_after].icloud [segment_it][iy-miniy][ix-minix]=1;
                    segment->edges.group_y[before_after].out_t[segment_it]=grasp_output_tile_pixel_information_out_t(results,it,ix,iy);
                    segment->edges.group_y[before_after].out_x [segment_it][iy-miniy][ix-minix]=grasp_output_tile_pixel_information_out_x(results,it,ix,iy);
                    segment->edges.group_y[before_after].out_y [segment_it][iy-miniy][ix-minix]=grasp_output_tile_pixel_information_out_y(results,it,ix,iy);
                    for (ipar = 0; ipar < npar; ipar++) { 
                        assert(ipar>=0 && ipar<_KPARS);
                        segment->edges.group_y[before_after].AP[segment_it][iy-miniy][ix-minix][ipar]=grasp_output_tile_retrieval_par_parameters(results,it,ix,iy,ipar);
                    }          
                }
            }
        }
        if(timeblock_has_data){
            segment_it++;
        }
    }   
    segment->edges.group_y[before_after].nt=segment_it;   
}

int grasp_input_edges_search_in_T(grasp_segment_t *segment, grasp_results_t *results, grasp_tile_dimensions_t *tile_dimensions, int it, int iy, int ix, int dir){   
    int iit; //index of loops
    int ipix; // index of pixel in result map
    bool found=false;
    
    iit=it+dir;

    while(iit>=0 && iit <tile_dimensions->tile_nt && !found){
        // Getting pixel position in result map
        ipix=index3D(iit,ix,iy, tile_dimensions->tile_nx,tile_dimensions->tile_ny);        
        if(results->tile_result_map[ipix]!=NULL){
            found=true;  
        }
        iit=iit+dir;
    }
    
    iit=iit-dir;
    
    if(found){
        return iit;
    }else{
        return -1;
    }
    
}
 

void grasp_input_edges_find_X_before(grasp_settings *settings, grasp_segment_t *segment, grasp_tile_dimensions_t *tile_dimensions, grasp_results_t *results, int icol, int irow, int itime){
    int maxit,minit,maxiy,miniy,maxix,minix; // range of loops
    
    minit=itime*settings->input.segment_nt;
    maxit=(itime*settings->input.segment_nt)+settings->input.segment_nt;
    if(maxit>tile_dimensions->tile_nt){
        maxit=tile_dimensions->tile_nt;
    }
    if(minit>maxit){
        minit=maxit;
    }
    
    miniy=(irow*settings->input.segment_ny);
    maxiy=(irow*settings->input.segment_ny)+settings->input.segment_ny;
    if(maxiy>tile_dimensions->tile_ny){
        maxiy=tile_dimensions->tile_ny;
    }  
    if(miniy>maxiy){
        miniy=maxiy;
    }
    
    minix=(icol*settings->input.segment_nx)-settings->retrieval.edges.nx;
    if(minix<0){
        minix=0;
    }
    maxix=(icol*settings->input.segment_nx); 
    if(minix>maxix){
        minix=maxix;
    }
    
    segment->edges.group_x[0].nx=settings->retrieval.edges.nx;
    segment->edges.group_x[0].ny=segment->sdata.ny;
        
    grasp_input_edges_loop_and_copy_data_x(settings, segment, tile_dimensions, results, minit,maxit,miniy,maxiy,minix,maxix,0);    
}

void grasp_input_edges_find_X_after(grasp_settings *settings, grasp_segment_t *segment, grasp_tile_dimensions_t *tile_dimensions, grasp_results_t *results, int icol, int irow, int itime){
    int maxit,minit,maxiy,miniy,maxix,minix; // range of loops
    
    minit=itime*settings->input.segment_nt;
    maxit=(itime*settings->input.segment_nt)+settings->input.segment_nt;
    if(maxit>tile_dimensions->tile_nt){
        maxit=tile_dimensions->tile_nt;
    }  
    if(minit>maxit){
        minit=maxit;
    }
    
    miniy=(irow*settings->input.segment_ny);
    maxiy=(irow*settings->input.segment_ny)+settings->input.segment_ny;
    if(maxiy>tile_dimensions->tile_ny){
        maxiy=tile_dimensions->tile_ny;
    }  
    if(miniy>maxiy){
        miniy=maxiy;
    }
    
    minix=(icol*settings->input.segment_nx)+settings->input.segment_nx;
    maxix=(icol*settings->input.segment_nx)+settings->input.segment_nx+settings->retrieval.edges.nx; 
    if(maxix>tile_dimensions->tile_nx){
        maxix=tile_dimensions->tile_nx;
    }  
    if(minix>maxix){
        minix=maxix;
    }
    
    segment->edges.group_x[1].nx=settings->retrieval.edges.nx;
    segment->edges.group_x[1].ny=segment->sdata.ny; 
    
    grasp_input_edges_loop_and_copy_data_x( settings, segment, tile_dimensions, results, minit,maxit,miniy,maxiy,minix,maxix,1);
}

void grasp_input_edges_find_Y_before(grasp_settings *settings, grasp_segment_t *segment, grasp_tile_dimensions_t *tile_dimensions, grasp_results_t *results, int icol, int irow, int itime){
    int maxit,minit,maxiy,miniy,maxix,minix; // range of loops
    
    minit=itime*settings->input.segment_nt;
    maxit=(itime*settings->input.segment_nt)+settings->input.segment_nt;
    if(maxit>tile_dimensions->tile_nt){
        maxit=tile_dimensions->tile_nt;
    }  
    if(minit>maxit){
        minit=maxit;
    }
    
    miniy=(irow*settings->input.segment_ny)-settings->retrieval.edges.ny;
    if(miniy<0){
        miniy=0;
    }
    maxiy=(irow*settings->input.segment_ny);
    if(miniy>maxiy){
        miniy=maxiy;
    }
    
    minix=(icol*settings->input.segment_nx);
    maxix=(icol*settings->input.segment_nx)+settings->input.segment_nx;
    if(maxix>tile_dimensions->tile_nx){
        maxix=tile_dimensions->tile_nx;
    } 
    if(minix>maxix){
        minix=maxix;
    }
    
    segment->edges.group_y[0].nx=segment->sdata.nx;
    segment->edges.group_y[0].ny=settings->retrieval.edges.ny;
    
    grasp_input_edges_loop_and_copy_data_y( settings, segment, tile_dimensions, results, minit,maxit,miniy,maxiy,minix,maxix,0);
}

void grasp_input_edges_find_Y_after(grasp_settings *settings, grasp_segment_t *segment, grasp_tile_dimensions_t *tile_dimensions, grasp_results_t *results, int icol, int irow, int itime){
    int maxit,minit,maxiy,miniy,maxix,minix; // range of loops
    
    minit=itime*settings->input.segment_nt;
    maxit=(itime*settings->input.segment_nt)+settings->input.segment_nt;
    if(maxit>tile_dimensions->tile_nt){
        maxit=tile_dimensions->tile_nt;
    } 
    if(minit>maxit){
        minit=maxit;
    }
    
    miniy=(irow*settings->input.segment_ny)+settings->input.segment_ny;
    maxiy=(irow*settings->input.segment_ny)+settings->input.segment_ny+settings->retrieval.edges.ny; 
    if(maxiy>tile_dimensions->tile_ny){
        maxiy=tile_dimensions->tile_ny;
    }       
    if(miniy>maxiy){
        miniy=maxiy;
    }
    
    minix=(icol*settings->input.segment_nx);
    maxix=(icol*settings->input.segment_nx)+settings->input.segment_nx;
    if(maxix>tile_dimensions->tile_nx){
        maxix=tile_dimensions->tile_nx;
    }   
    if(minix>maxix){
        minix=maxix;
    }
    
    segment->edges.group_y[1].nx=segment->sdata.nx; 
    segment->edges.group_y[1].ny=settings->retrieval.edges.ny;    
    
    grasp_input_edges_loop_and_copy_data_y( settings, segment, tile_dimensions, results, minit,maxit,miniy,maxiy,minix,maxix,1);
}

void grasp_input_edges_find_T_before(grasp_settings *settings, grasp_segment_t *segment, grasp_tile_dimensions_t *tile_dimensions, grasp_results_t *results, int icol, int irow, int itime){
    int maxiy,miniy,maxix,minix; // range of loops
    int it, iy, ix, i_nt, t_gotten, ipar, times[_KIEDGE], ntimes, ipix, npar;
    
    npar=settings->retrieval.KNSING;
        
    miniy=(irow*settings->input.segment_ny);
    maxiy=(irow*settings->input.segment_ny)+settings->input.segment_ny;
    if(maxiy>tile_dimensions->tile_ny){
        maxiy=tile_dimensions->tile_ny;
    }     
    if(miniy>maxiy){
        miniy=maxiy;
    }
    
    minix=(icol*settings->input.segment_nx);
    maxix=(icol*settings->input.segment_nx)+settings->input.segment_nx;
    if(maxix>tile_dimensions->tile_nx){
        maxix=tile_dimensions->tile_nx;
    } 
    if(minix>maxix){
        minix=maxix;
    }
        
    segment->edges.group_t[0].nx=segment->sdata.nx;
    segment->edges.group_t[0].ny=segment->sdata.ny;
    
    for (iy = miniy; iy < maxiy; iy++) {
        for (ix = minix; ix < maxix; ix++) {
            it=(itime*settings->input.segment_nt);
            ntimes=0;
            for (i_nt = 0; i_nt < settings->retrieval.edges.nt; i_nt++) {
                t_gotten=grasp_input_edges_search_in_T( segment, results, tile_dimensions, it, iy, ix, -1);
                if(t_gotten<0){
                    break;
                }
                it=t_gotten;
                assert(it>=0 && it <tile_dimensions->tile_nt);
                times[i_nt]=t_gotten;       
                ntimes++;
            }
            // We have to order the results like inverse order and with values in 0 always...
            for (i_nt = ntimes; i_nt >0; i_nt--) {
                // Copying the information
                assert(i_nt-1>=0);
                ipix=index3D(times[i_nt-1],ix,iy, tile_dimensions->tile_nx,tile_dimensions->tile_ny);
                assert(results->tile_result_map[ipix]!=NULL);     
                assert(ntimes-i_nt>=0);
                assert(ntimes-i_nt>=0 && ntimes-i_nt<_KIEDGE);
                assert(iy-miniy>=0 && iy-miniy<_KIY);
                assert(ix-minix>=0 && ix-minix<_KIX);    
                segment->edges.group_t[0].t[ntimes-i_nt]=grasp_output_tile_pixel_information_time(results,times[i_nt-1],ix,iy);
                segment->edges.group_t[0].x[ntimes-i_nt][iy-miniy][ix-minix]=grasp_output_tile_pixel_information_longitude(results,times[i_nt-1],ix,iy);
                segment->edges.group_t[0].y[ntimes-i_nt][iy-miniy][ix-minix]=grasp_output_tile_pixel_information_latitude(results,times[i_nt-1],ix,iy);
                segment->edges.group_t[0].icloud[ntimes-i_nt][iy-miniy][ix-minix]=1;
                segment->edges.group_t[0].out_t[ntimes-i_nt]=grasp_output_tile_pixel_information_out_t(results,times[i_nt-1],ix,iy);
                segment->edges.group_t[0].out_x[ntimes-i_nt][iy-miniy][ix-minix]=grasp_output_tile_pixel_information_out_x(results,times[i_nt-1],ix,iy);
                segment->edges.group_t[0].out_y[ntimes-i_nt][iy-miniy][ix-minix]=grasp_output_tile_pixel_information_out_y(results,times[i_nt-1],ix,iy);                
                for (ipar = 0; ipar < npar; ipar++) { 
                    assert(ipar>=0 && ipar<_KPARS);
                    segment->edges.group_t[0].AP[ntimes-i_nt][iy-miniy][ix-minix][ipar]=grasp_output_tile_retrieval_par_parameters(results,times[i_nt-1],ix,iy,ipar);
                }  
            } 
            if(segment->edges.group_t[0].nt<ntimes){
                segment->edges.group_t[0].nt=ntimes;
            }
        }
    }    
}

void grasp_input_edges_find_T_after(grasp_settings *settings, grasp_segment_t *segment, grasp_tile_dimensions_t *tile_dimensions, grasp_results_t *results, int icol, int irow, int itime){
    int maxiy,miniy,maxix,minix; // range of loops
    int it, iy, ix, i_nt, t_gotten, ipar, ipix, npar;
    
    npar=settings->retrieval.KNSING;
        
    miniy=(irow*settings->input.segment_ny);
    maxiy=(irow*settings->input.segment_ny)+settings->input.segment_ny;
    if(maxiy>tile_dimensions->tile_ny){
        maxiy=tile_dimensions->tile_ny;
    }   
    if(miniy>maxiy){
        miniy=maxiy;
    }
    
    minix=(icol*settings->input.segment_nx);
    maxix=(icol*settings->input.segment_nx)+settings->input.segment_nx;
    if(maxix>tile_dimensions->tile_nx){
        maxix=tile_dimensions->tile_nx;
    }  
    if(minix>maxix){
        minix=maxix;
    }    
       
    segment->edges.group_t[1].nt=0;
    segment->edges.group_t[1].nx=segment->sdata.nx;
    segment->edges.group_t[1].ny=segment->sdata.ny;   
    
    for (iy = miniy; iy < maxiy; iy++) {
        for (ix = minix; ix < maxix; ix++) {
            it=(itime*settings->input.segment_nt);
            for (i_nt = 0; i_nt < settings->retrieval.edges.nt; i_nt++) {
                t_gotten=grasp_input_edges_search_in_T(segment, results, tile_dimensions, it, iy, ix, +1);
                if(t_gotten<0){
                    break;
                }
                it=t_gotten;
                assert(it>=0 && it <tile_dimensions->tile_nt);
                ipix=index3D(it,ix,iy, tile_dimensions->tile_nx,tile_dimensions->tile_ny);
                assert(results->tile_result_map[ipix]!=NULL);  
                assert(i_nt>=0 && i_nt<_KIEDGE);
                assert(iy-miniy>=0 && iy-miniy<_KIY);
                assert(ix-minix>=0 && ix-minix<_KIX); 
                segment->edges.group_x[1].t[i_nt]=grasp_output_tile_pixel_information_time(results,it,ix,iy);
                segment->edges.group_x[1].x[i_nt][iy-miniy][ix-minix]=grasp_output_tile_pixel_information_longitude(results,it,ix,iy);
                segment->edges.group_x[1].y[i_nt][iy-miniy][ix-minix]=grasp_output_tile_pixel_information_latitude(results,it,ix,iy);
                segment->edges.group_x[1].icloud[i_nt][iy-miniy][ix-minix]=1;
                segment->edges.group_x[1].out_t[i_nt]=grasp_output_tile_pixel_information_out_t(results,it,ix,iy);
                segment->edges.group_x[1].out_x[i_nt][iy-miniy][ix-minix]=grasp_output_tile_pixel_information_out_x(results,it,ix,iy);
                segment->edges.group_x[1].out_y[i_nt][iy-miniy][ix-minix]=grasp_output_tile_pixel_information_out_y(results,it,ix,iy);                
                for (ipar = 0; ipar < npar; ipar++) { 
                    assert(ipar>=0 && ipar<_KPARS);
                    segment->edges.group_t[1].AP[i_nt][iy-miniy][ix-minix][ipar]=grasp_output_tile_retrieval_par_parameters(results,it,ix,iy,ipar);
                }       
            }
            if(i_nt>segment->edges.group_t[1].nt){
                segment->edges.group_t[1].nt=i_nt;
            }
        }
    }       
}


void grasp_input_edges_find(grasp_settings *settings, grasp_segment_t *segment, grasp_tile_dimensions_t *tile_dimensions, grasp_results_t *results, int icol, int irow, int itime) {    
    assert(icol >=0 && icol <=tile_dimensions->segment_ncols &&
           irow >=0 && irow <=tile_dimensions->segment_nrows &&
           itime>=0 && itime<=tile_dimensions->segment_ntimes);
    
    grasp_input_edges_initialization(&segment->edges);
    
    if(segment->sdata.npixels>0){
        // for x before: segment->edges.group_x[0]        
        grasp_input_edges_find_X_before(settings, segment, tile_dimensions, results, icol, irow, itime);
        
        // for x after: segment->edges.group_x[1]
        grasp_input_edges_find_X_after(settings, segment, tile_dimensions, results, icol, irow, itime);
        
        // for y before: segment->edges.group_y[0]
        grasp_input_edges_find_Y_before(settings, segment, tile_dimensions, results, icol, irow, itime);
        
        // for y after: segment->edges.group_y[1]
        grasp_input_edges_find_Y_after(settings, segment, tile_dimensions, results, icol, irow, itime);
        
        // For T before: segment->edges.group_t[0]
        grasp_input_edges_find_T_before(settings, segment, tile_dimensions, results, icol, irow, itime);
        
        // For T after: segment->edges.group_t[1]
        grasp_input_edges_find_T_after(settings, segment, tile_dimensions, results, icol, irow, itime);     
    }
}

/* provides the number of pixels in the group of pixels corresponding to 
 * the it index time (between 1 and NT).
 */
int grasp_input_get_num_pixels_at_it(const grasp_segment_t *segment, int it) {
  int ipixel;
  int npixels;
  int npixels_it;
  
  assert(segment != NULL);
  npixels = segment->sdata.npixels;

  assert(1 <= it && it <= segment->sdata.nt);

  npixels_it = 0;
  
  for (ipixel = 0 ; ipixel < npixels ; ipixel++) {
    const pixel_t *segment_pixel = &segment->sdata.pixel[ipixel];

    if (segment_pixel->it == it) {
      npixels_it++;
      /* one could leave the loop as soon as nx*ny has been reached (for increasing efficiency)
       * but it could mask bugs (that would be caught by the following assertions). I prefer to keep it
       * simple and stupid.
       */
    }
  }
  
  assert(npixels_it <= segment->sdata.nx*segment->sdata.ny);
  assert(npixels_it <= npixels);
  return npixels_it;
}

static bool is_valid_value(float value, float missing_value) {
  return ! ( isnan(value) || value == missing_value );
}

/* removes missing values from each pixel in the given segment (e.g.
 * measurements may be missing for some geometries in some channels).
 * The retrieval expects contiguous arrays only (without missing values)
 */
void grasp_input_clean_segment(grasp_segment_t *segment, float missing_value){
  int ipixel, iwl, ip;
  int ivm;
  int num_good_meas[_KIP];
  float value, *array;
  data_wl_t owl;
  sensor_data_t *sdata;
    
  sdata=&segment->sdata;
    
  for(ipixel=0;ipixel<sdata->npixels;ipixel++){
    for(iwl=0;iwl<sdata->pixel[ipixel].nwl;iwl++){
      memcpy(&owl,&sdata->pixel[ipixel].meas[iwl],sizeof(data_wl_t));
      memset(&sdata->pixel[ipixel].meas[iwl],0,sizeof(data_wl_t));
	
      grasp_segment_set_pixel_meas_wl(sdata,ipixel,iwl,owl.wl);
      grasp_segment_set_pixel_meas_ind_wl(sdata,ipixel,iwl,owl.ind_wl);
      grasp_segment_set_pixel_meas_sza(sdata,ipixel,iwl,owl.sza);
      grasp_segment_set_pixel_meas_nsurf(sdata,ipixel,iwl,owl.nsurf);
      grasp_segment_set_pixel_meas_gaspar(sdata,ipixel,iwl,owl.gaspar);
      grasp_segment_set_pixel_meas_nip(sdata,ipixel,iwl,owl.nip);
        
      for (ip=0;ip<sdata->pixel[ipixel].meas[iwl].nip;ip++) {                 
	num_good_meas[ip] = 0;
	
	memcpy(&sdata->pixel[ipixel].meas[iwl].groundpar[num_good_meas[ip]],owl.groundpar,sizeof(float)*_KSURF );
        grasp_segment_set_pixel_meas_ifcov(sdata,ipixel,iwl,ip,owl.ifcov[ip]);
        grasp_segment_set_pixel_meas_ifmp(sdata,ipixel,iwl,ip,owl.ifmp[ip]);
        grasp_segment_set_pixel_meas_meas_type(sdata,ipixel,iwl,ip,owl.meas_type[ip]);
        
        // In order to disable this assertion assert(ivalidmeas<sdata->pixel[ipixel].meas[iwl].nbvm[ip]); we set up initially nbvm as maximum and then we will set up as specific value
        grasp_segment_set_pixel_meas_nbvm(sdata,ipixel,iwl,ip,_NBVM);
        
	for(ivm=0 ; ivm<owl.nbvm[ip] ; ivm++) {
	  int ivm_valid;

	  switch(owl.meas_type[ip]){
	  case (MEAS_TYPE_TOD):
	  case (MEAS_TYPE_AOD):
	    {
	      array=owl.tau;
	      value = array[ivm];
	      if(is_valid_value(value, missing_value)) {
		ivm_valid = num_good_meas[ip];
		
		assert(0 <= ivm && ivm < _NBVM);
		assert(0 <= ivm_valid && ivm_valid < _NBVM);
		assert(is_valid_value(owl.thetav[ip][ivm], missing_value));
		assert(is_valid_value(owl.phi[ip][ivm], missing_value));
                grasp_segment_set_pixel_meas_thetav(sdata,ipixel,iwl,ip,ivm_valid,owl.thetav[ip][ivm]);
                grasp_segment_set_pixel_meas_phi(sdata,ipixel,iwl,ip,ivm_valid,owl.phi[ip][ivm]);               
                grasp_segment_set_pixel_meas_tau(sdata,ipixel,iwl,ivm_valid,value);
		num_good_meas[ip]++;
	      }
	    }
	    break;
	  case (MEAS_TYPE_P11):
	    {
	      array=owl.p11;
	      value = array[ivm];
	      if(is_valid_value(value, missing_value)) {
		ivm_valid = num_good_meas[ip];
		
		assert(0 <= ivm && ivm < _NBVM);
		assert(0 <= ivm_valid && ivm_valid < _NBVM);
		assert(is_valid_value(owl.thetav[ip][ivm], missing_value));
		assert(is_valid_value(owl.phi[ip][ivm], missing_value));
                grasp_segment_set_pixel_meas_thetav(sdata,ipixel,iwl,ip,ivm_valid,owl.thetav[ip][ivm]);
                grasp_segment_set_pixel_meas_phi(sdata,ipixel,iwl,ip,ivm_valid,owl.phi[ip][ivm]);               
                grasp_segment_set_pixel_meas_p11(sdata,ipixel,iwl,ivm_valid,value);
		num_good_meas[ip]++;
	      }
	    }
	    break;
	  case (MEAS_TYPE_P12):
	    {
	      array=owl.p12;
	      value = array[ivm];
	      if(is_valid_value(value, missing_value)) {
		ivm_valid = num_good_meas[ip];
		
		assert(0 <= ivm && ivm < _NBVM);
		assert(0 <= ivm_valid && ivm_valid < _NBVM);
		assert(is_valid_value(owl.thetav[ip][ivm], missing_value));
		assert(is_valid_value(owl.phi[ip][ivm], missing_value));
                grasp_segment_set_pixel_meas_thetav(sdata,ipixel,iwl,ip,ivm_valid,owl.thetav[ip][ivm]);
                grasp_segment_set_pixel_meas_phi(sdata,ipixel,iwl,ip,ivm_valid,owl.phi[ip][ivm]);               
                grasp_segment_set_pixel_meas_p12(sdata,ipixel,iwl,ivm_valid,value);
		num_good_meas[ip]++;
	      }
	    }
	    break;
	  case (MEAS_TYPE_P22):
	    {
	      array=owl.p22;
	      value = array[ivm];
	      if(is_valid_value(value, missing_value)) {
		ivm_valid = num_good_meas[ip];
		
		assert(0 <= ivm && ivm < _NBVM);
		assert(0 <= ivm_valid && ivm_valid < _NBVM);
		assert(is_valid_value(owl.thetav[ip][ivm], missing_value));
		assert(is_valid_value(owl.phi[ip][ivm], missing_value));
                grasp_segment_set_pixel_meas_thetav(sdata,ipixel,iwl,ip,ivm_valid,owl.thetav[ip][ivm]);
                grasp_segment_set_pixel_meas_phi(sdata,ipixel,iwl,ip,ivm_valid,owl.phi[ip][ivm]);               
                grasp_segment_set_pixel_meas_p22(sdata,ipixel,iwl,ivm_valid,value);
		num_good_meas[ip]++;
	      }
	    }
	    break;
	  case (MEAS_TYPE_P33):
	    {
	      array=owl.p33;
	      value = array[ivm];
	      if(is_valid_value(value, missing_value)) {
		ivm_valid = num_good_meas[ip];
		
		assert(0 <= ivm && ivm < _NBVM);
		assert(0 <= ivm_valid && ivm_valid < _NBVM);
		assert(is_valid_value(owl.thetav[ip][ivm], missing_value));
		assert(is_valid_value(owl.phi[ip][ivm], missing_value));
                grasp_segment_set_pixel_meas_thetav(sdata,ipixel,iwl,ip,ivm_valid,owl.thetav[ip][ivm]);
                grasp_segment_set_pixel_meas_phi(sdata,ipixel,iwl,ip,ivm_valid,owl.phi[ip][ivm]);               
                grasp_segment_set_pixel_meas_p33(sdata,ipixel,iwl,ivm_valid,value);
		num_good_meas[ip]++;
	      }
	    }
	    break;
	  case (MEAS_TYPE_P34):
	    {
	      array=owl.p34;
	      value = array[ivm];
	      if(is_valid_value(value, missing_value)) {
		ivm_valid = num_good_meas[ip];
		
		assert(0 <= ivm && ivm < _NBVM);
		assert(0 <= ivm_valid && ivm_valid < _NBVM);
		assert(is_valid_value(owl.thetav[ip][ivm], missing_value));
		assert(is_valid_value(owl.phi[ip][ivm], missing_value));
                grasp_segment_set_pixel_meas_thetav(sdata,ipixel,iwl,ip,ivm_valid,owl.thetav[ip][ivm]);
                grasp_segment_set_pixel_meas_phi(sdata,ipixel,iwl,ip,ivm_valid,owl.phi[ip][ivm]);               
                grasp_segment_set_pixel_meas_p34(sdata,ipixel,iwl,ivm_valid,value);
		num_good_meas[ip]++;
	      }
	    }
	    break;
	  case (MEAS_TYPE_P44):
	    {
	      array=owl.p44;
	      value = array[ivm];
	      if(is_valid_value(value, missing_value)) {
		ivm_valid = num_good_meas[ip];
		
		assert(0 <= ivm && ivm < _NBVM);
		assert(0 <= ivm_valid && ivm_valid < _NBVM);
		assert(is_valid_value(owl.thetav[ip][ivm], missing_value));
		assert(is_valid_value(owl.phi[ip][ivm], missing_value));
                grasp_segment_set_pixel_meas_thetav(sdata,ipixel,iwl,ip,ivm_valid,owl.thetav[ip][ivm]);
                grasp_segment_set_pixel_meas_phi(sdata,ipixel,iwl,ip,ivm_valid,owl.phi[ip][ivm]);               
                grasp_segment_set_pixel_meas_p44(sdata,ipixel,iwl,ivm_valid,value);
		num_good_meas[ip]++;
	      }
	    }
	    break;
	  case (MEAS_TYPE_LS):
	    {
	      array=owl.ls;
	      value = array[ivm];
	      if(is_valid_value(value, missing_value)) {
		ivm_valid = num_good_meas[ip];
		
		assert(0 <= ivm && ivm < _KVERTM);
		assert(0 <= ivm_valid && ivm_valid < _KVERTM);
		assert(is_valid_value(owl.mprof[ip][ivm], missing_value));
                grasp_segment_set_pixel_meas_mprof(sdata,ipixel,iwl,ip,ivm_valid,owl.mprof[ip][ivm]);
                grasp_segment_set_pixel_meas_ls(sdata,ipixel,iwl,ivm_valid,value);
		num_good_meas[ip]++;
	      }
	    }
	    break;
	  case (MEAS_TYPE_DP):
	    {
	      array=owl.dp;
	      value = array[ivm];
	      if(is_valid_value(value, missing_value)) {
		ivm_valid = num_good_meas[ip];
		
		assert(0 <= ivm && ivm < _KVERTM);
		assert(0 <= ivm_valid && ivm_valid < _KVERTM);
		assert(is_valid_value(owl.mprof[ip][ivm], missing_value));
                grasp_segment_set_pixel_meas_mprof(sdata,ipixel,iwl,ip,ivm_valid,owl.mprof[ip][ivm]);
                grasp_segment_set_pixel_meas_dp(sdata,ipixel,iwl,ivm_valid,value);
		num_good_meas[ip]++;
	      }
	    }
	    break;
	  case (MEAS_TYPE_RL):
	    {
	      array=owl.rl;
	      value = array[ivm];
	      if(is_valid_value(value, missing_value)) {
		ivm_valid = num_good_meas[ip];
		
		assert(0 <= ivm && ivm < _KVERTM);
		assert(0 <= ivm_valid && ivm_valid < _KVERTM);
		assert(is_valid_value(owl.mprof[ip][ivm], missing_value));
                grasp_segment_set_pixel_meas_mprof(sdata,ipixel,iwl,ip,ivm_valid,owl.mprof[ip][ivm]);
                grasp_segment_set_pixel_meas_rl(sdata,ipixel,iwl,ivm_valid,value);
		num_good_meas[ip]++;
	      }
	    }
	    break;
	  case (MEAS_TYPE_I):
	    {
	      array=owl.i;
	      value = array[ivm];
	      if(is_valid_value(value, missing_value)) {
		ivm_valid = num_good_meas[ip];
		  
		assert(0 <= ivm && ivm < _NBVM);
		assert(0 <= ivm_valid && ivm_valid < _NBVM);
		assert(is_valid_value(owl.thetav[ip][ivm], missing_value));
		assert(is_valid_value(owl.phi[ip][ivm], missing_value));
                grasp_segment_set_pixel_meas_thetav(sdata,ipixel,iwl,ip,ivm_valid,owl.thetav[ip][ivm]);
                grasp_segment_set_pixel_meas_phi(sdata,ipixel,iwl,ip,ivm_valid,owl.phi[ip][ivm]);               
                grasp_segment_set_pixel_meas_i(sdata,ipixel,iwl,ivm_valid,value);
		num_good_meas[ip]++;
	      }
	    }
	    break;
	  case (MEAS_TYPE_Q):
	    {
	      array=owl.q;
	      value = array[ivm];
	      if(is_valid_value(value, missing_value)) {
		ivm_valid = num_good_meas[ip];

		assert(0 <= ivm && ivm < _NBVM);
		assert(0 <= ivm_valid && ivm_valid < _NBVM);
		assert(is_valid_value(owl.thetav[ip][ivm], missing_value));
		assert(is_valid_value(owl.phi[ip][ivm], missing_value));
                grasp_segment_set_pixel_meas_thetav(sdata,ipixel,iwl,ip,ivm_valid,owl.thetav[ip][ivm]);
                grasp_segment_set_pixel_meas_phi(sdata,ipixel,iwl,ip,ivm_valid,owl.phi[ip][ivm]);               
                grasp_segment_set_pixel_meas_q(sdata,ipixel,iwl,ivm_valid,value);
		num_good_meas[ip]++;
	      }
	    }
	    break;
	  case (MEAS_TYPE_U):
	    {
	      array=owl.u;
	      value = array[ivm];
	      if(is_valid_value(value, missing_value)) {
		ivm_valid = num_good_meas[ip];
		
		assert(0 <= ivm && ivm < _NBVM);
		assert(0 <= ivm_valid && ivm_valid < _NBVM);
		assert(is_valid_value(owl.thetav[ip][ivm], missing_value));
		assert(is_valid_value(owl.phi[ip][ivm], missing_value));
                grasp_segment_set_pixel_meas_thetav(sdata,ipixel,iwl,ip,ivm_valid,owl.thetav[ip][ivm]);
                grasp_segment_set_pixel_meas_phi(sdata,ipixel,iwl,ip,ivm_valid,owl.phi[ip][ivm]);               
                grasp_segment_set_pixel_meas_u(sdata,ipixel,iwl,ivm_valid,value);
		num_good_meas[ip]++;
	      }
	    }
	    break;
	  case (MEAS_TYPE_P):
	    {
	      array=owl.p;
	      value = array[ivm];
	      if(is_valid_value(value, missing_value)) {
		ivm_valid = num_good_meas[ip];

		assert(0 <= ivm && ivm < _NBVM);
		assert(0 <= ivm_valid && ivm_valid < _NBVM);
		assert(is_valid_value(owl.thetav[ip][ivm], missing_value));
		assert(is_valid_value(owl.phi[ip][ivm], missing_value));
                grasp_segment_set_pixel_meas_thetav(sdata,ipixel,iwl,ip,ivm_valid,owl.thetav[ip][ivm]);
                grasp_segment_set_pixel_meas_phi(sdata,ipixel,iwl,ip,ivm_valid,owl.phi[ip][ivm]);               
                grasp_segment_set_pixel_meas_p(sdata,ipixel,iwl,ivm_valid,value);
		num_good_meas[ip]++;
	      }
	    }
	    break;
              
	  default:
	    fprintf(stderr, "%s:%d: unexpected value for meas_type: %d\n", __FILE__, __LINE__, owl.meas_type[ip]);
	    abort();
	  } /* switch */

	  ivm_valid = ivm;
	  assert(0 <= ivm && ivm < _KNBVM);
	  assert(0 <= ivm_valid && ivm_valid < _KNBVM);
	  assert(is_valid_value(owl.cmtrx[ip][ivm], missing_value));
          grasp_segment_set_pixel_meas_cmtrx(sdata,ipixel,iwl,ip,ivm_valid,owl.cmtrx[ip][ivm]);
	} /* for (ivm) */
	grasp_segment_set_pixel_meas_nbvm(sdata,ipixel,iwl,ip,num_good_meas[ip]);
      } /* ip */
    } /* iwl */
  } /* ipixel */
}

/*
 * auxiliary function used by grasp_input_dump_segment and grasp_input_print_segment
 */
static void convert_segment_pixel_to_sdata_pixel(SDATA_PIXEL *sdata_pixel, const pixel_t *segment_pixel) {
  int iwl;
  int nsurf;

  sdata_pixel->ix = segment_pixel->ix;
  sdata_pixel->iy = segment_pixel->iy;
  sdata_pixel->it = segment_pixel->it;
  sdata_pixel->unix_time = segment_pixel->t;
  sdata_pixel->index_of_day = 0;
  sdata_pixel->cloud_flag = segment_pixel->cloudy;
  sdata_pixel->irow = segment_pixel->irow;
  sdata_pixel->icol = segment_pixel->icol;
  sdata_pixel->latitude = segment_pixel->y;
  sdata_pixel->longitude = segment_pixel->x;
  sdata_pixel->hgr = segment_pixel->masl;
  sdata_pixel->hobs = segment_pixel->hobs;
  sdata_pixel->land_percent = segment_pixel->land_percent;
  sdata_pixel->num_wavelengths = segment_pixel->nwl;
  nsurf = sdata_pixel->nsurf = segment_pixel->meas[0].nsurf;
  sdata_pixel->ifgas = segment_pixel->ifgas;

  for (iwl = 0 ; iwl < segment_pixel->nwl ; iwl++) {
    int ip, isurf, ivm;
    int nip;
    const data_wl_t *meas = &segment_pixel->meas[iwl];
    
    sdata_pixel->wavelengths[iwl] = meas->wl;
    nip = sdata_pixel->num_meas_types[iwl] = meas->nip;
    sdata_pixel->thetas[iwl] = meas->sza;
    
    for (ip = 0 ; ip < nip ; ip++) {
      int nbvm;
      int meas_type;
      
      meas_type = sdata_pixel->meas_types[iwl][ip] = meas->meas_type[ip];
      
      if (meas_type == MEAS_TYPE_LS || meas_type == MEAS_TYPE_DP || meas_type == MEAS_TYPE_RL) {
	sdata_pixel->num_hvp = meas->nbvm[0];
      }

      sdata_pixel->num_valid_meas[iwl][ip] = meas->nbvm[ip];
      
      sdata_pixel->ifcov[iwl][ip] = meas->ifcov[ip];
      sdata_pixel->ifmp[iwl][ip] = meas->ifmp[ip];
      
      nbvm =  meas->nbvm[ip];
      for (ivm = 0 ; ivm < nbvm ; ivm++) {
	if (meas_type == MEAS_TYPE_LS || meas_type == MEAS_TYPE_DP || meas_type == MEAS_TYPE_RL) {
	  assert(ivm < _KVERTM);
	  sdata_pixel->hvp[ivm] = segment_pixel->hvp[ivm];
	}
	else {
	  assert(ivm < _NBVM);
	  sdata_pixel->thetav[iwl][ip][ivm] = meas->thetav[ip][ivm];
	}
	sdata_pixel->dphi[iwl][ip][ivm] = meas->phi[ip][ivm];
        if (meas->ifcov[ip] == 1) {
            sdata_pixel->cmtrx[iwl][ip][ivm] = meas->cmtrx[ip][ivm];
        }
        if (meas->ifmp[ip] == 1) { // one assigns matrix profiles only if the presence flag ifmp is set
	  sdata_pixel->mprof[iwl][ip][ivm] = meas->mprof[ip][ivm];
        }

	switch (meas_type) {
	case MEAS_TYPE_TOD:
	  sdata_pixel->tau[iwl][ivm] = meas->tau[ivm];
	  break;
	case MEAS_TYPE_AOD:
	  sdata_pixel->tau[iwl][ivm] = meas->tau[ivm];
	  break;          
	case MEAS_TYPE_P11:
	  sdata_pixel->p11[iwl][ivm] = meas->p11[ivm];
	  break;
	case MEAS_TYPE_P12:
	  sdata_pixel->p12[iwl][ivm] = meas->p12[ivm];
	  break;
	case MEAS_TYPE_P22:
	  sdata_pixel->p22[iwl][ivm] = meas->p22[ivm];
	  break;
	case MEAS_TYPE_P33:
	  sdata_pixel->p33[iwl][ivm] = meas->p33[ivm];
	  break;
	case MEAS_TYPE_P34:
	  sdata_pixel->p34[iwl][ivm] = meas->p34[ivm];
	  break;
	case MEAS_TYPE_P44:
	  sdata_pixel->p44[iwl][ivm] = meas->p44[ivm];
	  break;
	case MEAS_TYPE_LS:
	  sdata_pixel->LS[iwl][ivm] = meas->ls[ivm];
	  break;
	case MEAS_TYPE_DP:
	  sdata_pixel->DP[iwl][ivm] = meas->dp[ivm];
	  break;
	case MEAS_TYPE_RL:
	  sdata_pixel->RL[iwl][ivm] = meas->rl[ivm];
	  break;
	case MEAS_TYPE_I:
	  sdata_pixel->I[iwl][ivm] = meas->i[ivm];
	  break;
	case MEAS_TYPE_Q:
	  sdata_pixel->Q[iwl][ivm] = meas->q[ivm];
	  break;
	case MEAS_TYPE_U:
	  sdata_pixel->U[iwl][ivm] = meas->u[ivm];
	  break;
	case MEAS_TYPE_P:
	  sdata_pixel->P[iwl][ivm] = meas->p[ivm];
	  break;
	default:
	  fprintf(stderr, "%s:%d: fatal error: unexpected value for meas_type: %d\n", __FILE__, __LINE__, meas_type);
	  abort();
	}
      } /* ivm */

    } /* ip */

    for (isurf = 0 ; isurf < nsurf ; isurf++) {
      sdata_pixel->groundpar[iwl][isurf] = meas->groundpar[isurf];
    } /* isurf */

    sdata_pixel->gaspar[iwl] = meas->gaspar;
    
  } /* iwl */

}


void grasp_input_dump_iguess(FILE *output_stream, const grasp_settings *settings, const grasp_segment_t *segment){
    int ipixel;
    int ipar;

    assert(output_stream != NULL);
    assert(segment != NULL);   
    for (ipar=0;ipar<settings->retrieval.KNSING;ipar++){
        fprintf(output_stream,"%d",ipar+1);
        for (ipixel = 0 ; ipixel < segment->sdata.npixels ; ipixel++) {
            fprintf(output_stream,"\t%f",segment->iguess[ipixel][ipar]);
        }
        fprintf(output_stream,"\n");
    }
}

/* dumps a segment in the (quite unreadable) SDATA format */
void grasp_input_dump_segment(FILE *output_stream, const grasp_segment_t *segment) {
  SDATA_HEADER header;
  size_t ipixel;
  int previous_it;

  assert(output_stream != NULL);
  assert(segment != NULL);

  memset(&header, 0, sizeof(header));

#ifdef WARN_DRY
#warning "__SDATA_VERSION__ duplicated"
#endif   
  strcpy(header.format_version, "2.0");
  header.nx = segment->sdata.nx;
  header.ny = segment->sdata.ny;
  header.nrecords = segment->sdata.nt;

  sdata_dump_header(output_stream, &header);
  
  previous_it = 0;
  for (ipixel = 0 ; ipixel < segment->sdata.npixels ; ipixel++) {
    const pixel_t *segment_pixel = &segment->sdata.pixel[ipixel];
    SDATA_PIXEL sdata_pixel;
    SDATA_RECORD record;
    if (segment_pixel->it != previous_it) {      
      time_t timestamp;

      /* pixels are grouped by times. If the time increases, one displays a new record header (record == group of pixels).
       * The time index is supposed to increase of one unity at a time. If it's not the case,
       * the following assertion will stop the processing.
       */
      if (segment_pixel->it != previous_it + 1) {
        grasp_input_print_segment(stderr, "segment", segment);
        fprintf(stderr, "%s:%d: %s: previous_it: %d segment_pixel->it: %d\n", __FILE__, __LINE__, __FUNCTION__,
          previous_it, segment_pixel->it);
      }
      assert(segment_pixel->it == previous_it + 1);

      timestamp = segment_pixel->t;
      time_to_string_r(timestamp, "%FT%H:%M:%SZ" /* strftime specification*/,
		       sizeof(record.str_date_time), record.str_date_time);

      record.npixels = grasp_input_get_num_pixels_at_it(segment, segment_pixel->it);
      record.hobs = segment_pixel->hobs;
      record.nsurf = segment_pixel->meas[0].nsurf;
      record.ifgas = segment_pixel->ifgas;

      if (previous_it > 0) fprintf(output_stream, "\n");
      sdata_dump_record(output_stream, &record);
      previous_it = segment_pixel->it;
    }
    convert_segment_pixel_to_sdata_pixel(&sdata_pixel, segment_pixel);
    sdata_dump_pixel(output_stream, &sdata_pixel);
  }
}

/* prints a segment in a easy-to-read fashion (mainly for debugging) */
void grasp_input_print_segment(FILE *output_stream, const char *label, const grasp_segment_t *segment) {
  int ipixel;
  assert(output_stream != NULL);
  assert(segment != NULL);

  if (label != NULL && label[0] != '\0') {
    fprintf(output_stream, "%s ", label);
  }
  
  fprintf(output_stream, "segment->npixels: %d\n", segment->sdata.npixels);
  fprintf(output_stream, "segment->NX: %d\n", segment->sdata.nx);
  fprintf(output_stream, "segment->NY: %d\n", segment->sdata.ny);
  fprintf(output_stream, "segment->NT: %d\n", segment->sdata.nt);
  fprintf(output_stream, "\n");
  for (ipixel = 0 ; ipixel < segment->sdata.npixels ; ipixel++) {
    const pixel_t *segment_pixel = &segment->sdata.pixel[ipixel];
    SDATA_PIXEL sdata_pixel;
    char pixel_label[31 + 1];

    sprintf(pixel_label, "pixel %d", ipixel + 1);
    convert_segment_pixel_to_sdata_pixel(&sdata_pixel, segment_pixel);
    sdata_print_pixel(output_stream, pixel_label, &sdata_pixel);
    fprintf(output_stream, "\n");
  }
}


int grasp_input_expand_files(grasp_settings *settings, char ***files){
    int i,j,c;
    int totalfiles=0;    
    int nfiles_expanded[GRASP_INPUT_MAX_FILES];
    filepath_t *glob_expanded; 
    filepath_t *expanded_list[GRASP_INPUT_MAX_FILES];
    char **tmp_files;
    
    
    for (i = 0; i < settings->input.nfiles; i++) {
        nfiles_expanded[i]=find_files_from_pattern(settings->input.files[i],&glob_expanded);
        if(nfiles_expanded[i]<=0){
            switch(nfiles_expanded[i]){
                case -1:
                    fprintf(stderr, "WARNING: %s not matching. Pattern given to search files return 0 files.\n" , settings->input.files[i]);
                    break;
                case -2:
                    fprintf(stderr, "WARNING: %s glob aborted. Pattern given to search files is not valid.\n" , settings->input.files[i]);
                    break;
                case -3:
                    fprintf(stderr, "WARNING: %s causes fatal error. Pattern given to search files is not valid.\n" , settings->input.files[i]);
                    break;
                case -4:
                    fprintf(stderr, "WARNING: %s causes unexpected return value by glob library. Pattern given to search files is not valid.\n" , settings->input.files[i]);
                    break;
            }
        }else{
            totalfiles+=nfiles_expanded[i];
            expanded_list[i]=glob_expanded;
        }
    }

    tmp_files = (char **) trackmem_malloc(sizeof (char *)*totalfiles);
    assert(tmp_files!=NULL);
    
    c=0;
    for (i = 0; i < settings->input.nfiles; i++) {
        glob_expanded=expanded_list[i];
        for (j = 0; j < nfiles_expanded[i]; j++) {
            tmp_files[c] = (char *) trackmem_malloc(sizeof (char)*(strlen(glob_expanded[j])+1));
            assert(tmp_files[c]!=NULL);
            strcpy(tmp_files[c], glob_expanded[j]);
            c++;
        }   
        if(nfiles_expanded[i]>=0){
            trackmem_free(glob_expanded);
        }
    }


    *files=tmp_files;
    return c;
}

void grasp_input_deallocate_expanded_files(char **files, int nfiles){
    int i;
    
    for (i = 0; i < nfiles; i++) {
        trackmem_free(files[i]);        
    }
    
    trackmem_free(files);
}
