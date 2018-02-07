/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

#include "grasp_output_segment_classic.h"
#include "../../grasp_output_stream.h"
#include <assert.h>

grasp_output_segment_function_t grasp_output_segment_function_classic(){
    grasp_output_segment_function_t x;
    
    x.init=grasp_output_segment_function_classic_init;
    x.function=grasp_output_segment_function_classic_process;
    x.close=grasp_output_segment_function_classic_close;
    
    return x;
}

int grasp_output_segment_function_classic_init(grasp_settings *settings, grasp_tile_description_t *input_information){
    return 0;
}

int grasp_output_segment_function_classic_close(void){
    return 0;
}

int grasp_output_segment_function_classic_process(grasp_output_stream *stream, grasp_settings *settings, grasp_segment_t *segment, output_segment_general *output, grasp_tile_description_t *tile_description,int icol,int irow,int itime){    
    char *filename=NULL;
    if(segment->sdata.npixels<=0){
        fprintf(stderr, "The segment in position time=%d, x=%d and y=%d has no pixels. Results with classic segment output function can not be printed.\n", itime, icol, irow );
    }else{
        if(stream->screen==true){
            grasp_output_segment_print_classic_f90_screen(&settings->retrieval, &segment->sdata, output);
        }else{
            filename=grasp_output_stream_filename(stream,settings,segment,output, &tile_description->dimensions ,icol,irow,itime);
            
            grasp_output_segment_print_classic_f90_file(filename, &settings->retrieval, &segment->sdata, output);
            
            free(filename);
        }
    }
    
    return 0;
}

grasp_settings_parameter_array *grasp_output_segment_function_settings_classic(grasp_settings *settings){
    return NULL;
}
