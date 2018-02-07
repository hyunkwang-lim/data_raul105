/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

#include "grasp_output_segment_classic_plot.h"
#include <string.h>
#include <assert.h>
#include <stdlib.h>

grasp_output_segment_function_t grasp_output_segment_function_classic_plot(){
    grasp_output_segment_function_t x;
    
    x.init=grasp_output_segment_function_classic_plot_init;
    x.function=grasp_output_segment_function_classic_plot_process;
    x.close=grasp_output_segment_function_classic_plot_close;
    
    return x;
}

int grasp_output_segment_function_classic_plot_init(grasp_settings *settings, grasp_tile_description_t *input_information){
    return 0;
}

int grasp_output_segment_function_classic_plot_close(void){
    return 0;
}

int grasp_output_segment_function_classic_plot_process(grasp_output_stream *stream, grasp_settings *settings, grasp_segment_t *segment, output_segment_general *output, grasp_tile_description_t *tile_description,int icol,int irow,int itime){
    int i;
    char *filename;
    
    if(segment->sdata.npixels<=0){
        fprintf(stderr, "The segment in position time=%d, x=%d and y=%d has no pixels. Results with classic_plot segment output function can not be printed.\n", itime, icol, irow );
    }else{
        // Initialize stream and get filename
        filename=grasp_output_stream_filename(stream, settings, segment, output, &tile_description->dimensions, icol, irow, itime); 
        strcpy(settings->retrieval.plotting_output_file, filename);
        settings->retrieval.plotting_output_file[strlen(settings->retrieval.plotting_output_file)]=' ';

        // Check if it is writable
        
        if(grasp_output_stream_writable_file(stream)==false){
            printf("ERROR: classic_plot output segment function needs a valid file for being written. Current stream is not valid. Process aborted\n");
            exit(-1);
        }

        // Call plot function
        grasp_output_segment_print_classic_plot_f90(&settings->retrieval, &segment->sdata, output);


        // Clean plotting name. If you don't clean the name you can affect to next retrieval.
        for (i = 0; i < _GBL_FILE_PATH_LEN; i++) {
            settings->retrieval.plotting_output_file[i]=' ';
        }        
    }    
    
    return 0;
}


grasp_settings_parameter_array *grasp_output_segment_function_settings_classic_plot(grasp_settings *settings){
    return NULL;
}
