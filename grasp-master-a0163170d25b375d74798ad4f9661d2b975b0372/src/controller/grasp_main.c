/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

/* 
 * File:   main.c
 * Author: David Fuertes
 *
 * Created on 29 de septiembre de 2013, 22:27
 * Adapted for MPI by Fabrice Ducos, last change 8 october 2014
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "grasp_main.h"
#include "yamlsettings/yamlsettings.h"
#include "../input/grasp_input.h"
#include "../output/grasp_output.h"
#include "../settings/grasp_settings.h"
#include <grasp/utils.h>
#include "grasp_controller.h"
#include "grasp_mpi_engine.h"
#include "mo_grasp_controller.h"
#include "../global/grasp_runtime_information.h"

#ifdef USE_MPI
#include "mpi.h"
#endif

void main_sequential(int argc, char** argv) {
  grasp_settings *settings=NULL;
  grasp_processing_functions_t functions;
  grasp_tile_description_t tile_description;
  grasp_results_t results;
  grasp_output_stream *stream;
  benchmark_t benchmark;

  benchmark_start(&benchmark);
    
  // 1) Read settings and process controller options
  settings=grasp_controller_read_settings(argc,argv);
  
  // 2) Get tile information
  grasp_controller_initialize_inversion(settings, &tile_description, &functions, &results);
  
  // 3) Invert tile data
  grasp_controller_invert_tile(settings, &tile_description, &results, &functions); 
  
  // 4) Manage output
  grasp_controller_manage_tile(settings, &tile_description, &results, &functions);
  
  benchmark_stop(&benchmark);
  
  stream=grasp_controller_get_stream();
  if(grasp_output_tile_information_tile_npixels(&results)>0){
    gos_fprintf(stream, "Total Time: %d pixels processed in %f seconds (cpu time: %f). Average per pixel: %f (cpu time: %f)\n",        grasp_output_tile_information_tile_npixels(&results), benchmark.delta_ct,                                      benchmark.delta_ut,                                     benchmark.delta_ct/grasp_output_tile_information_tile_npixels(&results),                                       benchmark.delta_ut/grasp_output_tile_information_tile_npixels(&results));
    gos_fprintf(stream, "Algorithm Time: %d pixels processed in %f seconds (cpu time: %f). Average per pixel: %f (cpu time: %f)\n",    grasp_output_tile_information_tile_npixels(&results), grasp_controller_get_algorithm_ct(),                     grasp_controller_get_algorithm_ut(),                    grasp_controller_get_algorithm_ct()/grasp_output_tile_information_tile_npixels(&results),                      grasp_controller_get_algorithm_ut()/grasp_output_tile_information_tile_npixels(&results));
    gos_fprintf(stream, "Control Unit Time: %d pixels processed in %f seconds (cpu time: %f). Average per pixel: %f (cpu time: %f)\n", grasp_output_tile_information_tile_npixels(&results), benchmark.delta_ct-grasp_controller_get_algorithm_ct(),  benchmark.delta_ut-grasp_controller_get_algorithm_ut(), (benchmark.delta_ct-grasp_controller_get_algorithm_ct())/grasp_output_tile_information_tile_npixels(&results), (benchmark.delta_ct-grasp_controller_get_algorithm_ut())/grasp_output_tile_information_tile_npixels(&results));
  }
  if(grasp_controller_get_nerror_segment()>0){
    gos_fprintf(stream, "RETRIEVAL ERROR: During the retrieval process %d segments (%d pixels) were ignored because retrieval code returned an error\n", grasp_controller_get_nerror_segment(), grasp_controller_get_nerror_pixel());
  }
  
  // 5) Clean up 
  grasp_controller_clean_memory(settings,  &tile_description, &results, &functions);

}

#ifdef USE_MPI

void cleanup() {
  MPI_Finalize();
}

void main_mpi(int argc, char** argv) {
  grasp_settings *settings=NULL;
  grasp_processing_functions_t functions;
  grasp_tile_description_t tile_description;
  grasp_results_t results;
  grasp_output_stream *stream;
  benchmark_t benchmark;

  atexit(cleanup);
  MPI_Init(&argc, &argv);

  benchmark_start(&benchmark);
  
  // 1) Read settings and process controller options
  settings=grasp_controller_read_settings(argc,argv);
  
  /* in MPI mode, only the master has to load the tile */
  if (grasp_mpi_engine_is_master()) {
    // 2) Get Tile information
    grasp_controller_initialize_inversion(settings, &tile_description, &functions, &results);
      
    // 3) Invert tile data
    grasp_controller_invert_tile(settings, &tile_description, &results, &functions);

    // 4) Manage output
    grasp_controller_manage_tile(settings,&tile_description, &results, &functions);
  
  }
  else { /* no tile is loaded for the workers, but they have to start the processing too */
    grasp_controller_invert_tile(settings, NULL, NULL, NULL);
  } /* if (grasp_mpi_engine_is_master()) */
 
  benchmark_stop(&benchmark);

  if (grasp_mpi_engine_is_master()) {
        stream=grasp_controller_get_stream();
        if(grasp_output_tile_information_tile_npixels(&results)>0){
          gos_fprintf(stream,"Total Time: %d pixels processed in %f seconds (cpu time: %f). Average per pixel: %f (cpu time: %f)\n",        grasp_output_tile_information_tile_npixels(&results), benchmark.delta_ct,                                      benchmark.delta_ut,                                     benchmark.delta_ct/grasp_output_tile_information_tile_npixels(&results),                                       benchmark.delta_ut/grasp_output_tile_information_tile_npixels(&results));
          gos_fprintf(stream,"Algorithm Time (sum of each node): %d pixels processed in %f seconds (cpu time: %f). Average per pixel: %f (cpu time: %f)\n",    grasp_output_tile_information_tile_npixels(&results), grasp_controller_get_algorithm_ct(),                     grasp_controller_get_algorithm_ut(),                    grasp_controller_get_algorithm_ct()/grasp_output_tile_information_tile_npixels(&results),                      grasp_controller_get_algorithm_ut()/grasp_output_tile_information_tile_npixels(&results));
          //gos_fprintf(stream,"Control Unit Time: %d pixels processed in %f seconds (cpu time: %f). Average per pixel: %f (cpu time: %f)\n", grasp_output_tile_information_tile_npixels(&results), benchmark.delta_ct-grasp_controller_get_algorithm_ct(),  benchmark.delta_ut-grasp_controller_get_algorithm_ut(), (benchmark.delta_ct-grasp_controller_get_algorithm_ct())/grasp_output_tile_information_tile_npixels(&results), (benchmark.delta_ct-grasp_controller_get_algorithm_ut())/grasp_output_tile_information_tile_npixels(&results));
        }    
        if(grasp_controller_get_nerror_segment()>0){
          gos_fprintf(stream, "RETRIEVAL ERROR: During the retrieval process %d segments (%d pixels) were ignored because retrieval code returned an error\n", grasp_controller_get_nerror_segment(), grasp_controller_get_nerror_pixel());
        }
  } /* if (grasp_mpi_engine_is_master()) */

  // 5) Clean up
  grasp_controller_clean_memory(settings, &tile_description, &results, &functions);   
}
#endif /* USE_MPI */

int main(int argc, char** argv) {

#ifdef USE_MPI
  main_mpi(argc, argv);
#else
  main_sequential(argc, argv);
#endif

    return 0;
}


