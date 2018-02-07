/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

/**
 * @file grasp_controller.h
 * @author David Fuertes
 * @date 14 Oct 2013
 * @brief Functions needed by the controller to organize the workflow
 */

#ifndef GRASP_CONTROLLER_H
#define	GRASP_CONTROLLER_H

#ifdef	__cplusplus
extern "C" {
#endif

#include "yamlsettings/yamlsettings.h"    // used by retr_input
#include "../input/grasp_input.h"          // used by pixel_t
#include "../output/grasp_output.h"        // used by segment_retr_output    
#include "../output/grasp_output_stream.h"    
    
extern grasp_output_stream controller_stream;
    
/**
 * Return controller stream to be use in other submodules. This stream is opened
 * after reading the settings and it is close at the end of the workflow
 * @return memory direction of controller stream
 */
grasp_output_stream *grasp_controller_get_stream();

/**
 * Return real time in seconds used by retrieval algorithm (computed total time taken by all retrievals) 
 * @return real time in seconds used by retrieval algorithm
 */
double grasp_controller_get_algorithm_ut();

/**
 * Return user time in seconds used by retrieval algorithm (computed total time taken by all retrievals)
 * @return User time in seconds used by retrieval algorithm
 */
double grasp_controller_get_algorithm_ct();

/**
 * Return the number of segmenter that had to be ignored in the process because a retrieval error
 * @return number of segments
 */
int grasp_controller_get_nerror_segment();

/**
 * Return the number of pixels that were not processed in the process because a retrieval error
 * @return number of pixels
 */
int grasp_controller_get_nerror_pixel();

/**
 * Initilize, read and return settings from argc,argv parameters using settings module.
 * @param argc count of total input arguments
 * @param argv value of input arguments from main call (command line)
 * @return Settings structure allocated and filled 
 */
grasp_settings *grasp_controller_read_settings(int argc, char** argv);

/**
 * Function for initializing framework extension functions. It call init method from extensions.
 * @param settings current settings structure
 * @param tile_description current tile description
 * @param functions extension functions structure which will be initialize and updated
 * @return 0 if the process finished OK
 */
int grasp_controller_initialize_functions(grasp_settings *settings, grasp_tile_description_t *tile_description, grasp_processing_functions_t *functions);

/**
 * Initialize a tile from a settings. Return 0 is everything was ok
 * @param settings Current settings structure filled
 * @param tile_description Tile description which will be updated
 * @param functions extension functions structure which will be initialize and updated
 * @param results grasp result structure (output result tile) which will be initialized
 * @return 0 if the process finished OK
 */
int grasp_controller_initialize_inversion(grasp_settings *settings, grasp_tile_description_t *tile_description, grasp_processing_functions_t *functions, grasp_results_t *results);

/**
 * Iterate over segments in a tile for retrieve its
 * @param settings Current settings structure
 * @param tile_description Current description of the tile
 * @param results Result structure will be filled with results from retrievals
 * @param functions extension functions will be called during the process
 */
void grasp_controller_invert_tile(grasp_settings *settings, grasp_tile_description_t *tile_description, grasp_results_t *results, grasp_processing_functions_t *functions);

/**
 * Call inversion with settings, tile and output, knowing the output function of the segment and the number of inversion and its position in the tile (icol, irow, itime)
 * @param settings
 * @param segment
 * @param output
 * @param tile_description
 * @param results
 * @param functions
 * @param iinversion
 * @param id_inversion
 */
void grasp_controller_call_inversion(grasp_settings *settings, grasp_segment_t *segment, output_segment_general *output, grasp_tile_description_t *tile_description, grasp_results_t *results, grasp_processing_functions_t *functions, int iinversion, int id_inversion);

/**
 * Given a inversion number it return the position of the tile icol, irow and itime. The function return -1 if the number of inversion is out of bounds, 0 otherwise.
 * @param tile_description current tile description
 * @param iinversion number of the inversion. First, second...
 * @return 0 if the function worked OK
 */
int grasp_controller_get_next_inversion(grasp_tile_dimensions_t *tile_description, int iinversion);

/**
 * Called by grasp_controller_call_inversion. It performs a processing of a segment
 * @param settings Current input settings
 * @param segment Current segment to be processed
 * @param output Output obtained from retrieval
 * @param tile_description Description of current tile 
 * @param ninversion Number of current inversion used in print statements
 * @return error code. 0 if the process is ok.
 */
int grasp_controller_processor_unit(const grasp_settings *settings, const grasp_segment_t *segment, output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int ninversion);

/**
 * Called by grasp_controller_call_inversion after a segment has been processed. It will call output extension functions
 * @param settings Current input segment
 * @param segment Segment which has been processed
 * @param output Output obtained from the retrieval algorithm
 * @param tile_description Current tile description
 * @param results Result structure where output will be included
 * @param functions Extension functions
 */
void grasp_controller_post_process_segment(grasp_settings *settings, grasp_segment_t *segment, output_segment_general *output, grasp_tile_description_t *tile_description, grasp_results_t *results, grasp_processing_functions_t *functions);

/**
 * This function check if a segment is an range of segments that will be retrieved and true will be returned if the segment have to be retrieved, otherwise it will return false
 * @param ninversion Number of current inversion
 * @param settings Current settings
 * @param tile_description Description of current tile
 * @return true if it is invertible
 */
bool grasp_controller_segment_is_invertible(int ninversion, const grasp_settings *settings, const grasp_tile_dimensions_t *tile_description);

/**
 * Process controller options (call help, call debug, ...). It checks the input options and call controller functions to answer user demands
 * @param settings Settings from user yaml file
 */
void grasp_controller_process_options(grasp_settings *settings);

/**
 * It function call output extension function with result structure. After inverting all segments this function will work with output result structure
 * @param settings Current input settings
 * @param tile_description Current tile description
 * @param results Result structure with all output from retrievals
 * @param functions Output extension function to call after processing and generate output files
 */
void grasp_controller_manage_tile(grasp_settings *settings, grasp_tile_description_t *tile_description, grasp_results_t *results, grasp_processing_functions_t *functions);

/**
 * Clean functions, results, tile description and settings. The application is ready to finish.
 * @param settings Current input settings
 * @param tile_description Current tile description
 * @param results Output results obtained from all retrievals
 * @param functions Extension functions used during the process
 */
void grasp_controller_clean_memory(grasp_settings *settings, grasp_tile_description_t *tile_description, grasp_results_t *results, grasp_processing_functions_t *functions);



#ifdef	__cplusplus
}
#endif

#endif	/* GRASP_CONTROLLER_H */

