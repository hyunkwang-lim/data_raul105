/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

#ifndef GRASP_DATA_H
#define	GRASP_DATA_H

#include "mod_par_OS.inc"
#include "mod_par_inv.inc"
#include <time.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include "../settings/grasp_settings.h"
#include "../output/grasp_output_stream_t.h"
#include "../controller/grasp_controller_functions.h"

#define GRASP_INPUT_UNKNOWN_FILEINDEX -1


// Initialize grasp tile description
int grasp_input_initialize_tile_description(grasp_settings *settings, grasp_tile_description_t *tile_description);

// Main function for get segment data. 
// Return segment filled with segment data and return number of pixels retrieved. Return a number lower than 0 if there was a problem
int grasp_input_extract_segment(grasp_settings *settings, grasp_input_driver_t *driver, int ntransformers,grasp_input_transformer_t *transformers, grasp_segment_t *segment,grasp_results_t *results, grasp_tile_dimensions_t *tile_dimensions, int id_inversion);

// Initialization funtion from retrieval module. It is implemented in mod_edges.f90
extern void grasp_input_edges_initialization(segment_edges *edges);

// Set in segment.edges the information of available pixels around the segment
void grasp_input_edges_find(grasp_settings *settings, grasp_segment_t *segment, grasp_tile_dimensions_t *tile_dimensions, grasp_results_t *results, int icol, int irow, int itime);

/* provides the number of pixels in the group of pixels corresponding to 
 * the it index time (between 1 and NT).
 */
int grasp_input_get_num_pixels_at_it(const grasp_segment_t *segment, int it);

/**
 * Return the position in the tile of a inversion identified by id_inversion
 * @param tile_dimensions Current tile dimension information
 * @param id_inversion Identification of the inversion which user want to know the position in the tile
 * @param icol [output] Column of the segment in the tile
 * @param irow [output] Row of the segment in the tile
 * @param itime [output] index of time of the segment in the tile
 * @return 0 if the segment is in the tile, otherwise it will be -1
 */
int grasp_input_position_of_inversion(const grasp_tile_dimensions_t *tile_dimensions, int id_inversion, int *icol, int *irow, int *itime);

/*
 * This function remove missing_value from segment and update the corresponding related fields
 */
void grasp_input_clean_segment(grasp_segment_t *segment, float missing_value);

/**
 * Dump initial guess of each pixel in image.dat format
 * @param output_stream Stream where dump the information (screen, file, ...)
 * @param settings Settings file of current process needed to know number of characteristics that drives in forward model
 * @param segment Initial guess that will be dumped.
 */
void grasp_input_dump_iguess(FILE *output_stream, const grasp_settings *settings, const grasp_segment_t *segment);
/*
 * This function dumps a segment in the SDATA format.
 */
void grasp_input_dump_segment(FILE *output_stream, const grasp_segment_t *segment);

/*
 * This function prints a segment in a readable form, for debugging purpose.
 * The label is any string the caller wants to display in front of the segment
 * (may be set to NULL or an empty string if no label is needed).
 */
void grasp_input_print_segment(FILE *output_stream, const char *label, const grasp_segment_t *segment);

// Return all filenames that will be processed expading filenames array of configuration
int grasp_input_expand_files(grasp_settings *settings, char ***files);

// trackmem_free memory of input filenames used
void grasp_input_deallocate_expanded_files(char **files, int nfiles);

#endif	/* GRASP_DATA_H */
