/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

/**
 * @file  grasp_box.h
 * @author Fabrice Ducos <fabrice.ducos@univ-lille1.fr>
 *
 */

#ifndef GRASP_BOX_H_
#define GRASP_BOX_H_

#include <stdlib.h> /* size_t */
#include <stdio.h>  /* FILE * */
#include <stdbool.h>
#include <grasp/utils.h> /* filepath_t */

typedef struct grasp_box_vector_t_ {
  int x, y, z;
} grasp_box_vector_t;

typedef struct grasp_box_settings_t_ {
  grasp_box_vector_t orig;
  grasp_box_vector_t extent;
} grasp_box_settings_t;

typedef struct grasp_box_t_ grasp_box_t;

#if defined (__cplusplus)
extern "C" {
#endif

extern void grasp_box_set_debug_level(int level);
extern void grasp_box_rewind(grasp_box_t *box);
extern grasp_box_t *grasp_box_new(size_t record_size, const grasp_box_settings_t *settings);
extern void grasp_box_print(FILE *output_stream, const char *label, const grasp_box_t *box);
extern void grasp_box_delete(grasp_box_t *box);
extern bool grasp_box_attach_data_to_pixel(grasp_box_t *box, const grasp_box_vector_t *pixel, const size_t data_size, const void *data);
extern bool grasp_box_get_data_from_pixel(const grasp_box_t *box, const grasp_box_vector_t *pixel, void *data_copy);
extern bool grasp_box_set_data_on_pixel(grasp_box_t *box, const grasp_box_vector_t *pixel, const void *data_copy);
extern size_t grasp_box_get_num_records(const grasp_box_t* box);
extern bool grasp_box_get_next_pixel(grasp_box_t *box, void *data_copy);
extern void grasp_box_get_settings(const grasp_box_t *box, grasp_box_settings_t *settings_copy);

extern void grasp_box_add_parent_file(grasp_box_t *box, const char *parent_file);
extern size_t grasp_box_get_num_parent_files(const grasp_box_t *box);
extern void grasp_box_get_parent_files(const grasp_box_t *box, size_t num_parent_files, filepath_t *parent_files /* must be allocated by the caller */);

#if defined (__cplusplus)
}
#endif

#endif /* GRASP_BOX_H_ */
