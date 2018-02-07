/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

/**
 * @file  grasp_driver_settings.h
 * @author Fabrice Ducos <fabrice.ducos@univ-lille1.fr>
 *
 */

#ifndef DRIVER_SETTINGS_H_
#define DRIVER_SETTINGS_H_

#include <time.h>
#include <stdbool.h>

#include "grasp_box.h"

typedef struct grasp_driver_settings_t_ {
  /* private, do not access directly */
  //char file_pattern[1023 + 1]; // long enough to hold very long paths
  int nfiles;
  char **files;
  int xmin; /* coordinate of pixel */
  int xmax; /* coordinate of pixel */
  int ymin; /* coordinate of pixel */
  int ymax; /* coordinate of pixel */
  time_t gmt_time_min; /* time in seconds since the Unix Epoch */
  time_t gmt_time_max; /* time in seconds since the Unix Epoch */
  char str_gmt_time_min[20 + 1]; // format ISO 8601 YYYY-MM-DDThh:mm:ssZ
  char str_gmt_time_max[20 + 1]; // format ISO 8601 YYYY-MM-DDThh:mm:ssZ
  double missing_value; /* initialization value when data are not available */
} grasp_driver_settings_t;

#if defined (__cplusplus)
extern "C" {
#endif

extern void grasp_driver_settings_set_dimensions(size_t nrows, size_t ncols);
extern bool grasp_driver_settings_init_by_row_col(grasp_driver_settings_t *settings, 
					    int nfiles, char **files, 
					    double row_min, double col_min, 
					    double row_max, double col_max, 
					    const char *str_gmt_time_min, const char *str_gmt_time_max, 
					    double missing_value);

extern void grasp_driver_settings_clear(grasp_driver_settings_t *settings);
extern bool grasp_driver_settings_are_valid(const grasp_driver_settings_t *settings);
extern void grasp_driver_settings_convert(const grasp_driver_settings_t *user_settings, grasp_box_settings_t *box_settings);
extern void grasp_driver_settings_print(FILE *stream, const char *label, const grasp_driver_settings_t *settings);

#if defined (__cplusplus)
}
#endif

#endif /* DRIVER_SETTINGS_H_ */
