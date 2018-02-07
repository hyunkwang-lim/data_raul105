/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

/**
 * @file grasp_driver_settings.c 
 * @author Fabrice Ducos <fabrice.ducos@univ-lille1.fr>
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>

#include <grasp/utils/debug.h>
#include <grasp/utils/time_utils.h>
#include "grasp_box.h"

#define GRASPUTILS_BUILD_TIME_
#include "grasp_driver_settings.h"

static size_t NROWS = 0;
static size_t NCOLS = 0;

/* swaps a pair of values if necessary to make sure that
 * the first one is the min and the second one is the max
 */
static void sort_min_max(double *min, double *max) {
  if (*min > *max) {
    double tmp;
    
    tmp = *min;
    *min = *max;
    *max = tmp;
  }

  assert(*min <= *max);
}

void grasp_driver_settings_clear(grasp_driver_settings_t *settings) {
  memset(settings, 0, sizeof(grasp_driver_settings_t));
}

void grasp_driver_settings_print(FILE *stream, const char *label, const grasp_driver_settings_t *settings) {
  assert(stream != NULL);
  assert(settings != NULL);

  if (label != NULL && label[0] != '\0') {
    fprintf(stream, "%s:\n", label);
  }

  fprintf(stream, "  nfiles:  %d\n"      ,   settings->nfiles);
  fprintf(stream, "  xmin:          %d\n",   settings->xmin);
  fprintf(stream, "  xmax:          %d\n",   settings->xmax);
  fprintf(stream, "  ymin:          %d\n",   settings->ymin);
  fprintf(stream, "  ymax:          %d\n",   settings->ymax);
  fprintf(stream, "  gmt_time_min:  %s [%ld]\n", settings->str_gmt_time_min, (long) settings->gmt_time_min);
  fprintf(stream, "  gmt_time_max:  %s [%ld]\n", settings->str_gmt_time_max, (long) settings->gmt_time_max);
  fprintf(stream, "  missing_value: %g\n",   settings->missing_value);
  
}

static bool time_is_valid(const char *file, int line, time_t time) {
  bool retval;
  
  retval = (time != (time_t) 0 && time != (time_t) -1);
  if (! retval) {
#ifdef DEBUG
    const char *time_format = "%FT%H:%M:%SZ";
    fprintf(DEBUG_STREAM, "%s:%d: invalid time: %ld (%s)\n", file, line, time, time_to_string(time, time_format));
#endif
  }
  return retval;
}

#define LONGITUDE_IS_VALID(lon) longitude_is_valid(__FILE__, __LINE__, lon)
#define LATITUDE_IS_VALID(lat)  latitude_is_valid(__FILE__, __LINE__, lat)
#define TIME_IS_VALID(time)     time_is_valid(__FILE__, __LINE__, time)

void grasp_driver_settings_set_dimensions(size_t nrows, size_t ncols) {
  NROWS = nrows;
  NCOLS = ncols;
}

static void check_boundaries(double *row_min, double *row_max, double *col_min, double *col_max) {
  if (*row_min < 1) {
    fprintf(DEBUG_STREAM, "%s:%d: row_min (%f) out of the range [1..%d], forced to 1\n", __FILE__, __LINE__, *row_min, (int) NROWS);
    *row_min = 1;
  }
  
  if (*row_max > NROWS) {
    fprintf(DEBUG_STREAM, "%s:%d: row_min (%f) out of the range [1..%d], forced to %d\n", __FILE__, __LINE__, *row_max, (int) NROWS, (int) NROWS);
    *row_max = NROWS;
  }
  
  if (*col_min < 1) {
    fprintf(DEBUG_STREAM, "%s:%d: col_min (%f) out of the range [1..%d], forced to 1\n", __FILE__, __LINE__, *col_min, (int) NCOLS);
    *col_min = 1;
  }
  
  if (*col_max > NCOLS) {
    fprintf(DEBUG_STREAM, "%s:%d: col_max (%f) out of the range [1..%d], forced to %d\n", __FILE__, __LINE__, *col_max, (int) NCOLS, (int) NCOLS);
    *col_max = NCOLS;
  }

}

bool grasp_driver_settings_are_valid(const grasp_driver_settings_t *settings) {
  bool retval = true;
  
  /* one doesn't return at the first failure so one can display as many errors as possible */

  if (! TIME_IS_VALID(settings->gmt_time_min)) {
    fprintf(DEBUG_STREAM, "%s:%d: invalid gmt_time_min: %s [time %ld]\n", __FILE__, __LINE__, settings->str_gmt_time_min,
	    settings->gmt_time_min);
    retval = false;
  }
  
  if (! TIME_IS_VALID(settings->gmt_time_max)) {
    fprintf(DEBUG_STREAM, "%s:%d: invalid gmt_time_max: %s [time %ld]\n", __FILE__, __LINE__, settings->str_gmt_time_max,
	    settings->gmt_time_max);
    retval = false;
  }

   if (! (settings->gmt_time_min <= settings->gmt_time_max)) {
    fprintf(DEBUG_STREAM, "%s:%d: ymin (%s) is not less than ymax (%s)\n", __FILE__, __LINE__, 
	    settings->str_gmt_time_min, settings->str_gmt_time_max);
    retval = false;
  }
  
  if (! (1 <= settings->xmin && settings->xmin <= NCOLS)) {
    fprintf(DEBUG_STREAM, "%s:%d: invalid xmin: %d, out of the range [1..%d]\n", __FILE__, __LINE__, settings->xmin, (int) NCOLS);
    retval = false;
  }
  
  if (! (1 <= settings->xmax && settings->xmax <= NCOLS)) {
    fprintf(DEBUG_STREAM, "%s:%d: invalid xmax: %d, out of the range [1..%d]\n", __FILE__, __LINE__, settings->xmax, (int) NCOLS);
    retval = false;
  }
  
  if (! (1 <= settings->ymin && settings->ymin <= NROWS)) {
    fprintf(DEBUG_STREAM, "%s:%d: invalid ymin: %d, out of the range [1..%d]\n", __FILE__, __LINE__, settings->ymin, (int) NROWS);
    retval = false;
  }
  
  if (! (1 <= settings->ymax && settings->ymax <= NROWS)) {
    fprintf(DEBUG_STREAM, "%s:%d: invalid ymax: %d, out of the range [1..%d]\n", __FILE__, __LINE__, settings->ymax, (int) NROWS);
    retval = false;
  }
  
  if (! (settings->xmin <= settings->xmax)) {
    fprintf(DEBUG_STREAM, "%s:%d: xmin (%d) is not less than xmax (%d)\n", __FILE__, __LINE__, 
	    settings->xmin, settings->xmax);
    retval = false;
  }
  
  if (! (settings->ymin <= settings->ymax)) {
    fprintf(DEBUG_STREAM, "%s:%d: ymin (%d) is not less than ymax (%d)\n", __FILE__, __LINE__, 
	    settings->ymin, settings->ymax);
    retval = false;
  }
  
  return retval;
}


void grasp_driver_settings_convert(const grasp_driver_settings_t *user_settings, grasp_box_settings_t *box_settings) {
  int xmin, ymin, zmin;
  int xmax, ymax, zmax;

  assert(user_settings != NULL);
  assert(box_settings != NULL);

  xmin = user_settings->xmin;
  ymin = user_settings->ymin;
  zmin = user_settings->gmt_time_min / 86400; /* number of days since the 1st january of 1970 */

  xmax = user_settings->xmax;
  ymax = user_settings->ymax;
  zmax = user_settings->gmt_time_max / 86400; /* number of days since the 1st january of 1970 */

#ifdef DEBUG2
  fprintf(stderr, "%s: %d: ", __FILE__, __LINE__);
  fprintf(stderr, "xmin: %d xmax: %d ymin: %d ymax: %d zmin: %d zmax: %d\n", xmin, xmax, ymin, ymax, zmin, zmax);
#endif
  assert(xmin <= xmax);
  assert(ymin <= ymax);
  assert(zmin <= zmax);

  box_settings->orig.x = xmin;
  box_settings->orig.y = ymin;
  box_settings->orig.z = zmin;

  box_settings->extent.x = xmax - xmin + 1;
  box_settings->extent.y = ymax - ymin + 1;
  box_settings->extent.z = zmax - zmin + 1;
}

bool grasp_driver_settings_init_by_row_col(grasp_driver_settings_t *settings, int nfiles, char **files, double row_min, double col_min, 
				     double row_max, double col_max, const char *str_gmt_time_min, 
				     const char *str_gmt_time_max, double missing_value) {
  time_t gmt_time_min;
  time_t gmt_time_max;
  int err;

  grasp_driver_settings_clear(settings);
  
  /* make sure min and max values were not inverted by the caller */
  sort_min_max(&row_min, &row_max);
  sort_min_max(&col_min, &col_max);
  
  check_boundaries(&row_min, &row_max, &col_min, &col_max);

  err = convert_string_to_time(str_gmt_time_min, &gmt_time_min, TIMEFMT_ISO8601);
  if (err != 0) {
#ifdef DEBUG
    fprintf(DEBUG_STREAM, "%s:%d: failed at convert_string_to_time(str_gmt_time_min = \"%s\", &gmt_time_min, TIMEFMT_ISO8601)\n",
	    __FILE__, __LINE__, str_gmt_time_min);
#endif
    return false;
  }
  
  err = convert_string_to_time(str_gmt_time_max, &gmt_time_max, TIMEFMT_ISO8601);
  if (err != 0) {
#ifdef DEBUG
    fprintf(DEBUG_STREAM, "%s:%d: failed at convert_string_to_time(str_gmt_time_max = \"%s\", &gmt_time_max, TIMEFMT_ISO8601)\n",
	    __FILE__, __LINE__, str_gmt_time_max);
#endif
    return false;
  }

#ifdef DEBUG2
  {
    size_t ifile;
    for (ifile = 0 ; ifile < nfiles ; ifile++) {
      fprintf(stderr, "%s:%d: %s: file #%2.2d: ", __FILE__, __LINE__, __func__, (int) ifile + 1);
      fprintf(stderr, "%s\n", files[ifile]);
    }
  }
#endif


  settings->nfiles=nfiles;
  settings->files=files;
  strncpy(settings->str_gmt_time_min, str_gmt_time_min, sizeof(settings->str_gmt_time_min) - 1);
  strncpy(settings->str_gmt_time_max, str_gmt_time_max, sizeof(settings->str_gmt_time_max) - 1);
  
  settings->xmin = lround(col_min);
  settings->xmax = lround(col_max);
  settings->ymin = lround(row_min);
  settings->ymax = lround(row_max);
  settings->gmt_time_min = gmt_time_min;
  settings->gmt_time_max = gmt_time_max;
  settings->missing_value = missing_value;

  return grasp_driver_settings_are_valid(settings);
}


