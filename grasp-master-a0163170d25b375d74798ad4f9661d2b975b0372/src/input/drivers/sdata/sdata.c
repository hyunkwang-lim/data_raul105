/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

/**
 * @file  sdata.c
 * @author Fabrice Ducos <fabrice.ducos@univ-lille1.fr>
 *
 */

//#define DEBUG2

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdbool.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <errno.h>

#include <grasp/utils/error.h>
#include "../../util/grasp_driver_settings.h"

#include "mod_par_inv.inc"
#include "sdata.h"
#include "sdata-impl.h"
#include <grasp/utils/print_array.h>
#include "check_helpers.h"

#define CSTRING_MAXLEN 255

#ifndef NAN
#warning "NAN is not provided by <math.h>. I'll use my own definition."
#define NAN 0./0.
#endif

#ifdef WARN_DRY
#warning "__SDATA_VERSION__ duplicated"
#endif   
const int driver_major_version = 2;

inline double min(double a, double b) {
  return (a < b ? a : b);
}

inline double max(double a, double b) {
  return (a > b ? a : b);
}

const char separator[] = " ";
const double MISSING_VALUE = NAN;

#ifdef WARN_DRY
#warning "__MEAS_TYPE__ duplicated"
#endif
/* MEAS_TYPE_* defined in src/global/grasp_retrieval_meas_type.h */
const char *sdata_str_meas_types(int meas_type) {
    switch (meas_type) {
        case MEAS_TYPE_UNKNOWN: { return "UNKNOWN" ; break ; }
        case MEAS_TYPE_TOD: { return "TAU" ; break ; }
        case MEAS_TYPE_AOD: { return "EXT" ; break ; }
        case MEAS_TYPE_P11: { return "P11" ; break ; }
        case MEAS_TYPE_P12: { return "P12" ; break ; }
        case MEAS_TYPE_P22: { return "P22" ; break ; }
        case MEAS_TYPE_P33: { return "P33" ; break ; }
        case MEAS_TYPE_P34: { return "P34" ; break ; }
        case MEAS_TYPE_P44: { return "P44" ; break ; }
        case MEAS_TYPE_LS: { return "LS" ; break ; }
        case MEAS_TYPE_DP: { return "DP" ; break ; }
        case MEAS_TYPE_RL: { return "RL" ; break ; }
        case MEAS_TYPE_I: { return "I" ; break ; }
        case MEAS_TYPE_Q: { return "Q" ; break ; }
        case MEAS_TYPE_U: { return "U" ; break ; }
        case MEAS_TYPE_P: { return "P" ; break ; }
        default: {
            fprintf(sdata_get_error_stream(), "%s:%d: unexpected value for meas_type: %d\n", __FILE__, __LINE__,
                    meas_type);
            abort();
            break ; 
        }
        
        
    }
}


/* unfortunately, err_stream can't be statically initialized to stderr (initalizer not constant);
 * for this reason, err_stream must never be sent directly to stream-aware functions (such as fprintf).
 * Such functions must be used with sdata_get_error_stream that will take care automatically of the
 * right interpretation of err_stream.
 */
static FILE *err_stream = NULL; /* error messages will be sent to stderr by default (see sdata_get_error_stream) */
static FILE *debug_stream = NULL; /* no debugging info by default */

/*************************************************************************************************/
/* internal implementation of perror that accepts a stream argument
 * (instead of sending automatically to stderr) 
 */
static void my_perror(FILE *stream, const char *s) {
  if (s != NULL && s[0] != '\0') {
    fprintf(stream, "%s: ", s);
  }

  fprintf(stream, "%s\n", strerror(errno));
}

/*************************************************************************************************/

void sdata_set_error_stream(FILE *stream) {
  /* err_stream can't be NULL: error messages must always be displayed.
   * If the user really wants to get rid of them, he/she can always take the responsibility
   * to redirect them to /dev/null.
   */
  if (stream == NULL) {
    err_stream = stderr;
  }
  else {
    err_stream = stream;
  }
}

void sdata_set_debug_stream(FILE *stream) {
  /* NULL means no debugging */
  debug_stream = stream;
}

FILE *sdata_get_error_stream(void) {
  if (err_stream == NULL) {
    return stderr;
  }
  else {
    return err_stream;
  }
}

FILE *sdata_get_debug_stream(void) {
  return debug_stream;
}

/*************************************************************************************************/
static void init_valid_ranges(SDATA_VALID_RANGES *valid_ranges) {
  valid_ranges->ix_min = 1;
  valid_ranges->ix_max = 10;
  valid_ranges->iy_min = 1;
  valid_ranges->iy_max = 10;
  valid_ranges->cloud_flag_min = 0;
  valid_ranges->cloud_flag_max = 1;
  valid_ranges->longitude_min = -180.;
  valid_ranges->longitude_max =  180.;
  valid_ranges->latitude_min  =  -90.;
  valid_ranges->latitude_max  =   90.;
  valid_ranges->hgr_min    =  -500.;
  valid_ranges->hgr_max    = 10000.;
  valid_ranges->land_percent_min = 0.;
  valid_ranges->land_percent_max = 100.;
  valid_ranges->num_wavelengths_min = 1;
  valid_ranges->num_wavelengths_max = SDATA_MAX_NWL;
  valid_ranges->wavelength_min = 0.300;
  valid_ranges->wavelength_max = 2.500;
  valid_ranges->num_meas_type_min = 1;
  valid_ranges->num_meas_type_max = 3;
  valid_ranges->meas_type_min     = MEAS_TYPE_TOD;
  valid_ranges->meas_type_max     = MEAS_TYPE_P;
  valid_ranges->num_valid_meas_min = 0;
  valid_ranges->num_valid_meas_max = SDATA_MAX_KNBVM;
  valid_ranges->thetas_min = 0.;
  valid_ranges->thetas_max = 90.;
  valid_ranges->thetav_min = 0.;
  valid_ranges->thetav_max = 180.; /* should not be greater than 90 for POLDER, but can be greater for AERONET */
  valid_ranges->dphi_min = -720.;
  valid_ranges->dphi_max = 720.;
  valid_ranges->hvp_min = -500.; /* in meters */
  valid_ranges->hvp_max = 30000.; /* in meters */
  valid_ranges->groundpar_min = 0.;
  valid_ranges->groundpar_max = 0.;
  valid_ranges->gaspar_min = 0.;
  valid_ranges->gaspar_max = 1.;
  valid_ranges->meas_min = -9999.;
  valid_ranges->meas_max =  9999.;
  valid_ranges->ifcov_min = 0;
  valid_ranges->ifcov_max = 1;
  valid_ranges->ifmp_min  = 0;
  valid_ranges->ifmp_max  = 1;
  valid_ranges->cmtrx_min = 0.;
  valid_ranges->cmtrx_max = 10.;
  valid_ranges->mprof_min = 0.0;
  valid_ranges->mprof_max = 1.0e-3;
}

/*************************************************************************************************/

static void sdata_fill_timestamp(SDATA_TIMESTAMP *timestamp, const char *str_date_time,
				 time_t unix_time, long index_of_day) {
  assert(timestamp != NULL);
  assert(str_date_time != NULL);
  
  strncpy(timestamp->str_date_time, str_date_time, sizeof(timestamp->str_date_time) - 1);
  timestamp->str_date_time[sizeof(timestamp->str_date_time) - 1] = '\0';

  timestamp->unix_time = unix_time;
  timestamp->index_of_day = index_of_day;
}

/*************************************************************************************************/

static int sdata_get_sdata_limits(SDATA_HANDLE *handle, SDATA_TIME_LIMITS *time_limits, int *xmin, int *xmax, int *ymin, int *ymax, int *npixels_total) {
  SDATA_HEADER header;
  SDATA_RECORD record;
  int irecord;
  int err;
  int xmin_, xmax_;
  int ymin_, ymax_;
  int npixels_total_; /* total number of pixels in the file (sum of the pixels from all the records) */

  assert(handle != NULL);
  assert(time_limits != NULL);

  err = sdata_get_header(handle, &header);
  assert(err == SDATA_OK);
  
  irecord = 0;
  npixels_total_ = 0;

  /* yes, this is not an error. Min and max are initialized with reversed system limits
   * so the classic tests for finding actual limits can be applied.
   */
  xmin_ = INT_MAX;
  xmax_ = INT_MIN;
  ymin_ = INT_MAX;
  ymax_ = INT_MIN;

  while(true) {
    char str_date_time[20 + 1]; /* format YYYY-MM-DDThh:mm:ssZ */
    time_t unix_time;
    long index_of_day;
    int iret;
    bool bret;
    size_t ipixel;
    (void)iret; // This avoid -Wunused-but-set-variable warning
    (void)bret; // This avoid -Wunused-but-set-variable warning
    
    /* do not attempt to read more records than specified by the header */
    if (irecord == header.nrecords) break;
    
    err = sdata_read_next_record(handle, &record);
    
    if (err != SDATA_OK) {
      if (err == SDATA_EOF) break;
      
      fprintf(sdata_get_error_stream(), "%s:%d: ", __FILE__, __LINE__);
      fprintf(sdata_get_error_stream(), "sdata_read_next_record() failed in %s\n", __FUNCTION__);
      abort();
    }
    
    sdata_get_date_time(&record, str_date_time);
    iret = convert_string_to_time(str_date_time, &unix_time /* seconds since 1970-01-01T00-00-00Z */, TIMEFMT_ISO8601);
    assert(iret == 0);
    bret = get_index_of_day(str_date_time, &index_of_day);
    assert(bret == true);
    
    if (debug_stream) {
      fprintf(sdata_get_debug_stream(), "%s:%d: ", handle->file_name, handle->current_line);
      fprintf(sdata_get_debug_stream(), "str_date_time: %s (index_of_day: %ld) (unix_time: %ld)\n", str_date_time, index_of_day, (long) unix_time);
    }
    
    /* the first record can also be the last one ! */
    if (irecord == 0) {
      sdata_fill_timestamp(&time_limits->first_record, str_date_time, unix_time, index_of_day);
    }
    
    if (irecord == header.nrecords - 1) {
      sdata_fill_timestamp(&time_limits->last_record, str_date_time, unix_time, index_of_day);
    }

    for (ipixel = 0 ; ipixel < record.npixels ; ipixel++) {
      grasp_box_vector_t tpixel;
      const SDATA_PIXEL *spixel = sdata_get_pixel(&record, ipixel);
      tpixel.x = spixel->icol;
      tpixel.y = spixel->irow;
      
      if (tpixel.x < xmin_) xmin_ = tpixel.x;
      if (tpixel.x > xmax_) xmax_ = tpixel.x;
      if (tpixel.y < ymin_) ymin_ = tpixel.y;
      if (tpixel.y > ymax_) ymax_ = tpixel.y;
    }

    npixels_total_ += record.npixels;
    irecord++;
  } /* while (true) */
  assert(irecord == header.nrecords);
  assert(npixels_total_ <= header.nx * header.ny * header.nrecords);

  *xmin = xmin_;
  *xmax = xmax_;
  *ymin = ymin_;
  *ymax = ymax_;
  *npixels_total = npixels_total_;
  
  return SDATA_OK;

}

/*************************************************************************************************/

static int sdata_parse_version(const char *version, int *major, int *minor) {
  int nread;
  
  nread = sscanf(version, "%d.%d", major, minor);

  if (nread != 2) {
    return SDATA_ERROR;
  }

  return SDATA_OK;
}

static void sdata_rewind(SDATA_HANDLE *handle) {
  size_t nread;
  int major, minor;
  int err;

  rewind(handle->file_pointer);
  handle->current_line = 1;
  handle->current_record = 0;
  
  /* reads the header lines */
  nread = fscanf(handle->file_pointer, "SDATA version %s", handle->header.format_version);
  if (nread != 1) {
    fprintf(sdata_get_error_stream(), "%s: not a valid SDATA file (SDATA version missing on the first line)\n", handle->file_name);
    exit (EXIT_FAILURE);
  }

  err = sdata_parse_version(handle->header.format_version, &major, &minor);
  if (err == SDATA_ERROR) {
    fprintf(sdata_get_error_stream(), "%s: the version pattern could not be recognized (%s)\n", handle->file_name, handle->header.format_version);
    exit (EXIT_FAILURE);
  }
  
  if (major != driver_major_version) {
    fprintf(sdata_get_error_stream(), "%s: version %s not supported (driver's version is %d)\n", handle->file_name, handle->header.format_version, driver_major_version);
    exit (EXIT_FAILURE);
  }
  
  discard_rest_of_line(handle->file_pointer); /* discards comments at the end of the line */
  
  handle->current_line++;

  nread = fscanf(handle->file_pointer, "%d %d %d", 
		 &handle->header.nx, &handle->header.ny, &handle->header.nrecords);
  assert(nread == 3);

  assert(handle->header.nx <= SDATA_MAX_NX);
  assert(handle->header.ny <= SDATA_MAX_NY);

  init_valid_ranges(&handle->valid_ranges);

  discard_rest_of_line(handle->file_pointer); /* discards comments at the end of the line */
  handle->current_line++;
  discard_rest_of_line(handle->file_pointer); /* reads and discards an empty line */
  assert(handle->current_line == 3);
}

/*************************************************************************************************/

SDATA_HANDLE *sdata_open(const char *file_name) {
  SDATA_TIME_LIMITS time_limits;
  SDATA_HANDLE *handle;
  int err;
  (void)err; // This avoid -Wunused-but-set-variable warning
  
  assert(file_name != NULL);
#ifdef DEBUG2
  sdata_set_debug_stream(stderr);
#endif
  
  handle = trackmem_malloc(sizeof(*handle));
  assert(handle != NULL);
  
  strncpy(handle->file_name, file_name, CSTRING_MAXLEN);
  handle->file_name[CSTRING_MAXLEN] = '\0'; /* just for safety */
  handle->file_pointer = fopen(file_name, "r");
  if (handle->file_pointer == NULL) {
    if (debug_stream) {
      fprintf(sdata_get_debug_stream(), "%s:%d: sdata_open(\"%s\") failed: ", __FILE__, __LINE__, file_name);
      my_perror(sdata_get_debug_stream(), file_name);
    }
    return NULL;
  }

  if (debug_stream) {
    fprintf(sdata_get_debug_stream(), "%s:%d: ", __FILE__, __LINE__);
    fprintf(sdata_get_debug_stream(), "%s has been opened successfully\n", handle->file_name);
  }

  sdata_rewind(handle);

  if (debug_stream) {
    fprintf(sdata_get_debug_stream(), "%s:%d: ", __FILE__, __LINE__);
    fprintf(sdata_get_debug_stream(), "%s: a first pass is needed to get some metadata (count of records, spatial and time boundaries) from the file's contents. The data will be read but not stored (data allocation can be performed only at the end of this first pass).\n",
	    handle->file_name);
  }
  err = sdata_get_sdata_limits(handle, &time_limits, 
			       &handle->header.xmin, &handle->header.xmax, 
			       &handle->header.ymin, &handle->header.ymax, &handle->header.npixels);
  
  if(handle->header.nx>_KIX){
      fprintf(sdata_get_error_stream(), "grasp: %s: the dimension nx (%d) in the sdata file is too large. The current upper bound is %d. A workaround is to rebuild grasp with a different constant set (make CONSTANTS_SET=...)\n", handle->file_name, handle->header.nx, _KIX);
      exit (EXIT_FAILURE);
  }
  
  if(handle->header.ny>_KIY){
      fprintf(sdata_get_error_stream(), "grasp: %s: the dimension ny (%d) in the sdata file is too large. The current upper bound is %d. A workaround is to rebuild grasp with a different constant set (make CONSTANTS_SET=...)\n", handle->file_name, handle->header.ny, _KIY);
      exit (EXIT_FAILURE);
  }
  
  if(handle->header.nrecords>_KITIME){
      fprintf(sdata_get_error_stream(), "grasp: %s: the dimension ntime (%d) in the sdata file is too large. The current upper bound is %d. A workaround is to rebuild grasp with a different constant set (make CONSTANTS_SET=...)\n", handle->file_name, handle->header.nrecords, _KITIME);
      exit (EXIT_FAILURE);
  }
  
  if(handle->header.npixels>_KIMAGE){
      fprintf(sdata_get_error_stream(), "grasp: %s: the number of total pixels (%d) in the sdata file is too large. The current upper bound is %d. A workaround is to rebuild grasp with a different constant set (make CONSTANTS_SET=...)\n", handle->file_name, handle->header.npixels, _KIMAGE);
      exit (EXIT_FAILURE);
  }
  
  assert(err == SDATA_OK);
  sdata_rewind(handle);
  memcpy(&handle->header.time_limits, &time_limits, sizeof(SDATA_TIME_LIMITS));

  if (debug_stream) {
    fprintf(sdata_get_debug_stream(), "%s:%d: ", __FILE__, __LINE__);
    fprintf(sdata_get_debug_stream(), "%s: the header's metadata are now available. The file is ready for actual processing.\n",
	    handle->file_name);
  }
  
  return handle;
}

/*************************************************************************************************/
/* for debugging */
void sdata_print_handle(FILE *output_stream, const char *label, SDATA_HANDLE *handle) {
  assert(output_stream != NULL);
  assert(handle != NULL);

  if (label != NULL && label[0] != '\0') {
    fprintf(output_stream, "%s ", label);
  }
  fprintf(output_stream, "{\n");
  fprintf(output_stream, 
	  "  file_name:      %s\n"
	  "  file_pointer:   %p\n"
	  "  current_line:   %d\n"
	  "}\n", handle->file_name, handle->file_pointer, handle->current_line);
}

/*************************************************************************************************/

int sdata_get_header(const SDATA_HANDLE *handle, SDATA_HEADER *header) {

  assert(handle != NULL);
  assert(header != NULL);
  memcpy(header, &handle->header, sizeof(SDATA_HEADER));
  
  return SDATA_OK;
}

/*************************************************************************************************/

const char *sdata_get_file_name(const SDATA_HANDLE *handle) {
    assert(handle != NULL);
    
    return handle->file_name;
}



/*************************************************************************************************/

void sdata_dump_header(FILE *output_stream, const SDATA_HEADER *header) {
  const char *header_comments = ": NX NY NT";
  assert(output_stream != NULL);
  assert(header != NULL);

  fprintf(output_stream, "SDATA version %s\n", header->format_version);
  fprintf(output_stream, "%d%s%d%s%d%s%s\n\n", header->nx, separator, header->ny, separator, header->nrecords, separator, header_comments);
}

/*************************************************************************************************/

void sdata_print_header(FILE *output_stream, const char *label, const SDATA_HEADER *header) {
  assert(output_stream != NULL);
  assert(header != NULL);

  if (label != NULL && label[0] != '\0') {
    fprintf(output_stream, "%s\n", label);
  }

  fprintf(output_stream, "header->format_version: %s\n", header->format_version);
  fprintf(output_stream, "header->nx: %d\n", header->nx);
  fprintf(output_stream, "header->ny: %d\n", header->ny);
  fprintf(output_stream, "header->nrecords: %d\n", header->nrecords);
  fprintf(output_stream, "header->npixels: %d\n", header->npixels);
  fprintf(output_stream, "header->xmin: %d\n", header->xmin);
  fprintf(output_stream, "header->xmax: %d\n", header->xmax);
  fprintf(output_stream, "header->ymin: %d\n", header->ymin);
  fprintf(output_stream, "header->ymax: %d\n", header->ymax);
  fprintf(output_stream, "header->time_limits.first_record.str_date_time: %s\n", header->time_limits.first_record.str_date_time);
  fprintf(output_stream, "header->time_limits.first_record.unix_time: %lu\n", (long unsigned) header->time_limits.first_record.unix_time);
  fprintf(output_stream, "header->time_limits.first_record.index_of_day: %ld\n", header->time_limits.first_record.index_of_day);
  fprintf(output_stream, "header->time_limits.last_record.str_date_time:  %s\n", header->time_limits.last_record.str_date_time);
  fprintf(output_stream, "header->time_limits.last_record.unix_time:  %lu\n", (long unsigned) header->time_limits.last_record.unix_time);
  fprintf(output_stream, "header->time_limits.last_record.index_of_day:  %ld\n", header->time_limits.last_record.index_of_day);
}

/*************************************************************************************************/

static void print_array_int_2D(FILE *output_stream, const char *label, const SDATA_PIXEL *pixel, const int data[SDATA_MAX_NWL][SDATA_MAX_NIP], const char *format) {
  int iwl;
  
  assert(output_stream != NULL);
  assert(label != NULL);
  assert(pixel != NULL);
  assert(data != NULL);

  for (iwl = 0 ; iwl < pixel->num_wavelengths ; iwl++) {
    int nip = pixel->num_meas_types[iwl];
    fprintf(output_stream, "%s[%d {%.3f}] (%d): ", label, iwl, pixel->wavelengths[iwl], nip);
    println_array(output_stream, "", data[iwl], nip, 0, format, " ");
  }
}

/*************************************************************************************************/
/* NOT USED FOR THE MOMENT (disabled by a directive put to avoid warnings) */
#if 0
static void print_array_double_2D(FILE *output_stream, const char *label, const SDATA_PIXEL *pixel, const double data[SDATA_MAX_NWL][SDATA_MAX_NIP], const char *format) {
  int iwl;

  assert(output_stream != NULL);
  assert(label != NULL);
  assert(pixel != NULL);
  assert(data != NULL);

  for (iwl = 0 ; iwl < pixel->num_wavelengths ; iwl++) {
    int nip = pixel->num_meas_types[iwl];
    fprintf(output_stream, "%s[%.3f] (%d): ", label, pixel->wavelengths[iwl], nip);
    println_array(output_stream, "", data[iwl], nip, 0, format, " ");
  }
}
#endif
/*************************************************************************************************/

static void print_cmtrx(FILE *output_stream, const char *label, const SDATA_PIXEL *pixel, const double data[SDATA_MAX_NWL][SDATA_MAX_NIP][SDATA_MAX_KNBVM], const char *format) {
  int iwl;
  int ip;

  assert(output_stream != NULL);
  assert(label != NULL);
  assert(pixel != NULL);
  assert(data != NULL);

  for (iwl = 0 ; iwl < pixel->num_wavelengths ; iwl++) {
    int nip = pixel->num_meas_types[iwl];
    for (ip = 0 ; ip < nip ; ip++) {
      int meas_type = pixel->meas_types[iwl][ip];
      const char *str_meas_type = sdata_str_meas_types(meas_type);
      fprintf(output_stream, "%s[%.3f%s] (%d): ", label, pixel->wavelengths[iwl], str_meas_type, pixel->num_valid_meas[iwl][ip]);
      println_array(output_stream, "", data[iwl][ip], pixel->num_valid_meas[iwl][ip], 0, format, " ");
    }
  }
}

/*************************************************************************************************/

static void print_array_double_3D(FILE *output_stream, const char *label, const SDATA_PIXEL *pixel, const double data[SDATA_MAX_NWL][SDATA_MAX_NIP][SDATA_MAX_KNBVM], const char *format) {
  int iwl;
  int ip;

  assert(output_stream != NULL);
  assert(label != NULL);
  assert(pixel != NULL);
  assert(data != NULL);

  for (iwl = 0 ; iwl < pixel->num_wavelengths ; iwl++) {
    int nip = pixel->num_meas_types[iwl];
    for (ip = 0 ; ip < nip ; ip++) {
      int meas_type = pixel->meas_types[iwl][ip];
      const char *str_meas_type = sdata_str_meas_types(meas_type);
      fprintf(output_stream, "%s[%.3f%s] (%d): ", label, pixel->wavelengths[iwl], str_meas_type, pixel->num_valid_meas[iwl][ip]);
      println_array(output_stream, "", data[iwl][ip], pixel->num_valid_meas[iwl][ip], 0, format, " ");
    }
  }
}

/*************************************************************************************************/

void sdata_print_pixel(FILE *output_stream, const char *label, const SDATA_PIXEL *pixel) {
  char str_date_time[255 + 1];
  int iwl;
  int ip;

  assert(output_stream != NULL);
  assert(pixel != NULL);

  if (label != NULL && label[0] != '\0') {
    fprintf(output_stream, "%s\n", label);
  }

  strncpy(str_date_time, time_to_string(pixel->unix_time, NULL), 255);
  fprintf(output_stream, "pixel->ix: %d (relative to the segment)\n", pixel->ix);
  fprintf(output_stream, "pixel->iy: %d (relative to the segment)\n", pixel->iy);
  fprintf(output_stream, "pixel->it: %d (index of observation, i.e. the index of record in the SDATA)\n", pixel->it);
  fprintf(output_stream, "pixel->unix_time: %ld (%s)\n", (long) pixel->unix_time, str_date_time);
  fprintf(output_stream, "pixel->index_of_day: %ld (number of days since the Unix Epoch)\n", (long) pixel->index_of_day);
  fprintf(output_stream, "pixel->cloud_flag: %d (0 = cloudy, 1 = clear)\n", pixel->cloud_flag);
  fprintf(output_stream, "pixel->irow: %d (relative to the satellite grid)\n", pixel->irow);
  fprintf(output_stream, "pixel->icol: %d (relative to the satellite grid)\n", pixel->icol);

  fprintf(output_stream, "pixel->latitude: %.3f\n", pixel->latitude);
  fprintf(output_stream, "pixel->longitude: %.3f\n", pixel->longitude);
  fprintf(output_stream, "pixel->hgr: %.1f\n", pixel->hgr);
  fprintf(output_stream, "pixel->hobs: %.1f\n", pixel->hobs);
  fprintf(output_stream, "pixel->land_percent: %.1f\n", pixel->land_percent);
  
  fprintf(output_stream, "pixel->num_wavelengths: %d\n", pixel->num_wavelengths);
  
  println_array(output_stream, "pixel->wavelengths", pixel->wavelengths, pixel->num_wavelengths, 0, "%g", separator);
  println_array(output_stream, "pixel->num_meas_types", pixel->num_meas_types, pixel->num_wavelengths, 0, "%d", separator);
  
  print_array_int_2D(output_stream, "pixel->num_valid_meas", pixel, pixel->num_valid_meas, "%d");
  print_array_int_2D(output_stream, "pixel->ifcov", pixel, pixel->ifcov, "%d");
  print_array_int_2D(output_stream, "pixel->ifmp", pixel, pixel->ifmp, "%d");

  println_array(output_stream, "pixel->thetas", pixel->thetas, pixel->num_wavelengths, 0, "%.3f", separator);
  
  for (iwl = 0 ; iwl < pixel->num_wavelengths ; iwl++) {
    for (ip = 0 ; ip < pixel->num_meas_types[iwl] ; ip++) {
      int meas_type = pixel->meas_types[iwl][ip];
      const char *str_meas_type = sdata_str_meas_types(meas_type);

      switch (meas_type) {
      case MEAS_TYPE_TOD:
      case MEAS_TYPE_AOD:
      case MEAS_TYPE_P11:
      case MEAS_TYPE_P12:
      case MEAS_TYPE_P22:
      case MEAS_TYPE_P33:
      case MEAS_TYPE_P34:
      case MEAS_TYPE_P44:
      case MEAS_TYPE_I:
      case MEAS_TYPE_Q:
      case MEAS_TYPE_U:
      case MEAS_TYPE_P:
	{
	  char label[255 + 1];
	  snprintf(label, sizeof(label) - 1, "pixel->thetav[%.3f%s]", pixel->wavelengths[iwl], str_meas_type);
	  print_array(output_stream, label, pixel->thetav[iwl][ip], pixel->num_valid_meas[iwl][ip], 0, "%g", separator);
	  fprintf(output_stream, "\n");
	}
	break;
      case MEAS_TYPE_LS:
      case MEAS_TYPE_DP:
      case MEAS_TYPE_RL:
	{
	  char label[255 + 1];
	  snprintf(label, sizeof(label) - 1, "pixel->hvp[%.3f%s]", pixel->wavelengths[iwl], str_meas_type);
	  print_array(output_stream, label, pixel->hvp, pixel->num_hvp, 0, "%g", separator);
	  fprintf(output_stream, "\n");
	}
	break;
      default:
	fprintf(sdata_get_error_stream(), "%s:%d: fatal error: unexpected value for meas_type: %d\n", __FILE__, __LINE__, meas_type);
	abort();
      } /* switch */
    } /* for (ip) */
  } /* for (iwl) */


  print_array_double_3D(output_stream, "pixel->dphi", pixel, pixel->dphi, "%.3f");
  print_cmtrx(output_stream, "pixel->cmtrx", pixel, pixel->cmtrx, "%.3g");

  for (iwl = 0 ; iwl < pixel->num_wavelengths ; iwl++) {
    int nip = pixel->num_meas_types[iwl];
    for (ip = 0 ; ip < nip ; ip++) {
      int meas_type = pixel->meas_types[iwl][ip];
      const char *str_meas_type = sdata_str_meas_types(meas_type);
      switch (meas_type) {
      case MEAS_TYPE_TOD:
      case MEAS_TYPE_AOD:
	{
	  fprintf(output_stream, "pixel->tau[%.3f%s] (%d): ", pixel->wavelengths[iwl], str_meas_type, pixel->num_valid_meas[iwl][ip]);
	  println_array(output_stream, "", pixel->tau[iwl], pixel->num_valid_meas[iwl][ip], 0, "%.3e", separator);
	}
	break;
      case MEAS_TYPE_P11:
	{
	  fprintf(output_stream, "pixel->p11[%.3f%s] (%d): ", pixel->wavelengths[iwl], str_meas_type, pixel->num_valid_meas[iwl][ip]);
	  println_array(output_stream, "", pixel->p11[iwl], pixel->num_valid_meas[iwl][ip], 0, "%.3e", separator);
	}
	break;
	
      case MEAS_TYPE_P12:
	{
	  fprintf(output_stream, "pixel->p12[%.3f%s] (%d): ", pixel->wavelengths[iwl], str_meas_type, pixel->num_valid_meas[iwl][ip]);
	  println_array(output_stream, "", pixel->p12[iwl], pixel->num_valid_meas[iwl][ip], 0, "%.3e", separator);
	}
	break;
	
      case MEAS_TYPE_P22:
	{
	  fprintf(output_stream, "pixel->p22[%.3f%s] (%d): ", pixel->wavelengths[iwl], str_meas_type, pixel->num_valid_meas[iwl][ip]);
	  println_array(output_stream, "", pixel->p22[iwl], pixel->num_valid_meas[iwl][ip], 0, "%.3e", separator);
	}
	break;
	
      case MEAS_TYPE_P33:
	{
	  fprintf(output_stream, "pixel->p33[%.3f%s] (%d): ", pixel->wavelengths[iwl], str_meas_type, pixel->num_valid_meas[iwl][ip]);
	  println_array(output_stream, "", pixel->p33[iwl], pixel->num_valid_meas[iwl][ip], 0, "%.3e", separator);
	}
	break;
	
      case MEAS_TYPE_P34:
	{
	  fprintf(output_stream, "pixel->p34[%.3f%s] (%d): ", pixel->wavelengths[iwl], str_meas_type, pixel->num_valid_meas[iwl][ip]);
	  println_array(output_stream, "", pixel->p34[iwl], pixel->num_valid_meas[iwl][ip], 0, "%.3e", separator);
	}
	break;
	
      case MEAS_TYPE_P44:
	{
	  fprintf(output_stream, "pixel->p44[%.3f%s] (%d): ", pixel->wavelengths[iwl], str_meas_type, pixel->num_valid_meas[iwl][ip]);
	  println_array(output_stream, "", pixel->p44[iwl], pixel->num_valid_meas[iwl][ip], 0, "%.3e", separator);
	}
	break;
	
      case MEAS_TYPE_LS:
	{
	  fprintf(output_stream, "pixel->LS[%.3f%s] (%d): ", pixel->wavelengths[iwl], str_meas_type, pixel->num_valid_meas[iwl][ip]);
	  println_array(output_stream, "", pixel->LS[iwl], pixel->num_valid_meas[iwl][ip], 0, "%.3e", separator);
	}
	break;

      case MEAS_TYPE_DP:
	{
	  fprintf(output_stream, "pixel->DP[%.3f%s] (%d): ", pixel->wavelengths[iwl], str_meas_type, pixel->num_valid_meas[iwl][ip]);
	  println_array(output_stream, "", pixel->DP[iwl], pixel->num_valid_meas[iwl][ip], 0, "%.3e", separator);
	}
	break;

      case MEAS_TYPE_RL:
	{
	  fprintf(output_stream, "pixel->RL[%.3f%s] (%d): ", pixel->wavelengths[iwl], str_meas_type, pixel->num_valid_meas[iwl][ip]);
	  println_array(output_stream, "", pixel->RL[iwl], pixel->num_valid_meas[iwl][ip], 0, "%.3e", separator);
	}
	break;
	
      case MEAS_TYPE_I:
	{
	  fprintf(output_stream, "pixel->I[%.3f%s] (%d): ", pixel->wavelengths[iwl], str_meas_type, pixel->num_valid_meas[iwl][ip]);
	  println_array(output_stream, "", pixel->I[iwl], pixel->num_valid_meas[iwl][ip], 0, "%.3e", separator);
	}
	break;
      case MEAS_TYPE_Q:
	{
	  fprintf(output_stream, "pixel->Q[%.3f%s] (%d): ", pixel->wavelengths[iwl], str_meas_type, pixel->num_valid_meas[iwl][ip]);
	  println_array(output_stream, "", pixel->Q[iwl], pixel->num_valid_meas[iwl][ip], 0, "%.3e", separator);
	}
	break;
      case MEAS_TYPE_U:
	{
	  fprintf(output_stream, "pixel->U[%.3f%s] (%d): ", pixel->wavelengths[iwl], str_meas_type, pixel->num_valid_meas[iwl][ip]);
	  println_array(output_stream, "", pixel->U[iwl], pixel->num_valid_meas[iwl][ip], 0, "%.3e", separator);
	}
	break;
      case MEAS_TYPE_P:
	{
	  fprintf(output_stream, "pixel->P[%.3f%s] (%d): ", pixel->wavelengths[iwl], str_meas_type, pixel->num_valid_meas[iwl][ip]);
	  println_array(output_stream, "", pixel->P[iwl], pixel->num_valid_meas[iwl][ip], 0, "%.3e", separator);
	}
	break;
      default: 
	{
	  fprintf(sdata_get_error_stream(), "%s:%d: unexpected (or not supported yet) value for meas_type (%d)", __FILE__, __LINE__, meas_type);
	  abort();
	}
      } /* switch (meas_type) */
    }
  }

#ifdef SHOW_WARNINGS
#warning "TODO: display of groundpar, gaspar and mprof not yet implemented"
#endif
}

/*************************************************************************************************/

void sdata_read_array_double_1D(SDATA_HANDLE *handle, const char *label, size_t dimension, double *array, double val_min, double val_max) {
  int nread;
  int i;
  bool success;
 (void)nread; // This avoid -Wunused-but-set-variable warning

  for (i = 0 ; i < dimension ; i++) {
    nread = fscanf(handle->file_pointer, "%lf", &array[i]);
    assert(nread == 1);
  }
  
  success = check_array_double(sdata_get_debug_stream(), label, dimension, array, val_min, val_max, handle->file_name, handle->current_line);
  if (success == false) {
    fprintf(sdata_get_error_stream(), "grasp: error at %s:%d: %s out of bounds [%g .. %g]: ", handle->file_name, handle->current_line, label, val_min, val_max);
    println_array(sdata_get_error_stream(), "", array, dimension, 0, "%g", ", ");
    fprintf(sdata_get_error_stream(), "grasp: for easy SDATA debugging, you can use the sdata_dump utility provided with grasp: sdata_dump -D %s (sdata_dump can be found in the grasp/bin directory)\n", handle->file_name);
    exit (EXIT_FAILURE);
  }
  
}

/*************************************************************************************************/

void sdata_read_array_int_1D(SDATA_HANDLE *handle, const char *label, size_t dimension, int *array, int val_min, int val_max) {
  int nread;
  int i;
  bool success;
  (void)nread; // This avoid -Wunused-but-set-variable warning
  
  for (i = 0 ; i < dimension ; i++) {
    nread = fscanf(handle->file_pointer, "%d", &array[i]);
    assert(nread == 1);
  }
  success = check_array_int(sdata_get_debug_stream(), label, dimension, array, val_min, val_max, handle->file_name, handle->current_line);
  if (success == false) {
    fprintf(sdata_get_error_stream(), "grasp: error at %s:%d: %s out of bounds [%d .. %d]: ", handle->file_name, handle->current_line, label, val_min, val_max);
    println_array(sdata_get_error_stream(), "", array, dimension, 0, "%d", ", ");
    fprintf(sdata_get_error_stream(), "grasp: for easy SDATA debugging, you can use the sdata_dump utility provided with grasp: sdata_dump -D %s (sdata_dump can be found in the grasp/bin directory)\n", handle->file_name);
    exit (EXIT_FAILURE);
  }

}

/*************************************************************************************************/

void sdata_read_array_int_2D(SDATA_HANDLE *handle, const char *label, int num_wavelengths, double wavelengths[SDATA_MAX_NWL], int num_meas_types[SDATA_MAX_NWL], int array[SDATA_MAX_NWL][SDATA_MAX_NIP], int val_min, int val_max) {
  int iwl, ip;

  assert(num_wavelengths <= SDATA_MAX_NWL);

  for (iwl = 0 ; iwl < num_wavelengths ; iwl++) {
     char label2[CSTRING_MAXLEN + 1];
     bool success;

     snprintf(label2, CSTRING_MAXLEN, "%s[%.3g]", label, wavelengths[iwl]);
     
     assert(num_meas_types[iwl] <= SDATA_MAX_NIP);
     for (ip = 0 ; ip < num_meas_types[iwl] ; ip++) {
       int nread;
       nread = fscanf(handle->file_pointer, "%d", &array[iwl][ip]);
       
       if (nread != 1) {
	 fprintf(sdata_get_error_stream(), "%s:%d: ", __FILE__, __LINE__);
	 fprintf(sdata_get_error_stream(), "%s: fatal error: failed to read %s (unexpected value for nread: %d)\n", handle->file_name, label2, nread);
	 exit (EXIT_FAILURE);
       }
     }
     
     success = check_array_int(sdata_get_debug_stream(), label2, num_meas_types[iwl], array[iwl], val_min, val_max, handle->file_name, handle->current_line);
     if (success == false) {
       fprintf(sdata_get_error_stream(), "grasp: error at %s:%d: %s out of bounds [%d .. %d]: ", handle->file_name, handle->current_line, label2, val_min, val_max);
       println_array(sdata_get_error_stream(), "", array[iwl], num_meas_types[iwl], 0, "%d", ", ");
       fprintf(sdata_get_error_stream(), "grasp: for easy SDATA debugging, you can use the sdata_dump utility provided with grasp: sdata_dump -D %s (sdata_dump can be found in the grasp/bin directory)\n", handle->file_name);
       exit (EXIT_FAILURE);
     }

  }
}

/*************************************************************************************************/

void sdata_read_groundpar(SDATA_HANDLE *handle, const char *label, int num_wavelengths, double wavelengths[SDATA_MAX_NWL], int nsurf, double array[SDATA_MAX_NWL][SDATA_MAX_NSURF], double val_min, double val_max) {
  int iwl, isurf;

  assert(nsurf <= SDATA_MAX_NSURF);

  for (iwl = 0 ; iwl < num_wavelengths ; iwl++) {
    char label2[CSTRING_MAXLEN + 1];
    bool success;

    snprintf(label2, CSTRING_MAXLEN, "%s[%.3g]", label, wavelengths[iwl]);
    
    for (isurf = 0 ; isurf < nsurf ; isurf++) {
      int nread;
      (void)nread; // This avoid -Wunused-but-set-variable warning
      nread = fscanf(handle->file_pointer, "%lf", &array[iwl][isurf]);
      assert(nread == 1);
    }
  
    success = check_array_double(sdata_get_debug_stream(), label2, nsurf, array[iwl], val_min, val_max, handle->file_name, handle->current_line);
    if (success == false) {
      fprintf(sdata_get_error_stream(), "grasp: error at %s:%d: %s out of bounds [%g .. %g]: ", handle->file_name, handle->current_line, label2, val_min, val_max);
      println_array(sdata_get_error_stream(), "", array[iwl], nsurf, 0, "%g", ", ");
      fprintf(sdata_get_error_stream(), "grasp: for easy SDATA debugging, you can use the sdata_dump utility provided with grasp: sdata_dump -D %s (sdata_dump can be found in the grasp/bin directory)\n", handle->file_name);
      exit (EXIT_FAILURE);
    }

  }
}

/*************************************************************************************************/

void sdata_read_meas_types_2D(SDATA_HANDLE *handle, const char *label, int num_wavelengths, double wavelengths[SDATA_MAX_NWL], int num_meas_types[SDATA_MAX_NWL], int array[SDATA_MAX_NWL][SDATA_MAX_NIP], int val_min, int val_max) {
  int iwl, ip;

  assert(num_wavelengths <= SDATA_MAX_NWL);

  for (iwl = 0 ; iwl < num_wavelengths ; iwl++) {
     char label2[CSTRING_MAXLEN + 1];
     bool success;

     snprintf(label2, CSTRING_MAXLEN, "%s[%.3g]", label, wavelengths[iwl]);
     
     assert(num_meas_types[iwl] <= SDATA_MAX_NIP);
     for (ip = 0 ; ip < num_meas_types[iwl] ; ip++) {
       int nread;
       (void)nread; // This avoid -Wunused-but-set-variable warning
       nread = fscanf(handle->file_pointer, "%d", &array[iwl][ip]);
       assert(nread == 1);
     }
     success = check_meas_types(sdata_get_debug_stream(), label2, num_meas_types[iwl], array[iwl], val_min, val_max, handle->file_name, handle->current_line);
     if (success == false) {
       fprintf(sdata_get_error_stream(), "grasp: error at %s:%d: %s out of bounds [%d .. %d]: ", handle->file_name, handle->current_line, label2, val_min, val_max);
       println_array(sdata_get_error_stream(), "", array[iwl], num_meas_types[iwl], 0, "%d", ", ");
       fprintf(sdata_get_error_stream(), "grasp: for easy SDATA debugging, you can use the sdata_dump utility provided with grasp: sdata_dump -D %s (sdata_dump can be found in the grasp/bin directory)\n", handle->file_name);
       exit (EXIT_FAILURE);
     }
  }
}

/*************************************************************************************************/
void sdata_read_geometry(SDATA_HANDLE *handle, const char *label, SDATA_PIXEL *pixel, const SDATA_VALID_RANGES *valid_ranges) {
  int iwl;

  assert(pixel->num_wavelengths <= SDATA_MAX_NWL);
  for (iwl = 0 ; iwl < pixel->num_wavelengths ; iwl++) {
    int ip;
    
    assert(pixel->num_meas_types[iwl] <= SDATA_MAX_NIP);
    for (ip = 0 ; ip < pixel->num_meas_types[iwl] ; ip++) {
      int meas_type = pixel->meas_types[iwl][ip];
      bool success;
      int ivm;
      
      switch (meas_type) {
      case MEAS_TYPE_TOD:
      case MEAS_TYPE_AOD:
      case MEAS_TYPE_P11:
      case MEAS_TYPE_P12:
      case MEAS_TYPE_P22:
      case MEAS_TYPE_P33:
      case MEAS_TYPE_P34:
      case MEAS_TYPE_P44:
      case MEAS_TYPE_I:
      case MEAS_TYPE_Q:
      case MEAS_TYPE_U:
      case MEAS_TYPE_P:
	{
	  char label_thetav[CSTRING_MAXLEN + 1];
	  snprintf(label_thetav, CSTRING_MAXLEN, "pixel->thetav[%s%.0f]", sdata_str_meas_types(meas_type), 1000. * pixel->wavelengths[iwl]);
	     
	  assert(pixel->num_valid_meas[iwl][ip] <= SDATA_MAX_NBVM);
	  for (ivm = 0 ; ivm < pixel->num_valid_meas[iwl][ip] ; ivm++) {
	    int nread;
            (void)nread; // This avoid -Wunused-but-set-variable warning
	    nread = fscanf(handle->file_pointer, "%lf", &pixel->thetav[iwl][ip][ivm]);
	    assert(nread == 1);
	  }
	  success = check_array_double(sdata_get_debug_stream(), label_thetav, pixel->num_valid_meas[iwl][ip], 
				       pixel->thetav[iwl][ip], valid_ranges->thetav_min, valid_ranges->thetav_max, handle->file_name, handle->current_line);
	  if (success == false) {
	    fprintf(sdata_get_error_stream(), "grasp: error at %s:%d: %s out of bounds [%g .. %g]: ", handle->file_name, handle->current_line, label_thetav, valid_ranges->thetav_min, valid_ranges->thetav_max);
	    println_array(sdata_get_error_stream(), "", pixel->thetav[iwl][ip], pixel->num_valid_meas[iwl][ip], 0, "%g", ", ");
	    fprintf(sdata_get_error_stream(), "grasp: for easy SDATA debugging, you can use the sdata_dump utility provided with grasp: sdata_dump -D %s (sdata_dump can be found in the grasp/bin directory)\n", handle->file_name);
	    exit (EXIT_FAILURE);
	  }
	}
	break;
      case MEAS_TYPE_LS:
      case MEAS_TYPE_DP:
      case MEAS_TYPE_RL:
	{
	  char label_hvp[CSTRING_MAXLEN + 1];
	  snprintf(label_hvp, CSTRING_MAXLEN, "pixel->hvp[%s%.0f]", sdata_str_meas_types(meas_type), 1000. * pixel->wavelengths[iwl]);
	     
	  assert(pixel->num_valid_meas[iwl][ip] <= SDATA_MAX_KNBVM);
	  pixel->num_hvp = pixel->num_valid_meas[iwl][ip]; /* we are in a loop on iwl and ip, but the indices can be ignored (values are duplicate for all iwl and ip) */
	  
	  for (ivm = 0 ; ivm < pixel->num_valid_meas[iwl][ip] ; ivm++) {
	    int nread;
	    double hvp;
            (void)nread; // This avoid -Wunused-but-set-variable warning
	    
            assert(pixel->num_valid_meas[iwl][ip] == pixel->num_hvp); /* the SDATA files have duplicate values of altitudes for each iwl and for each ip.
								       * Don't ask me why it has been designed like that.
								       */
	    
	    nread = fscanf(handle->file_pointer, "%lf", &hvp);
	    assert(nread == 1);
	       
#ifdef SHOW_WARNINGS
#warning "SDATA files may contain duplicate values of heights of vertical profiles. We'll keep only the first set of values."
#endif
	    pixel->hvp[ivm] = hvp; /* we are in a loop on iwl, ip and ivm, but the values for iwl and ip can be ignored (values are duplicate for all iwl and ip) */
	  }
	  
	  success = check_array_double(sdata_get_debug_stream(), label_hvp, pixel->num_valid_meas[iwl][ip], 
				       pixel->hvp, valid_ranges->hvp_min, valid_ranges->hvp_max, handle->file_name, handle->current_line);
	  
	  if (success == false) {
	    fprintf(sdata_get_error_stream(), "grasp: error at %s:%d: %s out of bounds [%g .. %g]: ", handle->file_name, handle->current_line, label_hvp, valid_ranges->hvp_min, valid_ranges->hvp_max);
	    println_array(sdata_get_error_stream(), "", pixel->hvp, pixel->num_valid_meas[iwl][ip], 0, "%g", ", ");
	    fprintf(sdata_get_error_stream(), "grasp: for easy SDATA debugging, you can use the sdata_dump utility provided with grasp: sdata_dump -D %s (sdata_dump can be found in the grasp/bin directory)\n", handle->file_name);
	    exit (EXIT_FAILURE);
	  }
	}
	break;
      default: 
	{
	  fprintf(sdata_get_error_stream(), "%s:%d: %s: unexpected (or not supported yet) value for meas_type (%d) at line %d\n", __FILE__, __LINE__, handle->file_name, meas_type, handle->current_line);
	  abort();
	}
      } /* switch (meas_type) */
	 
    }
  }

  assert(pixel->num_wavelengths <= SDATA_MAX_NWL);
  for (iwl = 0 ; iwl < pixel->num_wavelengths ; iwl++) {
    int ip;
       
    assert(pixel->num_meas_types[iwl] <= SDATA_MAX_NIP);
    for (ip = 0 ; ip < pixel->num_meas_types[iwl] ; ip++) {
      int meas_type = pixel->meas_types[iwl][ip];
      bool success;
      int ivm;
	 
      switch (meas_type) {
      case MEAS_TYPE_TOD:
      case MEAS_TYPE_AOD:
      case MEAS_TYPE_P11:
      case MEAS_TYPE_P12:
      case MEAS_TYPE_P22:
      case MEAS_TYPE_P33:
      case MEAS_TYPE_P34:
      case MEAS_TYPE_P44:
      case MEAS_TYPE_LS:
      case MEAS_TYPE_DP:
      case MEAS_TYPE_RL:
      case MEAS_TYPE_I:
      case MEAS_TYPE_Q:
      case MEAS_TYPE_U:
      case MEAS_TYPE_P:
	{
	  char label_dphi[CSTRING_MAXLEN + 1];
	  snprintf(label_dphi, CSTRING_MAXLEN, "pixel->dphi[%s%.0f]", sdata_str_meas_types(meas_type), 1000. * pixel->wavelengths[iwl]);
	     
	  assert(pixel->num_valid_meas[iwl][ip] <= SDATA_MAX_KNBVM);
	  for (ivm = 0 ; ivm < pixel->num_valid_meas[iwl][ip] ; ivm++) {
	    int nread;
            (void)nread; // This avoid -Wunused-but-set-variable warning
	    nread = fscanf(handle->file_pointer, "%lf", &pixel->dphi[iwl][ip][ivm]);
	    assert(nread == 1);
	  }
	  success = check_array_double(sdata_get_debug_stream(), label_dphi, pixel->num_valid_meas[iwl][ip], 
				       pixel->dphi[iwl][ip], valid_ranges->dphi_min, valid_ranges->dphi_max, handle->file_name, handle->current_line);
	  if (success == false) {
	    fprintf(sdata_get_error_stream(), "grasp: error at %s:%d: %s out of bounds [%g .. %g]: ", handle->file_name, handle->current_line, label_dphi, valid_ranges->dphi_min, valid_ranges->dphi_max);
	    println_array(sdata_get_error_stream(), "", pixel->dphi[iwl][ip], pixel->num_valid_meas[iwl][ip], 0, "%g", ", ");
	    fprintf(sdata_get_error_stream(), "grasp: for easy SDATA debugging, you can use the sdata_dump utility provided with grasp: sdata_dump -D %s (sdata_dump can be found in the grasp/bin directory)\n", handle->file_name);
	    exit (EXIT_FAILURE);
	  }
	     
	}
	break;
      default: 
	{
	  fprintf(sdata_get_error_stream(), "%s:%d: unexpected (or not supported yet) value for meas_type (%d)", __FILE__, __LINE__, meas_type);
	  abort();
	}
      } /* switch (meas_type) */
	 
    }
  }
  
} /* sdata_read_geometry */

/*************************************************************************************************/

void sdata_read_measurements(SDATA_HANDLE *handle, const char *label, int num_wavelengths, double wavelengths[SDATA_MAX_NWL], int num_meas_types[SDATA_MAX_NWL], int meas_types[SDATA_MAX_NWL][SDATA_MAX_NIP], int num_valid_meas[SDATA_MAX_NWL][SDATA_MAX_NIP], double array[SDATA_MAX_NWL][SDATA_MAX_NIP][SDATA_MAX_KNBVM], double val_min, double val_max) {
  int iwl, ip, ivm;

  assert(num_wavelengths <= SDATA_MAX_NWL);

  for (iwl = 0 ; iwl < num_wavelengths ; iwl++) {
    assert(num_meas_types[iwl] <= SDATA_MAX_NIP);
    for (ip = 0 ; ip < num_meas_types[iwl] ; ip++) {
      char label2[CSTRING_MAXLEN + 1];
      int meas_type = meas_types[iwl][ip];
      
      snprintf(label2, CSTRING_MAXLEN, "%s[%s%.0f]", label, sdata_str_meas_types(meas_type), 1000. * wavelengths[iwl]);
      
      assert(num_valid_meas[iwl][ip] <= SDATA_MAX_KNBVM);
      for (ivm = 0 ; ivm < num_valid_meas[iwl][ip] ; ivm++) {
	int nread;
        (void)nread; // This avoid -Wunused-but-set-variable warning
	nread = fscanf(handle->file_pointer, "%lf", &array[iwl][ip][ivm]);
	assert(nread == 1);
       }
      assert( check_array_double(sdata_get_debug_stream(), label2, num_valid_meas[iwl][ip], array[iwl][ip], val_min, val_max, handle->file_name, handle->current_line) );
    }
   }
}

/*************************************************************************************************/

void sdata_read_cmtrx(SDATA_HANDLE *handle, const char *label, int num_wavelengths, double wavelengths[SDATA_MAX_NWL], int num_meas_types[SDATA_MAX_NWL], int meas_types[SDATA_MAX_NWL][SDATA_MAX_NIP], int num_valid_meas[SDATA_MAX_NWL][SDATA_MAX_NIP], double cmtrx[SDATA_MAX_NWL][SDATA_MAX_NIP][SDATA_MAX_KNBVM], double val_min, double val_max, int ifcov[SDATA_MAX_NWL][SDATA_MAX_NIP]) {
  int iwl, ip, ivm;

  assert(num_wavelengths <= SDATA_MAX_NWL);

  for (iwl = 0 ; iwl < num_wavelengths ; iwl++) {
    assert(num_meas_types[iwl] <= SDATA_MAX_NIP);
    for (ip = 0 ; ip < num_meas_types[iwl] ; ip++) {
      char label2[CSTRING_MAXLEN + 1];
      int meas_type = meas_types[iwl][ip];
      bool success;
      
      snprintf(label2, CSTRING_MAXLEN, "%s[%s%.0f]", label, sdata_str_meas_types(meas_type), 1000. * wavelengths[iwl]);
      
      assert(num_valid_meas[iwl][ip] <= SDATA_MAX_KNBVM);
      for (ivm = 0 ; ivm < num_valid_meas[iwl][ip] * ifcov[iwl][ip] ; ivm++) {
	int nread;
        (void)nread; // This avoid -Wunused-but-set-variable warning
	nread = fscanf(handle->file_pointer, "%lf", &cmtrx[iwl][ip][ivm]);
	assert(nread == 1);
      }
      success = check_array_double(sdata_get_debug_stream(), label2, num_valid_meas[iwl][ip] * ifcov[iwl][ip], cmtrx[iwl][ip], val_min, val_max, handle->file_name, handle->current_line);
      if (success == false) {
	fprintf(sdata_get_error_stream(), "grasp: error at %s:%d: %s out of bounds [%g .. %g]: ", handle->file_name, handle->current_line, label2, val_min, val_max);
	println_array(sdata_get_error_stream(), "", cmtrx[iwl][ip], num_valid_meas[iwl][ip] * ifcov[iwl][ip], 0, "%g", ", ");
	fprintf(sdata_get_error_stream(), "grasp: for easy SDATA debugging, you can use the sdata_dump utility provided with grasp: sdata_dump -D %s (sdata_dump can be found in the grasp/bin directory)\n", handle->file_name);
	exit (EXIT_FAILURE);
      }
      
    }
   }
}

/*************************************************************************************************/
void sdata_read_mprof(SDATA_HANDLE *handle, const char *label, int num_wavelengths, double wavelengths[SDATA_MAX_NWL], int num_meas_types[SDATA_MAX_NWL], int meas_types[SDATA_MAX_NWL][SDATA_MAX_NIP], int num_valid_meas[SDATA_MAX_NWL][SDATA_MAX_NIP], double mprof[SDATA_MAX_NWL][SDATA_MAX_NIP][SDATA_MAX_KVERTM], double val_min, double val_max, int ifmp[SDATA_MAX_NWL][SDATA_MAX_NIP]) {
  int iwl, ip, ivm;

  assert(num_wavelengths <= SDATA_MAX_NWL);

  for (iwl = 0 ; iwl < num_wavelengths ; iwl++) {
    assert(num_meas_types[iwl] <= SDATA_MAX_NIP);
    for (ip = 0 ; ip < num_meas_types[iwl] ; ip++) {
      char label2[CSTRING_MAXLEN + 1];
      int meas_type = meas_types[iwl][ip];
      bool success;
      
      snprintf(label2, CSTRING_MAXLEN, "%s[%s%.0f]", label, sdata_str_meas_types(meas_type), 1000. * wavelengths[iwl]);
      
      assert(num_valid_meas[iwl][ip] * ifmp[iwl][ip] <= SDATA_MAX_KVERTM);
      for (ivm = 0 ; ivm < num_valid_meas[iwl][ip] * ifmp[iwl][ip] ; ivm++) {
	int nread;
        (void)nread; // This avoid -Wunused-but-set-variable warning
	nread = fscanf(handle->file_pointer, "%lf", &mprof[iwl][ip][ivm]);
	assert(nread == 1);
      }
      success = check_array_double(sdata_get_debug_stream(), label2, num_valid_meas[iwl][ip] * ifmp[iwl][ip], mprof[iwl][ip], val_min, val_max, handle->file_name, handle->current_line);
      if (success == false) {
	fprintf(sdata_get_error_stream(), "grasp: error at %s:%d: %s out of bounds [%g .. %g]: ", handle->file_name, handle->current_line, label2, val_min, val_max);
	println_array(sdata_get_error_stream(), "", mprof[iwl][ip], num_valid_meas[iwl][ip] * ifmp[iwl][ip], 0, "%g", ", ");
	fprintf(sdata_get_error_stream(), "grasp: for easy SDATA debugging, you can use the sdata_dump utility provided with grasp: sdata_dump -D %s (sdata_dump can be found in the grasp/bin directory)\n", handle->file_name);
	exit (EXIT_FAILURE);
      }
    }
   }
}

/*************************************************************************************************/

SDATA_RECORD *sdata_alloc_record(SDATA_HANDLE *handle) {
  SDATA_RECORD *record;

  assert(handle != NULL);

  record = trackmem_malloc(sizeof(*record));
  assert(record != NULL);

  record->handle = handle;

  return record;
}

/*************************************************************************************************/

void sdata_free_record(SDATA_RECORD *record) {
  trackmem_free(record);
}

/*************************************************************************************************/

void sdata_get_date_time(SDATA_RECORD *record, char str_date_time[20 + 1]) { /* format YYYY-MM-DDThh:mm:ssZ */
  assert(record != NULL);
  assert(str_date_time != NULL);
  strncpy(str_date_time, record->str_date_time, sizeof(record->str_date_time));
  str_date_time[20] = '\0'; // paranoid statement
}

/*************************************************************************************************/

size_t sdata_get_num_pixels(SDATA_RECORD *record) {
  assert(record != NULL);
  return record->npixels;
}

/*************************************************************************************************/

const SDATA_PIXEL *sdata_get_pixel(const SDATA_RECORD *record, size_t ipixel) {
  assert(record != NULL);
  assert(0 <= ipixel && ipixel < record->npixels);
  return &record->pixels[ipixel];
}

/*************************************************************************************************/

void compute_polarized_component(double Q[SDATA_MAX_NWL][SDATA_MAX_NBVM], double U[SDATA_MAX_NWL][SDATA_MAX_NBVM], double P[SDATA_MAX_NWL][SDATA_MAX_NBVM]) {
  int iwl; /* index of wavelength */
  int ivm; /* index of valid measurement (e.g. index of direction for PARASOL) */

  /* even if only some of the SDATA_MAX_NWL indices and some of the SDATA_MAX_NBVM indices are valid,
   * it's harmless (and simpler) to copy all of them (invalid indices should never be accessed).
   */
  for (iwl = 0 ; iwl < SDATA_MAX_NWL ; iwl++) {
    for (ivm = 0 ; ivm < SDATA_MAX_NBVM ; ivm++) {
      double Q2 = Q[iwl][ivm] * Q[iwl][ivm];
      double U2 = U[iwl][ivm] * U[iwl][ivm];
      P[iwl][ivm] = sqrt(Q2 + U2);
    } /* */
  }
}

/*************************************************************************************************/

void copy_data_to_measurement_fields(SDATA_PIXEL *pixel, double data[SDATA_MAX_NWL][SDATA_MAX_NIP][SDATA_MAX_KNBVM]) {
  int iwl;
  int ip;

  for (iwl = 0 ; iwl < pixel->num_wavelengths ; iwl++) {
    int nip = pixel->num_meas_types[iwl];
    for (ip = 0 ; ip < nip ; ip++) {
      int meas_type = pixel->meas_types[iwl][ip];
      switch (meas_type) {
      case MEAS_TYPE_TOD:
      case MEAS_TYPE_AOD:
	{
	  memcpy(pixel->tau[iwl], data[iwl][ip], pixel->num_valid_meas[iwl][ip] * sizeof(pixel->tau[iwl][0]));
	}
	break;
	
      case MEAS_TYPE_P11:
	{
	  memcpy(pixel->p11[iwl], data[iwl][ip], pixel->num_valid_meas[iwl][ip] * sizeof(pixel->p11[iwl][0]));
	}
	break;
	
      case MEAS_TYPE_P12:
	{
	   memcpy(pixel->p12[iwl], data[iwl][ip], pixel->num_valid_meas[iwl][ip] * sizeof(pixel->p12[iwl][0]));
	}
	break;
	
      case MEAS_TYPE_P22:
	{
	  memcpy(pixel->p22[iwl], data[iwl][ip], pixel->num_valid_meas[iwl][ip] * sizeof(pixel->p22[iwl][0]));
	}
	break;
	
      case MEAS_TYPE_P33:
	{
	  memcpy(pixel->p33[iwl], data[iwl][ip], pixel->num_valid_meas[iwl][ip] * sizeof(pixel->p33[iwl][0]));
	}
	break;
	
      case MEAS_TYPE_P34:
	{
	  memcpy(pixel->p34[iwl], data[iwl][ip], pixel->num_valid_meas[iwl][ip] * sizeof(pixel->p34[iwl][0]));
	}
	break;
	
      case MEAS_TYPE_P44:
	{
	  memcpy(pixel->p44[iwl], data[iwl][ip], pixel->num_valid_meas[iwl][ip] * sizeof(pixel->p44[iwl][0]));
	}
	break;

      case MEAS_TYPE_LS:
	{
	  memcpy(pixel->LS[iwl], data[iwl][ip], pixel->num_valid_meas[iwl][ip] * sizeof(pixel->LS[iwl][0]));
	}
	break;

      case MEAS_TYPE_DP:
	{
	  memcpy(pixel->DP[iwl], data[iwl][ip], pixel->num_valid_meas[iwl][ip] * sizeof(pixel->DP[iwl][0]));
	}
	break;

      case MEAS_TYPE_RL:
	{
	  memcpy(pixel->RL[iwl], data[iwl][ip], pixel->num_valid_meas[iwl][ip] * sizeof(pixel->RL[iwl][0]));
	}
	break;
	
      case MEAS_TYPE_I:
	{
	  memcpy(pixel->I[iwl], data[iwl][ip], pixel->num_valid_meas[iwl][ip] * sizeof(pixel->I[iwl][0]));
	}
	break;
      case MEAS_TYPE_Q:
	{
	  memcpy(pixel->Q[iwl], data[iwl][ip], pixel->num_valid_meas[iwl][ip] * sizeof(pixel->Q[iwl][0]));
	}
	break;
      case MEAS_TYPE_U:
	{
	  memcpy(pixel->U[iwl], data[iwl][ip], pixel->num_valid_meas[iwl][ip] * sizeof(pixel->U[iwl][0]));
	}
	break;
      case MEAS_TYPE_P:
	{
	  memcpy(pixel->P[iwl], data[iwl][ip], pixel->num_valid_meas[iwl][ip] * sizeof(pixel->P[iwl][0]));
	}
	break;
      default: 
	{
	  fprintf(sdata_get_error_stream(), "%s:%d: unexpected (or not supported yet) value for meas_type (%d)", __FILE__, __LINE__, meas_type);
	  abort();
	}
      } /* switch (meas_type) */
    }
  }

}

/*************************************************************************************************/
#ifdef SHOW_WARNINGS
#warning "NOT USED YET"
#endif

bool check_heights_of_vertical_profile(FILE *log_stream, const char *label, int nelements, double values[], double val_min, double val_max, const char *filename, int line) {
  int i;
  
  if (log_stream) {
    fprintf(log_stream, "%s:%d: checking %s: ", filename, line, label);
    print_array(log_stream, "", values, nelements, 0, "%g", ", ");
    fprintf(log_stream, " --> ");
  }
  
  for (i = 0 ; i < nelements ; i++) {
    double value = values[i];
    
    if (! (val_min <= value && value <= val_max)) {
      if (log_stream) {
	fprintf(log_stream, "out of bounds [%g..%g]\n", val_min, val_max);
      }
      return false;
    }

    if (i > 0 && (values[i] >= values[i - 1])) {
      if (log_stream) {
	fprintf(log_stream, "not decreasing\n");
      }
      return false;
    }
  }
  
  if (log_stream) {
    fprintf(log_stream, "ok\n");
  }
  return true;

}

/*************************************************************************************************/

int sdata_read_next_record(SDATA_HANDLE *handle, SDATA_RECORD *next_record) {
  int nread, nread_expected;
  int ipixel; /* index on pixels */
  SDATA_VALID_RANGES *valid_ranges;
  long index_of_day;
  int err;
  bool bret;
  double measurements[SDATA_MAX_NWL][SDATA_MAX_NIP][SDATA_MAX_KNBVM]; /* measurements (e.g. Stokes Parameters I,Q,U) */
  (void)bret; // This avoid -Wunused-but-set-variable warning
  
  assert(handle != NULL);
  assert(next_record != NULL);
  
  assert(handle->current_line >= 2);

  if (handle->current_record == handle->header.nrecords) {
    if (debug_stream) {
      fprintf(sdata_get_debug_stream(), "%s:%d: ", __FILE__, __LINE__);
      fprintf(sdata_get_debug_stream(), "%s: the last record has been read\n", handle->file_name);
    }
    return SDATA_EOF;
  }

  handle->current_line++;

  valid_ranges = &handle->valid_ranges;

  memset(next_record, 0, sizeof(SDATA_RECORD));;

  nread_expected = 5;
  nread = fscanf(handle->file_pointer, "%d %s %lf %d %d",
                 &next_record->npixels, next_record->str_date_time,
                 &next_record->hobs, &next_record->nsurf, &next_record->ifgas);

  if (nread != nread_expected) {
    if (feof(handle->file_pointer)) {
      if (debug_stream) {
	fprintf(sdata_get_debug_stream(), "%s:%d: ", __FILE__, __LINE__);
	fprintf(sdata_get_debug_stream(), "grasp: %s: end of file reached\n", handle->file_name);
      }
      return SDATA_EOF;
    }
    else {
#ifdef DEBUG
      fprintf(sdata_get_error_stream(), "%s:%d: ", __FILE__, __LINE__);
#endif
      fprintf(sdata_get_error_stream(), "grasp: error at %s:%d: malformed file\n", handle->file_name, handle->current_line);
#ifdef DEBUG
      fprintf(sdata_get_error_stream(), "grasp: nread = %d (%d was expected)\n",  nread, nread_expected);
#endif
      exit(EXIT_FAILURE);
    }
  }
  
  handle->current_record++;

  err = convert_string_to_time(next_record->str_date_time, &next_record->unix_time, TIMEFMT_ISO8601);
  if (err != 0) {
    fprintf(sdata_get_error_stream(), "grasp: error at %s:%d: malformed time specification: %s (ISO8601 format expected: YYYY-MM-DDThh:mm:ssZ)\n",
	    handle->file_name, handle->current_line, next_record->str_date_time);
    exit (EXIT_FAILURE);
  }

  discard_rest_of_line(handle->file_pointer); /* ignore comments after the useful data */

  if (debug_stream) {
    fprintf(sdata_get_debug_stream(), "%s:%d: debug: npixels: %d str_date_time: %s hobs: %f nsurf: %d ifgas: %d\n",
	    handle->file_name, handle->current_line, next_record->npixels, next_record->str_date_time,
	    next_record->hobs, next_record->nsurf, next_record->ifgas);
  }

  assert(next_record->npixels <= SDATA_MAX_NX*SDATA_MAX_NY);


  bret = get_index_of_day(next_record->str_date_time, &index_of_day);
  assert(bret == true);

  for (ipixel = 0 ; ipixel < next_record->npixels ; ipixel++) {
    SDATA_PIXEL *pixel;
    bool success;
    
    pixel = &next_record->pixels[ipixel];
    pixel->nsurf = next_record->nsurf;
    pixel->ifgas  = next_record->ifgas;
    pixel->unix_time = next_record->unix_time;
    pixel->index_of_day = index_of_day;
    pixel->it = handle->current_record;
    pixel->hobs = next_record->hobs;

    handle->current_line++;
    nread = fscanf(handle->file_pointer, "%d %d %d %d %d %lf %lf %lf %lf %d",
		   &pixel->ix, &pixel->iy, &pixel->cloud_flag, 
		   &pixel->irow, &pixel->icol, &pixel->longitude, &pixel->latitude,
		   &pixel->hgr, &pixel->land_percent, &pixel->num_wavelengths);
    
    if (nread != 10) {
      fprintf(sdata_get_error_stream(), "grasp: error at %s:%d: malformed line (failed to read the 10 first fields of the record)\n", handle->file_name, handle->current_line);
      exit (EXIT_FAILURE);
    }

    assert(pixel->hgr <= pixel->hobs); /* make sure that the height of observation is higher than the height of ground */
    
    success = check_scalar_int(sdata_get_debug_stream(), "pixel->ix", pixel->ix, valid_ranges->ix_min, valid_ranges->ix_max, handle->file_name, handle->current_line);
    if (success == false) {
      fprintf(sdata_get_error_stream(), "grasp: error at %s:%d: pixel->ix out of bounds [%d .. %d]: %d\n", handle->file_name, handle->current_line, valid_ranges->ix_min, valid_ranges->ix_max, pixel->ix);
      exit (EXIT_FAILURE);
    }

    success = check_scalar_int(sdata_get_debug_stream(), "pixel->iy", pixel->iy, valid_ranges->iy_min, valid_ranges->iy_max, handle->file_name, handle->current_line);
    if (success == false) {
      fprintf(sdata_get_error_stream(), "grasp: error at %s:%d: pixel->iy out of bounds [%d .. %d]: %d\n", handle->file_name, handle->current_line, valid_ranges->iy_min, valid_ranges->iy_max, pixel->iy);
      exit (EXIT_FAILURE);
    }

    success = check_scalar_int(sdata_get_debug_stream(), "pixel->cloud_flag", pixel->cloud_flag, valid_ranges->cloud_flag_min, valid_ranges->cloud_flag_max, handle->file_name, handle->current_line);
    if (success == false) {
      fprintf(sdata_get_error_stream(), "grasp: error at %s:%d: pixel->cloud_flag out of bounds [%d .. %d]: %d\n", handle->file_name, handle->current_line, valid_ranges->cloud_flag_min, valid_ranges->cloud_flag_max, 
	      pixel->cloud_flag);
      exit (EXIT_FAILURE);
    }

    success = check_scalar_double(sdata_get_debug_stream(), "pixel->longitude", pixel->longitude, valid_ranges->longitude_min, valid_ranges->longitude_max, handle->file_name, handle->current_line);
    if (success == false) {
      fprintf(sdata_get_error_stream(), "grasp: error at %s:%d: pixel->longitude out of bounds [%g .. %g]: %g\n", handle->file_name, handle->current_line, 
	      valid_ranges->longitude_min, valid_ranges->longitude_max, pixel->longitude);
      exit (EXIT_FAILURE);
    }

    success = check_scalar_double(sdata_get_debug_stream(), "pixel->latitude", pixel->latitude, valid_ranges->latitude_min, valid_ranges->latitude_max, handle->file_name, handle->current_line);
    if (success == false) {
      fprintf(sdata_get_error_stream(), "grasp: error at %s:%d: pixel->latitude out of bounds [%g .. %g]: %g\n", handle->file_name, handle->current_line, 
	      valid_ranges->latitude_min, valid_ranges->latitude_max, pixel->latitude);
      exit (EXIT_FAILURE);
    }

    success = check_scalar_double(sdata_get_debug_stream(), "pixel->hgr", pixel->hgr, valid_ranges->hgr_min, valid_ranges->hgr_max, handle->file_name, handle->current_line);
    if (success == false) {
      fprintf(sdata_get_error_stream(), "grasp: error at %s:%d: pixel->hgr (height of ground) out of bounds [%g .. %g]: %g\n", handle->file_name, handle->current_line, 
	      valid_ranges->hgr_min, valid_ranges->hgr_max, pixel->hgr);
      exit (EXIT_FAILURE);
    }

    success = check_scalar_double(sdata_get_debug_stream(), "pixel->land_percent", pixel->land_percent, valid_ranges->land_percent_min, valid_ranges->land_percent_max, handle->file_name, handle->current_line);
    if (success == false) {
      fprintf(sdata_get_error_stream(), "grasp: error at %s:%d: pixel->land_percent out of bounds [%g .. %g]: %g\n", handle->file_name, handle->current_line, 
	      valid_ranges->land_percent_min, valid_ranges->land_percent_max, pixel->land_percent);
      exit (EXIT_FAILURE);
    }

    success = check_scalar_int(sdata_get_debug_stream(), "pixel->num_wavelengths", pixel->num_wavelengths, valid_ranges->num_wavelengths_min, valid_ranges->num_wavelengths_max, handle->file_name, handle->current_line);
    if (success == false) {
      fprintf(sdata_get_error_stream(), "grasp: error at %s:%d: pixel->num_wavelengths out of bounds [%d .. %d]: %d\n", handle->file_name, handle->current_line, 
	      valid_ranges->num_wavelengths_min, valid_ranges->num_wavelengths_max, pixel->num_wavelengths);
      exit (EXIT_FAILURE);
    }
    
    sdata_read_array_double_1D(handle, "pixel->wavelengths", pixel->num_wavelengths, pixel->wavelengths, valid_ranges->wavelength_min, valid_ranges->wavelength_max);
    
    sdata_read_array_int_1D(handle, "pixel->num_meas_types", pixel->num_wavelengths, pixel->num_meas_types, valid_ranges->num_meas_type_min, valid_ranges->num_meas_type_max);
    
    sdata_read_meas_types_2D(handle, "pixel->meas_types", pixel->num_wavelengths, pixel->wavelengths, pixel->num_meas_types, pixel->meas_types, valid_ranges->meas_type_min, valid_ranges->meas_type_max);

    sdata_read_array_int_2D(handle, "pixel->num_valid_meas", pixel->num_wavelengths, pixel->wavelengths, pixel->num_meas_types, pixel->num_valid_meas, valid_ranges->num_valid_meas_min, valid_ranges->num_valid_meas_max);
    
    sdata_read_array_double_1D(handle, "pixel->thetas", pixel->num_wavelengths, pixel->thetas, valid_ranges->thetas_min, valid_ranges->thetas_max);
    
    sdata_read_geometry(handle, "pixel", pixel, valid_ranges);
    
    sdata_read_measurements(handle, "measurements", pixel->num_wavelengths, pixel->wavelengths, pixel->num_meas_types, pixel->meas_types, pixel->num_valid_meas, measurements, valid_ranges->meas_min, valid_ranges->meas_max);
    copy_data_to_measurement_fields(pixel, measurements);
    
    sdata_read_groundpar(handle, "pixel->groundpar", pixel->num_wavelengths, pixel->wavelengths, pixel->nsurf, pixel->groundpar, valid_ranges->groundpar_min, valid_ranges->groundpar_max);
    
    if (pixel->ifgas != 0) {
      sdata_read_array_double_1D(handle, "pixel->gaspar", pixel->num_wavelengths, pixel->gaspar, valid_ranges->gaspar_min, valid_ranges->gaspar_max);
    }
    
    sdata_read_array_int_2D(handle, "pixel->ifcov", pixel->num_wavelengths, pixel->wavelengths, pixel->num_meas_types, pixel->ifcov, valid_ranges->ifcov_min, valid_ranges->ifcov_max);
    
    sdata_read_cmtrx(handle, "pixel->cmtrx", pixel->num_wavelengths, pixel->wavelengths, pixel->num_meas_types, pixel->meas_types, pixel->num_valid_meas, pixel->cmtrx, valid_ranges->cmtrx_min, valid_ranges->cmtrx_max, pixel->ifcov);
    
    sdata_read_array_int_2D(handle, "pixel->ifmp", pixel->num_wavelengths, pixel->wavelengths, pixel->num_meas_types, pixel->ifmp, valid_ranges->ifmp_min, valid_ranges->ifmp_max);
    
    sdata_read_mprof(handle, "pixel->mprof", pixel->num_wavelengths, pixel->wavelengths, pixel->num_meas_types, pixel->meas_types, pixel->num_valid_meas, pixel->mprof, valid_ranges->mprof_min, valid_ranges->mprof_max, pixel->ifmp);
    
    discard_rest_of_line(handle->file_pointer); /* ignores possible data after mprof (e.g. MNOISEI) */

#ifdef DEBUG2
    sdata_print_pixel(stderr, "pixel", pixel);
#endif
  } /* ipixel */
  
  handle->current_line++;
  discard_rest_of_line(handle->file_pointer); /* reads an empty line at the end of the current record */

  return SDATA_OK;
}

/*************************************************************************************************/

int sdata_close(SDATA_HANDLE *handle) {
  int err;

  assert(handle != NULL);

  err = fclose(handle->file_pointer);
  
  if (err != 0) {
    if (debug_stream) {
      fprintf(sdata_get_debug_stream(), "%s:%d: ", __FILE__, __LINE__);
      fprintf(sdata_get_debug_stream(), "sdata_close failed: ");
      my_perror(sdata_get_debug_stream(), handle->file_name);
    }
    trackmem_free(handle);
    return SDATA_ERROR;
  }

  if (debug_stream) {
    fprintf(sdata_get_debug_stream(), "%s:%d: ", __FILE__, __LINE__);
    fprintf(sdata_get_debug_stream(), "%s has been closed successfully\n", handle->file_name);
  }

  trackmem_free(handle);
  return SDATA_OK;
}

/*************************************************************************************************/

void sdata_dump_record(FILE *output_stream, SDATA_RECORD *record) {
  const char *record_comments = ": NPIXELS TIMESTAMP HOBS NSURF IFGAS";
  const char *str_format = "%d%s%s%s%g%s%d%s%d%s%s\n";
  assert(record != NULL);

  fprintf(output_stream, str_format,
	  record->npixels, separator, record->str_date_time, separator,
	  record->hobs, separator, record->nsurf, separator, record->ifgas, separator, record_comments);
}

/*************************************************************************************************/

void sdata_dump_pixel(FILE *output_stream, const SDATA_PIXEL *pixel) {
  const char *str_format = "%d%s%d%s%d%s%d%s%d%s%.3f%s%.3f%s%.3f%s%.0f%s%d";
  const char *sep = separator;
  int iwl;
  int ip;

  assert(output_stream != NULL);
  assert(pixel != NULL);
  assert(pixel->num_wavelengths <= SDATA_MAX_NWL);
  
  fprintf(output_stream, str_format, 
	  pixel->ix, sep, pixel->iy, sep, pixel->cloud_flag, sep, pixel->irow, sep, pixel->icol, sep,
	  pixel->longitude, sep, pixel->latitude, sep, pixel->hgr, sep, pixel->land_percent, sep,
	  pixel->num_wavelengths);

  fprintf(output_stream, separator);
  print_array(output_stream, "", pixel->wavelengths, pixel->num_wavelengths, 0, "%.4f", separator);
  
  fprintf(output_stream, separator);
  print_array(output_stream, "", pixel->num_meas_types, pixel->num_wavelengths, 0, "%d", separator);
  
  for (iwl = 0 ; iwl < pixel->num_wavelengths ; iwl++) {
    fprintf(output_stream, separator);
    print_array(output_stream, "", pixel->meas_types[iwl], pixel->num_meas_types[iwl], 0, "%d", separator);
  }

  for (iwl = 0 ; iwl < pixel->num_wavelengths ; iwl++) {
    fprintf(output_stream, separator);
    print_array(output_stream, "", pixel->num_valid_meas[iwl], pixel->num_meas_types[iwl], 0, "%d", separator);
  }
  
  fprintf(output_stream, separator);
  print_array(output_stream, "", pixel->thetas, pixel->num_wavelengths, 0, "%f", separator);

  for (iwl = 0 ; iwl < pixel->num_wavelengths ; iwl++) {
    for (ip = 0 ; ip < pixel->num_meas_types[iwl] ; ip++) {
      int meas_type = pixel->meas_types[iwl][ip];
      
      switch (meas_type) {
      case MEAS_TYPE_TOD:
      case MEAS_TYPE_AOD:
      case MEAS_TYPE_P11:
      case MEAS_TYPE_P12:
      case MEAS_TYPE_P22:
      case MEAS_TYPE_P33:
      case MEAS_TYPE_P34:
      case MEAS_TYPE_P44:
      case MEAS_TYPE_I:
      case MEAS_TYPE_Q:
      case MEAS_TYPE_U:
      case MEAS_TYPE_P:
	{
	  fprintf(output_stream, separator);
	  print_array(output_stream, "", pixel->thetav[iwl][ip], pixel->num_valid_meas[iwl][ip], 0, "%g", separator);
	}
	break;
      case MEAS_TYPE_LS:
      case MEAS_TYPE_DP:
      case MEAS_TYPE_RL:
	{
	  fprintf(output_stream, separator);
	  print_array(output_stream, "", pixel->hvp, pixel->num_hvp, 0, "%g", separator);
	}
	break;
      default:
	fprintf(sdata_get_error_stream(), "%s:%d: fatal error: unexpected value for meas_type: %d\n", __FILE__, __LINE__, meas_type);
	abort();
      } /* switch */
    } /* for (ip) */
  } /* for (iwl) */

  for (iwl = 0 ; iwl < pixel->num_wavelengths ; iwl++) {
    for (ip = 0 ; ip < pixel->num_meas_types[iwl] ; ip++) {
      fprintf(output_stream, separator);
      print_array(output_stream, "", pixel->dphi[iwl][ip], pixel->num_valid_meas[iwl][ip], 0, "%g", separator);
    }
  }

  for (iwl = 0 ; iwl < pixel->num_wavelengths ; iwl++) {
    for (ip = 0 ; ip < pixel->num_meas_types[iwl] ; ip++) {
      int meas_type;
      int nbvm;
      fprintf(output_stream, separator);
      
      meas_type = pixel->meas_types[iwl][ip];
      nbvm = pixel->num_valid_meas[iwl][ip];

      switch (meas_type) {
      case MEAS_TYPE_TOD:
      case MEAS_TYPE_AOD:
	print_array(output_stream, "", pixel->tau[iwl], nbvm, 0, "%g", separator);
	break;
      case MEAS_TYPE_P11:
	print_array(output_stream, "", pixel->p11[iwl], nbvm, 0, "%g", separator);
	break;
      case MEAS_TYPE_P12:
	print_array(output_stream, "", pixel->p12[iwl], nbvm, 0, "%g", separator);
	break;
      case MEAS_TYPE_P22:
	print_array(output_stream, "", pixel->p22[iwl], nbvm, 0, "%g", separator);
	break;
      case MEAS_TYPE_P33:
	print_array(output_stream, "", pixel->p33[iwl], nbvm, 0, "%g", separator);
	break;
      case MEAS_TYPE_P34:
	print_array(output_stream, "", pixel->p34[iwl], nbvm, 0, "%g", separator);
	break;
      case MEAS_TYPE_P44:
	print_array(output_stream, "", pixel->p44[iwl], nbvm, 0, "%g", separator);
	break;
      case MEAS_TYPE_LS:
	print_array(output_stream, "", pixel->LS[iwl], nbvm, 0, "%g", separator);
	break;
      case MEAS_TYPE_DP:
	print_array(output_stream, "", pixel->DP[iwl], nbvm, 0, "%g", separator);
	break;
      case MEAS_TYPE_RL:
	print_array(output_stream, "", pixel->RL[iwl], nbvm, 0, "%g", separator);
	break;
      case MEAS_TYPE_I:
	print_array(output_stream, "", pixel->I[iwl], nbvm, 0, "%g", separator);
	break;
      case MEAS_TYPE_Q:
	print_array(output_stream, "", pixel->Q[iwl], nbvm, 0, "%g", separator);
	break;
      case MEAS_TYPE_U:
	print_array(output_stream, "", pixel->U[iwl], nbvm, 0, "%g", separator);
	break;
      case MEAS_TYPE_P:
	print_array(output_stream, "", pixel->P[iwl], nbvm, 0, "%g", separator);
	break;
      default:
	fprintf(sdata_get_error_stream(), "%s:%d: fatal error: unexpected value for meas_type: %d\n", __FILE__, __LINE__, meas_type);
	abort();
      }
      
    }
  }

  /* THIS MAY BE REPLACED */
  
  for (iwl = 0 ; iwl < pixel->num_wavelengths ; iwl++) {
    fprintf(output_stream, separator);
    print_array(output_stream, "", pixel->groundpar[iwl], pixel->nsurf, 0, "%g", separator);
  }
  
  /* BY THIS (AS LONG AS THIS IS NOT CONFIRMED, WE'LL KEEP THE FORMER IMPLEMENTATION) */
#if 0
  for (iwl = 0 ; iwl < pixel->num_wavelengths ; iwl++) {
    fprintf(output_stream, separator);
    fprintf(output_stream, "%.1f%s" , pixel->groundpar[iwl][0] /* one anticipates one groundpar value per wavelength */, separator);
  }
#endif
  
  if(pixel->ifgas>0){
    for (iwl = 0 ; iwl < pixel->num_wavelengths ; iwl++) {
      fprintf(output_stream, separator);
      fprintf(output_stream, "%f%s" , pixel->gaspar[iwl], separator);
    }
  }
  
  for (iwl = 0 ; iwl < pixel->num_wavelengths ; iwl++) {
    fprintf(output_stream, separator);
    print_array(output_stream, "", pixel->ifcov[iwl], pixel->num_meas_types[iwl], 0, "%d", separator);
  }

  for (iwl = 0 ; iwl < pixel->num_wavelengths ; iwl++) {
    for (ip = 0 ; ip < pixel->num_meas_types[iwl] * pixel->ifcov[iwl][ip] ; ip++) {
      fprintf(output_stream, separator);
      print_array(output_stream, "", pixel->cmtrx[iwl][ip], pixel->num_valid_meas[iwl][ip] * pixel->ifcov[iwl][ip], 0, "%g", separator);
    }
  }

  for (iwl = 0 ; iwl < pixel->num_wavelengths ; iwl++) {
    fprintf(output_stream, separator);
    print_array(output_stream, "", pixel->ifmp[iwl], pixel->num_meas_types[iwl], 0, "%d", separator);
  }
  
  for (iwl = 0 ; iwl < pixel->num_wavelengths ; iwl++) {
    for (ip = 0 ; ip < pixel->num_meas_types[iwl]  * pixel->ifmp[iwl][ip] ; ip++) {
      fprintf(output_stream, separator);
      print_array(output_stream, "", pixel->mprof[iwl][ip], pixel->num_valid_meas[iwl][ip] * pixel->ifmp[iwl][ip], 0, "%g", separator);
    }
  }
  fprintf(output_stream, "\n");

}

/*************************************************************************************************/

int sdata_dump_handle(FILE *output_stream, SDATA_HANDLE *handle) {
  SDATA_HEADER header;
  SDATA_RECORD record;
  int err;

  assert(handle != NULL);
  err = sdata_get_header(handle, &header);
  assert(err == SDATA_OK);
  sdata_dump_header(output_stream, &header);
  
  while(true) {
    size_t npixels;
    size_t ipixel;

    err = sdata_read_next_record(handle, &record);
    if (err != SDATA_OK) break;
    sdata_dump_record(output_stream, &record);
    npixels = sdata_get_num_pixels(&record);
    
    for (ipixel = 0 ; ipixel < npixels ; ipixel++) {
      const SDATA_PIXEL *pixel = sdata_get_pixel(&record, ipixel);
      assert(pixel != NULL);
      sdata_dump_pixel(output_stream, pixel);
      fprintf(output_stream, "\n");
    }
  }
  if (err != SDATA_EOF) {
    if (debug_stream) {
      fprintf(sdata_get_debug_stream(), "%s:%d: ", __FILE__, __LINE__);
      fprintf(sdata_get_debug_stream(), "%s: unexpected end of loop (end of file was not reached)\n", handle->file_name);
    }
    return SDATA_ERROR;
  }
  
  return SDATA_OK;
}

/*************************************************************************************************/

int sdata_dump_file(FILE *output_stream, const char *sdata_file) {
  int err, err_close;
  SDATA_HANDLE *handle;
  
  assert(output_stream != NULL);
  handle = sdata_open(sdata_file);
  if (handle == NULL) {
    if (debug_stream) {
      fprintf(sdata_get_debug_stream(), "%s:%d: sdata_dump_file(output_stream = %p, \"%s\") failed: ", __FILE__, __LINE__, output_stream, sdata_file);
      my_perror(sdata_get_debug_stream(), sdata_file);
    }
    return SDATA_ERROR;
  }

  err = sdata_dump_handle(output_stream, handle);
  
  err_close = sdata_close(handle);
  if ((err == SDATA_OK) && (err_close == SDATA_OK)) {
    err = SDATA_OK;
  }
  else {
    err = SDATA_ERROR;
  }

  return err;
}


static bool sdata_load_box_aux(SDATA_HANDLE *handle, const grasp_driver_settings_t *settings, grasp_box_t *box) {
  size_t irecord;
  int xmin, xmax;
  int ymin, ymax;
  int zmin, zmax; /* time in days since the reference day (1st jan 1970) */
  int err;
  
  assert(handle != NULL);
  assert(box != NULL);

  xmin = settings->xmin;
  xmax = settings->xmax;
  ymin = settings->ymin;
  ymax = settings->ymax;
  zmin = settings->gmt_time_min / 86400;
  zmax = (settings->gmt_time_max + 86400) / 86400;

  irecord = 0;
  while (true) {
    SDATA_RECORD record; 
    size_t npixels;
    size_t ipixel;
    
    
    err = sdata_read_next_record(handle, &record);
    if (err != SDATA_OK) break; /* error or end-of-file */
    
    npixels = sdata_get_num_pixels(&record);

    if (debug_stream) {
      fprintf(sdata_get_debug_stream(), "%s:%d: ", __FILE__, __LINE__);
      fprintf(sdata_get_debug_stream(), "str_date_time: %s npixels: %d\n", record.str_date_time, (int) npixels);
    }

    for (ipixel = 0 ; ipixel < npixels ; ipixel++) {
      grasp_box_vector_t tpixel;
      const SDATA_PIXEL *spixel = sdata_get_pixel(&record, ipixel);
      assert(spixel != NULL);
      
      tpixel.x = spixel->icol;
      tpixel.y = spixel->irow;
      tpixel.z = record.unix_time / 86400;
      
      if (debug_stream) {
	fprintf(sdata_get_debug_stream(), "%s:%d: ", __FILE__, __LINE__);
	fprintf(sdata_get_debug_stream(), "tpixel.x: %d tpixel.y: %d tpixel.z: %d (unix_time: %ld)\n", tpixel.x, tpixel.y, tpixel.z, (long) record.unix_time);
      }

      if ((xmin <= tpixel.x && tpixel.x <= xmax)
	  && (ymin <= tpixel.y && tpixel.y <= ymax)
	  && (zmin <= tpixel.z && tpixel.z <= zmax)) {
	grasp_box_attach_data_to_pixel(box, &tpixel, sizeof(SDATA_PIXEL), spixel);
	irecord++;

	if (debug_stream) {
	  char cstr_time[255 + 1];
	  char cstr_time_min[255 + 1];
	  char cstr_time_max[255 + 1];
	  time_t tmin, tmax; /* time in seconds since the reference day */
	  time_t tmeasure;
	  
	  tmin = 86400*zmin;
	  tmax = 86400*zmax;
	  tmeasure = record.unix_time;
	  
	  strncpy(cstr_time, time_to_string(tmeasure, NULL), 255);
	  strncpy(cstr_time_min, time_to_string(tmin, NULL), 255);
	  strncpy(cstr_time_max, time_to_string(tmax, NULL), 255);
	  
	  fprintf(sdata_get_debug_stream(), "%s:%d: pixel[row: %d col: %d time: %s (%lu) -- %s (%lu) .. %s (%lu)] registered\n", 
		  __FILE__, __LINE__, tpixel.y, tpixel.x, cstr_time, tmeasure,
		  cstr_time_min, tmin, cstr_time_max, tmax);
	}
      } /* if (xmin <= ...) */
      else {
	if (debug_stream) {
	  fprintf(sdata_get_debug_stream(), "tpixel.x: %d tpixel.y: %d tpixel.z: %d (unix_time: %ld) not in the ranges x: [%d..%d], y: [%d..%d], z: [%d..%d]. Skipped.\n", tpixel.x, tpixel.y, tpixel.z, (long) record.unix_time, xmin, xmax, ymin, ymax, zmin, zmax);
	}
      }
    } /* for (ipixel) */
    
  } /* while (true) */
  
  if (err != SDATA_EOF) {
    if (debug_stream) {
      fprintf(sdata_get_debug_stream(), "%s:%d: ", __FILE__, __LINE__);
      fprintf(sdata_get_debug_stream(), "%s: unexpected end of loop (end of file was not reached)\n", handle->file_name);
    }
    return false;
  }

  if (debug_stream) {
    grasp_box_print(sdata_get_debug_stream(), "box", box);
  }
  return true;
} /* sdata_load_box_aux */

/*************************************************************************************************/

grasp_box_t *sdata_load_box(SDATA_HANDLE *handle, const grasp_driver_settings_t *user_settings) {
  SDATA_HEADER header;
  grasp_box_t *box;
  grasp_box_settings_t box_settings;
  const grasp_driver_settings_t *user_settings_;
  grasp_driver_settings_t default_settings;
  int err;
  (void)err; // This avoid -Wunused-but-set-variable warning
  
  err = sdata_get_header(handle, &header);
  assert(err == SDATA_OK);
  if (debug_stream) {
    fprintf(sdata_get_debug_stream(), "%s:%d: header[%s]\n", __FILE__, __LINE__, handle->file_name);
    sdata_print_header(sdata_get_debug_stream(), NULL, &header);
  }
  
  default_settings.xmin = header.xmin;
  default_settings.xmax = header.xmax;
  default_settings.ymin = header.ymin;
  default_settings.ymax = header.ymax;
  default_settings.gmt_time_min = header.time_limits.first_record.unix_time;
  default_settings.gmt_time_max = header.time_limits.last_record.unix_time;
  default_settings.missing_value = MISSING_VALUE;
  
  if (user_settings == NULL) {
    user_settings_ = &default_settings;
  }
  else {
    user_settings_ = user_settings;
  }
  
  if (! (user_settings_->xmin <= user_settings_->xmax)) {
    fprintf(sdata_get_error_stream(), "%s:%d: ", __FILE__, __LINE__);
    fprintf(sdata_get_error_stream(), "__FUNCTION__: invalid x-range: [xmin..xmax] = [%d..%d]\n", user_settings_->xmin, user_settings_->xmax);
    abort();
    return NULL;
  }

  if (! (user_settings_->ymin <= user_settings_->ymax)) {
    fprintf(sdata_get_error_stream(), "%s:%d: ", __FILE__, __LINE__);
    fprintf(sdata_get_error_stream(), "__FUNCTION__: invalid x-range: [ymin..ymax] = [%d..%d]\n", user_settings_->ymin, user_settings_->ymax);
    abort();
    return NULL;
  }

  if (! (user_settings_->gmt_time_min <= user_settings_->gmt_time_max)) {
    char cstr_time_min[255 + 1];
    char cstr_time_max[255 + 1];
    time_t tmin, tmax; /* time in seconds since the reference day */
	  
    tmin = user_settings_->gmt_time_min;
    tmax = user_settings_->gmt_time_max;
    
    strncpy(cstr_time_min, time_to_string(tmin, NULL), 255);
    strncpy(cstr_time_max, time_to_string(tmax, NULL), 255);
    fprintf(sdata_get_error_stream(), "%s:%d: ", __FILE__, __LINE__);
    fprintf(sdata_get_error_stream(), "%s: invalid time range: [tmin..tmax] = [%s (%lu s).. %s (%lu s)]\n", 
	    __FUNCTION__, cstr_time_min, (long unsigned) tmin, cstr_time_max, (long unsigned) tmax);
    abort();
    return NULL;
  }
  
  grasp_driver_settings_convert(user_settings_, &box_settings);

  box = grasp_box_new(sizeof(SDATA_PIXEL), &box_settings);
  
  if (box == NULL) {
    if (debug_stream) {
      fprintf(sdata_get_debug_stream(), "%s:%d: ", __FILE__, __LINE__);
      fprintf(sdata_get_debug_stream(), "grasp_box_new(sizeof(SDATA_PIXEL), &box_settings) failed\n");
    }
    return NULL;
  }

  sdata_load_box_aux(handle, user_settings_, box);
  
  return box;

}

/*************************************************************************************************/
#ifndef UNIT_TEST
#define UNIT_TEST 0
#endif

#if (UNIT_TEST != 0)

int main(int argc, char *argv[]) {
  const char *sdata_file;

  if (argc != 2) {
    fprintf(stderr, "usage: %s <sdata_file>\n", argv[0]);
    exit (EXIT_FAILURE);
  }

  sdata_file = argv[1];
  sdata_dump_file(stdout, sdata_file);
  return EXIT_SUCCESS;
}
#endif
