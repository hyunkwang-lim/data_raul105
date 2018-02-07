/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

/**
 * @file  sdata-impl.h
 * @author Fabrice Ducos <fabrice.ducos@univ-lille1.fr>
 *
 */

#ifndef SDATA_IMPL_H
#define SDATA_IMPL_H

#include <time.h>
#include "sdata.h"

extern int sdata_dump_handle(FILE *output_stream, SDATA_HANDLE *handle);
extern int sdata_read_next_record(SDATA_HANDLE *handle, SDATA_RECORD *record);
extern void sdata_get_date_time(SDATA_RECORD *record, char str_date_time[20 + 1]); /* format YYYY-MM-DDThh:mm:ssZ */
extern size_t sdata_get_num_pixels(SDATA_RECORD *record);
extern const SDATA_PIXEL *sdata_get_pixel(const SDATA_RECORD *record, size_t ipixel);


typedef struct SDATA_VALID_RANGES_ {
  int    ix_min, ix_max;
  int    iy_min, iy_max;
  int    cloud_flag_min, cloud_flag_max;
  double longitude_min, longitude_max;
  double latitude_min, latitude_max;
  double hgr_min, hgr_max;
  double land_percent_min, land_percent_max;
  int    num_wavelengths_min, num_wavelengths_max;
  double wavelength_min, wavelength_max;
  int    num_meas_type_min, num_meas_type_max;
  int    meas_type_min, meas_type_max;
  int    num_valid_meas_min, num_valid_meas_max;
  double thetas_min, thetas_max;
  double thetav_min, thetav_max;
  double dphi_min,   dphi_max;
  double hvp_min, hvp_max;
  double groundpar_min, groundpar_max;
  double gaspar_min, gaspar_max;
  double meas_min,  meas_max;
  int    ifcov_min, ifcov_max;
  int    ifmp_min,  ifmp_max;
  double cmtrx_min, cmtrx_max;
  double mprof_min, mprof_max;
} SDATA_VALID_RANGES;

struct SDATA_HANDLE_ {
  char file_name[255 + 1];
  FILE *file_pointer;
  SDATA_HEADER header;
  SDATA_VALID_RANGES valid_ranges;
  int current_line;
  int current_record;
};

struct SDATA_RECORD_ {
  int npixels;
  double sdata_time; /* number of days since a reference day (e.g. 1st january of the current year), + fraction of the day */
  time_t unix_time; /* time of acquisition, in seconds since 1970-01-01T00:00:00Zday of year */
  char str_date_time[20 + 1]; /* date of the acquisition, format YYYY-MM-DDThh:mm:ssZ (redundant with unix_time, but provided for convenience and control) */
  double hobs;
  int nsurf;
  int ifgas;
  SDATA_HANDLE *handle;
  SDATA_PIXEL pixels[SDATA_MAX_NX*SDATA_MAX_NY];
};

#endif /* SDATA_IMPL_H */
