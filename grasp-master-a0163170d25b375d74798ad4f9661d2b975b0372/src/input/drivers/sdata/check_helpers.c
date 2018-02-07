/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

/**
 * @file  check_helpers.c
 * @author Fabrice Ducos <fabrice.ducos@univ-lille1.fr>
 *
 */

#include <grasp/utils/print_array.h>
#include "check_helpers.h"
#include "sdata-impl.h"

/*
#define DEBUG2
#ifdef DEBUG2
#define DEFAULT_VERBOSE true
#else
#define DEFAULT_VERBOSE false
#endif

#warning "TODO: set verbose from settings"
const bool verbose = DEFAULT_VERBOSE;
*/

bool check_scalar_int(FILE *log_stream, const char *label, int value, int val_min, int val_max, const char *filename, int line) {

  if (log_stream) {
    fprintf(log_stream, "%s:%d: checking %s: %d --> ", filename, line, label, value);
  }

  if (! (val_min <= value && value <= val_max)) {
    if (log_stream) {
      fprintf(log_stream, "out of bounds [%d..%d]\n", val_min, val_max);
    }
    return false;
  }

  if (log_stream) {
    fprintf(log_stream, "ok\n");
  }
  return true;
}

bool check_scalar_double(FILE *log_stream, const char *label, double value, double val_min, double val_max, const char *filename, int line) {

  if (log_stream) {
    fprintf(log_stream, "%s:%d: checking %s: %g --> ", filename, line, label, value);
  }

  if (! (val_min <= value && value <= val_max)) {
    if (log_stream) {
      fprintf(log_stream, "out of bounds [%g..%g]\n", val_min, val_max);
    }
    return false;
  }

  if (log_stream) {
    fprintf(log_stream, "ok\n");
  }
  return true;
}

bool check_array_int(FILE *log_stream, const char *label, int nelements, int values[], int val_min, int val_max, const char *filename, int line) {
  int i;

  if (log_stream) {
    fprintf(log_stream, "%s:%d: checking %s: ", filename, line, label);
    print_array(log_stream, "", values, nelements, 0, "%d", ", ");
    fprintf(log_stream, " --> ");
  }

  for (i = 0 ; i < nelements ; i++) {
    int value = values[i];
    
    if (! (val_min <= value && value <= val_max)) {
      if (log_stream) {
	fprintf(log_stream, "out of bounds [%d..%d]\n", val_min, val_max);
      }
      return false;
    }
  }

  if (log_stream) {
    fprintf(log_stream, "ok\n");
  }
  return true;
}

bool check_array_double(FILE *log_stream, const char *label, int nelements, double values[], double val_min, double val_max, const char *filename, int line) {
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
  }

  if (log_stream) {
    fprintf(log_stream, "ok\n");
  }
  return true;
}

bool check_meas_types(FILE *log_stream, const char *label, int nelements, int values[], int val_min, int val_max, const char *filename, int line) {
  int i;
  int value;

  if (log_stream) {
    fprintf(log_stream, "%s:%d: checking %s: ", filename, line, label);
    
    for (i = 0 ; i < nelements - 1 ; i++) {
      value = values[i];
      fprintf(log_stream, "%s, ", sdata_str_meas_types(value));
    }
    value = values[nelements - 1];
    fprintf(log_stream, "%s --> ", sdata_str_meas_types(value));
  }

  for (i = 0 ; i < nelements ; i++) {
    double value = values[i];
    
    if (! (val_min <= value && value <= val_max)) {
      if (log_stream) {
	fprintf(log_stream, "out of bounds [%s..%s]\n", sdata_str_meas_types(val_min), sdata_str_meas_types(val_max));
      }
      return false;
    }
  }

  if (log_stream) {
    fprintf(log_stream, "ok\n");
  }
  return true;
}
