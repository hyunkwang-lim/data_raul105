/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

/**
 * @file  check_helpers.h
 * @author Fabrice Ducos <fabrice.ducos@univ-lille1.fr>
 *
 */

#ifndef CHECK_HELPERS_
#define CHECK_HELPERS_

#include <stdio.h>
#include <stdbool.h>

extern bool check_scalar_int(FILE *log_stream, const char *label, int value, int val_min, int val_max, const char *filename, int line);
extern bool check_scalar_double(FILE *log_stream, const char *label, double value, double val_min, double val_max, const char *filename, int line);
extern bool check_array_int(FILE *log_stream, const char *label, int nelements, int values[], int val_min, int val_max, const char *filename, int line);
extern bool check_array_double(FILE *log_stream, const char *label, int nelements, double values[], double val_min, double val_max, const char *filename, int line);
extern bool check_meas_types(FILE *log_stream, const char *label, int nelements, int values[], int val_min, int val_max, const char *filename, int line);

#endif /* CHECK_HELPERS_ */
