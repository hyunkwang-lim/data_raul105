/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

/**
 * @file  time_utils.h
 * @author Fabrice Ducos <fabrice.ducos@univ-lille1.fr>
 *
 * embedded version for GRASP retrieval
 */

#ifndef TIME_UTILS_H
#define TIME_UTILS_H

#include <stdbool.h>
#include <time.h>

/* values accepted by the argument time_format of convert_string_to_time() 
 * In the current version, only a subset of ISO8601 is implemented (times ending with Z suffix, i.e. GMT times)
 */
enum {
  TIMEFMT_ISO8601, /* format "YYYY-MM-DDThh:mm:ssZ" */
  TIMEFMT_COMPACT  /* format "YYYYMMDDhhmmss" */
};

/* returns 0 in case of success, -1 in case of failure 
 * author: F. Ducos
 */
extern int convert_string_to_time(const char *str_time, time_t *time /* seconds since 1970-01-01T00-00-00Z */, int time_format);


/* This routine converts a unix timestamp into a readable string.
 * The result is stored in a static buffer.
 * Because of this, the function is easy and natural to use in non multithreaded applications, 
 * but is not reentrant and especially not thread-safe.
 * Make sure the output won't be larger than TIME_TO_STRING_BUFSIZE characters (255 by default)
 *
 * never use two calls of this function in the same function call
 * always save or display the output pointee (not simply the pointer) before calling
 * this function again
 *
 * For advanced usages (especially in multithreaded environments), see time_to_string_r
 *
 * author: F. Ducos
 */
extern const char *time_to_string(time_t numeric_time, const char *time_format);

/* a "reentrant" (actually thread-safe) version of time_to_string
 * author: F. Ducos
 */
extern const char *time_to_string_r(time_t numeric_time, const char *time_format /* look for strftime specification*/,
			     size_t cstr_size, char *cstr_time);

/* provides the number of days since the Unix Epoch (1st jan 1970 00:00:00)
 * returns true in case of success, false in case of failure (wrong date specification)
 *
 * author: F. Ducos
 */
extern bool get_index_of_day(const char *cstr_time /* format "YYYY-MM-DD" or "YYYY-MM-DDThh:mm:ssZ" */, long *index_of_day);


/* Change de format of a string date. example: 20:12:2013 -> 2013-12-02
 * author: D. Fuertes
 */
char *dateDDMMAAAA2AAAAMMDD(const char *original_date);

/* Return julian date of specific date
 * author: D. Fuertes
 */
float getJulianDate(time_t x);


#endif /* TIME_UTILS_H */
