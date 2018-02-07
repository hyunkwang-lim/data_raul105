/**
 * @file: time_utils.c 
 * @Author: Fabrice Ducos <fabrice.ducos@univ-lille1.fr>
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

#include "utils/time_utils.h"

/* must be compatible with the enum in time_utils.h */
static const char *str_time_formats[] = {
  "%4d-%2d-%2dT%2d:%2d:%2dZ", /* TIMEFMT_ISO8601 */
  "%4d%2d%2d%2d%2d%2d"        /* TIMEFMT_COMPACT */
};

/* a portable version of GNU timgm, from the man page of timegm */
static time_t my_timegm (struct tm *tm) {
  time_t ret;
  char *tz;

  tz = getenv("TZ");
  setenv("TZ", "UTC", 1);
  tzset();
  ret = mktime(tm);
  if (tz)
    setenv("TZ", tz, 1);
  else
    unsetenv("TZ");
  tzset();
  return ret;
}

int convert_string_to_time(const char *str_time, time_t *time /* seconds since 1970-01-01T00:00:00Z */, int time_format) {
  struct tm obj_tm;
  const char *str_time_format;

  assert(str_time != NULL);
  assert(str_time[0] != '\0');
  assert(time != NULL);

  assert((time_format == TIMEFMT_ISO8601) || (time_format == TIMEFMT_COMPACT));
  str_time_format = str_time_formats[time_format];

  if (sscanf(str_time, str_time_format,
             &obj_tm.tm_year, &obj_tm.tm_mon, &obj_tm.tm_mday,
             &obj_tm.tm_hour, &obj_tm.tm_min, &obj_tm.tm_sec) != 6) {
#ifdef DEBUG
    fprintf(stderr, "\"%s\" couldn't be converted to a valid time\n", str_time);
#endif
    return -1;
  }

  obj_tm.tm_year -= 1900;
  obj_tm.tm_mon--;
  obj_tm.tm_isdst = 0; /* day saving time flag, not relevant for UTC times */

  *time = my_timegm(&obj_tm);

  if (*time == -1) {
#ifdef DEBUG
    fprintf(stderr, "\"%s\" couldn't be converted to a valid time\n", str_time);
#endif
    return -1;
  }

  return 0;
}

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
 */
const char *time_to_string(time_t numeric_time, const char *time_format /* look for strftime specification */) {
#ifndef TIME_TO_STRING_BUFSIZE
#define TIME_TO_STRING_BUFSIZE 255
#endif
  static char cstr_time[TIME_TO_STRING_BUFSIZE + 1];
  const char *time_format_;

  if (time_format == NULL) {
    time_format_ = "%FT%H:%M:%SZ";
  }
  else {
    time_format_ = time_format;
  }

  memset(cstr_time, 0, sizeof(cstr_time));
  strftime(cstr_time, sizeof(cstr_time) - 1, time_format_, gmtime(&numeric_time));
  
  return cstr_time;
}

/* a "reentrant" (actually thread-safe) version of time_to_string */
const char *time_to_string_r(time_t numeric_time, const char *time_format /* look for strftime specification*/,
			     size_t cstr_size, char *cstr_time) {
  const char *time_format_;

  if (time_format == NULL) {
    time_format_ = "%FT%H:%M:%SZ";
  }
  else {
    time_format_ = time_format;
  }

  memset(cstr_time, 0, cstr_size);
  strftime(cstr_time, cstr_size, time_format_, gmtime(&numeric_time));

  return cstr_time;
}

bool get_index_of_day(const char *cstr_time /* format "YYYY-MM-DD" or "YYYY-MM-DDThh:mm:ssZ" */, long *index_of_day) {
  /* NOTE 1: actually cstr_time can be anything starting with "YYYY-MM-DD", but this
   * behaviour is not documented, and future versions of the routine may be stricter
   * on what it accepts as input.
   *
   * NOTE 2: the routine is known to work even for dates older than the 1st january 1970, at least
   * on some implementations. In that case it will return negative indices. This feature has not been
   * tested extensively yet so one shouldn't rely on it without careful tests (it depends on the
   * time_t definition of the implementation). 
   * On a MacOSX 10.7 platform 64 bits, it was possible to reach the 14th of december 1901
   * (index -24855, that is, 24855 days before the Epoch).
   */
  char cstr_day[10 + 1];
  char cstr_time_0[21 + 1];
  time_t time;
  int err;
  
  /* keeps only the substring "YYYY-MM-DD" */
  strncpy(cstr_day, cstr_time, sizeof(cstr_day) - 1);
  cstr_day[10] = '\0';

  snprintf(cstr_time_0, 21, "%sT00:00:00Z", cstr_day);
  err = convert_string_to_time(cstr_time_0, &time, TIMEFMT_ISO8601);
  if (err != 0) {
#ifdef DEBUG
    fprintf(stderr, "%s:%d: failed to convert \"%s\" to a day index\n", __FILE__, __LINE__, cstr_time);
#endif
    *index_of_day = 0;
    return false;
  }

  *index_of_day = time / 86400;
  return true;
}


char *dateDDMMAAAA2AAAAMMDD(const char *original_date){
    char *result=NULL;
    
    if(strlen(original_date)!=10){
        result=(char *)malloc(sizeof(char)*1);
        result[0]='\0';
        return result;
    }
    
    result=(char *)malloc(sizeof(char)*11);
    
    result[0]=original_date[6];
    result[1]=original_date[7];
    result[2]=original_date[8];
    result[3]=original_date[9];
    result[4]='-';
    result[5]=original_date[3];
    result[6]=original_date[4];
    result[7]='-';
    result[8]=original_date[0];
    result[9]=original_date[1];
    result[10]='\0';
    
    return result;
}


