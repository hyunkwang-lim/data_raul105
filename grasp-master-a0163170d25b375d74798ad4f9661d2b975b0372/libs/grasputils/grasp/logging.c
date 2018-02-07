/**
 * @file:  logging.c
 * @Author: Fabrice Ducos <fabrice.ducos@univ-lille1.fr>
 *
 */

#include <string.h>
#include <errno.h>
#include <assert.h>

#include "utils/logging.h"

static int debug_level = 0;
static FILE *debug_stream = NULL;
static FILE *err_stream   = NULL;

void set_error_stream(FILE *stream) {
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

void set_debug_stream(FILE *stream) {
  /* NULL means no debugging */
  debug_stream = stream;
  
  if (stream == NULL) {
    debug_level = 0;
  }
  else {
    if (debug_level == 0) {
      debug_level = 1;
    }
  }
}

FILE *get_error_stream(void) {
  if (err_stream == NULL) {
    return stderr;
  }
  else {
    return err_stream;
  }
}

FILE *get_debug_stream(void) {
  return debug_stream;
}

bool debug_is_enabled(void) {
  return (debug_level > 0);
}

/* implementation of perror that accepts a stream argument
 * (instead of sending automatically to stderr) 
 */
void perror2(FILE *stream, const char *s) {
  if (s != NULL && s[0] != '\0') {
    fprintf(stream, "%s: ", s);
  }

  fprintf(stream, "%s\n", strerror(errno));
}

void set_debug_level(int level) {
  if (level <= 0) {
    debug_level = 0;
    debug_stream = NULL;
  }
  else {
    debug_level = level;
    if (debug_stream == NULL) {
      debug_stream = stderr;
    }
  }
  
}

int get_debug_level(void) {
  return debug_level;
}
