/**
 * @file: error.c 
 * @Author: Fabrice Ducos <fabrice.ducos@univ-lille1.fr>
 *
 */

#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <assert.h>
#include "utils/error.h"

#define ERROR_PREFIX_LEN 255
#define LAST_ERROR_LEN 1024

static char error_prefix[ERROR_PREFIX_LEN + 1] = "";
static char last_error_message[LAST_ERROR_LEN + 1] = "";

void set_error_prefix(const char *prefix) {
  strncpy(error_prefix, prefix, ERROR_PREFIX_LEN + 1);
}

void set_last_error_(const char *file, int line, int sys_errnum, const char *format, ...) {
  va_list arg;

  assert(format != NULL);

  va_start(arg, format);
  set_last_error_va_(file, line, sys_errnum, format, arg);
  va_end(arg);
}

void set_last_error_va_(const char *file, int line, int sys_errnum, const char *format, va_list arg) {
  size_t message_len;

  assert(format != NULL);

  message_len = 0;
  if (error_prefix[0] != '\0') {
    snprintf(last_error_message, LAST_ERROR_LEN, "%s: ", error_prefix);
    message_len = strlen(last_error_message);
  }
  
  vsnprintf(last_error_message + message_len, LAST_ERROR_LEN - message_len, format, arg);
  message_len = strlen(last_error_message);

  if (sys_errnum) {
    snprintf(last_error_message + message_len, LAST_ERROR_LEN - message_len, ": %s", strerror(sys_errnum));
    message_len = strlen(last_error_message);
  }
  last_error_message[LAST_ERROR_LEN] = '\0'; /* a paranoid statement to prevent possible buffer overflows */

  if (strlen(last_error_message) > LAST_ERROR_LEN - 4) {
    /* the error message is likely to be truncated, put an ellipsis to emphasize that fact */
    strncpy(&last_error_message[LAST_ERROR_LEN - 4], " ...", 4);
  }

#ifdef DEBUG
  if (error_prefix[0] != '\0') {
    fprintf(stderr, "%s:", error_prefix);
  }
  fprintf(stderr, "%s:%d: %s\n", file, line, get_last_error());
#endif
}

const char *get_last_error(void) {
  return last_error_message;
}

void print_last_error(void) {
  fprintf(stderr, "%s", last_error_message);
}

void println_last_error(void) {
  fprintf(stderr, "%s\n", last_error_message);
}

void reset_last_error(void) {
  last_error_message[0] = '\0';
}
