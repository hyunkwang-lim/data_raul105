/**
 * @file:  trackmem.c
 * @Author: Fabrice Ducos <fabrice.ducos@univ-lille1.fr>
 *
 */

#include <stdlib.h>
#include <stdio.h>

#include "utils/trackmem.h"

static FILE *g_debug_stream_ = NULL;
static bool g_abort_enabled_ = false;

void trackmem_set_debug_stream(FILE *debug_stream) {
  g_debug_stream_ = debug_stream;
}

void trackmem_abort_on_failed_allocations(bool abort_enabled) {
  g_abort_enabled_ = abort_enabled;
}

void *trackmem_malloc_(const char *file, size_t line, const char *arg_size, size_t size) {
  void *ptr = malloc(size);
  if (g_debug_stream_ != NULL) {
    if (ptr != NULL) {
      fprintf(g_debug_stream_, "%s:%zu: malloc(%s [%zu]) --> %p\n", file, line, arg_size, size, ptr);
    }
    else {
      fprintf(g_debug_stream_, "%s:%zu: malloc(%s [%zu]) --> (nil)\n", file, line, arg_size, size);
    }
  }

  if (ptr == NULL) {
    if (g_abort_enabled_ == true) {
      fprintf(stderr, "%s:%zu: malloc(%s [%zu]) has failed. Failed allocations have been declared fatal. Aborting.\n", file, line, arg_size, size);
      abort();
    }
  }
  
  return ptr;
}

void trackmem_free_(const char *file, size_t line, const char *arg_ptr, void *ptr) {
  if (g_debug_stream_ != NULL) {
    if (ptr != NULL) {
      fprintf(g_debug_stream_, "%s:%zu: free(%s [%p])\n", file, line, arg_ptr, ptr);
    }
    else {
      fprintf(g_debug_stream_, "%s:%zu: free(%s [(nil)])\n", file, line, arg_ptr);
    }
  }
  free(ptr);
}

