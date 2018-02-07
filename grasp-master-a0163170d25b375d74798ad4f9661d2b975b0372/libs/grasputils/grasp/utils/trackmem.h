/**
 * @file:  
 * @Author: Fabrice Ducos <fabrice.ducos@univ-lille1.fr>
 *
 * Trackable allocation routines
 */

/*
Example of code:
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include "trackmem.h"

int main(void) {
  void *ptr;

  trackmem_set_debug_stream(stderr);
  trackmem_abort_on_failed_allocations(false); // change between true and false for testing 
 
  ptr = trackmem_malloc(10);
  trackmem_free(ptr);

  ptr = trackmem_malloc(-1);
  trackmem_free(ptr);

  return EXIT_SUCCESS;
}
*/

#ifndef TRACKMEM_H_
#define TRACKMEM_H_

#include <stdio.h>
#include <stdbool.h>

#define trackmem_malloc(size) trackmem_malloc_(__FILE__, __LINE__, #size, size)
#define trackmem_free(ptr) do { trackmem_free_(__FILE__, __LINE__, #ptr, ptr); ptr = NULL ; } while(0)

/* if a debug stream is defined (not NULL), a message will be displayed 
 * at each trackmem_malloc or trackmem_free call.
 * A usual debug_stream is stderr.
 *
 * Default value for debug_stream: NULL
 */
void trackmem_set_debug_stream(FILE *debug_stream);

/* if abort_enabled is true, trackmem_malloc will stop with an error message
 * instead of returning NULL.
 *
 * Default value for abort_enabled: false
 */
void trackmem_abort_on_failed_allocations(bool abort_enabled);

/* the following routines are not supposed to be called directly */
extern void *trackmem_malloc_(const char *file, size_t line, const char *arg_size, size_t size);
extern void trackmem_free_(const char *file, size_t line, const char *arg_ptr, void *ptr);

#endif /* TRACKMEM_H_ */
