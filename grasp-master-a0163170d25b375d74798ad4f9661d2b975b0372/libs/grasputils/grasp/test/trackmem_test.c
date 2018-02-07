/**
 * @file: trackmem_test.c 
 * @Author: Fabrice Ducos <fabrice.ducos@univ-lille1.fr>
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include "../utils/trackmem.h"

int main(void) {
  void *ptr;

  trackmem_set_debug_stream(stderr);
  trackmem_abort_on_failed_allocations(false); /* change between true and false for testing */
 
  ptr = trackmem_malloc(10);
  trackmem_free(ptr);

  ptr = trackmem_malloc(-1);
  trackmem_free(ptr);

  return EXIT_SUCCESS;
}
