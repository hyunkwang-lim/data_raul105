/**
 * @file: allocator.c 
 * @Author: Fabrice Ducos <fabrice.ducos@univ-lille1.fr>
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <errno.h>
#include "utils/allocator.h"
#include "utils/trackmem.h"

#include "utils/debug.h"

struct allocator_t_ {
  size_t num_alloc_records; /* number of allocated records */
  size_t num_valid_records; /* number of records with real data */
  size_t record_size;
  void  *records;
  size_t iterator; /* to iterate on the valid records */
};

static void allocator_release(allocator_t *allocator) {
  assert(allocator != NULL);

  PRINT_DEBUG("%lu records (%lu bytes) released at %p\n",
              (long unsigned int) allocator->num_alloc_records,
              (long unsigned int) allocator->num_alloc_records * allocator->record_size,
              allocator->records);

  trackmem_free(allocator->records);
  allocator->records  = NULL;
  allocator->num_alloc_records = 0;
  allocator->num_valid_records = 0;
  allocator->iterator = 0;
}

void allocator_rewind(allocator_t *allocator) {
  assert(allocator != NULL);

  allocator->iterator = 0;
}

static void allocator_init(allocator_t *allocator, size_t record_size) {
  assert(allocator != NULL);
  assert(record_size > 0);
  
  allocator->num_alloc_records = 0;
  allocator->num_valid_records = 0;
  allocator->records = NULL;
  allocator->record_size = record_size;
  allocator->iterator = 0;
}

allocator_t *allocator_new(size_t record_size) {
  allocator_t *allocator;

  if (record_size == 0) {
    PRINT_DEBUG("attempt to allocate a 0-byte record with allocator_new");
    errno = EINVAL;
    return NULL;
  }

  allocator = trackmem_malloc(sizeof(allocator_t));
  if (allocator == NULL) {
    return NULL;
  }
  
  allocator_init(allocator, record_size);
  return allocator;
}

void allocator_delete(allocator_t *allocator) {
  allocator_release(allocator);
  trackmem_free(allocator);
}

void allocator_print(FILE *output_stream, const char *label, const allocator_t *allocator) {
  size_t alloc_size = allocator->num_alloc_records * allocator->record_size;

  fprintf(output_stream, "%s: { num_alloc_records: %lu num_valid_records: %lu iterator: %ld record_size: %lu memory: between %p and %p (%lu bytes) }\n", label, 
	  (long unsigned int) allocator->num_alloc_records, 
	  (long unsigned int) allocator->num_valid_records,
	  (long unsigned int) allocator->iterator,
	  (long unsigned int) allocator->record_size, allocator->records, allocator->records + alloc_size - 1, 
	  (long unsigned int) alloc_size);
}

bool allocator_increase_capacity(allocator_t *allocator, size_t n) {
  void *ptr;
  size_t num_alloc_records_new;
  
  assert(allocator != NULL);
  assert(n > 0);
  
  num_alloc_records_new = allocator->num_alloc_records + n;

  ptr = realloc(allocator->records, num_alloc_records_new * allocator->record_size);
  if (ptr == NULL) {
    PRINT_DEBUG("failed to realloc(allocator->records = %p, %lu*record_size = %lu bytes)\n",
		allocator->records, (long unsigned int) num_alloc_records_new, 
		(long unsigned int) num_alloc_records_new * allocator->record_size);
    
    return false;
  }

  allocator->records  = ptr;
  allocator->num_alloc_records = num_alloc_records_new;
  
  PRINT_DEBUG("%lu new records allocated, for a total of %lu records (%lu bytes) at %p\n",
	      (long unsigned int) n,
	      (long unsigned int) allocator->num_alloc_records, 
	      (long unsigned int) allocator->num_alloc_records * allocator->record_size, 
	      allocator->records);

  return true;
}

size_t allocator_get_num_alloc_records(const allocator_t *allocator) {
  assert(allocator != NULL);
  
  return allocator->num_alloc_records;
  
}

size_t allocator_get_num_valid_records(const allocator_t *allocator) {
  assert(allocator != NULL);
  
  return allocator->num_valid_records;
  
}

size_t allocator_get_allocated_memory(const allocator_t *allocator) {
  return allocator->num_alloc_records * allocator->record_size;
}

size_t allocator_get_record_size(const allocator_t *allocator) {
  assert(allocator != NULL);

  return allocator->record_size;
}

static void *get_record_offset(const allocator_t *allocator, size_t irecord) {
  unsigned char *offset;

  assert(allocator != NULL);
  
  assert(irecord < allocator->num_alloc_records);
  offset = allocator->records;
  offset += irecord * allocator->record_size;

  return offset;
}

bool allocator_get_record(const allocator_t *allocator, size_t irecord, void *record_copy) {
  void *offset;
  
  assert(allocator != NULL);
  
  offset = get_record_offset(allocator, irecord);
  assert(offset != NULL);

  if (record_copy != NULL) {
    memcpy(record_copy, offset, allocator->record_size);
  }

  return true;
}

bool allocator_get_next_record(allocator_t *allocator, void *record_copy) {
  bool retval;
  
  if (allocator->iterator == allocator->num_valid_records) {

#ifdef DEBUG
    fprintf(stderr, "%s:%d: ", __FILE__, __LINE__);
    fprintf(stderr, "%s: iterator's end for the allocator [%p] has been reached (num_valid_records: %ld)", __FUNCTION__, allocator,
	    allocator->num_valid_records);
#endif
    return false;
  }
  retval = allocator_get_record(allocator, allocator->iterator, record_copy);
  assert(retval == true);
  allocator->iterator++;
  return retval;
}

bool allocator_set_record(allocator_t *allocator, size_t irecord, const void *record_copy) {
  void *offset;
  
  /* this routine allows to overwrite old records or to set new ones (in the pool
   * of allocated records).
   */

  assert(allocator != NULL);
  assert(record_copy != NULL);
  
  if (! (irecord < allocator->num_alloc_records)) {
    fprintf(stderr, "%s:%d: ", __FILE__, __LINE__);
    fprintf(stderr, "%s: bad value for the argument irecord: %zu (was expected in the range [0..%zu]). Please fix the caller.\n",
	    __FILE__, irecord, allocator->num_alloc_records - 1);
    abort();
  }

  offset = get_record_offset(allocator, irecord);
  assert(offset != NULL);
  
  if (irecord == allocator->num_valid_records) {
    allocator->num_valid_records++;
  }

  /* TODO: this should be replaced by a call to allocator_increase_capacity() */
  assert(allocator->num_valid_records <= allocator->num_alloc_records);

  memcpy(offset, record_copy, allocator->record_size);
  
  return true;
}


#ifndef UNIT_TEST
#define UNIT_TEST 0
#endif

#if (UNIT_TEST != 0)

static bool test(void) {
  allocator_t *allocator;
  typedef char record_t[10];
  FILE *log_stream = stderr;

  fprintf(log_stream, "allocator_new(0) ");
  allocator = allocator_new(0);
  if (allocator != NULL) {
    allocator_print(log_stream, "passed, but shouldn't have", allocator);
    return false;
  }
  else {
    perror("failed (the failure was expected, the test passed)");
  }

  fprintf(log_stream, "allocator_new(sizeof(record_t)) ");
  allocator = allocator_new(sizeof(record_t));
  if (allocator != NULL) {
    allocator_print(log_stream, "passed", allocator);
  }
  else {
    perror("failed");
    return false;
  }
  
  fprintf(log_stream, "allocator_increase_capacity(allocator, 10) ");
  if (allocator_increase_capacity(allocator, 10)) {
    allocator_print(log_stream, "passed", allocator);
  }
  else {
    perror("failed");
    return false;
  }
  
  fprintf(log_stream, "allocator_increase_capacity(allocator, 10) ");
  if (allocator_increase_capacity(allocator, 10)) {
    allocator_print(log_stream, "passed", allocator);
  }
  else {
    perror("failed");
    return false;
  }
  
  allocator_delete(allocator);
  allocator_print(log_stream, "state after allocator_delete(allocator): ", allocator);
  
  return true;
}

int main(int argc, char *argv[]) {
  return (test() ? EXIT_SUCCESS : EXIT_FAILURE);
}
#endif
