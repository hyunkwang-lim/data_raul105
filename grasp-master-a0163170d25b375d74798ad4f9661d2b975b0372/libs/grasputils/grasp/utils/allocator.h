/**
 * @file: allocator.h 
 * @Author: Fabrice Ducos <fabrice.ducos@univ-lille1.fr>
 *
 */

#ifndef ALLOCATOR_H_
#define ALLOCATOR_H_

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

typedef struct allocator_t_ allocator_t;

extern allocator_t *allocator_new(size_t record_size);
extern void allocator_rewind(allocator_t *allocator);
extern  void allocator_delete(allocator_t *allocator);
extern void allocator_print(FILE *output_stream, const char *label, const allocator_t *allocator);
extern bool allocator_increase_capacity(allocator_t *allocator, size_t n /* number of records to add */);
extern size_t allocator_get_num_alloc_records(const allocator_t *allocator);
extern size_t allocator_get_num_valid_records(const allocator_t *allocator);
extern size_t allocator_get_record_size(const allocator_t *allocator);
extern size_t allocator_get_allocated_memory(const allocator_t *allocator);
extern bool allocator_get_record(const allocator_t *allocator, size_t irecord, void *record_copy);
extern bool allocator_get_next_record(allocator_t *allocator, void *record_copy);
extern bool allocator_set_record(allocator_t *allocator, size_t irecord, const void *record_copy);

#endif /* ALLOCATOR_H_ */
