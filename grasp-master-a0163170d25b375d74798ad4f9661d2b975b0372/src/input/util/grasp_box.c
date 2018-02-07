/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

/**
 * @file:  grasp_box.c
 * @author: Fabrice Ducos <fabrice.ducos@univ-lille1.fr>
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <glob.h>
#include <errno.h>  /* errno */
#include <string.h> /* memset */
#include <stdint.h> /* C99 only */
#include <stdbool.h> /* C99 only */
#include <assert.h>

#include <grasp/utils.h>
#include "grasp_reindexer3D.h"
#include "grasp_box.h"

static int debug_level = 0;

struct grasp_box_t_ {
  allocator_t *allocator;
  grasp_reindexer3D_t *reindexer;
  size_t num_records;
  
  size_t num_parent_files; /* number of elements of parent_files (optional) */
  size_t num_parent_files_alloc; /* number of allocated strings for files (num_parent_files <= num_parent_files_alloc) */
  filepath_t *parent_files; /* list of parent files used for filling the box (optional) */

  grasp_reindexer3D_props_t props; /* properties (origin and dimensions) of the box */
};

void grasp_box_set_debug_level(int level) {
  debug_level = level;
}

size_t grasp_box_get_num_parent_files(const grasp_box_t *box) {
  assert(box != NULL);
  return box->num_parent_files;
}

void grasp_box_get_parent_files(const grasp_box_t *box, size_t num_parent_files, filepath_t *parent_files /* must be allocated by the caller */) {
  size_t ifile;
  assert(box != NULL);

  assert(num_parent_files == box->num_parent_files);
  for (ifile = 0 ; ifile < num_parent_files ; ifile++) {
    strncpy(parent_files[ifile], box->parent_files[ifile], sizeof(filepath_t) - 1);
  }
}

void grasp_box_add_parent_file(grasp_box_t *box, const char *parent_file) {
  size_t nfiles;

  assert(box != NULL);
  
  nfiles = box->num_parent_files;

  assert(nfiles <= box->num_parent_files_alloc);
  if (nfiles == box->num_parent_files_alloc) {
    box->num_parent_files_alloc += 128;
    filepath_t *ptr = realloc(box->parent_files, (box->num_parent_files_alloc) * sizeof(filepath_t));
    if (ptr != NULL) {
      box->parent_files = ptr;
    }
  }
  assert(box->parent_files != NULL);
  memset(box->parent_files[nfiles], 0, sizeof(filepath_t));
  strncpy(box->parent_files[nfiles], parent_file, sizeof(filepath_t) - 1);
  box->num_parent_files++;
}

void grasp_box_rewind(grasp_box_t *box) {

  assert(box != NULL);
  allocator_rewind(box->allocator);
}

grasp_box_t *grasp_box_new(size_t record_size, const grasp_box_settings_t *settings) {
  grasp_reindexer3D_props_t props;
  grasp_box_t *box;

  assert(settings != NULL);
 
#ifdef DEBUG
  if (debug_level == 0) { debug_level = 1; }
#endif
 
  props.x0 = settings->orig.x;
  props.y0 = settings->orig.y;
  props.z0 = settings->orig.z;
  props.nx = settings->extent.x;
  props.ny = settings->extent.y;
  props.nz = settings->extent.z;
  
  box = trackmem_malloc(sizeof(grasp_box_t));
  assert(box!=NULL);
  assert(box != NULL);
  
  box->num_records = 0;
  box->num_parent_files = 0;
  box->num_parent_files_alloc = 0;
  box->props = props;
  box->parent_files = NULL;
  box->allocator = allocator_new(record_size); /* will manage the memory allocation */
  assert(box->allocator != NULL);
  box->reindexer = grasp_reindexer3D_new(&props); /* will manage the memory indexation */
  assert(box->reindexer != NULL);
  grasp_box_rewind(box);
  return box;
}

void grasp_box_print(FILE *output_stream, const char *label, const grasp_box_t *box) {
  assert(output_stream != NULL);
  assert(box != NULL);

  size_t alloc_mem = allocator_get_allocated_memory(box->allocator);
  size_t reind_mem = grasp_reindexer3D_get_allocated_memory(box->reindexer);
  size_t total_mem = alloc_mem + reind_mem;

  if (label != NULL && label[0] != '\0') {
    fprintf(output_stream, "%s: ", label);
  }
  fprintf(output_stream, "num_records: %lu ", (long unsigned int) box->num_records);
  allocator_print(output_stream, "allocator", box->allocator);
  fprintf(output_stream, "memory: allocator[%lu bytes] + reindexer3D[%lu bytes] = total[%lu bytes]\n",
	  (unsigned long) alloc_mem, (unsigned long) reind_mem, (unsigned long) total_mem);
}

void grasp_box_delete(grasp_box_t *box) {
  assert(box != NULL);
  
  assert(box->num_records <= allocator_get_num_alloc_records(box->allocator));
  allocator_delete(box->allocator);
  grasp_reindexer3D_delete(box->reindexer);
  trackmem_free(box->parent_files);
  box->parent_files = NULL;
  box->num_parent_files = 0;
  box->num_parent_files_alloc = 0;
  box->num_records = 0;
  trackmem_free(box);
}

bool grasp_box_attach_data_to_pixel(grasp_box_t *box, const grasp_box_vector_t *pixel, const size_t data_size, const void *data) {
  size_t num_records_per_alloc = 128; /* number of records allocated at each allocation call */
  size_t irecord;
  size_t record_size;
  int x0, y0, z0; /* origin of the box */
  int nx, ny, nz; /* dimensions of the box */
  int x, y, z; /* coordinates of the pixel to attach */
  bool success; 
  (void)record_size; // This avoid -Wunused-but-set-variable warning
  
  assert(box != NULL);
  assert(pixel != NULL);
  
  /* may allow NULL data later */
  assert(data != NULL);

  x0 = box->props.x0;
  y0 = box->props.y0;
  z0 = box->props.z0;
  nx = box->props.nx;
  ny = box->props.ny;
  nz = box->props.nz;
  x  = pixel->x;
  y  = pixel->y;
  z  = pixel->z;

  if (! ((x0 <= x && x < x0 + nx) && (y0 <= y && y < y0 + ny) && (z0 <= z && z < z0 + nz)) ) {
    if (debug_level > 0) {
      fprintf(stderr, "%s:%d: %s: attempt to attach a pixel (%d,%d,%d), outside of the boundaries of the box [%d..%d), [%d..%d), [%d..%d). Failed.\n",
	      __FILE__, __LINE__, __func__, x, y, z, x0, x0 + (int) nx, y0, y0 + (int) ny, z0, z0 + (int) nz);
    }
    return false;
  } 

  irecord = box->num_records;
  record_size = allocator_get_record_size(box->allocator);
  
  assert(data_size == record_size);
  
  if (irecord % num_records_per_alloc == 0) {
    allocator_increase_capacity(box->allocator, num_records_per_alloc);
  }
  
  allocator_set_record(box->allocator, irecord, data);
  success = grasp_reindexer3D_set_linear_index(box->reindexer, x, y, z, irecord);
  if (success == false) {
    if (debug_level > 0) {
      println_last_error();
    }
    return false;
  }

  box->num_records++;
  
  assert(box->num_records == allocator_get_num_valid_records(box->allocator));

  return true;
}

bool grasp_box_get_data_from_pixel(const grasp_box_t *box, const grasp_box_vector_t *pixel, void *data_copy) {
  int irecord;

  assert(box != NULL);
  assert(pixel != NULL);
  assert(data_copy != NULL);

  irecord = grasp_reindexer3D_get_linear_index(box->reindexer, pixel->x, pixel->y, pixel->z);
  if (irecord == GRASP_INVALID_INDEX) {
    return false;
  }
  allocator_get_record(box->allocator, irecord, data_copy);
  return true;
}

bool grasp_box_set_data_on_pixel(grasp_box_t *box, const grasp_box_vector_t *pixel, const void *data_copy) {
  int irecord;

  assert(box != NULL);
  assert(pixel != NULL);
  assert(data_copy != NULL);

  irecord = grasp_reindexer3D_get_linear_index(box->reindexer, pixel->x, pixel->y, pixel->z);
  if (irecord == GRASP_INVALID_INDEX) {
    return false;
  }

  allocator_set_record(box->allocator, irecord, data_copy);
  return true;
}

bool grasp_box_get_next_pixel(grasp_box_t *box, void *data_copy) {
  bool retval;
  
  assert(box != NULL);
  assert(data_copy != NULL);
  retval = allocator_get_next_record(box->allocator, data_copy);
  if (retval == false) {

#ifdef DEBUG
    fprintf(stderr, "%s:%d: ", __FILE__, __LINE__);
    fprintf(stderr, "%s: the end of the box [%p] has been reached (no more pixels)\n", __FUNCTION__, box);
#endif
  }
  return retval;
}

size_t grasp_box_get_num_records(const grasp_box_t *box) {
  assert(box != NULL);
  return box->num_records;
}

void grasp_box_get_settings(const grasp_box_t *box, grasp_box_settings_t *settings_copy) {

  assert(box != NULL);
  assert(settings_copy != NULL);
  
  settings_copy->orig.x = box->props.x0;
  settings_copy->orig.y = box->props.y0;
  settings_copy->orig.z = box->props.z0;
  settings_copy->extent.x = box->props.nx;
  settings_copy->extent.y = box->props.ny;
  settings_copy->extent.z = box->props.nz;
  
}

/*************************************************************************************************/
#ifndef UNIT_TEST
#define UNIT_TEST 0
#endif

#if (UNIT_TEST != 0)

static bool test(void) {
  grasp_box_t *box;
  const size_t NX = 2, NY = 2, NZ = 2;
  int ix, iy, iz;
  grasp_box_vector_t orig = { 0, 0, 0 };
  grasp_box_vector_t extent = { NX, NY, NZ };
  grasp_box_settings_t settings;
  int values_in[3];
  int count = 0;
  
  settings.orig = orig;
  settings.extent = extent;

  box = grasp_box_new(sizeof(values_in), &settings);
  
  for (ix = 0 ; ix < NX ; ix++) {
    for (iy = 0 ; iy < NY ; iy++) {
      for (iz = 0 ; iz < NZ ; iz++) {
	grasp_box_vector_t pixel;
	bool ret;
	
	pixel.x = ix;
	pixel.y = iy;
	pixel.z = iz;
	values_in[0] = count;
	values_in[1] = 2*count;
	values_in[2] = 3*count;
	ret = grasp_box_attach_data_to_pixel(box, &pixel, sizeof(values_in), values_in);
	if (ret == false) {
	  fprintf(stderr, "grasp_box: unit test failed while trying to attach data to the pixel (%d %d %d)\n",
		  ix, iy, iz);
	  goto test_failed;
	}
	count++;
      }
    }
  }

  count = 0;
  for (ix = 0 ; ix < NX ; ix++) {
    for (iy = 0 ; iy < NY ; iy++) {
      for (iz = 0 ; iz < NZ ; iz++) {
	grasp_box_vector_t pixel;
	int values_out[3];
	bool ret;

	pixel.x = ix;
	pixel.y = iy;
	pixel.z = iz;
	
	ret = grasp_box_get_data_from_pixel(box, &pixel, values_out);
	if (ret == false) {
	  fprintf(stderr, "grasp_box: unit test failed while trying to retrieve data from the pixel (%d %d %d)\n",
		  ix, iy, iz);
	  goto test_failed;
	}

	if ((values_out[0] != count) || (values_out[1] != 2*count) || (values_out[2] != 3*count)) {
	  fprintf(stderr, "grasp_box: unit test failed at (%d %d %d): retrieved values: %d %d %d (expected values: %d %d %d)\n", 
		 ix, iy, iz, values_out[0], values_out[1], values_out[2], count, 2*count, 3*count);
	  goto test_failed;
	}
	count++;
      }
    }
  }

  fprintf(stderr, "grasp_box: unit test passed\n");
  grasp_box_delete(box);
  return true;

 test_failed:
  grasp_box_delete(box);
  return false;
}

int main(int argc, char *argv[]) {
  
  return (test() ? EXIT_SUCCESS : EXIT_FAILURE);
}
#endif
