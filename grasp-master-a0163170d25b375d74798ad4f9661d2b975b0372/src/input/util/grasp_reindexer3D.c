/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

/**
 * @file:  grasp_reindexer3D.c
 * @author: Fabrice Ducos <fabrice.ducos@univ-lille1.fr>
 *
 */

#undef DEBUG /* remove this for local debugging */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* memset */
#include <errno.h>  /* errno */
#include <assert.h>
#include <stdbool.h>
#include <grasp/utils.h>
#include "grasp_reindexer3D.h"

typedef int index_type;

struct grasp_reindexer3D_t_ {
  grasp_reindexer3D_props_t props;
  int ***linear_indices;
  size_t allocated_memory;
};

static int ***alloc3d(const size_t n1, const size_t n2, const size_t n3) {
  int *buffer;
  int **ptr2;
  int ***ptr1;
  int i;

  buffer = trackmem_malloc(n1 * n2 * n3 * sizeof(int));
  assert(buffer!=NULL);
  if (buffer == NULL) {
    errno = ENOMEM;
    return NULL;
  }

  ptr2 = trackmem_malloc(n1 * n2 * sizeof(int *));
  assert(ptr2!=NULL);
  if (ptr2 == NULL) {
    trackmem_free(buffer);
    errno = ENOMEM;
    return NULL;
  }

  ptr1 = trackmem_malloc(n1 * sizeof(int **));
  assert(ptr1!=NULL);
  if (ptr1 == NULL) {
    trackmem_free(buffer);
    trackmem_free(ptr2);
    errno = ENOMEM;
    return NULL;
  }

  for (i = 0 ; i < n1*n2 ; i++) {
    if (i < n1) { ptr1[i] = &ptr2[i*n2]; }
    ptr2[i] = &buffer[i*n3];
  }

  return ptr1;
}

static void free3d(int ***ptr1) {
  int *buffer;
  int **ptr2;

  if (ptr1 == NULL) return;

  ptr2 = ptr1[0];
  assert(ptr2 != NULL);
  buffer = ptr2[0];

  trackmem_free(buffer);
  trackmem_free(ptr2);
  trackmem_free(ptr1);
}

static bool reindexer3D_init(grasp_reindexer3D_t *reindexer, const grasp_reindexer3D_props_t *props) {
  int x, y, z;
  size_t nx, ny, nz;

  assert(reindexer != NULL);
  assert(props != NULL);

  nx = props->nx;
  ny = props->ny;
  nz = props->nz;

  memcpy(&reindexer->props, props, sizeof(*props));
  
  reindexer->linear_indices = alloc3d(nx, ny, nz);
  reindexer->allocated_memory = (nx*ny*nz + nx*ny + nx)*sizeof(int);

  if (reindexer->linear_indices != NULL) {
    for (x = 0 ; x < nx ; x++) {
      for (y = 0 ; y < ny ; y++) {
	for (z = 0 ; z < nz ; z++) {
	  reindexer->linear_indices[x][y][z] = GRASP_INVALID_INDEX;
	}
      }
    }
  }

  PRINT_DEBUG("leaving grasp_reindexer3D_init\n");
  return (reindexer->linear_indices != NULL);
}

static void reindexer3D_release(grasp_reindexer3D_t *reindexer) {
  assert(reindexer != NULL);
  free3d(reindexer->linear_indices);

  memset(reindexer, 0, sizeof(grasp_reindexer3D_t));
  PRINT_DEBUG("leaving grasp_reindexer3D_release\n");
}

grasp_reindexer3D_t *grasp_reindexer3D_new(const grasp_reindexer3D_props_t *props) {
  grasp_reindexer3D_t *reindexer;

  reindexer = trackmem_malloc(sizeof(grasp_reindexer3D_t));
  assert(reindexer!=NULL);
  if (reindexer == NULL) {
    PRINT_DEBUG("trackmem_malloc failed\n");
    return NULL;
  }

  if (! reindexer3D_init(reindexer, props)) {
    PRINT_DEBUG("grasp_reindexer3D_init failed\n");
  }

  PRINT_DEBUG("leaving grasp_reindexer3D_new\n");
  return reindexer;
}

void grasp_reindexer3D_delete(grasp_reindexer3D_t *reindexer) {
  reindexer3D_release(reindexer);
  trackmem_free(reindexer);
  PRINT_DEBUG("leaving grasp_reindexer3D_delete\n");
}

int grasp_reindexer3D_get_linear_index(const grasp_reindexer3D_t *reindexer, int x, int y, int z) {
  assert(reindexer != NULL);
  assert(reindexer->linear_indices != NULL);
  int x_ = x - reindexer->props.x0;
  int y_ = y - reindexer->props.y0;
  int z_ = z - reindexer->props.z0;

  reset_last_error();
  
  if (! (0 <= x_ && x_ < reindexer->props.nx)) {
    SET_LAST_ERROR(0, "%s: 1st arg (%d) out of bounds (should be in the range [%d..%d])", 
		   __FUNCTION__, x, reindexer->props.x0, reindexer->props.x0 + reindexer->props.nx - 1);
    return GRASP_INVALID_INDEX;
  }
  if (! (0 <= y_ && y_ < reindexer->props.ny)) {
    SET_LAST_ERROR(0, "%s: 2nd arg (%d) out of bounds (should be in the range [%d..%d])", 
		   __FUNCTION__, y, reindexer->props.y0, reindexer->props.y0 + reindexer->props.ny - 1);
    return GRASP_INVALID_INDEX;
  }
  if (! (0 <= z_ && z_ < reindexer->props.nz)) {
    SET_LAST_ERROR(0, "%s: 3rd arg (%d) out of bounds (should be in the range [%d..%d])", 
		   __FUNCTION__, z, reindexer->props.z0, reindexer->props.z0 + reindexer->props.nz - 1);
    return GRASP_INVALID_INDEX;
  }
  
  PRINT_DEBUG("leaving grasp_reindexer3D_get_linear_index(reindexer: %p, x: %d, y: %d, z: %d)\n",
	      reindexer, x, y, z);
  return reindexer->linear_indices[x_][y_][z_];
}

bool grasp_reindexer3D_set_linear_index(grasp_reindexer3D_t *reindexer, int x, int y, int z, int linear_index) {
  assert(reindexer != NULL);
  assert(reindexer->linear_indices != NULL);
  int x_ = x - reindexer->props.x0;
  int y_ = y - reindexer->props.y0;
  int z_ = z - reindexer->props.z0;
  
  reset_last_error();

  if (! (0 <= x_ && x_ < reindexer->props.nx)) {
    SET_LAST_ERROR(0, "%s: 1st arg (%d) out of bounds (should be in the range [%d..%d])", 
		   __FUNCTION__, x, reindexer->props.x0, reindexer->props.x0 + reindexer->props.nx - 1);
    return false;
  }
  if (! (0 <= y_ && y_ < reindexer->props.ny)) {
    SET_LAST_ERROR(0, "%s: 2nd arg (%d) out of bounds (should be in the range [%d..%d])", 
		   __FUNCTION__, y, reindexer->props.y0, reindexer->props.y0 + reindexer->props.ny - 1);
    return false;
  }
  if (! (0 <= z_ && z_ < reindexer->props.nz)) {
    SET_LAST_ERROR(0, "%s: 3rd arg (%d) out of bounds (should be in the range [%d..%d])", 
		   __FUNCTION__, z, reindexer->props.z0, reindexer->props.z0 + reindexer->props.nz - 1);
    return false;
  }
  
  reindexer->linear_indices[x_][y_][z_] = linear_index;
  PRINT_DEBUG("leaving grasp_reindexer3D_set_linear_index(reindexer: %p, x: %d, y: %d, z: %d, linear_index: %d)\n",
	      reindexer, x, y, z, linear_index);
  return true;
}

void grasp_reindexer3D_get_properties(const grasp_reindexer3D_t *reindexer, grasp_reindexer3D_props_t *props) {
  memcpy(props, &reindexer->props, sizeof(*props));
}

size_t grasp_reindexer3D_get_allocated_memory(const grasp_reindexer3D_t *reindexer) {
  return reindexer->allocated_memory;
}

#ifndef UNIT_TEST
#define UNIT_TEST 0
#endif

#if (UNIT_TEST != 0)

static bool test(void) {
  grasp_reindexer3D_t *reindexer;
  grasp_reindexer3D_props_t props = {
    1, 1, 1, /* x0, y0, z0 */
    4, 4, 1  /* nx, ny, nz */
  };
  int x, y, z;
  int i;

  double values_in[] = {
     0,  1,  2,  3, 
     4,  5,  6,  7, 
     8,  9, 10, 11, 
     12, 13, 14, 15
  };
  
  set_error_prefix("error");
  reindexer = grasp_reindexer3D_new(&props);
  assert(reindexer != NULL);
  
  /* transposition */
  i = 0;
  for (x = props.x0 ; x < props.x0 + props.nx ; x++) {
    for (y = props.y0 ; y < props.y0 + props.ny ; y++) {
      for (z = props.z0 ; z < props.z0 + props.nz; z++) {
	if (i < sizeof(values_in)/sizeof(*values_in)) {
	  bool success = grasp_reindexer3D_set_linear_index(reindexer, y, x, z, i);
	  if (! success) {
	    println_last_error();
	  }
	}
	else {
	  grasp_reindexer3D_set_linear_index(reindexer, y, x, z, GRASP_INVALID_INDEX);
	}
	i++;
      }
    }
  }
  
  for (x = props.x0 ; x < props.x0 + props.nx ; x++) {
    for (y = props.y0 ; y < props.y0 + props.ny ; y++) {
      for (z = props.z0 ; z < props.z0 + props.nz ; z++) {
	int i = reindexer3D_get_linear_index(reindexer, x, y, z);
	if (i >= 0) {
	  printf("x = %d y = %d z = %d value = %f\n", x, y ,z, values_in[i]);
	}
	else {
	  printf("x = %d y = %d z = %d value = NONE\n", x, y ,z);
	}
      }
    }
  }

  grasp_reindexer3D_delete(reindexer);

  return true;
}


int main(int argc, char *argv[]) {
  return (test() ? EXIT_SUCCESS : EXIT_FAILURE);
}
#endif
