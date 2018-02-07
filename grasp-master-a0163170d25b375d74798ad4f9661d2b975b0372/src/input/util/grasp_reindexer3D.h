/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

/**
 * @file  grasp_reindexer3D.h
 * @author Fabrice Ducos <fabrice.ducos@univ-lille1.fr>
 *
 */

#ifndef GRASP_REINDEXER3D_H_
#define GRASP_REINDEXER3D_H_

static const int GRASP_INVALID_INDEX = -1;

/* This module provides an interface to map linear indices to triplets of coordinates (it could
 * be easily generalized to other sets of coordinates).
 * 
 * The coordinates don't need to start at 0, or even to be positive. That way,
 * one can express indices directly in natural coordinates for the problem being dealt with.
 */

typedef struct grasp_reindexer3D_t_ grasp_reindexer3D_t;

typedef struct grasp_reindexer3D_props_t_ {
  int x0;
  int y0;
  int z0;
  size_t nx;
  size_t ny;
  size_t nz;
} grasp_reindexer3D_props_t;

extern grasp_reindexer3D_t *grasp_reindexer3D_new(const grasp_reindexer3D_props_t *props);
extern void grasp_reindexer3D_delete(grasp_reindexer3D_t *reindexer);
extern void grasp_reindexer3D_get_properties(const grasp_reindexer3D_t *reindexer, grasp_reindexer3D_props_t *props);
extern int  grasp_reindexer3D_get_linear_index(const grasp_reindexer3D_t *reindexer, int x, int y, int z);
extern bool grasp_reindexer3D_set_linear_index(grasp_reindexer3D_t *reindexer, int x, int y, int z, int linear_index);
extern size_t grasp_reindexer3D_get_allocated_memory(const grasp_reindexer3D_t *reindexer);

#endif /* GRASP_REINDEXER3D_H_ */
