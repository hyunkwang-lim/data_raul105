/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include <grasp/utils.h>
#include "sdata.h"
#include "../../util/grasp_box.h"

void test_dump_file(const char *file_name) {
  sdata_dump_file(stdout, file_name);
}

void test_get_header(const char *file_name) {
  SDATA_HANDLE *handle;
  SDATA_HEADER header;
  int err;
  (void)err; // This avoid -Wunused-but-set-variable warning
  
  handle = sdata_open(file_name);
  if (handle == NULL) {
    perror(file_name);
  }
  
  err = sdata_get_header(handle, &header);
  assert(err == 0);
  
  sdata_print_header(stdout, NULL, &header);
  
  err = sdata_close(handle);
  assert(err == 0);
}

void test_direct_access_to_pixel(const char *file_name)  {
  SDATA_HANDLE *handle;
  grasp_box_t *box;
  int err;
  bool bret; /* boolean return value */
  (void)err; // This avoid -Wunused-but-set-variable warning
  
  grasp_box_vector_t pixel_position = { 1377, 3286, 13884 }; /* for use with the file "Parasol.sdat2" provided in the test folder */
  SDATA_PIXEL pixel_data;
  
  assert(file_name != NULL);
  
  handle = sdata_open(file_name);
  
  if (handle == NULL) {
#ifdef DEBUG
    fprintf(stderr, "%s:%d: ", __FILE__, __LINE__);
    fprintf(stderr, "sdata_open(\"%s\") failed: ", file_name);
#endif
    perror(file_name);
    exit (EXIT_FAILURE);
  }
  
  box = sdata_load_box(handle, NULL /* settings inferred from the SDATA file itself */);
  assert(box != NULL);
  
  bret = grasp_box_get_data_from_pixel(box, &pixel_position, &pixel_data);
  if (bret == false) {
    fprintf(stderr, "%s: pixel [x: %d, y: %d, z: %d] not found (z: number of days since the Unix Epoch)\n", 
	    file_name, pixel_position.x, pixel_position.y, pixel_position.z);
    exit (EXIT_FAILURE);
  }
  sdata_print_pixel(stdout, NULL, &pixel_data);
  
  err = sdata_close(handle);
  assert(err == SDATA_OK);
  grasp_box_delete(box);
}

void test_iterate_on_pixels(const char *file_name) {
  SDATA_HANDLE *handle;
  grasp_box_t *box;
  int err;
  bool last_pixel_not_met;
  SDATA_PIXEL pixel_data;
  (void)err; // This avoid -Wunused-but-set-variable warning
  
  assert(file_name != NULL);
  
  handle = sdata_open(file_name);
  
  if (handle == NULL) {
#ifdef DEBUG
    fprintf(stderr, "%s:%d: ", __FILE__, __LINE__);
    fprintf(stderr, "sdata_open(\"%s\") failed: ", file_name);
    perror(file_name);
#endif
    exit (EXIT_FAILURE);
  }
  
  box = sdata_load_box(handle, NULL /* settings inferred from the SDATA file itself */);
  assert(box != NULL);

  printf("num_pixels: %ld\n", (unsigned long) grasp_box_get_num_records(box));

  do {
    last_pixel_not_met = grasp_box_get_next_pixel(box, &pixel_data);
    if (last_pixel_not_met) {
      sdata_print_pixel(stdout, NULL, &pixel_data);
      printf("\n");
    }
  } while (last_pixel_not_met);

  err = sdata_close(handle);
  assert(err == SDATA_OK);
  grasp_box_delete(box);

}

int main(int argc, char *argv[]) {
  const char *sdata_file;

  if (argc != 2 && argc != 3) {
    fprintf(stderr, "usage: %s <sdata_file> [-d]\n", argv[0]);
    exit (EXIT_FAILURE);
  }

  sdata_file = argv[1];

  if (argc == 3) {
    if (strcmp(argv[2], "-d") == 0) {
      sdata_set_debug_stream(stderr);
    }
    else {
      fprintf(stderr, "%s: invalid option: %s\n", argv[0], argv[2]);
      exit (EXIT_FAILURE);
    }
  }

//#define TEST_GET_HEADER
#ifdef TEST_GET_HEADER
  test_get_header(sdata_file);
#endif

#define TEST_DUMP_FILE
#ifdef TEST_DUMP_FILE
  test_dump_file(sdata_file);
#endif

  //#define TEST_DIRECT_ACCESS_TO_PIXEL  
#ifdef TEST_DIRECT_ACCESS_TO_PIXEL
  test_direct_access_to_pixel(sdata_file);
#endif
  
//#define TEST_ITERATE_ON_PIXELS
#ifdef TEST_ITERATE_ON_PIXELS
  test_iterate_on_pixels(sdata_file);
#endif

  return EXIT_SUCCESS;
}
