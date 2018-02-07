/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

/**
 * @file sdata_dump.c 
 * @author Fabrice Ducos <fabrice.ducos@univ-lille1.fr>
 *
 */


/* #define _GNU_SOURCE */ /* modify basename's behavior and enables getopt_long */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h> /* getopt */
#include <libgen.h> /* basename, dirname */
#include <errno.h>  /* errno */
#include <string.h> /* memset */
#include <stdint.h> /* C99 only */
#include <stdbool.h> /* C99 only */
#include <assert.h>

#ifdef _GNU_SOURCE
#include <getopt.h> /* getopt_long, getopt_long_only */
#endif

#include "sdata.h"

#define CSTRING_MAXLEN 255 /* MAXLEN, not MAXSIZE, so the terminal NULL character is not to be counted;
			    * for instance, a cstring should be allocated with:
			    * char cstring[CSTRING_MAXLEN + 1];
			    * for safety, functions like snprintf should be used like this:
			    * ret = snprintf(cstring, sizeof(cstring), format, ...)
			    */

char *g_appname;
char *g_appdir;
char g_version[CSTRING_MAXLEN + 1] = "0.1.0";
char g_output_file[CSTRING_MAXLEN + 1] = "";
uint32_t g_flags = 0; /* C99 only */
int g_verbosity_level = 0;

/* powers of two only */
enum {
  OPT_DEBUG = 1,
  OPT_UNFOLD = 2,
  OPT_HEADER = 4
};


void usage() {
  fprintf(stderr, "%s %s\n\nusage: %s [OPTIONS] <sdata_file>\n", g_appname, g_version, g_appname);
  fprintf(stderr, 
	  "OPTIONS:\n"
	  "  -D                  debug mode\n"
	  "  -H                  print the SDATA header\n"
	  "  -R                  displays sdata in raw, but normalized format (default)\n"
	  "  -U                  unfold sdata\n");
  exit (EXIT_FAILURE);
}

	
void parse_options(int *argc, char **argv[]) {
  int option;
  const char *optstring = "DHRU";

#ifdef _GNU_SOURCE
  int *longindex = NULL;

  struct option longopts[] = {
    /* name       has_arg                      flag,        val  */
    {  "debug",   no_argument,                 NULL,        'D' },
    {  "header",  no_argument,                 NULL,        'H' },
    {  "raw",     no_argument,                 NULL,        'R' },
    {  "unfold",  required_argument,           NULL,        'U' },
    {  NULL,      no_argument,                 NULL,         0  }
  };

  while ((option = getopt_long(*argc, *argv, optstring, longopts, longindex)) != EOF)
#else
  while ((option = getopt(*argc, *argv, optstring)) != EOF )
#endif
    {
      switch(option) {
      case 'D' :
	g_flags |= OPT_DEBUG;
	break;
      case 'H' :
	g_flags |= OPT_HEADER;
	break;
      case 'R' :
	g_flags &= ~OPT_UNFOLD;
	break;
      case 'U' :
	g_flags |= OPT_UNFOLD;
	break;
      case ':' :
	fprintf(stderr, "%s: option %c: an argument was expected", g_appname, optopt);
	exit (EXIT_FAILURE);
	break;
      case '?' :
	fprintf(stderr, "%s: option %c: unknown\n\n", g_appname, optopt);
	exit (EXIT_FAILURE);
	break;
      default : /* should never arrive here in a stable version */
	fprintf(stderr,"%s: parse_options: unexpected value for option: %d\n", g_appname, option);
	abort();
	break;
      } /* switch */    
    } /* while */
  
  *argc -= optind - 1;
  *argv += optind - 1;
}

void unfold_sdata(const char *file_name) {
  SDATA_HANDLE *handle;
  grasp_box_t *box;
  int err;
  bool last_pixel_not_met;
  SDATA_PIXEL pixel_data;

  assert(file_name != NULL);
  
  handle = sdata_open(file_name);
  
  if (handle == NULL) {
    if (g_flags & OPT_DEBUG) {
      fprintf(stderr, "%s:%d: ", __FILE__, __LINE__);
      fprintf(stderr, "sdata_open(\"%s\") failed: ", file_name);
    }
    perror(file_name);
    exit (EXIT_FAILURE);
  }
 
  box = sdata_load_box(handle, NULL /* settings inferred from the SDATA file itself */);
  if (box == NULL) {
    if (g_flags & OPT_DEBUG) {
      fprintf(stderr, "%s:%d: ", __FILE__, __LINE__);
      fprintf(stderr, "sdata_load_box(handle [\"%s\"], NULL) failed: ", file_name);
      abort();
    }
    fprintf(stderr, "failed to load a box from %s\n", file_name);
    exit (EXIT_FAILURE);
  }

#ifdef DEBUG
  printf("num_pixels: %ld\n", (unsigned long) grasp_box_get_num_records(box));
#endif

  do {
    last_pixel_not_met = grasp_box_get_next_pixel(box, &pixel_data);
    if (last_pixel_not_met) {
      sdata_print_pixel(stdout, NULL, &pixel_data);
      printf("\n");
    }
  } while (last_pixel_not_met);

  err = sdata_close(handle);
  if (err != SDATA_OK) {
    if (g_flags & OPT_DEBUG) {
      fprintf(stderr, "%s:%d: sdata_close(handle [\"%s\"]) failed\n", __FILE__, __LINE__, file_name);
    }
    
    perror(file_name);
    exit (EXIT_FAILURE);
  }
  grasp_box_delete(box);

}

void print_header(const char *file_name) {
  SDATA_HEADER header;
  SDATA_HANDLE *handle;
  int err;
  
  handle = sdata_open(file_name);
  
  if (handle == NULL) {
    if (g_flags & OPT_DEBUG) {
      fprintf(stderr, "%s:%d: ", __FILE__, __LINE__);
      fprintf(stderr, "sdata_open(\"%s\") failed: ", file_name);
    }
    perror(file_name);
    exit (EXIT_FAILURE);
  }
  
  err = sdata_get_header(handle, &header);
  if (err != 0) {
    fprintf(stderr, "%s:%d: ", __FILE__, __LINE__);
    perror(file_name);
    exit (EXIT_FAILURE);
  }
  sdata_print_header(stdout, NULL, &header);
  
  err = sdata_close(handle);
  if (err != SDATA_OK) {
    if (g_flags & OPT_DEBUG) {
      fprintf(stderr, "%s:%d: sdata_close(handle [\"%s\"]) failed\n", __FILE__, __LINE__, file_name);
    }
    
    perror(file_name);
    exit (EXIT_FAILURE);
  }

}

int main(int argc, char *argv[]) {
  char long_appname[CSTRING_MAXLEN + 1];
  char appname[CSTRING_MAXLEN + 1];
  char appdir[CSTRING_MAXLEN + 1];
  const char *input_file;
  
  memset(long_appname, 0, sizeof(long_appname));
  strncpy(long_appname, argv[0], sizeof(long_appname) - 1);
  strncpy(appname, long_appname, sizeof(appname) - 1);
  strncpy(appdir, long_appname, sizeof(appdir) - 1);
  g_appname = basename(appname);
  g_appdir  = dirname(appdir);

  parse_options(&argc, &argv);
  if (argc != 2) {
    usage();
  }

  input_file = argv[1];

  if (g_flags & OPT_DEBUG) {
    sdata_set_debug_stream(stderr);
  }

  if (g_flags & OPT_HEADER) {
    print_header(input_file);
    exit (EXIT_SUCCESS);
  }

  if (g_flags & OPT_UNFOLD) {
    unfold_sdata(input_file);
  }
  else {
    int err;
    err = sdata_dump_file(stdout, input_file);
    if (err == SDATA_ERROR) {
      if (errno != 0) {
	fprintf(stderr, "%s: ", appname);
	perror(input_file);
	exit (EXIT_FAILURE);
      }
      else {
	fprintf(stderr, "%s: failed to read %s (try the -D flag for more info)\n", appname, input_file);
	exit (EXIT_FAILURE);
      }
    }
  }
  
  return EXIT_SUCCESS;
}
