/**
 * @file:  cmdline_utils.c
 * @Author: Fabrice Ducos <fabrice.ducos@univ-lille1.fr>
 *
 */

#include "utils/cmdline_utils.h"

void print_command_line(FILE *stream, int argc, char *argv[]) {
  int i;

  for (i = 0 ; i < argc - 1 ; i++) {
    fprintf(stream, "%s ", argv[i]);
  }
  fprintf(stream, "%s\n", argv[argc - 1]);
}

