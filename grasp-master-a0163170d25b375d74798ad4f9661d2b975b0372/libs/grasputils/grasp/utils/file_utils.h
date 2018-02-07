/**
 * @file:  file_utils.c
 * @Author: Fabrice Ducos <fabrice.ducos@univ-lille1.fr>
 *
 */

#ifndef FILE_UTILS_
#define FILE_UTILS_

#include <stdlib.h>
#include <stdbool.h>

/* may rely on PATH_MAX later, but only if it can be done portably */
#ifndef FILEPATH_LEN_
#define FILEPATH_LEN_ 1023
#endif

typedef char filepath_t[FILEPATH_LEN_ + 1];
// Expand glob pattern and return paths in matches.
// Return value is number of files finded or:
// -1 No match
// -2 glob problem
// -3 not memory enogh for glob pattern
// -4 fatal error in glob library
extern int find_files_from_pattern(const char *pattern, filepath_t **matches);

// Checks if a file exists and is readable
extern bool is_readable_file(const char *fname);

// Return path from a string 
char *pathoffile(const char *file);

// Return filename from a string (all what is after /)
char *name_of_file(const char *file);

// Return filename from a string (all what is after /) removing everything after last point (.)
char *name_of_file_without_extension(const char *file);

// Return true if the path is absolute
bool isabsolute(const char *file);

#endif
