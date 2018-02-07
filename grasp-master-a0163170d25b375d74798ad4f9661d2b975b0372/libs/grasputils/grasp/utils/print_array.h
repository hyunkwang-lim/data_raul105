/**
 * @file:  print_array.h
 * @Author: Fabrice Ducos <fabrice.ducos@univ-lille1.fr>
 *
 */

#ifndef PRINT_ARRAY_H
#define PRINT_ARRAY_H

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

/* a debugging facility to display one-dimensional arrays of any type (multidimensional arrays may be added later).
 * <stream>          : type FILE *
 *                     is the stream (or file pointer) on which to write.
 *                     in most cases : stdout or stderr
 * <array_label>     : type const char *
 *                     a label to name the array during the display. It is optional and can be set to NULL
 *                     or "". When a label is to be displayed, it is followed by a colon and a space,
 *                     then by the array itself.
 *                     NOTE: the current implementation causes some compilers (e.g. gcc) to emit a
 *                     warning when <array_label> is set to NULL. This is harmless (the implementation
 *                     behaves correctly with a NULL <array_label>) but to avoid the warning 
 *                     use an empty string "" instead of NULL, it will have the same effect.
 * <array>           : any type of array with content displayable by printf (int *, double *, char * ...)
 *                     the address of the array
 * <nelements>       : type size_t
 *                     the number of valid elements of the array
 * <nmax_to_display> : type size_t
 *                     the maximal number of elements to display
 *                     (only the first and last elements will be displayed if nelements is too large, 
 *                     and an ellipsis will replace the supernumerary elements)
 *                     all the elements will be displayed if it is set to 0
 *                     if <nmax_to_display> is odd, the first elements will count one
 * <fmt_string>      : type const char *
 *                     a format string from the printf specification. It must be consistent
 *                     with the type of array to be displayed. The consistency control
 *                     is compiler-dependent.
 * <separator>       : type const char *
 *                     a short string that will be displayed between the elements
 *                     there is no limit in the number of characters in <separator>, but
 *                     most of the time it should be kept short (one or two characters)
 *                     for readability.
 *                     Examples of separators: spaces, tabulations, a comma, ...
 *
 * The difference between print_array and println_array is that the latter
 * adds a end-of-line character.
 *
 * Examples of usage:
 * const int k[] = { -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5 };
 * const double x[] = { 0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10. };
 * const char *str[] = { "hello", "the", "world" };
 *
 * size_t nk = sizeof(k)/sizeof(*k);
 * size_t nx = sizeof(x)/sizeof(*x);
 * size_t nstr = sizeof(str)/sizeof(*str);
 *
 * println_array(stdout, "k", k, nk, 0, "%d", ", ");
 * println_array(stdout, "x", x, nx, 6, "%g", ", ");
 * println_array(stdout, "", str, nstr, 0, "%s", " ");
 */

/* C does not support templates, so we will be using a macro here, to avoid code duplication between types */
#define print_array(stream, array_label, array, nelements, nmax_to_display, fmt_string, separator) do { \
    int i;								\
    int counter;							\
    									\
    assert(array != NULL);						\
    assert(fmt_string != NULL);						\
    if ((void *) array_label != NULL) {					\
      /* cast needed to allow compilation when array_label is a NULL pointer */	\
      if (*((char *) array_label) != '\0') {				\
	/* this statement causes a warning with some compilers (e.g. gcc) */ \
	/* when array_label is NULL. This is harmless since array_label was */ \
	/* checked earlier. Though, it is recommended to use empty strings */ \
	/* instead of NULL when no label is needed */			\
	fprintf(stream, "%s (%d): ", (const char *) array_label, (int) nelements); \
      }									\
    }									\
    if (nmax_to_display == 0 || nelements <= nmax_to_display) {		\
      /* displays all the elements */					\
      for (i = 0 ; i < nelements ; i++) {				\
	fprintf(stream, fmt_string, array[i]);				\
	if (i < nelements - 1) { fprintf(stream, "%s", separator); }	\
      }									\
    }									\
    else {								\
      counter = nmax_to_display;					\
      /* displays at most nmax_to_display elements */			\
      for (i = 0 ; i < (nmax_to_display + 1) / 2 ; i++) {		\
	fprintf(stream, fmt_string, array[i]);				\
	if (i < (nmax_to_display + 1) / 2 - 1) { fprintf(stream, "%s", separator); } \
	counter--;							\
      }									\
      fprintf(stream, " ... ");						\
      for (i = nelements - counter ; i < nelements ; i++) {		\
	fprintf(stream, fmt_string, array[i]);				\
	if (i < nelements - 1) { fprintf(stream, "%s", separator); }	\
      }									\
    }									\
  } while (0)

#define println_array(stream, array_label, array, nelements, nmax_to_display, fmt_string, separator) do { \
    print_array(stream, array_label, array, nelements, nmax_to_display, fmt_string, separator); \
    fprintf(stream, "\n");						\
  } while(0)

#endif /* PRINT_ARRAY_H */
