lib_LTLIBRARIES = libgrasputils.la
pkginclude_HEADERS = utils.h
nobase_pkginclude_HEADERS = utils/error.h utils/file_utils.h utils/time_utils.h \
                            utils/cmdline_utils.h utils/allocator.h utils/benchmark.h \
                            utils/logging.h utils/print_array.h utils/debug.h  \
                            utils/string_utils.h utils/fortran_utils.h utils/array_utils.h

libgrasputils_la_SOURCES = \
  error.c \
  file_utils.c \
  time_utils.c \
  cmdline_utils.c \
  allocator.c \
  logging.c \
  fortran_utils.c \
  array_utils.c \
  benchmark.c \
  string_utils.c

libgrasputils_la_LDFLAGS = -version-info @VERSION_INFO_CURRENT@:@VERSION_INFO_REVISION@:@VERSION_INFO_AGE@
