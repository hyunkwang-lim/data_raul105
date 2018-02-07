# - Find the libcsv library
# This module defines the following variables:
#  CSV_INCLUDE_DIRS - include directories for CSV
#  CSV_LIBRARIES - libraries to link against CSV
#  CSV_FOUND - true if CSV has been found and can be used

find_path(CSV_INCLUDE_DIR csv.h)
find_library(CSV_LIBRARY NAMES csv)

set(CSV_INCLUDE_DIRS ${CSV_INCLUDE_DIR})
set(CSV_LIBRARIES ${CSV_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CSV
	FOUND_VAR CSV_FOUND
	REQUIRED_VARS CSV_LIBRARY CSV_INCLUDE_DIR)

mark_as_advanced(CSV_INCLUDE_DIR CSV_LIBRARY)
