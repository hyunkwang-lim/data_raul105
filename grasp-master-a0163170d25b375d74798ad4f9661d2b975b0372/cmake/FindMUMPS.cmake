# - Find the MUMPS direct solver for large sparse linear systems
# This module defines the following variables:
#  MUMPS_INCLUDE_DIRS - include directories for MUMPS
#  MUMPS_LIBRARIES - libraries to link against MUMPS
#  MUMPS_FOUND - true if MUMPS has been found and can be used

find_path(MUMPS_INCLUDE_DIR mumps_compat.h)
find_library(MUMPS_LIBRARY NAMES dmumps)

set(MUMPS_INCLUDE_DIRS ${MUMPS_INCLUDE_DIR})
set(MUMPS_LIBRARIES ${MUMPS_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MUMPS
	FOUND_VAR MUMPS_FOUND
	REQUIRED_VARS MUMPS_LIBRARY MUMPS_INCLUDE_DIR)

mark_as_advanced(MUMPS_INCLUDE_DIR MUMPS_LIBRARY)
