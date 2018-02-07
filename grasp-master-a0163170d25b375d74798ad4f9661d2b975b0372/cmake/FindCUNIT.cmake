# - Find the CUNIT C unit testing framework
# This module defines the following variables:
#  CUNIT_INCLUDE_DIRS - include directories for CUNIT
#  CUNIT_LIBRARIES - libraries to link against CUNIT
#  CUNIT_FOUND - true if CUNIT has been found and can be used

find_path(CUNIT_INCLUDE_DIR CUnit/CUnit.h)
find_library(CUNIT_LIBRARY NAMES cunit)

set(CUNIT_INCLUDE_DIRS ${CUNIT_INCLUDE_DIR})
set(CUNIT_LIBRARIES ${CUNIT_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CUNIT
	FOUND_VAR CUNIT_FOUND
	REQUIRED_VARS CUNIT_LIBRARY CUNIT_INCLUDE_DIR)

mark_as_advanced(CUNIT_INCLUDE_DIR CUNIT_LIBRARY)
