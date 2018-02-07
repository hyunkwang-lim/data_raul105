# - Find the GribApi
# This module defines the following variables:
#  GribApi_INCLUDE_DIRS - include directories for GribApi
#  GribApi_LIBRARIES - libraries to link against GribApi
#  GribApi_FOUND - true if GribApi has been found and can be used

find_path(GribApi_INCLUDE_DIR grib_api.h)
find_library(GribApi_LIBRARY NAMES grib_api)

set(GribApi_INCLUDE_DIRS ${GribApi_INCLUDE_DIR})
set(GribApi_LIBRARIES ${GribApi_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GribApi
	FOUND_VAR GribApi_FOUND
	REQUIRED_VARS GribApi_LIBRARY GribApi_INCLUDE_DIR)

mark_as_advanced(GribApi_INCLUDE_DIR GribApi_LIBRARY)
