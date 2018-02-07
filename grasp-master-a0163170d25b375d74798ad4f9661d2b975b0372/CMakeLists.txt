# Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
# Licensed under the GRASP Open Source License V1.0 (see LICENSE file)

project(grasp)

macro(subdirlist result curdir)
  file(GLOB children RELATIVE ${curdir} *)
  set(dirlist "")
  foreach(child ${children})
    if(IS_DIRECTORY ${curdir}/${child})
        list(APPEND dirlist ${child})
    endif()
  endforeach()
  set(${result} ${dirlist})
endmacro()

# Macro for removing a custom flag from the c/c++ flags for build types Release and RelWithDebinfo
macro(remove_c_cxx_flag_for_release_and_relwithdebinfo flag)
  string(REPLACE "${flag}" "" CMAKE_C_FLAGS_RELEASE          "${CMAKE_C_FLAGS_RELEASE}")
  string(REPLACE "${flag}" "" CMAKE_CXX_FLAGS_RELEASE        "${CMAKE_CXX_FLAGS_RELEASE}")
  string(REPLACE "${flag}" "" CMAKE_C_FLAGS_RELWITHDEBINFO   "${CMAKE_C_FLAGS_RELWITHDEBINFO}")
  string(REPLACE "${flag}" "" CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
endmacro()

set(USE_MPI off CACHE BOOL "Build the controller with MPI support.")
set(DEBUG_MPI off CACHE BOOL "Build the MPI engine in advanced debug mode.")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

cmake_minimum_required(VERSION 2.8.0)
enable_language(Fortran)

#Set build type if not set
if (NOT CMAKE_BUILD_TYPE)
  message(STATUS "No build type selected, default to Release")
  set(CMAKE_BUILD_TYPE "Release")
endif()

# Custom build types
set(CMAKE_C_FLAGS_FAST "${CMAKE_CFLAGS_Release}")
set(CMAKE_CXX_FLAGS_FAST "${CMAKE_CXXFLAGS_Release}")
set(CMAKE_EXE_LINKER_FLAGS_FAST "${CMAKE_EXE_LINKER_FLAGS_Release}")
set(CMAKE_SHARED_LINKER_FLAGS_FAST "${CMAKE_SHARED_LINKER_FLAGS_Release}")

# Set configuration for Fortran compilers
if("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
  set(CMAKE_Fortran_FLAGS " -fpp ")
  set(CMAKE_Fortran_FLAGS_DEBUG " -DDEBUG -O0 -g -check bounds -traceback ")
  set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO " -O2 -g -check bounds -traceback ")
  set(CMAKE_Fortran_FLAGS_RELEASE " -O3 ") # Depending on the version of ifort 01 does not affect to results. 
  set(CMAKE_Fortran_FLAGS_FAST " -Ofast ")
else() # compiler like gfortran
  set(CMAKE_Fortran_FLAGS " -cpp -lgfortran -ffree-line-length-none")
  set(CMAKE_Fortran_FLAGS_DEBUG " -DDEBUG -O0 -g -fbacktrace -fbounds-check -finit-real=snan -ffpe-trap=invalid,zero,overflow,underflow ")
  set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO " -O2 -g -fbacktrace -fbounds-check ")
  set(CMAKE_Fortran_FLAGS_RELEASE " -O3 ")
  set(CMAKE_Fortran_FLAGS_FAST " -Ofast ")
endif()

# Tune compiler warnings for the various build configurations
#set(CMAKE_C_FLAGS_RELEASE          "${CMAKE_C_FLAGS_RELEASE}")
#set(CMAKE_CXX_FLAGS_RELEASE        "${CMAKE_CXX_FLAGS_RELEASE}")

set(CMAKE_C_FLAGS_RELWITHDEBINFO   "${CMAKE_C_FLAGS_RELWITHDEBINFO}   ")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} ")

set(CMAKE_C_FLAGS_DEBUG            "${CMAKE_C_FLAGS_DEBUG}            -Wall -Wno-unused-function")
set(CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_CXX_FLAGS_DEBUG}          -Wall -Wno-unused-function")

# For now, we don't want -DNDEBUG set for Release and RelWithDebInfo
remove_c_cxx_flag_for_release_and_relwithdebinfo("-DNDEBUG")

if (APPLE)
  set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} -Wl,-no_compact_unwind")
  set(CMAKE_C_LINK_FLAGS "${CMAKE_C_LINK_FLAGS} -Wl,-no_compact_unwind")
endif (APPLE)  
 
# Resource path
set(RESOURCES_PREFIX share/grasp/)
add_definitions(-DRESOURCES_PREFIX="${CMAKE_INSTALL_PREFIX}/${RESOURCES_PREFIX}")

# Set search directories for grasp-deps
set(ADDITIONAL_DEPENDENCIES_PATH "/usr/local/grasp-deps" CACHE PATH "Search prefix for additional grasp dependencies")

set(CMAKE_PREFIX_PATH
  "${ADDITIONAL_DEPENDENCIES_PATH}"
  "${CMAKE_PREFIX_PATH}")

include_directories(${ADDITIONAL_DEPENDENCIES_PATH}/include)

if(CMAKE_BUILD_TYPE STREQUAL Debug)
  add_definitions(-DSHOW_WARNINGS)
  add_definitions(-DDEBUG)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fcheck=bounds")
endif()

set(PROFILING off CACHE BOOL "Enable profiling")
if(PROFILING)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pg")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -pg")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
endif()

find_library(MATH_LIB NAMES m)

enable_testing()

set(CONSTANTS_INCLUDE_DIR ${PROJECT_BINARY_DIR}/src/retrieval)
set(CONSTANTS_SRC_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/src/retrieval)
include_directories(${CONSTANTS_INCLUDE_DIR} ${CONSTANTS_SRC_INCLUDE_DIR})

set(RETRIEVAL_BIN ${CMAKE_CURRENT_BINARY_DIR}/src/retrieval)
set(INPUT_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/input)
set(INPUT_BIN ${CMAKE_CURRENT_BINARY_DIR}/src/input)
set(OUTPUT_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/output)
set(OUTPUT_BIN ${CMAKE_CURRENT_BINARY_DIR}/src/output)


set(ENABLE_CODE_COVERAGE "off"	  CACHE BOOL    "Enable code coverage report")
if(ENABLE_CODE_COVERAGE)
	if(NOT(CMAKE_BUILD_TYPE MATCHES "Debug"))
		message(FATAL_ERROR "Code coverage requires debug mode (CMAKE_BUILD_TYPE)")
	endif()
	message("Code coverage enabled")
	include(CodeCoverage)
	set(CMAKE_CXX_FLAGS "-g -O0 -fprofile-arcs -ftest-coverage")
	set(CMAKE_C_FLAGS "-g -O0 -fprofile-arcs -ftest-coverage")
	set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -O0 -fprofile-arcs -ftest-coverage")
	configure_file(src/test/grasp-coverage.sh ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/grasp-coverage.sh COPYONLY)
	setup_target_for_coverage(coverage ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/grasp-coverage.sh coverage)
endif()

add_subdirectory(libs)
#add_subdirectory(resources)
add_subdirectory(src)
