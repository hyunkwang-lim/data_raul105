#
# CMake build script for the ViennaCL sparse solver.
#
# Author: Michael Aspetsberger

project(ViennaCLSolver)

find_path(VIENNACL_DIR NAMES compressed_matrix.hpp PATH_SUFFIXES viennacl)
if(NOT VIENNACL_DIR)
	message(FATAL_ERROR "ViennaCL not found")
endif()

# set file sets
set(SOURCE
	viennacl_solve.cpp
)
set(HEADERS
	ViennaCLSolver.h
	ViennaCLSolver.inl
)

add_library(${PROJECT_NAME} ${SOURCE} ${HEADERS})
