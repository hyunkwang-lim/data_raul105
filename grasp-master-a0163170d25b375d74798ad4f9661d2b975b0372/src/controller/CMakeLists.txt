# Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
# Licensed under the GRASP Open Source License V1.0 (see LICENSE file)

include_directories(${CMAKE_CURRENT_BINARY_DIR} ${INPUT_BIN} ${INPUT_SRC}
	${OUTPUT_SRC} ${OUTPUT_BIN})

add_definitions(-DSPHEROID -DOSH -DGARRLIC)

if(${USE_MPI})
	find_package(MPI REQUIRED)
	add_definitions(-DUSE_MPI)
        include_directories(SYSTEM ${MPI_INCLUDE_PATH})
endif()

if(${DEBUG_MPI})
        add_definitions(-DDEBUG_MPI)
endif()

add_library(grasp_controller
  grasp_mpi_engine.c grasp_mpi_engine.h
  grasp_controller.c
)
target_link_libraries(grasp_controller
	grasp_retrieval_lib
	grasp_settings
	grasp_input
	grasp_global
	grasp_output
	grasputils
)
if(${USE_MPI})
        target_link_libraries(grasp_controller ${MPI_C_LIBRARIES})
endif()

add_executable(grasp_app grasp_main.c)
target_link_libraries(grasp_app grasp_controller)
set_property(TARGET grasp_app PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)
install(TARGETS grasp_app
	DESTINATION bin)

configure_file(grasp.in ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/grasp @ONLY)
install(PROGRAMS ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/grasp
  DESTINATION bin)

