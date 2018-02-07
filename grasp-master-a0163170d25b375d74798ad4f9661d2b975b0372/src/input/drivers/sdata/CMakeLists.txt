# Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
# Licensed under the GRASP Open Source License V1.0 (see LICENSE file)

add_library(input_driver_sdata
	sdata.c
	sdata.h
	sdata-impl.h
	grasp_input_driver_sdata.c
	grasp_input_driver_sdata.h
	check_helpers.c
	check_helpers.h
)
target_link_libraries(input_driver_sdata
	grasp_input_util
        grasp_input
	grasputils
	${MATH_LIB}
)

add_executable(input_driver_sdata_test test.c)
target_link_libraries(input_driver_sdata_test input_driver_sdata)

add_test(NAME input_driver_sdata
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
	COMMAND input_driver_sdata_test ${CMAKE_CURRENT_SOURCE_DIR}/testsdata.dat)

add_executable(sdata_dump sdata_dump.c)
target_link_libraries(sdata_dump grasp_input input_driver_sdata)
install(TARGETS sdata_dump DESTINATION bin)
