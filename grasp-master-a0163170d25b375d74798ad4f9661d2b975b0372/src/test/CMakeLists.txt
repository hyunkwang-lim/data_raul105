# Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
# Licensed under the GRASP Open Source License V1.0 (see LICENSE file)

find_package(CUNIT REQUIRED)

include_directories(${CMAKE_CURRENT_BINARY_DIR} ${INPUT_BIN} ${INPUT_SRC}
  ${OUTPUT_SRC} ${OUTPUT_BIN} ${CUNIT_INCLUDE_DIRS})

add_definitions(-DSPHEROID -DOSH -DGARRLIC)

add_executable(grasp_test_aod
	grasp_test_aod.c
)
target_link_libraries(grasp_test_aod
	grasp_controller
	grasp_settings
	grasp_input
	grasp_output
        ${CUNIT_LIBRARIES}
	grasputils
)

add_executable(grasp_test_sdata
	grasp_test_sdata.c
)
target_link_libraries(grasp_test_sdata
	grasp_controller
	grasp_settings
	grasp_input
	grasp_output
        ${CUNIT_LIBRARIES}
	grasputils
)

add_executable(grasp_test_output
	grasp_test_output.c
)
target_link_libraries(grasp_test_output
	grasp_controller
	grasp_settings
	grasp_output
        ${CUNIT_LIBRARIES}
)

add_test(NAME aod
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
	COMMAND grasp_test_aod)
add_test(NAME sdata
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
	COMMAND grasp_test_sdata)
add_test(NAME output
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
	COMMAND grasp_test_output)
