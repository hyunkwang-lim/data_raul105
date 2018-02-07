# Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
# Licensed under the GRASP Open Source License V1.0 (see LICENSE file)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(retrieval)
add_subdirectory(controller)
add_subdirectory(input)
add_subdirectory(output)
add_subdirectory(global)
add_subdirectory(settings)
add_subdirectory(test)
