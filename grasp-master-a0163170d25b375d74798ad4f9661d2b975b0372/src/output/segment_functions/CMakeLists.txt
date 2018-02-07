# Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
# Licensed under the GRASP Open Source License V1.0 (see LICENSE file)

subdirlist(SUBDIRS ${CMAKE_CURRENT_LIST_DIR})

foreach(SUBDIR ${SUBDIRS})
    if(EXISTS "${CMAKE_CURRENT_LIST_DIR}/${SUBDIR}/CMakeLists.txt")
        add_subdirectory(${SUBDIR})
        LIST(APPEND segment_functions ${SUBDIR})
        set_property(GLOBAL APPEND PROPERTY OUTPUT_SEGMENT_FUNCTIONS ${SUBDIR})

        execute_process(WORKING_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/${SUBDIR}/"
            COMMAND /usr/bin/env git rev-parse --show-toplevel
            RESULT_VARIABLE result
            OUTPUT_STRIP_TRAILING_WHITESPACE
            OUTPUT_VARIABLE git_root)
        if(NOT ${git_root} STREQUAL ${CMAKE_SOURCE_DIR})
            execute_process(WORKING_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/${SUBDIR}/"
                    COMMAND /usr/bin/env git describe --exact-match HEAD
                    RESULT_VARIABLE result
                    ERROR_VARIABLE ignoring
                    OUTPUT_VARIABLE EXTENSION_VERSION)
            string(STRIP "${EXTENSION_VERSION}" EXTENSION_VERSION) 
            if( (NOT ${result} EQUAL 0) OR ("${EXTENSION_VERSION}" STREQUAL "undefined") ) #no tag
                    execute_process(WORKING_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/${SUBDIR}/"
                    COMMAND /usr/bin/env git rev-parse --short HEAD
                    RESULT_VARIABLE result
                    OUTPUT_VARIABLE EXTENSION_VERSION)

                    string(STRIP "${EXTENSION_VERSION}" EXTENSION_VERSION)  
                    if (NOT ${result} EQUAL 0) #no commit
                        set(EXTENSION_VERSION "unknown")
                    endif ()
            endif()
            set_property(GLOBAL APPEND PROPERTY OUTPUT_SEGMENT_FUNCTION_VERSIONS "${EXTENSION_VERSION}")
        else()
            set_property(GLOBAL APPEND PROPERTY OUTPUT_SEGMENT_FUNCTION_VERSIONS "-")
        endif()
    endif()
endforeach()

message(STATUS "Output segment functions detected: " "${segment_functions}" )

