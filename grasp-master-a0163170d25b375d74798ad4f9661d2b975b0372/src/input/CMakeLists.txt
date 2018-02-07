# Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
# Licensed under the GRASP Open Source License V1.0 (see LICENSE file)

include_directories(${RETRIEVAL_BIN} ${CMAKE_CURRENT_SOURCE_DIR}
	${CMAKE_CURRENT_BINARY_DIR} ${OUTPUT_SRC} ${OUTPUT_BIN}
	${CMAKE_CURRENT_BINARY_DIR}/drivers/sdata)   # FIXME ugly hack to fix relative include

add_subdirectory(util)
add_subdirectory(drivers)
add_subdirectory(transformers)

set(GRASP_INPUT_DRIVER_SETTINGS_DRIVER_INCLUDE "")
set(GRASP_INPUT_DRIVER_SETTINGS_STRUCT "")
set(GRASP_INPUT_LOAD_DRIVER_INCLUDE "")
set(GRASP_INPUT_LOAD_DRIVER_FUNCTION "")
set(GRASP_INPUT_LOAD_DRIVER_PARAMETER "")
set(GRASP_INPUT_LOAD_DRIVER_NAMES "")
get_property(INPUT_DRIVERS GLOBAL PROPERTY INPUT_DRIVERS)
list(LENGTH INPUT_DRIVERS GRASP_INPUT_LOAD_DRIVER_NDRIVERS)
set(index 0)
set(nsettings 0)
foreach(d ${INPUT_DRIVERS})
	list(APPEND DRIVER_LIBS input_driver_${d})

	if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/drivers/${d}/grasp_input_driver_settings_${d}.h")
		set(GRASP_INPUT_DRIVER_SETTINGS_DRIVER_INCLUDE "${GRASP_INPUT_DRIVER_SETTINGS_DRIVER_INCLUDE}\n#include \"drivers/${d}/grasp_input_driver_settings_${d}.h\"")
		set(GRASP_INPUT_DRIVER_SETTINGS_STRUCT "${GRASP_INPUT_DRIVER_SETTINGS_STRUCT}\n\tgrasp_input_driver_settings_${d}_t ${d};")
                math(EXPR nsettings "${nsettings}+1")
        endif()
	set(GRASP_INPUT_LOAD_DRIVER_INCLUDE "${GRASP_INPUT_LOAD_DRIVER_INCLUDE}\n#include \"drivers/${d}/grasp_input_driver_${d}.h\"")
        set(GRASP_INPUT_LOAD_DRIVER_FUNCTION "${GRASP_INPUT_LOAD_DRIVER_FUNCTION}\n\t}else if(strcmp(name,\"${d}\")==0){ return grasp_input_driver_${d}();")
	set(GRASP_INPUT_LOAD_DRIVER_PARAMETER "${GRASP_INPUT_LOAD_DRIVER_PARAMETER}\n\tresult[${index}]=grasp_input_driver_settings_${d}(settings);")
        set(GRASP_INPUT_LOAD_DRIVER_NAMES "${GRASP_INPUT_LOAD_DRIVER_NAMES}            case ${index}: return \"${d}\";\n")
	MATH(EXPR index "${index}+1")
endforeach()
if(nsettings EQUAL 0)
    set(GRASP_INPUT_DRIVER_SETTINGS_STRUCT "${GRASP_INPUT_DRIVER_SETTINGS_STRUCT}\n\t\tvoid *noelements;")
endif()

set(GRASP_INPUT_TRANSFORMER_SETTINGS_TRANSFORMER_INCLUDE "")
set(GRASP_INPUT_TRANSFORMER_SETTINGS_STRUCT "")
set(GRASP_INPUT_LOAD_TRANSFORMER_INCLUDE "")
set(GRASP_INPUT_LOAD_TRANSFORMER_FUNCTION "")
set(GRASP_INPUT_LOAD_TRANSFORMER_PARAMETER "")
set(GRASP_INPUT_LOAD_TRANSFORMER_NAMES "")
get_property(INPUT_TRANSFORMERS GLOBAL PROPERTY INPUT_TRANSFORMERS)
list(LENGTH INPUT_TRANSFORMERS GRASP_INPUT_LOAD_TRANSFORMER_NTRANSFORMERS)
set(index 0)
set(nsettings 0)
foreach(d ${INPUT_TRANSFORMERS})
	list(APPEND TRANSFORMER_LIBS input_transformer_${d})

	if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/transformers/${d}/grasp_input_transformer_settings_${d}.h")
		set(GRASP_INPUT_TRANSFORMER_SETTINGS_TRANSFORMER_INCLUDE "${GRASP_INPUT_TRANSFORMER_SETTINGS_TRANSFORMER_INCLUDE}\n#include \"transformers/${d}/grasp_input_transformer_settings_${d}.h\"")
		set(GRASP_INPUT_TRANSFORMER_SETTINGS_STRUCT "${GRASP_INPUT_TRANSFORMER_SETTINGS_STRUCT}\n\tgrasp_input_transformer_settings_${d}_t ${d};")
                math(EXPR nsettings "${nsettings}+1")
        endif()
	set(GRASP_INPUT_LOAD_TRANSFORMER_INCLUDE "${GRASP_INPUT_LOAD_TRANSFORMER_INCLUDE}\n#include \"transformers/${d}/grasp_input_transformer_${d}.h\"")
        set(GRASP_INPUT_LOAD_TRANSFORMER_FUNCTION "${GRASP_INPUT_LOAD_TRANSFORMER_FUNCTION}\n\t}else if(strcmp(name,\"${d}\")==0){ return grasp_input_transformer_${d}();")
	set(GRASP_INPUT_LOAD_TRANSFORMER_PARAMETER "${GRASP_INPUT_LOAD_TRANSFORMER_PARAMETER}\n\tresult[${index}]=grasp_input_transformer_settings_${d}(settings);")
        set(GRASP_INPUT_LOAD_TRANSFORMER_NAMES "${GRASP_INPUT_LOAD_TRANSFORMER_NAMES}            case ${index}: return \"${d}\";\n")
	MATH(EXPR index "${index}+1")
endforeach()
if(nsettings EQUAL 0)
    set(GRASP_INPUT_TRANSFORMER_SETTINGS_STRUCT "${GRASP_INPUT_TRANSFORMER_SETTINGS_STRUCT}\n\t\tvoid *noelements;")
endif()

configure_file(grasp_input_functions_settings_template.h
	${CMAKE_CURRENT_BINARY_DIR}/grasp_input_functions_settings.h)
configure_file(grasp_input_load_functions_template.c
	${CMAKE_CURRENT_BINARY_DIR}/grasp_input_load_functions.c)

# Full input library with all drivers.
add_library(grasp_input
        grasp_input.c
        grasp_input_segment.c
	${CMAKE_CURRENT_BINARY_DIR}/grasp_input_load_functions.c
)
target_link_libraries(grasp_input
  grasp_input_util
  ${DRIVER_LIBS} ${TRANSFORMER_LIBS})
