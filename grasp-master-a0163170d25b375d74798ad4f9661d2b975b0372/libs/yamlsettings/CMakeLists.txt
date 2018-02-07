find_library(YAML_LIB NAMES yaml)
if (NOT YAML_LIB)
  message(FATAL_ERROR "YAML library not found")
endif()
find_path(LIBYAML_INCLUDE_DIR NAMES yaml.h)
if (NOT LIBYAML_INCLUDE_DIR)
  message(FATAL_ERROR "YAML header files not found")
endif()

include(FindPkgConfig)
pkg_check_modules(GLIB REQUIRED glib-2.0)
find_library(GLIB_LIB
	NAMES glib-2.0
	HINTS ${PC_GLIB_LIBDIR} ${PC_GLIB_LIBRARY_DIRS}
)
if (NOT GLIB_LIB)
  message(FATAL_ERROR "GLIB not found")
endif()

add_library(yamlsettings
	yamlsettings.c
	yamlsettings_input_yaml.c
	yamlsettings_assign_data.c
	yamlsettings_process_yaml.c
	yamlsettings_data_types.c
	yamlsettings_validators.c
	yamlsettings_dictionary.c
	yamlsettings_error.c
        yamlsettings_defaults.c
)
target_include_directories(yamlsettings PUBLIC ${GLIB_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/.. ${LIBYAML_INCLUDE_DIR})
target_link_libraries(yamlsettings ${GLIB_LIB} ${YAML_LIB} grasputils)

add_subdirectory(test)
