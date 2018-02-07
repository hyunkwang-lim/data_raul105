project(grasputils)

enable_language(Fortran)

add_library(grasputils
  grasp/allocator.c
  grasp/array_utils.c
  grasp/benchmark.c
  grasp/cmdline_utils.c
  grasp/error.c
  grasp/file_utils.c
  grasp/fortran_utils.c
  grasp/logging.c
  grasp/string_utils.c
  grasp/time_utils.c
  grasp/trackmem.c
  grasp_f90/mod_c_utils.f90
  grasp_f90/get_sizeof_time_t.c
  )
target_include_directories(grasputils PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
