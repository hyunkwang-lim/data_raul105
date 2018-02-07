# Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
# Licensed under the GRASP Open Source License V1.0 (see LICENSE file)

project(retrieval)
cmake_minimum_required(VERSION 2.8.0)

enable_language(Fortran)

# Take care of retrieval modules
# Parameter defining which modules are used
# Example: -DENABLE_MODULES=module_a,module_b
set(ENABLE_RETRIEVAL_MODULES "" CACHE LIST "List of retrieval modules to use (comma separated)")

# Initialize variables used by modules
set(RETRIEVAL_MODULES_SRCS "")
set(RETRIEVAL_MODULES_LIBS "")
set(RETRIEVAL_MODULES_DEFINITIONS "")
set(INVERSION_SUBSYSTEM "default")

# Add all submodules
foreach(module ${ENABLE_RETRIEVAL_MODULES})
	message(STATUS "Using module ${module}")
	add_subdirectory(${module})
endforeach()
message(STATUS "Using inversion subsystem \"${INVERSION_SUBSYSTEM}\"")

# Apply modified variable to retrieval
add_definitions(${RETRIEVAL_MODULES_DEFINITIONS})


# Calculate size of time type
include(CheckTypeSize)
check_type_size ("time_t" SIZEOF_TIME_T)
math(EXPR SIZEOF_TIME_T "${SIZEOF_TIME_T}*8")
message(STATUS "SIZEOF_TIME_T: ${SIZEOF_TIME_T}")

# Explicitly request C preprocessing
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DSIZEOF_TIME_T=${SIZEOF_TIME_T}")

add_definitions(-DSPHEROID -DOSH -DGARRLIC)
#add_definitions(-DSPHEROID)

set(SPARSE_LIB "SuperLU" CACHE STRING "Sparse Solver Library")
message(STATUS "Sparse solver: ${SPARSE_LIB}")
if(SPARSE_LIB STREQUAL SuperLU)
	add_definitions(-DSUPERLU_43)

	find_library(SUPERLU_LIB NAMES superlu superlu_4.3)
        if (NOT SUPERLU_LIB)
          message(FATAL_ERROR "Could not find SuperLU")
        endif()
        find_path(SUPERLU_INCLUDES NAMES slu_ddefs.h PATH_SUFFIXES superlu)
	find_package(BLAS REQUIRED)

        include_directories(${SUPERLU_INCLUDES})
        add_library(superlu_bridge external_interfaces/c_fortran_dgssv.c)
        target_link_libraries(superlu_bridge ${SUPERLU_LIB} ${BLAS_LIBRARIES})
        set(SOLVER_LIB superlu_bridge)
elseif(SPARSE_LIB STREQUAL SuperLU_MT)
	add_definitions(-DSUPERLU_MT)

	find_library(SUPERLU_MT_LIB NAMES superlu_mt)
        if (NOT SUPERLU_MT_LIB)
          message(FATAL_ERROR "Could not find SuperLU_MT")
        endif()
        find_path(SUPERLU_MT_INCLUDES NAMES pdsp_defs.h PATH_SUFFIXES superlu_mt)
	find_package(BLAS REQUIRED)

        include_directories(${SUPERLU_MT_INCLUDES})
        add_library(superlu_mt_bridge external_interfaces/c_bridge_pdgssv.c)
        target_link_libraries(superlu_mt_bridge ${SUPERLU_MT_LIB} ${BLAS_LIBRARIES})
	set(SOLVER_LIB superlu_mt_bridge)
elseif(SPARSE_LIB STREQUAL ViennaCL)
	add_definitions(-DVIENNA_CL)

        add_subdirectory(../../libs/src/ViennaCL viennacl)
        set(SOLVER_LIB ViennaCLSolver)
elseif(SPARSE_LIB STREQUAL MUMPS)
	add_definitions(-DMUMPS)

	find_package(MPI REQUIRED)
	find_package(MUMPS REQUIRED)
        find_package(BLAS REQUIRED)

        include_directories(${MPI_Fortran_INCLUDE_PATH} ${MUMPS_INCLUDE_DIR})
        set(SOLVER_LIB ${MPI_Fortran_LIBRARIES} ${MUMPS_LIBRARIES} ${BLAS_LIBRARIES})
else()
	message(SEND_ERROR "Unknown solver ${SPARSE_LIB}")
endif()
include_directories(${SOLVER_INCLUDES})


set(CONSTANTS_SET "generic" CACHE STRING "Used constants set")

if((NOT EXISTS ${CMAKE_SOURCE_DIR}/src/retrieval/constants_set/${CONSTANTS_SET}/mod_par_DLS.inc) OR 
   (NOT EXISTS ${CMAKE_SOURCE_DIR}/src/retrieval/constants_set/${CONSTANTS_SET}/mod_par_DLS_bin.inc) OR 
   (NOT EXISTS ${CMAKE_SOURCE_DIR}/src/retrieval/constants_set/${CONSTANTS_SET}/mod_par_inv.inc) OR 
   (NOT EXISTS ${CMAKE_SOURCE_DIR}/src/retrieval/constants_set/${CONSTANTS_SET}/mod_par_OS.inc) )
    message(FATAL_ERROR "ERROR: Constants set " ${CONSTANTS_SET} " does not exist")
else() 
    message(STATUS "Constants set used: " ${CONSTANTS_SET})
endif()

configure_file(constants_set/${CONSTANTS_SET}/mod_par_DLS.inc
  ${CMAKE_CURRENT_BINARY_DIR}/mod_par_DLS.inc COPYONLY)
configure_file(constants_set/${CONSTANTS_SET}/mod_par_DLS_bin.inc
  ${CMAKE_CURRENT_BINARY_DIR}/mod_par_DLS_bin.inc COPYONLY)
configure_file(constants_set/${CONSTANTS_SET}/mod_par_inv.inc
  ${CMAKE_CURRENT_BINARY_DIR}/mod_par_inv.inc COPYONLY)
configure_file(constants_set/${CONSTANTS_SET}/mod_par_OS.inc
  ${CMAKE_CURRENT_BINARY_DIR}/mod_par_OS.inc COPYONLY)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

if(NOT INVERSION_SUBSYSTEM OR INVERSION_SUBSYSTEM STREQUAL "default")
    set(RETRIEVAL_MODULES_SRCS
        inversion/inversion_subsystem_default.f90)
endif()

if( EXISTS ${CMAKE_SOURCE_DIR}/src/retrieval/models/phmx_models.f90 ) 
    add_subdirectory(${CMAKE_SOURCE_DIR}/src/retrieval/models/)
    set(GRASP_MODELS_LIB "grasp_models")
    set_property(GLOBAL PROPERTY GRASP_MODELS "yes")
    add_definitions(-DGRASP_MODELS="on")
endif()

set(SRC
    constants_set/mod_globals.f90

    inversion/getoptions.f90
    inversion/inversion.f90
    inversion/mod_par_inv.f90
    inversion/mod_alloc_arrays.f90
    inversion/mod_edges.f90
    inversion/iguess_apriori.f90
    inversion/random_noise.f90
    inversion/mod_covariance_matrix.f90
    inversion/smoothness.f90
    inversion/fisher_matrix.f90
    inversion/matrix_inversion.f90
    inversion/residual.f90
    inversion/ERR_estimates.f90

    interfaces/forw_model.f90
    interfaces/forw_set_pixel_fit.f90
    interfaces/validators.f90
    interfaces/mod_fisher_matrix_ccs.f90
    interfaces/mod_retr_settings_derived_type.f90
    interfaces/mod_retr_settings_par_type_particles.f90
    interfaces/mod_retr_settings_par_type_surface.f90
    interfaces/read_input_settings.f90
    interfaces/set_retr_settings.f90
    interfaces/mod_retr_general_output_derived_type.f90
    interfaces/set_retr_general_output.f90
    interfaces/mod_sdata_derived_type.f90
    interfaces/mod_sdata_meas_type.f90
    interfaces/mod_sdata.f90
    interfaces/mod_index_cloud.f90
    interfaces/unpack_par_vect_ap.f90
    interfaces/mod_vertical_distr_derived_type.f90
    interfaces/vertical_distr.f90

    utilities/mod_stop_report.f90
    utilities/print_GOUT.f90
    utilities/print_utils.f90
    utilities/mod_inversion_utils.f90
    utilities/subset_index.f90
    utilities/mod_intrpl_linear.f90
    utilities/mod_intrpl_spline.f90
    utilities/mod_c_utils.f90
    utilities/time_utils.c
    utilities/mod_time_utils.f90
    utilities/inssor.f90

    forward_model/mixture_ref_index_chem.f90
    forward_model/ammno3.f
    forward_model/mod_molecular_scattering.f90

    forward_model/phase_matrix/spheroid_package/forw_spheroid_package.f90
    forward_model/phase_matrix/spheroid_package/mod_par_DLS.f90
    forward_model/phase_matrix/spheroid_package/mod_par_DLS_bin.f90
    forward_model/phase_matrix/spheroid_package/mod_type_DLS.f90
    forward_model/phase_matrix/spheroid_package/mod_alloc_kernels.f90
    forward_model/phase_matrix/spheroid_package/DLS_optchr_bin.f90
    forward_model/phase_matrix/spheroid_package/DLS_read_input_bin.f90
    forward_model/phase_matrix/spheroid_package/DLS_intrpl_orgn_bin.f90
    forward_model/phase_matrix/spheroid_package/DLS_intrpl_bin.f90
    forward_model/phase_matrix/spheroid_package/DLS_fixget_bin.f90
    forward_model/phase_matrix/spheroid_package/sizedstr.f90
    forward_model/phase_matrix/spheroid_package/phase_func.f90

    forward_model/radiative_transfer/RTSOS/forw_SOS_RT.f90
    forward_model/radiative_transfer/RTSOS/mod_par_OS.f90
    forward_model/radiative_transfer/RTSOS/mod_SOS_RT.f90
    forward_model/surface/Mod_BRM.f90

    forward_model/lidar_signal/GARRLiC/forw_lidar_garrlic.f90
    forward_model/lidar_signal/GARRLiC/lidar_signal.f90

    forward_model/bbflux/mod_abs_kd.f90
    forward_model/bbflux/mod_bbflux.f90
    forward_model/bbflux/BIB/sub_abs_kd.f

#    models/phmx_models.f90

    ../controller/mo_grasp_controller.f90

    ${RETRIEVAL_MODULES_SRCS}
)


add_library(grasp_retrieval_lib ${SRC})
target_link_libraries(grasp_retrieval_lib
	${SOLVER_LIB}
	${RETRIEVAL_MODULES_LIBS}
        ${GRASP_MODELS_LIB}
)

add_executable(grasp_retrieval inversion/main.f90)
target_link_libraries(grasp_retrieval
	grasp_retrieval_lib
)

install(TARGETS grasp_retrieval
	DESTINATION bin)

install(DIRECTORY internal_files/
  DESTINATION ${RESOURCES_PREFIX}/kernels)

