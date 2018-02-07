/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

/* 
 * File:   grasp_compilation_information.h
 * Author: fuertes
 *
 * Created on June 14, 2014, 9:46 PM
 */

#ifndef GRASP_COMPILATION_INFORMATION_H
#define	GRASP_COMPILATION_INFORMATION_H

#ifdef	__cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdbool.h>
    
#ifndef GRASP_VERSION
#define GRASP_VERSION "[undefined]"
#endif
    
#ifndef GRASP_BRANCH_NAME   
#define GRASP_BRANCH_NAME "[undefined]"
#endif
    
#ifndef GRASP_COMMIT_REF     
#define GRASP_COMMIT_REF "[undefined]"
#endif
    
#ifndef GRASP_COMMIT_DATE     
#define GRASP_COMMIT_DATE "[undefined]"
#endif    
    
#ifndef GRASP_COMPILATION_DATE     
#define GRASP_COMPILATION_DATE "[undefined]"
#endif
        
#ifndef GRASP_CONSTANTS_SET
#define GRASP_CONSTANTS_SET "[undefined]"
#endif
    
#ifndef RESOURCES_PREFIX
#define RESOURCES_PREFIX "/usr/local/share/grasp/"
#endif    
    
#ifndef GRASP_CC_VERSION
#define GRASP_CC_VERSION "[undefined]"
#endif

#ifndef GRASP_FC_VERSION
#define GRASP_FC_VERSION "[undefined]"
#endif  
    
#ifndef GRASP_BUILD_TYPE
#define GRASP_BUILD_TYPE "[undefined]"
#endif 
    
#ifndef GRASP_BUILD_SYSTEM
#define GRASP_BUILD_SYSTEM "[undefined]"
#endif     
    
#ifndef GRASP_SPARSE_SOLVER
#define GRASP_SPARSE_SOLVER "[undefined]"
#endif     
    
#ifndef GRASP_INPUT_DRIVERS
#define GRASP_INPUT_DRIVERS "no_input_drivers"
#endif      
    
#ifndef GRASP_INPUT_TRANSFORMERS
#define GRASP_INPUT_TRANSFORMERS "no_input_transformers"
#endif          
    
#ifndef GRASP_OUTPUT_SEGMENT_FUNCTIONS
#define GRASP_OUTPUT_SEGMENT_FUNCTIONS "no_output_segment_functions"
#endif         
    
#ifndef GRASP_OUTPUT_TILE_FUNCTIONS
#define GRASP_OUTPUT_TILE_FUNCTIONS "no_output_tile_functions"
#endif         
    
#ifndef GRASP_OUTPUT_CURRENT_FUNCTIONS
#define GRASP_OUTPUT_CURRENT_FUNCTIONS "no_output_current_functions"
#endif        
   
#ifndef GRASP_INPUT_DRIVER_VERSIONS
#define GRASP_INPUT_DRIVER_VERSIONS "no-defined"
#endif      
    
#ifndef GRASP_INPUT_TRANSFORMER_VERSIONS
#define GRASP_INPUT_TRANSFORMER_VERSIONS "no-defined"
#endif          
    
#ifndef GRASP_OUTPUT_SEGMENT_FUNCTION_VERSIONS
#define GRASP_OUTPUT_SEGMENT_FUNCTION_VERSIONS "no-defined"
#endif         
    
#ifndef GRASP_OUTPUT_TILE_FUNCTION_VERSIONS
#define GRASP_OUTPUT_TILE_FUNCTION_VERSIONS "no-defined"
#endif         
    
#ifndef GRASP_OUTPUT_CURRENT_FUNCTION_VERSIONS
#define GRASP_OUTPUT_CURRENT_FUNCTION_VERSIONS "no-defined"
#endif   
    
#ifndef GRASP_MODELS
#define GRASP_MODELS "no"
#endif       
         
#ifdef USE_MPI
#define GRASP_MPI "yes"
#else
#define GRASP_MPI "no"    
#endif          
    
// Function that print in screen compilation information. f is output stream. argv0 is argv[0] in main function of application
void grasp_compilation_information_print(FILE *f);

// This function will print GRASP code version information in a line
void grasp_version_print(FILE *f);

// Return GRASP_INPUT_DRIVERS allocated
char *grasp_compilation_information_input_drivers();

// Return GRASP_INPUT_TRANSFORMERS allocated
char *grasp_compilation_information_input_transformers();

// Return GRASP_OUTPUT_SEGMENT_FUNCTIONS allocated
char *grasp_compilation_information_output_segment_functions();

// Return GRASP_OUTPUT_TILE_FUNCTIONS allocated
char *grasp_compilation_information_output_tile_functions();

// Return GRASP_OUTPUT_CURRENT_FUNCTIONS allocated
char *grasp_compilation_information_output_current_functions();

// Return GRASP_INPUT_DRIVERS with version information allocated
char *grasp_compilation_information_input_drivers_versioned();

// Return version of a specific driver allocated
char *grasp_compilation_information_input_drivers_version(const char *driver);

// Return GRASP_INPUT_TRANSFORMERS with version information allocated
char *grasp_compilation_information_input_transformers_versioned();

// Return version of a specific transformer allocated
char *grasp_compilation_information_input_transformers_version(const char *transformer);

// Return GRASP_OUTPUT_SEGMENT_FUNCTIONS with version information allocated
char *grasp_compilation_information_output_segment_functions_versioned();

// Return version of a specific segment function allocated
char *grasp_compilation_information_output_segment_functions_version(const char *segment_function);

// Return GRASP_OUTPUT_TILE_FUNCTIONS with version information allocated
char *grasp_compilation_information_output_tile_functions_versioned();

// Return version of a specific tile function allocated
char *grasp_compilation_information_output_tile_functions_version(const char *tile_function);

// Return GRASP_OUTPUT_CURRENT_FUNCTIONS with version information allocated
char *grasp_compilation_information_output_current_functions_versioned();

// Return version of a specific current function allocated
char *grasp_compilation_information_output_current_functions_version(const char *current_function);

// Return GRASP_VERSION allocated
char *grasp_compilation_information_version();

// Return an allocated string with GRASP_VERSION or commit name if GRASP_VERSION is 'undefined'
char *grasp_compilation_information_version_named();

// Return GRASP_BRANCH_NAME allocated
char *grasp_compilation_information_branch_name();

// Return GRASP_COMMIT_REF allocated
char *grasp_compilation_information_commit_ref();

// Return GRASP_CONSTANTS_SET allocated
char *grasp_compilation_information_constants_set();

// Return if GRASP_MODELS module is present or not
bool grasp_compilation_information_models_present();

#ifdef	__cplusplus
}
#endif

#endif	/* GRASP_COMPILATION_INFORMATION_H */

