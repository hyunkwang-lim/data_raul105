/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <grasp/utils.h>
#include "grasp_compilation_information.h"
#include "mod_par_inv.inc"
#include "../global/grasp_runtime_information.h"

void grasp_compilation_information_print(FILE *f){
    char *tmp;
    
    fprintf(f, " GRASP core version: %s (commit: %s ; branch_name: %s)\n",GRASP_VERSION, GRASP_COMMIT_REF, GRASP_BRANCH_NAME);
    fprintf(f, " Compiled on %s commit of %s\n", GRASP_COMPILATION_DATE, GRASP_COMMIT_DATE);
    fprintf(f, " With C compiler: %s \n", GRASP_CC_VERSION);
    fprintf(f, " With FORTRAN compiler: %s \n", GRASP_FC_VERSION);
    fprintf(f, " Using %s constant set and build type %s ", GRASP_CONSTANTS_SET, GRASP_BUILD_TYPE);
    if(strcmp(GRASP_MPI,"yes")==0){
        fprintf(f, "(using mpi)");
    }
    if(strcmp(GRASP_MODELS,"yes")==0){
        fprintf(f, " (models module compiled)");
    }
    printf("\n");
    fprintf(f, " Maximum segment size: nx=%d ; ny=%d ; nt=%d\n", _KIX, _KIY, _KITIME);
    tmp=grasp_compilation_information_input_drivers_versioned();
    fprintf(f, " Input drivers loaded: %s\n", tmp);
    free(tmp);
    tmp=grasp_compilation_information_input_transformers_versioned();
    fprintf(f, " Input transformers loaded: %s\n", tmp);
    free(tmp);
    tmp=grasp_compilation_information_output_segment_functions_versioned();
    fprintf(f, " Output segment functions loaded: %s\n", tmp);
    free(tmp);
    tmp=grasp_compilation_information_output_tile_functions_versioned();
    fprintf(f, " Output tile functions loaded: %s\n", tmp);
    free(tmp);
    tmp=grasp_compilation_information_output_current_functions_versioned();
    fprintf(f, " Output current functions loaded: %s\n", tmp);
    free(tmp);
    fprintf(f, " Path to resources: %s\n", RESOURCES_PREFIX);
    fprintf(f, " Sparse solver used: %s\n", GRASP_SPARSE_SOLVER);
    fprintf(f, " Build System: %s\n", GRASP_BUILD_SYSTEM);
    
    // Set GRASP_EXEC global variable with the absolute path to current executable
    fprintf(f, " Executable path: ");
    if(grasp_exec_file[0]=='/'){
        fprintf(f, "%s\n", grasp_exec_file);
    }else{
        fprintf(f, "undefined absolute path (%s)\n", grasp_exec_file);
    }
}

void grasp_version_print(FILE *f){
    if(strcmp(GRASP_VERSION, "[undefined]")==0 || strcmp(GRASP_VERSION, "[no version]")==0 ){
        fprintf(f, "GRASP core is not in a specific version. Reference commit: %s ; branch_name: %s\n", GRASP_COMMIT_REF, GRASP_BRANCH_NAME);        
    }else{
        fprintf(f, "GRASP core version: %s\n", GRASP_VERSION);
    }
}

char *grasp_compilation_information_input_drivers(){
    char *r;
    
    r = (char *) trackmem_malloc(sizeof (char)*(strlen(GRASP_INPUT_DRIVERS)+1));
    assert(r!=NULL);
    
    strcpy(r, GRASP_INPUT_DRIVERS);
        
    return r;
}

char *grasp_compilation_information_input_transformers(){
    char *r;
    
    r = (char *) trackmem_malloc(sizeof (char)*(strlen(GRASP_INPUT_TRANSFORMERS)+1));
    assert(r!=NULL);
    strcpy(r, GRASP_INPUT_TRANSFORMERS);
        
    return r;    
}

char *grasp_compilation_information_output_segment_functions(){
    char *r;
    
    r = (char *) trackmem_malloc(sizeof (char)*(strlen(GRASP_OUTPUT_SEGMENT_FUNCTIONS)+1));
    assert(r!=NULL);
    
    strcpy(r, GRASP_OUTPUT_SEGMENT_FUNCTIONS);
        
    return r;    
}

char *grasp_compilation_information_output_tile_functions(){
    char *r;
    
    r = (char *) trackmem_malloc(sizeof (char)*(strlen(GRASP_OUTPUT_TILE_FUNCTIONS)+1));
    assert(r!=NULL);
    
    strcpy(r, GRASP_OUTPUT_TILE_FUNCTIONS);
        
    return r;    
}

char *grasp_compilation_information_output_current_functions(){
    char *r;
    
    r = (char *) trackmem_malloc(sizeof (char)*(strlen(GRASP_OUTPUT_CURRENT_FUNCTIONS)+1));
    assert(r!=NULL);
    
    strcpy(r, GRASP_OUTPUT_CURRENT_FUNCTIONS);
        
    return r;    
}


char *grasp_compilation_information_input_drivers_versioned(){
    char extensions[]=GRASP_INPUT_DRIVERS;
    char versions[]=GRASP_INPUT_DRIVER_VERSIONS;
    char *extension, *version;
    char *rest_extensions, *rest_versions;
    char *result;
    char *ptr_extensions, *ptr_versions;
    
    result = (char *) trackmem_malloc(sizeof (char)*(strlen(GRASP_INPUT_DRIVERS)*3));
    assert(result!=NULL);
    
    strcpy(result, "");
    
    ptr_extensions=extensions;
    ptr_versions=versions;

    while((extension = strtok_r(ptr_extensions, " ", &rest_extensions))) {
        version = strtok_r(ptr_versions, " ", &rest_versions);
        
        strcat(result,extension);
        if(strcmp(version,"-")!=0){
            strcat(result,"(");
            strcat(result,version);
            strcat(result,")");
        }
        strcat(result, " ");
        
        ptr_extensions = rest_extensions;
        ptr_versions = rest_versions;
    }

    return result;
}

char *grasp_compilation_information_input_drivers_version(const char *driver){
    char extensions[]=GRASP_INPUT_DRIVERS;
    char versions[]=GRASP_INPUT_DRIVER_VERSIONS;
    char *extension, *version;
    char *rest_extensions, *rest_versions;
    char *result;
    char *ptr_extensions, *ptr_versions;
    
    result = (char *) trackmem_malloc(sizeof(char)*25);
    assert(result!=NULL);
    
    strcpy(result, "undefined");
    
    ptr_extensions=extensions;
    ptr_versions=versions;

    while((extension = strtok_r(ptr_extensions, " ", &rest_extensions))) {
        version = strtok_r(ptr_versions, " ", &rest_versions);
        
        if(strcmp(version,"-")!=0 && strcmp(extension,driver)==0){
            strcpy(result,version);
        }

        ptr_extensions = rest_extensions;
        ptr_versions = rest_versions;
    }

    return result;
}

char *grasp_compilation_information_input_transformers_versioned(){
    char extensions[]=GRASP_INPUT_TRANSFORMERS;
    char versions[]=GRASP_INPUT_TRANSFORMER_VERSIONS;
    char *extension, *version;
    char *rest_extensions, *rest_versions;
    char *result;
    char *ptr_extensions, *ptr_versions;
    
    result = (char *) trackmem_malloc(sizeof (char)*(strlen(GRASP_INPUT_TRANSFORMERS)*3));
    assert(result!=NULL);
    
    strcpy(result, "");
    
    ptr_extensions=extensions;
    ptr_versions=versions;

    while((extension = strtok_r(ptr_extensions, " ", &rest_extensions))) {
        version = strtok_r(ptr_versions, " ", &rest_versions);
        
        strcat(result,extension);
        if(strcmp(version,"-")!=0){
            strcat(result,"(");
            strcat(result,version);
            strcat(result,")");
        }
        strcat(result, " ");
        
        ptr_extensions = rest_extensions;
        ptr_versions = rest_versions;
    }

    return result;
}

char *grasp_compilation_information_input_transformers_version(const char *transformer){
    char extensions[]=GRASP_INPUT_TRANSFORMERS;
    char versions[]=GRASP_INPUT_TRANSFORMER_VERSIONS;
    char *extension, *version;
    char *rest_extensions, *rest_versions;
    char *result;
    char *ptr_extensions, *ptr_versions;
    
    result = (char *) trackmem_malloc(sizeof(char)*25);
    assert(result!=NULL);
    
    strcpy(result, "undefined");
    
    ptr_extensions=extensions;
    ptr_versions=versions;

    while((extension = strtok_r(ptr_extensions, " ", &rest_extensions))) {
        version = strtok_r(ptr_versions, " ", &rest_versions);
        
        if(strcmp(version,"-")!=0 && strcmp(extension,transformer)==0){
            strcpy(result,version);
        }

        ptr_extensions = rest_extensions;
        ptr_versions = rest_versions;
    }

    return result;
}

char *grasp_compilation_information_output_segment_functions_versioned(){
    char extensions[]=GRASP_OUTPUT_SEGMENT_FUNCTIONS;
    char versions[]=GRASP_OUTPUT_SEGMENT_FUNCTION_VERSIONS;
    char *extension, *version;
    char *rest_extensions, *rest_versions;
    char *result;
    char *ptr_extensions, *ptr_versions;
    
    result = (char *) trackmem_malloc(sizeof (char)*(strlen(GRASP_OUTPUT_SEGMENT_FUNCTIONS)*3));
    assert(result!=NULL);
    
    strcpy(result, "");
    
    ptr_extensions=extensions;
    ptr_versions=versions;

    while((extension = strtok_r(ptr_extensions, " ", &rest_extensions))) {
        version = strtok_r(ptr_versions, " ", &rest_versions);
        
        strcat(result,extension);
        if(strcmp(version,"-")!=0){
            strcat(result,"(");
            strcat(result,version);
            strcat(result,")");
        }
        strcat(result, " ");
        
        ptr_extensions = rest_extensions;
        ptr_versions = rest_versions;
    }

    return result;
}

char *grasp_compilation_information_output_segment_functions_version(const char *segment_function){
    char extensions[]=GRASP_OUTPUT_SEGMENT_FUNCTIONS;
    char versions[]=GRASP_OUTPUT_SEGMENT_FUNCTION_VERSIONS;
    char *extension, *version;
    char *rest_extensions, *rest_versions;
    char *result;
    char *ptr_extensions, *ptr_versions;
    
    result = (char *) trackmem_malloc(sizeof(char)*25);
    assert(result!=NULL);
    
    strcpy(result, "undefined");
    
    ptr_extensions=extensions;
    ptr_versions=versions;

    while((extension = strtok_r(ptr_extensions, " ", &rest_extensions))) {
        version = strtok_r(ptr_versions, " ", &rest_versions);
        
        if(strcmp(version,"-")!=0 && strcmp(extension,segment_function)==0){
            strcpy(result,version);
        }

        ptr_extensions = rest_extensions;
        ptr_versions = rest_versions;
    }

    return result;
}

char *grasp_compilation_information_output_tile_functions_versioned(){
    char extensions[]=GRASP_OUTPUT_TILE_FUNCTIONS;
    char versions[]=GRASP_OUTPUT_TILE_FUNCTION_VERSIONS;
    char *extension, *version;
    char *rest_extensions, *rest_versions;
    char *result;
    char *ptr_extensions, *ptr_versions;
    
    result = (char *) trackmem_malloc(sizeof (char)*(strlen(GRASP_OUTPUT_TILE_FUNCTIONS)*3));
    assert(result!=NULL);
    
    strcpy(result, "");
    
    ptr_extensions=extensions;
    ptr_versions=versions;

    while((extension = strtok_r(ptr_extensions, " ", &rest_extensions))) {
        version = strtok_r(ptr_versions, " ", &rest_versions);
        
        strcat(result,extension);
        if(strcmp(version,"-")!=0){
            strcat(result,"(");
            strcat(result,version);
            strcat(result,")");
        }
        strcat(result, " ");
        
        ptr_extensions = rest_extensions;
        ptr_versions = rest_versions;
    }

    return result;
}

char *grasp_compilation_information_output_tile_functions_version(const char *tile_function){
    char extensions[]=GRASP_OUTPUT_TILE_FUNCTIONS;
    char versions[]=GRASP_OUTPUT_TILE_FUNCTION_VERSIONS;
    char *extension, *version;
    char *rest_extensions, *rest_versions;
    char *result;
    char *ptr_extensions, *ptr_versions;
    
    result = (char *) trackmem_malloc(sizeof(char)*25);
    assert(result!=NULL);
    
    strcpy(result, "undefined");
    
    ptr_extensions=extensions;
    ptr_versions=versions;

    while((extension = strtok_r(ptr_extensions, " ", &rest_extensions))) {
        version = strtok_r(ptr_versions, " ", &rest_versions);
        
        if(strcmp(version,"-")!=0 && strcmp(extension,tile_function)==0){
            strcpy(result,version);
        }

        ptr_extensions = rest_extensions;
        ptr_versions = rest_versions;
    }

    return result;
}

char *grasp_compilation_information_output_current_functions_versioned(){
    char extensions[]=GRASP_OUTPUT_CURRENT_FUNCTIONS;
    char versions[]=GRASP_OUTPUT_CURRENT_FUNCTION_VERSIONS;
    char *extension, *version;
    char *rest_extensions, *rest_versions;
    char *result;
    char *ptr_extensions, *ptr_versions;
    
    result = (char *) trackmem_malloc(sizeof (char)*(strlen(GRASP_OUTPUT_CURRENT_FUNCTIONS)*3));
    assert(result!=NULL);
    
    strcpy(result, "");
    
    ptr_extensions=extensions;
    ptr_versions=versions;

    while((extension = strtok_r(ptr_extensions, " ", &rest_extensions))) {
        version = strtok_r(ptr_versions, " ", &rest_versions);
        
        strcat(result,extension);
        if(strcmp(version,"-")!=0){
            strcat(result,"(");
            strcat(result,version);
            strcat(result,")");
        }
        strcat(result, " ");
        
        ptr_extensions = rest_extensions;
        ptr_versions = rest_versions;
    }

    return result;    
}

char *grasp_compilation_information_output_current_functions_version(const char *current_function){
    char extensions[]=GRASP_OUTPUT_CURRENT_FUNCTIONS;
    char versions[]=GRASP_OUTPUT_CURRENT_FUNCTION_VERSIONS;
    char *extension, *version;
    char *rest_extensions, *rest_versions;
    char *result;
    char *ptr_extensions, *ptr_versions;
    
    result = (char *) trackmem_malloc(sizeof(char)*25);
    assert(result!=NULL);
    
    strcpy(result, "undefined");
    
    ptr_extensions=extensions;
    ptr_versions=versions;

    while((extension = strtok_r(ptr_extensions, " ", &rest_extensions))) {
        version = strtok_r(ptr_versions, " ", &rest_versions);
        
        if(strcmp(version,"-")!=0 && strcmp(extension,current_function)==0){
            strcpy(result,version);
        }

        ptr_extensions = rest_extensions;
        ptr_versions = rest_versions;
    }

    return result;
}

char *grasp_compilation_information_version(){
    char *r;
    
    r = (char *) trackmem_malloc(sizeof (char)*(strlen(GRASP_VERSION)+1));
    assert(r!=NULL);
    
    strcpy(r, GRASP_VERSION);
        
    return r;        
}

char *grasp_compilation_information_version_named(){
    char *result;
    
    result=grasp_compilation_information_version();

    if(strcmp(result,"undefined")==0){
        free(result);
        result=grasp_compilation_information_commit_ref();
    }

    return result;
}

char *grasp_compilation_information_branch_name(){
    char *r;
    
    r = (char *) trackmem_malloc(sizeof (char)*(strlen(GRASP_BRANCH_NAME)+1));
    assert(r!=NULL);
    
    strcpy(r, GRASP_BRANCH_NAME);
        
    return r;        
}

char *grasp_compilation_information_commit_ref(){
    char *r;
    
    r = (char *) trackmem_malloc(sizeof (char)*(strlen(GRASP_COMMIT_REF)+1));
    assert(r!=NULL);
    
    strcpy(r, GRASP_COMMIT_REF);
        
    return r;        
}

char *grasp_compilation_information_constants_set(){
    char *r;
    
    r = (char *) trackmem_malloc(sizeof (char)*(strlen(GRASP_CONSTANTS_SET)+1));
    assert(r!=NULL);
    
    strcpy(r, GRASP_CONSTANTS_SET);
        
    return r;        
}

bool grasp_compilation_information_models_present(){
    if (strcmp(GRASP_MODELS, "yes")==0){
        return true;
    }else{
        return false;
    }
}