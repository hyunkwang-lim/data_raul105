/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

#include <stdio.h>
#include <string.h>
#include <grasp/utils.h>
#include "../retrieval/constants_set/mod_globals.inc"

char grasp_exec_file[_GBL_FILE_PATH_LEN];

char grasp_current_path[_GBL_FILE_PATH_LEN];

char grasp_main_settings_file[_GBL_FILE_PATH_LEN];


void grasp_runtime_initialize(char **argv){
    strcpy(grasp_exec_file, "");
    strcpy(grasp_current_path, "");
    strcpy(grasp_main_settings_file, "");
}

void grasp_runtime_set(char *exec_file, char *current_path, char *main_settings_file){
    char *f;
    strcpy(grasp_exec_file, exec_file);
    
    f=pathoffile(current_path);
    strcpy(grasp_current_path, f);
    trackmem_free(f);
    
    strcpy(grasp_main_settings_file, main_settings_file);    
}

void grasp_runtime_debug(){
    printf("exec file: %s\n", grasp_exec_file);
    printf("current path: %s\n", grasp_current_path);
    printf("main settings file: %s\n", grasp_main_settings_file);
}

char *grasp_runtime_settings_file(){
    char *result;
    
    result=(char *)malloc(sizeof(char)*_GBL_FILE_PATH_LEN);
    
    strcpy(result, grasp_main_settings_file);
    
    return result;
}
    
char *grasp_runtime_settings_file_filename(){
    return name_of_file(grasp_main_settings_file);  
}
    
