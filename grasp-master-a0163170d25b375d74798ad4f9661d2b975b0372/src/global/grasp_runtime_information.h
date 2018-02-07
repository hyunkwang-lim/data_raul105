/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

/*
 * File:   grasp_runtime_information.h
 * Author: david
 *
 * Created on 2 de julio de 2014, 14:08
 */

#ifndef GRASP_RUNTIME_INFORMATION_H
#define	GRASP_RUNTIME_INFORMATION_H

#ifdef	__cplusplus
extern "C" {
#endif

#include "../retrieval/constants_set/mod_globals.inc"

    extern char grasp_exec_file[_GBL_FILE_PATH_LEN];

    extern char grasp_current_path[_GBL_FILE_PATH_LEN];

    extern char grasp_main_settings_file[_GBL_FILE_PATH_LEN];

    void grasp_runtime_initialize(char *argv[]);
    
    void grasp_runtime_set(char *exec_file, char *current_path, char *main_settings_file);

    void grasp_runtime_debug(); 
    
    /**
     * Get current settings file (complete path)
     * @return Allocated string with settings file (including complete path)
     */
    char *grasp_runtime_settings_file();

    /**
     * Get filename of current settings file
     * @return Allocated string with filename of current settings file
     */    
    char *grasp_runtime_settings_file_filename();


#ifdef	__cplusplus
}
#endif

#endif	/* GRASP_RUNTIME_INFORMATION_H */

