/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

/**
 * @file grasp_controller_settings.h
 * @author David Fuertes
 * @date 4 Dec 2013
 * @brief Controller settings to specify controller behavior
 *
 * Input YAML settings file contains a "controller" section where the user can
 * specify different options for the controller. Following file contains the structure
 * which will store these values
 */

#ifndef GRASP_CONTROLLER_SETTINGS_H
#define	GRASP_CONTROLLER_SETTINGS_H

#ifdef	__cplusplus
extern "C" {
#endif

#include "../retrieval/constants_set/mod_globals.inc" 

/**
 * Structure which contains settings to customize controller behavior
 */
typedef struct controller_settings_ {
    int  segment_range[2]; /**< @brief Range of segments that will be processed. If only one value is defined a specific retrieval will be done */ 
    int nsegment_range; /**< @brief Number of values used for defining segment_range */ 
    bool perform_retrieval; /**< @brief True if retrieval algorithm will be called. If you set false to this parameter only retrieval algorithm will be executed */ 
    bool compilation_information; /**< @brief True if compilation information will be printed */ 
    char track_mem_stream[_GBL_FILE_PATH_LEN]; /**< @brief Stream where all mallocs can dump the information */ 
    char stream_pattern[_GBL_FILE_PATH_LEN]; /**< @brief Stream where all controller information will be printed */ 
    int maximum_job_time; /**< @brief Maximum time of node job. If the process take more time than this (in seconds) it will be killed */ 
    int polling_time; /**< @brief Waiting tima in master node after chacks that everything is done */ 
} controller_settings_t;

#ifdef	__cplusplus
}
#endif

#endif	/* GRASP_CONTROLLER_SETTINGS_H */

