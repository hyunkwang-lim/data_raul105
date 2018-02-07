/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

/* 
 * File:   grasp_input_bridge_sdata_settings.h
 * Author: david
 *
 * Created on 3 de septiembre de 2014, 12:42
 */

#ifndef GRASP_INPUT_BRIDGE_SDATA_SETTINGS_H
#define	GRASP_INPUT_BRIDGE_SDATA_SETTINGS_H

#ifdef	__cplusplus
extern "C" {
#endif

    typedef struct grasp_input_driver_settings_sdata_t_ {
        // Activate debug mode in sdata driver
        bool debug;
    }grasp_input_driver_settings_sdata_t;


#ifdef	__cplusplus
}
#endif

#endif	/* GRASP_INPUT_BRIDGE_SDATA_SETTINGS_H */

