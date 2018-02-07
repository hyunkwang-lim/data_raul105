/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

/* 
 * File:   grasp_output_stream_t.h
 * Author: fuertes
 *
 * Created on June 27, 2014, 4:36 PM
 */

#ifndef GRASP_OUTPUT_STREAM_T_H
#define	GRASP_OUTPUT_STREAM_T_H

#ifdef	__cplusplus
extern "C" {
#endif

#include <stdbool.h>
#include <stdio.h>
    
typedef struct GRASP_STREAM_{
    char *filename; // Patern of filenames that it will generate
    bool none;      // If this is not going to print anything
    bool screen;    // If it stream is going to print in screen
    FILE *file;     // File opened by this stram
    bool writable;  // If it stream is writable
    bool open;      // If the file of this stream is opened
} grasp_output_stream;


#ifdef	__cplusplus
}
#endif

#endif	/* GRASP_OUTPUT_STREAM_T_H */

