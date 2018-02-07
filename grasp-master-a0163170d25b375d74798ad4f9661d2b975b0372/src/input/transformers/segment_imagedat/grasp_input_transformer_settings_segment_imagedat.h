/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */


#ifndef GRASP_INPUT_TRANSFORMER_SETTINGS_SEGMENT_IMAGEDAT_H
#define GRASP_INPUT_TRANSFORMER_SETTINGS_SEGMENT_IMAGEDAT_H

#ifdef __cplusplus
extern "C" {
#endif

    typedef struct grasp_input_transformer_settings_segment_imagedat_t_ {
        // Input dat file path
        char input_file[_GBL_FILE_PATH_LEN];
    }grasp_input_transformer_settings_segment_imagedat_t;


#ifdef __cplusplus
}
#endif

#endif /* GRASP_INPUT_TRANSFORMER_SETTINGS_SEGMENT_IMAGEDAT_H */

