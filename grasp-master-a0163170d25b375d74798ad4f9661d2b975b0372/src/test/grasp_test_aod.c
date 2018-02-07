/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

/*
 * File:   newcunittest.c
 * Author: david
 *
 * Created on 21-jul-2012, 21:37:59
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <CUnit/CUnit.h>
#include <CUnit/Basic.h>
#include <CUnit/Console.h>
#include <CUnit/Automated.h>
#include <grasp/utils.h>
#include <time.h>
#include <stdbool.h>
#include "../controller/grasp_main.h"
#include "yamlsettings/yamlsettings.h"
#include "../input/grasp_input.h"
#include "../output/grasp_output.h"
#include "../settings/grasp_settings.h"
#include "../controller/grasp_controller.h"
#include "../controller/mo_grasp_controller.h"
//#include "../output/mo_grasp_output.h"
/*
 * CUnit Test Suite
 */

int init_suite(void) {
    return 0;
}

int clean_suite(void) {
    return 0;
}

void grasp_aod_test() {
    /*
    int argc=2;
    char *argv[]={"","aod/settings_aeronet.yml"};
    grasp_settings *settings=NULL;
    grasp_processing_functions_t functions;
    grasp_tile_description_t tile_description;
    grasp_results_t results;
  
    
    // 1) Read settings and process controller options
    settings=grasp_controller_read_settings(argc,argv);
    
    // 2) Get Tile information
    grasp_controller_initialize_inversion(settings, &tile_description, &functions, &results);
    
    // 3) Invert tile data
    settings->retrieval.DLSF.keyEL = CONSTANT_KEYEL;
    grasp_controller_invert_tile(settings, &tile_description, &results, &functions);
    
    // 4) Manage output
    grasp_controller_manage_tile(settings, &tile_description, &results, &functions);
    
    // 5) Make test
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[0] , 0.008728, (( 0.008728 +results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[0])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[1] , 0.011616, (( 0.011616 +results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[1])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[2] , 0.015264, (( 0.015264 +results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[2])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[3] , 0.019219, (( 0.019219 +results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[3])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[4] , 0.022074, (( 0.022074 +results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[4])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[5] , 0.021656, (( 0.021656 +results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[5])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[6] , 0.017600, (( 0.017600 +results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[6])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[7] , 0.012325, (( 0.012325 +results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[7])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[8] , 0.008297, (( 0.008297 +results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[8])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[9] , 0.006069, (( 0.006069 +results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[9])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[10] , 0.005131, (( 0.005131 +results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[10])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[11] , 0.004926, (( 0.004926 +results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[11])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[12] , 0.005120, (( 0.005120 +results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[12])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[13] , 0.005548, (( 0.005548 +results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[13])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[14] , 0.006157, (( 0.006157 +results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[14])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[15] , 0.006942, (( 0.006942 +results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[15])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[16] , 0.007887, (( 0.007887 +results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[16])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[17] , 0.008986, (( 0.008986 +results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[17])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[18] , 0.010242, (( 0.010242 +results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[18])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[19] , 0.011666, (( 0.011666 +results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[19])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[20] , 0.013278, (( 0.013278 +results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[20])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[21] , 0.015108, (( 0.015108 +results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[21])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[22] , 1.449975, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[23] , 1.450011, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[24] , 1.450018, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[25] , 1.450021, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[26] , 1.449986, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[27] , 1.449993, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[28] , 1.449992, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[29] , 1.450005, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[30] , 0.020000, (( 0.020000 +results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[30])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[31] , 0.020000, (( 0.020000 +results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[31])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[32] , 0.020000, (( 0.020000 +results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[32])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[33] , 0.020000, (( 0.020000 +results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[33])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[34] , 0.020000, (( 0.020000 +results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[34])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[35] , 0.020000, (( 0.020000 +results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[35])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[36] , 0.020000, (( 0.020000 +results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[36])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[37] , 0.020000, (( 0.020000 +results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[37])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[38] , 0.999900, 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[39] , 2000.000244, 200.0);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[40] , 0.069888, 0.01);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[0] , 0.015435, (( 0.015435 +results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[0])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[1] , 0.021245, (( 0.021245 +results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[1])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[2] , 0.028798, (( 0.028798 +results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[2])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[3] , 0.037094, (( 0.037094 +results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[3])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[4] , 0.042774, (( 0.042774 +results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[4])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[5] , 0.040592, (( 0.040592 +results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[5])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[6] , 0.031228, (( 0.031228 +results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[6])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[7] , 0.020779, (( 0.020779 +results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[7])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[8] , 0.013405, (( 0.013405 +results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[8])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[9] , 0.009343, (( 0.009343 +results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[9])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[10] , 0.007252, (( 0.007252 +results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[10])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[11] , 0.006046, (( 0.006046 +results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[11])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[12] , 0.005261, (( 0.005261 +results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[12])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[13] , 0.004725, (( 0.004725 +results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[13])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[14] , 0.004351, (( 0.004351 +results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[14])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[15] , 0.004086, (( 0.004086 +results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[15])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[16] , 0.003897, (( 0.003897 +results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[16])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[17] , 0.003760, (( 0.003760 +results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[17])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[18] , 0.003656, (( 0.003656 +results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[18])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[19] , 0.003573, (( 0.003573 +results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[19])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[20] , 0.003502, (( 0.003502 +results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[20])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[21] , 0.003435, (( 0.003435 +results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[21])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[22] , 1.449945, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[23] , 1.450037, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[24] , 1.450035, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[25] , 1.450039, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[26] , 1.449937, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[27] , 1.449991, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[28] , 1.450015, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[29] , 1.450002, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[30] , 0.020000, (( 0.020000 +results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[30])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[31] , 0.020000, (( 0.020000 +results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[31])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[32] , 0.020000, (( 0.020000 +results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[32])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[33] , 0.020000, (( 0.020000 +results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[33])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[34] , 0.020000, (( 0.020000 +results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[34])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[35] , 0.020000, (( 0.020000 +results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[35])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[36] , 0.020000, (( 0.020000 +results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[36])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[37] , 0.020000, (( 0.020000 +results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[37])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[38] , 0.999900, 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[39] , 2000.000244, 200.0);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[40] , 0.044740, 0.01);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[0] , 0.006252, (( 0.006252 +results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[0])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[1] , 0.008990, (( 0.008990 +results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[1])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[2] , 0.012721, (( 0.012721 +results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[2])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[3] , 0.017037, (( 0.017037 +results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[3])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[4] , 0.020219, (( 0.020219 +results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[4])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[5] , 0.019355, (( 0.019355 +results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[5])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[6] , 0.014702, (( 0.014702 +results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[6])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[7] , 0.009616, (( 0.009616 +results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[7])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[8] , 0.006294, (( 0.006294 +results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[8])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[9] , 0.004771, (( 0.004771 +results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[9])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[10] , 0.004389, (( 0.004389 +results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[10])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[11] , 0.004632, (( 0.004632 +results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[11])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[12] , 0.005256, (( 0.005256 +results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[12])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[13] , 0.006192, (( 0.006192 +results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[13])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[14] , 0.007438, (( 0.007438 +results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[14])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[15] , 0.009018, (( 0.009018 +results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[15])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[16] , 0.010957, (( 0.010957 +results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[16])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[17] , 0.013289, (( 0.013289 +results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[17])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[18] , 0.016064, (( 0.016064 +results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[18])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[19] , 0.019355, (( 0.019355 +results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[19])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[20] , 0.023264, (( 0.023264 +results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[20])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[21] , 0.027934, (( 0.027934 +results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[21])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[22] , 1.449937, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[23] , 1.450046, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[24] , 1.450039, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[25] , 1.450039, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[26] , 1.449935, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[27] , 1.449993, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[28] , 1.450008, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[29] , 1.450004, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[30] , 0.020000, (( 0.020000 +results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[30])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[31] , 0.020000, (( 0.020000 +results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[31])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[32] , 0.020000, (( 0.020000 +results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[32])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[33] , 0.020000, (( 0.020000 +results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[33])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[34] , 0.020000, (( 0.020000 +results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[34])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[35] , 0.020000, (( 0.020000 +results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[35])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[36] , 0.020000, (( 0.020000 +results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[36])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[37] , 0.020000, (( 0.020000 +results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[37])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[38] , 0.999900, 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[39] , 2000.000244, 200.0);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(2,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[40] , 0.082027, 0.01);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[0] , 0.014571, (( 0.014571 +results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[0])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[1] , 0.021259, (( 0.021259 +results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[1])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[2] , 0.030494, (( 0.030494 +results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[2])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[3] , 0.041294, (( 0.041294 +results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[3])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[4] , 0.049208, (( 0.049208 +results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[4])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[5] , 0.046470, (( 0.046470 +results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[5])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[6] , 0.034616, (( 0.034616 +results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[6])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[7] , 0.022314, (( 0.022314 +results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[7])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[8] , 0.014268, (( 0.014268 +results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[8])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[9] , 0.010202, (( 0.010202 +results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[9])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[10] , 0.008296, (( 0.008296 +results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[10])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[11] , 0.007217, (( 0.007217 +results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[11])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[12] , 0.006479, (( 0.006479 +results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[12])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[13] , 0.005964, (( 0.005964 +results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[13])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[14] , 0.005596, (( 0.005596 +results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[14])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[15] , 0.005324, (( 0.005324 +results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[15])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[16] , 0.005118, (( 0.005118 +results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[16])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[17] , 0.004957, (( 0.004957 +results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[17])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[18] , 0.004829, (( 0.004829 +results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[18])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[19] , 0.004721, (( 0.004721 +results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[19])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[20] , 0.004624, (( 0.004624 +results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[20])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[21] , 0.004533, (( 0.004533 +results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[21])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[22] , 1.449912, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[23] , 1.450072, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[24] , 1.450047, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[25] , 1.450056, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[26] , 1.449895, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[27] , 1.449985, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[28] , 1.450033, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[29] , 1.450001, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[30] , 0.020000, (( 0.020000 +results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[30])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[31] , 0.020000, (( 0.020000 +results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[31])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[32] , 0.020000, (( 0.020000 +results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[32])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[33] , 0.020000, (( 0.020000 +results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[33])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[34] , 0.020000, (( 0.020000 +results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[34])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[35] , 0.020000, (( 0.020000 +results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[35])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[36] , 0.020000, (( 0.020000 +results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[36])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[37] , 0.020000, (( 0.020000 +results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[37])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[38] , 0.999900, 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[39] , 2000.000244, 200.0);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(3,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[40] , 0.045427, 0.01);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[0] , 0.017170, (( 0.017170 +results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[0])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[1] , 0.022621, (( 0.022621 +results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[1])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[2] , 0.029411, (( 0.029411 +results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[2])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[3] , 0.036574, (( 0.036574 +results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[3])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[4] , 0.041248, (( 0.041248 +results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[4])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[5] , 0.039132, (( 0.039132 +results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[5])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[6] , 0.030681, (( 0.030681 +results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[6])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[7] , 0.020953, (( 0.020953 +results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[7])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[8] , 0.013793, (( 0.013793 +results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[8])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[9] , 0.009678, (( 0.009678 +results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[9])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[10] , 0.007497, (( 0.007497 +results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[10])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[11] , 0.006251, (( 0.006251 +results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[11])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[12] , 0.005473, (( 0.005473 +results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[12])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[13] , 0.004974, (( 0.004974 +results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[13])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[14] , 0.004659, (( 0.004659 +results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[14])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[15] , 0.004470, (( 0.004470 +results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[15])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[16] , 0.004369, (( 0.004369 +results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[16])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[17] , 0.004328, (( 0.004328 +results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[17])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[18] , 0.004328, (( 0.004328 +results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[18])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[19] , 0.004352, (( 0.004352 +results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[19])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[20] , 0.004389, (( 0.004389 +results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[20])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[21] , 0.004432, (( 0.004432 +results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[21])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[22] , 1.449952, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[23] , 1.450035, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[24] , 1.450028, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[25] , 1.450030, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[26] , 1.449950, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[27] , 1.449997, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[28] , 1.450007, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[29] , 1.450002, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[30] , 0.020000, (( 0.020000 +results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[30])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[31] , 0.020000, (( 0.020000 +results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[31])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[32] , 0.020000, (( 0.020000 +results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[32])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[33] , 0.020000, (( 0.020000 +results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[33])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[34] , 0.020000, (( 0.020000 +results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[34])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[35] , 0.020000, (( 0.020000 +results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[35])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[36] , 0.020000, (( 0.020000 +results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[36])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[37] , 0.020000, (( 0.020000 +results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[37])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[38] , 0.999900, 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[39] , 2000.000244, 200.0);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(4,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[40] , 0.048270, 0.01);
    //output[0].pixels[0].opt.ssa 
    //output[0].pixels[0].opt.ext
    //output[0].pixels[0].opt.salbedo
    //AOD, SSD, surface ALbedo
    //CU_ASSERT_DOUBLE_EQUAL(output[0].pixels[0].opt.ssat[0],0.87668 , 0.001);
    //CU_ASSERT_DOUBLE_EQUAL(output[0].pixels[1].opt.ssat[0],0.88189, 0.001); 

    
    // 6) Clean up 
    grasp_controller_clean_memory(settings, &tile_description, &results, &functions); 
    */
}

int main(int argc, char *argv[]) {
    CU_pSuite grasp_Suite = NULL;

    /* Initialize the CUnit test registry */
    if (CUE_SUCCESS != CU_initialize_registry())
        return CU_get_error();

    /* Add a suite to the registry */
    grasp_Suite = CU_add_suite("grasp_settings_test", init_suite, clean_suite);
    if (NULL == grasp_Suite) {
        CU_cleanup_registry();
        return CU_get_error();
    }
   if (
       (NULL == CU_add_test(grasp_Suite, "grasp_aod_test", grasp_aod_test))
       ) {
      CU_cleanup_registry();
      return CU_get_error();
   }        

    /* Run all tests using the CUnit Basic interface */
    //CU_basic_set_mode(CU_BRM_VERBOSE);
    CU_basic_run_tests();
    CU_cleanup_registry();
    return CU_get_error();
}
