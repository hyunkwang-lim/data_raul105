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


void grasp_polder_test() {
    /*
    int argc=2;
    char *argv[]={"","sdata/settings_sdata.yml"};
    grasp_settings *settings=NULL;
    grasp_processing_functions_t functions;
    grasp_tile_description_t tile_description;
    grasp_results_t results;
    
    // 1) Read settings and process controller options
    settings=grasp_controller_read_settings(argc,argv);
    
    // 2) Get Tile information
    grasp_controller_initialize_inversion(settings, &tile_description, &functions, &results);
    
    // 3) Invert tile data
    grasp_controller_invert_tile(settings, &tile_description, &results, &functions);
    
    // 4) Manage output
    grasp_controller_manage_tile(settings, &tile_description, &results, &functions);
    
    // 5) Make test
    
   CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[0] , 0.005223, (( 0.005223 +results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[0])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[1] , 0.005162, (( 0.005162 +results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[1])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[2] , 0.005264, (( 0.005264 +results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[2])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[3] , 0.006092, (( 0.006092 +results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[3])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[4] , 0.008215, (( 0.008215 +results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[4])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[5] , 0.019901, (( 0.019901 +results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[5])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[6] , 0.033468, (( 0.033468 +results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[6])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[7] , 0.033142, (( 0.033142 +results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[7])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[8] , 0.022947, (( 0.022947 +results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[8])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[9] , 1.469859, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[10] , 1.452866, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[11] , 1.436932, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[12] , 1.430197, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[13] , 1.426319, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[14] , 1.426342, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[15] , 0.005791, (( 0.005791 +results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[15])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[16] , 0.006190, (( 0.006190 +results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[16])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[17] , 0.006770, (( 0.006770 +results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[17])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[18] , 0.007072, (( 0.007072 +results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[18])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[19] , 0.007404, (( 0.007404 +results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[19])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[20] , 0.007556, (( 0.007556 +results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[20])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[21] , 0.049775, 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[22] , 4635.174316, 200.0);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[23] , 0.104850, 0.01);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[0] , 0.005224, (( 0.005224 +results.tile_result_map[index3D(0,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[0])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[1] , 0.005176, (( 0.005176 +results.tile_result_map[index3D(0,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[1])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[2] , 0.005292, (( 0.005292 +results.tile_result_map[index3D(0,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[2])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[3] , 0.006117, (( 0.006117 +results.tile_result_map[index3D(0,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[3])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[4] , 0.008223, (( 0.008223 +results.tile_result_map[index3D(0,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[4])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[5] , 0.019860, (( 0.019860 +results.tile_result_map[index3D(0,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[5])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[6] , 0.033309, (( 0.033309 +results.tile_result_map[index3D(0,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[6])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[7] , 0.032985, (( 0.032985 +results.tile_result_map[index3D(0,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[7])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[8] , 0.022881, (( 0.022881 +results.tile_result_map[index3D(0,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[8])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[9] , 1.490043, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[10] , 1.474743, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[11] , 1.459549, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[12] , 1.452929, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[13] , 1.454262, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[14] , 1.455382, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[15] , 0.005813, (( 0.005813 +results.tile_result_map[index3D(0,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[15])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[16] , 0.006173, (( 0.006173 +results.tile_result_map[index3D(0,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[16])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[17] , 0.006668, (( 0.006668 +results.tile_result_map[index3D(0,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[17])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[18] , 0.007012, (( 0.007012 +results.tile_result_map[index3D(0,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[18])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[19] , 0.007374, (( 0.007374 +results.tile_result_map[index3D(0,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[19])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[20] , 0.007506, (( 0.007506 +results.tile_result_map[index3D(0,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[20])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[21] , 0.049685, 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[22] , 4632.973633, 200.0);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[23] , 0.093840, 0.01);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[0] , 0.005235, (( 0.005235 +results.tile_result_map[index3D(0,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[0])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[1] , 0.005165, (( 0.005165 +results.tile_result_map[index3D(0,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[1])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[2] , 0.005254, (( 0.005254 +results.tile_result_map[index3D(0,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[2])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[3] , 0.006070, (( 0.006070 +results.tile_result_map[index3D(0,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[3])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[4] , 0.008234, (( 0.008234 +results.tile_result_map[index3D(0,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[4])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[5] , 0.019832, (( 0.019832 +results.tile_result_map[index3D(0,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[5])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[6] , 0.033505, (( 0.033505 +results.tile_result_map[index3D(0,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[6])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[7] , 0.033284, (( 0.033284 +results.tile_result_map[index3D(0,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[7])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[8] , 0.023031, (( 0.023031 +results.tile_result_map[index3D(0,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[8])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[9] , 1.448108, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[10] , 1.429603, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[11] , 1.412230, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[12] , 1.405000, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[13] , 1.408818, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[14] , 1.410823, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[15] , 0.005739, (( 0.005739 +results.tile_result_map[index3D(0,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[15])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[16] , 0.006182, (( 0.006182 +results.tile_result_map[index3D(0,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[16])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[17] , 0.006788, (( 0.006788 +results.tile_result_map[index3D(0,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[17])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[18] , 0.007152, (( 0.007152 +results.tile_result_map[index3D(0,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[18])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[19] , 0.007483, (( 0.007483 +results.tile_result_map[index3D(0,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[19])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[20] , 0.007568, (( 0.007568 +results.tile_result_map[index3D(0,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[20])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[21] , 0.049864, 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[22] , 4570.128906, 200.0);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[23] , 0.116986, 0.01);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[0] , 0.005233, (( 0.005233 +results.tile_result_map[index3D(0,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[0])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[1] , 0.005179, (( 0.005179 +results.tile_result_map[index3D(0,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[1])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[2] , 0.005287, (( 0.005287 +results.tile_result_map[index3D(0,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[2])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[3] , 0.006095, (( 0.006095 +results.tile_result_map[index3D(0,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[3])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[4] , 0.008210, (( 0.008210 +results.tile_result_map[index3D(0,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[4])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[5] , 0.019782, (( 0.019782 +results.tile_result_map[index3D(0,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[5])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[6] , 0.033256, (( 0.033256 +results.tile_result_map[index3D(0,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[6])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[7] , 0.033063, (( 0.033063 +results.tile_result_map[index3D(0,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[7])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[8] , 0.022937, (( 0.022937 +results.tile_result_map[index3D(0,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[8])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[9] , 1.461516, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[10] , 1.445101, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[11] , 1.427854, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[12] , 1.420931, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[13] , 1.421720, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[14] , 1.422720, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[15] , 0.005754, (( 0.005754 +results.tile_result_map[index3D(0,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[15])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[16] , 0.006172, (( 0.006172 +results.tile_result_map[index3D(0,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[16])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[17] , 0.006738, (( 0.006738 +results.tile_result_map[index3D(0,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[17])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[18] , 0.007101, (( 0.007101 +results.tile_result_map[index3D(0,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[18])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[19] , 0.007451, (( 0.007451 +results.tile_result_map[index3D(0,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[19])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[20] , 0.007554, (( 0.007554 +results.tile_result_map[index3D(0,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[20])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[21] , 0.049682, 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[22] , 4590.514648, 200.0);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(0,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[23] , 0.098389, 0.01);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[0] , 0.003374, (( 0.003374 +results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[0])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[1] , 0.003457, (( 0.003457 +results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[1])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[2] , 0.003893, (( 0.003893 +results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[2])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[3] , 0.005331, (( 0.005331 +results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[3])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[4] , 0.008919, (( 0.008919 +results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[4])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[5] , 0.035390, (( 0.035390 +results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[5])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[6] , 0.052441, (( 0.052441 +results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[6])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[7] , 0.036436, (( 0.036436 +results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[7])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[8] , 0.021449, (( 0.021449 +results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[8])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[9] , 1.436614, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[10] , 1.443463, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[11] , 1.456501, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[12] , 1.471510, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[13] , 1.494439, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[14] , 1.504884, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[15] , 0.001774, (( 0.001774 +results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[15])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[16] , 0.001722, (( 0.001722 +results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[16])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[17] , 0.001667, (( 0.001667 +results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[17])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[18] , 0.001643, (( 0.001643 +results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[18])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[19] , 0.001575, (( 0.001575 +results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[19])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[20] , 0.001564, (( 0.001564 +results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[20])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[21] , 0.044921, 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[22] , 4015.826416, 200.0);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[23] , 0.142797, 0.01);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[0] , 0.003364, (( 0.003364 +results.tile_result_map[index3D(1,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[0])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[1] , 0.003447, (( 0.003447 +results.tile_result_map[index3D(1,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[1])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[2] , 0.003904, (( 0.003904 +results.tile_result_map[index3D(1,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[2])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[3] , 0.005390, (( 0.005390 +results.tile_result_map[index3D(1,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[3])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[4] , 0.009012, (( 0.009012 +results.tile_result_map[index3D(1,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[4])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[5] , 0.035430, (( 0.035430 +results.tile_result_map[index3D(1,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[5])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[6] , 0.052120, (( 0.052120 +results.tile_result_map[index3D(1,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[6])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[7] , 0.036344, (( 0.036344 +results.tile_result_map[index3D(1,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[7])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[8] , 0.021407, (( 0.021407 +results.tile_result_map[index3D(1,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[8])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[9] , 1.442435, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[10] , 1.449716, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[11] , 1.458018, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[12] , 1.471177, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[13] , 1.497234, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[14] , 1.505905, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[15] , 0.001773, (( 0.001773 +results.tile_result_map[index3D(1,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[15])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[16] , 0.001722, (( 0.001722 +results.tile_result_map[index3D(1,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[16])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[17] , 0.001678, (( 0.001678 +results.tile_result_map[index3D(1,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[17])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[18] , 0.001648, (( 0.001648 +results.tile_result_map[index3D(1,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[18])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[19] , 0.001577, (( 0.001577 +results.tile_result_map[index3D(1,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[19])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[20] , 0.001565, (( 0.001565 +results.tile_result_map[index3D(1,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[20])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[21] , 0.045018, 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[22] , 4076.826416, 200.0);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,0,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[23] , 0.127732, 0.01);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[0] , 0.003365, (( 0.003365 +results.tile_result_map[index3D(1,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[0])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[1] , 0.003440, (( 0.003440 +results.tile_result_map[index3D(1,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[1])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[2] , 0.003880, (( 0.003880 +results.tile_result_map[index3D(1,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[2])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[3] , 0.005305, (( 0.005305 +results.tile_result_map[index3D(1,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[3])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[4] , 0.008909, (( 0.008909 +results.tile_result_map[index3D(1,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[4])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[5] , 0.035603, (( 0.035603 +results.tile_result_map[index3D(1,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[5])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[6] , 0.052768, (( 0.052768 +results.tile_result_map[index3D(1,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[6])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[7] , 0.036630, (( 0.036630 +results.tile_result_map[index3D(1,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[7])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[8] , 0.021533, (( 0.021533 +results.tile_result_map[index3D(1,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[8])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[9] , 1.430248, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[10] , 1.438237, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[11] , 1.452102, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[12] , 1.467352, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[13] , 1.493489, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[14] , 1.499396, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[15] , 0.001787, (( 0.001787 +results.tile_result_map[index3D(1,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[15])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[16] , 0.001731, (( 0.001731 +results.tile_result_map[index3D(1,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[16])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[17] , 0.001670, (( 0.001670 +results.tile_result_map[index3D(1,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[17])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[18] , 0.001638, (( 0.001638 +results.tile_result_map[index3D(1,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[18])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[19] , 0.001570, (( 0.001570 +results.tile_result_map[index3D(1,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[19])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[20] , 0.001562, (( 0.001562 +results.tile_result_map[index3D(1,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[20])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[21] , 0.045057, 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[22] , 4120.938965, 200.0);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,1,0,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[23] , 0.140343, 0.01);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[0] , 0.003354, (( 0.003354 +results.tile_result_map[index3D(1,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[0])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[1] , 0.003432, (( 0.003432 +results.tile_result_map[index3D(1,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[1])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[2] , 0.003877, (( 0.003877 +results.tile_result_map[index3D(1,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[2])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[3] , 0.005330, (( 0.005330 +results.tile_result_map[index3D(1,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[3])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[4] , 0.008967, (( 0.008967 +results.tile_result_map[index3D(1,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[4])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[5] , 0.035598, (( 0.035598 +results.tile_result_map[index3D(1,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[5])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[6] , 0.052659, (( 0.052659 +results.tile_result_map[index3D(1,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[6])/2) * 0.010000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[7] , 0.036671, (( 0.036671 +results.tile_result_map[index3D(1,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[7])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[8] , 0.021540, (( 0.021540 +results.tile_result_map[index3D(1,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[8])/2) * 0.030000);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[9] , 1.408391, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[10] , 1.417077, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[11] , 1.429990, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[12] , 1.447318, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[13] , 1.471122, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[14] , 1.478174, 0.005);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[15] , 0.001795, (( 0.001795 +results.tile_result_map[index3D(1,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[15])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[16] , 0.001737, (( 0.001737 +results.tile_result_map[index3D(1,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[16])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[17] , 0.001683, (( 0.001683 +results.tile_result_map[index3D(1,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[17])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[18] , 0.001648, (( 0.001648 +results.tile_result_map[index3D(1,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[18])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[19] , 0.001573, (( 0.001573 +results.tile_result_map[index3D(1,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[19])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[20] , 0.001563, (( 0.001563 +results.tile_result_map[index3D(1,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[20])/2) * 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[21] , 0.045022, 0.05);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[22] , 4223.780762, 200.0);
CU_ASSERT_DOUBLE_EQUAL(results.tile_result_map[index3D(1,1,1,tile_description.dimensions.tile_nx,tile_description.dimensions.tile_ny)]->parameters[23] , 0.146381, 0.01); 
    
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
            (NULL == CU_add_test(grasp_Suite, "grasp_polder_test", grasp_polder_test))
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
