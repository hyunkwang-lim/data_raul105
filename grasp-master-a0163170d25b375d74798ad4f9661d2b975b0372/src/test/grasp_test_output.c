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
#include <time.h>
#include <stdbool.h>
#include <assert.h>
#include "../output/grasp_output_stream.h"
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

void grasp_test_output_stream() {
    int argc=4;
    char *argv[]={"","sdata/settings_sdata.yml","input.time.from=2015-03-23T05:00:00Z","input.time.to=2015-03-23T06:00:00Z"};
    grasp_settings *settings=NULL;
    grasp_processing_functions_t functions;
    grasp_tile_description_t tile_description;
    grasp_results_t results;
    grasp_segment_t *segment;
    output_segment_general *output;
    
    char *good_pattern1="file.txt";    
    char *good_pattern2="file_{auto}.txt";
    char *good_pattern3="file_{auto(3)}.txt";    
    char *good_pattern4="file_{auto}_repeated(3).txt";
    char *good_pattern5="{yml}file_{auto}.txt";
    char *good_pattern6="";
    char *good_pattern7="{yml}file_{auto(2)}.txt";
    char *bad_pattern1="file_{noexist}.txt";
    char *bad_pattern2="file_{yml}.txt";
    char *bad_pattern3="file_{}.txt";
    char *bad_pattern4="file_{.txt";
    char *bad_pattern5="file_{{auto}}.txt";
    char *bad_pattern6="file_{auto{}}.txt";
    char *bad_pattern7="file_{auto}{.txt";
    char *bad_pattern8="file_{auto}.txt{";
    
    // 1) We need a settings structure to do examples. We will load one
    settings=grasp_controller_read_settings(argc,argv);
    
    // 2) Get Tile information
    grasp_controller_initialize_inversion(settings, &tile_description, &functions, &results);
    
    // 3) Prepare environment
    segment = (grasp_segment_t *) malloc(sizeof (grasp_segment_t)*1);
    assert( segment!= NULL);
    output = (output_segment_general *) malloc(sizeof (output_segment_general)*1);    
    assert(output);
    grasp_input_extract_segment(settings, &(functions.driver), functions.ntransformers, functions.transformers, segment, &results, &tile_description.dimensions, 0);   
    // 4) Make test  
    CU_ASSERT(true == grasp_output_stream_filename_validation(good_pattern1));
    CU_ASSERT(true == grasp_output_stream_filename_validation(good_pattern2));
    CU_ASSERT(true == grasp_output_stream_filename_validation(good_pattern3));
    CU_ASSERT(true == grasp_output_stream_filename_validation(good_pattern4));
    CU_ASSERT(true == grasp_output_stream_filename_validation(good_pattern5));
    CU_ASSERT(true == grasp_output_stream_filename_validation(good_pattern6));
    CU_ASSERT(true == grasp_output_stream_filename_validation(good_pattern7));
    CU_ASSERT(false == grasp_output_stream_filename_validation(bad_pattern1));
    CU_ASSERT(false == grasp_output_stream_filename_validation(bad_pattern2));
    CU_ASSERT(false == grasp_output_stream_filename_validation(bad_pattern3));
    CU_ASSERT(false == grasp_output_stream_filename_validation(bad_pattern4));
    CU_ASSERT(false == grasp_output_stream_filename_validation(bad_pattern5));
    CU_ASSERT(false == grasp_output_stream_filename_validation(bad_pattern6));
    CU_ASSERT(false == grasp_output_stream_filename_validation(bad_pattern7));
    CU_ASSERT(false == grasp_output_stream_filename_validation(bad_pattern8));
    
    CU_ASSERT_STRING_EQUAL("sdata/file_1x1x1.txt",grasp_output_stream_filename_generator(good_pattern5, settings, segment, output, &tile_description.dimensions , 1,1,1));
    CU_ASSERT_STRING_EQUAL("sdata/file_010101.txt",grasp_output_stream_filename_generator(good_pattern7, settings, segment, output, &tile_description.dimensions , 1,1,1));
    CU_ASSERT_STRING_EQUAL("",grasp_output_stream_filename_generator(bad_pattern1, settings, segment, output, &tile_description.dimensions , 1,1,1));
    // HERE WE SHOULD TEST ALL WILLCARDS USING grasp_output_stream_filename_generator . IT IS NOT DONE YET.
    
    // Tests to grasp_output_stream_filename_generate_by_* functions
    CU_ASSERT_STRING_EQUAL("1x1x1", grasp_output_stream_filename_generate_by_auto("auto", NULL, settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("001001001", grasp_output_stream_filename_generate_by_auto("auto", "3", settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_auto("xxx", NULL, settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_auto("auto", "a", settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_auto("auto", NULL, settings, segment, output, &tile_description.dimensions, -1,1,1)); 
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_auto("auto", NULL, settings, segment, output, &tile_description.dimensions, 1,-1,1)); 
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_auto("auto", NULL, settings, segment, output, &tile_description.dimensions, 1,1,-1)); 
    
    CU_ASSERT_STRING_EQUAL("1", grasp_output_stream_filename_generate_by_icol("icol", NULL, settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("001", grasp_output_stream_filename_generate_by_icol("icol", "3", settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_icol("xxx", NULL, settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_icol("icol", "a", settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_icol("icol", NULL, settings, segment, output, &tile_description.dimensions, -1,1,1)); 
    
    CU_ASSERT_STRING_EQUAL("1", grasp_output_stream_filename_generate_by_irow("irow", NULL, settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("001", grasp_output_stream_filename_generate_by_irow("irow", "3", settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_irow("xxx", NULL, settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_irow("irow", "a", settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_irow("irow", NULL, settings, segment, output, &tile_description.dimensions, 1,-1,1)); 
    
    CU_ASSERT_STRING_EQUAL("1", grasp_output_stream_filename_generate_by_itime("itime", NULL, settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("001", grasp_output_stream_filename_generate_by_itime("itime", "3", settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_itime("xxx", NULL, settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_itime("itime", "a", settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_itime("itime", NULL, settings, segment, output, &tile_description.dimensions, 1,1,-1)); 
    
    CU_ASSERT_STRING_EQUAL("2", grasp_output_stream_filename_generate_by_segment_nx("segment_nx", NULL, settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("002", grasp_output_stream_filename_generate_by_segment_nx("segment_nx", "3", settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_segment_nx("xxx", NULL, settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_segment_nx("segment_nx", "a", settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_segment_nx("segment_nx", NULL, NULL, segment, output, &tile_description.dimensions, 1,1,1)); 

    CU_ASSERT_STRING_EQUAL("2", grasp_output_stream_filename_generate_by_segment_ny("segment_ny", NULL, settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("002", grasp_output_stream_filename_generate_by_segment_ny("segment_ny", "3", settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_segment_ny("xxx", NULL, settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_segment_ny("segment_ny", "a", settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_segment_ny("segment_ny", NULL, NULL, segment, output, &tile_description.dimensions, 1,1,1)); 
    
    CU_ASSERT_STRING_EQUAL("2", grasp_output_stream_filename_generate_by_segment_nt("segment_nt", NULL, settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("002", grasp_output_stream_filename_generate_by_segment_nt("segment_nt", "3", settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_segment_nt("xxx", NULL, settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_segment_nt("segment_nt", "a", settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_segment_nt("segment_nt", NULL, NULL, segment, output, &tile_description.dimensions, 1,1,1)); 
    
    CU_ASSERT_STRING_EQUAL("2015-03-23T05:00:00Z", grasp_output_stream_filename_generate_by_tile_from("tile_from", NULL, settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("2015-03-23", grasp_output_stream_filename_generate_by_tile_from("tile_from", "%Y-%m-%d", settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_tile_from("xxx", NULL, settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_tile_from("tile_from", NULL, NULL, segment, output, &tile_description.dimensions, 1,1,1)); 
    
    CU_ASSERT_STRING_EQUAL("2015-03-23T06:00:00Z", grasp_output_stream_filename_generate_by_tile_to("tile_to", NULL, settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("2015-03-23", grasp_output_stream_filename_generate_by_tile_to("tile_to", "%Y-%m-%d", settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_tile_to("xxx", NULL, settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_tile_to("tile_to", NULL, NULL, segment, output, &tile_description.dimensions, 1,1,1)); 
    
    settings->input.coordinates_grid.col=3286;
    CU_ASSERT_STRING_EQUAL("3286", grasp_output_stream_filename_generate_by_tile_corner_column("tile_corner_column", NULL, settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("3286", grasp_output_stream_filename_generate_by_tile_corner_column("tile_corner_column", "3", settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("03286", grasp_output_stream_filename_generate_by_tile_corner_column("tile_corner_column", "5", settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("03286", grasp_output_stream_filename_generate_by_tile_coordinate_x("tile_coordinate_x", "5", settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_tile_coordinate_x("xxx", NULL, settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_tile_coordinate_x("tile_coordinate_x", NULL, NULL, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_tile_corner_column("xxx", NULL, settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_tile_corner_column("tile_corner_column", "a", settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_tile_corner_column("tile_corner_column", NULL, NULL, segment, output, &tile_description.dimensions, 1,1,1)); 
    settings->input.grid_offset.col=1;
    CU_ASSERT_STRING_EQUAL("3285", grasp_output_stream_filename_generate_by_tile_corner_column("tile_corner_column", NULL, settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("3285", grasp_output_stream_filename_generate_by_tile_coordinate_x("tile_coordinate_x", NULL, settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    settings->input.grid_offset.col=0;
    settings->input.coordinates_grid.col=-1;
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_tile_corner_column("tile_corner_column", NULL, NULL, segment, output, &tile_description.dimensions, 1,1,1)); 
    
    settings->input.coordinates_grid.row=1376;
    CU_ASSERT_STRING_EQUAL("1376", grasp_output_stream_filename_generate_by_tile_corner_row("tile_corner_row", NULL, settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("01376", grasp_output_stream_filename_generate_by_tile_corner_row("tile_corner_row", "5", settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("01376", grasp_output_stream_filename_generate_by_tile_coordinate_y("tile_coordinate_y", "5", settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_tile_coordinate_y("tile_coordinate_y", NULL, NULL, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_tile_corner_row("xxx", NULL, settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_tile_corner_row("tile_corner_row", "a", settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_tile_corner_row("tile_corner_row", NULL, NULL, segment, output, &tile_description.dimensions, 1,1,1)); 
    settings->input.grid_offset.row=1;
    CU_ASSERT_STRING_EQUAL("1375", grasp_output_stream_filename_generate_by_tile_corner_row("tile_corner_row", NULL, settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    settings->input.grid_offset.row=0;
    settings->input.coordinates_grid.row=-1;
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_tile_corner_row("tile_corner_row", NULL, settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    
    settings->input.coordinates.lon=-4.2;
    CU_ASSERT_STRING_EQUAL("-4.200000", grasp_output_stream_filename_generate_by_tile_center_longitude("tile_center_longitude", NULL, settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("-4.200", grasp_output_stream_filename_generate_by_tile_center_longitude("tile_center_longitude", "0.3", settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("-4.200", grasp_output_stream_filename_generate_by_tile_coordinate_x("tile_coordinate_x", "3", settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_tile_coordinate_x("tile_coordinate_", "3", settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_tile_coordinate_x("tile_coordinate_x", "3", NULL, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_tile_center_longitude("xxx", NULL, settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_tile_center_longitude("tile_center_longitude", "a", settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_tile_center_longitude("tile_center_longitude", NULL, NULL, segment, output, &tile_description.dimensions, 1,1,1)); 
    settings->input.coordinates.lon=0.0;
    settings->input.coordinates_grid.col=3;
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_tile_center_longitude("tile_center_longitude", NULL, NULL, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("03", grasp_output_stream_filename_generate_by_tile_coordinate_x("tile_coordinate_x", "2", settings, segment, output, &tile_description.dimensions, 1,1,1));     
    settings->input.coordinates_grid.col=-1;
    
    settings->input.coordinates.lat=40.5;
    CU_ASSERT_STRING_EQUAL("40.500000", grasp_output_stream_filename_generate_by_tile_center_latitude("tile_center_latitude", NULL, settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("40.500", grasp_output_stream_filename_generate_by_tile_center_latitude("tile_center_latitude", "0.3", settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("40.500", grasp_output_stream_filename_generate_by_tile_coordinate_y("tile_coordinate_y", "3", settings, segment, output, &tile_description.dimensions, 1,1,1));     
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_tile_coordinate_y("tile_coordinate_", "3", settings, segment, output, &tile_description.dimensions, 1,1,1));     
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_tile_coordinate_y("tile_coordinate_", "3", NULL, segment, output, &tile_description.dimensions, 1,1,1));     
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_tile_center_latitude("xxx", NULL, settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_tile_center_latitude("tile_center_latitude", "a", settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_tile_center_latitude("tile_center_latitude", NULL, NULL, segment, output, &tile_description.dimensions, 1,1,1)); 
    settings->input.coordinates.lat=0.0;
    settings->input.coordinates_grid.row=3;
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_tile_center_latitude("tile_center_latitude", NULL, NULL, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("03", grasp_output_stream_filename_generate_by_tile_coordinate_y("tile_coordinate_y", "2", settings, segment, output, &tile_description.dimensions, 1,1,1));         
    settings->input.coordinates_grid.row=-1;   
    
    
    CU_ASSERT_STRING_EQUAL("2", grasp_output_stream_filename_generate_by_tile_width("tile_width", NULL, settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("00002", grasp_output_stream_filename_generate_by_tile_width("tile_width", "5", settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_tile_width("xxx", NULL, settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_tile_width("tile_width", "a", settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_tile_width("tile_width", NULL, NULL, segment, output, &tile_description.dimensions, 1,1,1)); 
    settings->input.area_width=-1;
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_tile_width("tile_width", NULL, settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    settings->input.area_width=2;
    
    CU_ASSERT_STRING_EQUAL("2", grasp_output_stream_filename_generate_by_tile_height("tile_height", NULL, settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("00002", grasp_output_stream_filename_generate_by_tile_height("tile_height", "5", settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_tile_height("xxx", NULL, settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_tile_height("tile_height", "a", settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_tile_height("tile_height", NULL, NULL, segment, output, &tile_description.dimensions, 1,1,1)); 
    settings->input.area_height=-1;
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_tile_height("tile_height", NULL, settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    settings->input.area_height=2;
            
    settings->input.coordinates_grid.col=3286;
    CU_ASSERT_STRING_EQUAL("3288", grasp_output_stream_filename_generate_by_segment_corner_column("segment_corner_column", NULL, settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("03288", grasp_output_stream_filename_generate_by_segment_corner_column("segment_corner_column", "5", settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    settings->input.grid_offset.col=1;
    CU_ASSERT_STRING_EQUAL("3287", grasp_output_stream_filename_generate_by_segment_corner_column("segment_corner_column", NULL, settings, segment, output, &tile_description.dimensions, 1,1,1));     
    settings->input.grid_offset.col=0;
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_segment_corner_column("xxx", NULL, settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_segment_corner_column("segment_corner_column", "a", settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_segment_corner_column("segment_corner_column", NULL, NULL, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_segment_corner_column("segment_corner_column", NULL, NULL, segment, output, &tile_description.dimensions, -1,1,1));     
    settings->input.coordinates_grid.col=-1;
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_segment_corner_column("segment_corner_column", NULL, settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    
    settings->input.coordinates_grid.row=1376;
    CU_ASSERT_STRING_EQUAL("1378", grasp_output_stream_filename_generate_by_segment_corner_row("segment_corner_row", NULL, settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("01378", grasp_output_stream_filename_generate_by_segment_corner_row("segment_corner_row", "5", settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    settings->input.grid_offset.row=1;
    CU_ASSERT_STRING_EQUAL("1377", grasp_output_stream_filename_generate_by_segment_corner_row("segment_corner_row", NULL, settings, segment, output, &tile_description.dimensions, 1,1,1));     
    settings->input.grid_offset.row=0;
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_segment_corner_row("xxx", NULL, settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_segment_corner_row("segment_corner_row", "a", settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_segment_corner_row("segment_corner_row", NULL, NULL, segment, output, &tile_description.dimensions, 1,1,1)); 
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_segment_corner_row("segment_corner_row", NULL, NULL, segment, output, &tile_description.dimensions, -1,1,1));     
    settings->input.coordinates_grid.row=-1;
    CU_ASSERT_STRING_EQUAL("", grasp_output_stream_filename_generate_by_segment_corner_row("segment_corner_row", NULL, settings, segment, output, &tile_description.dimensions, 1,1,1)); 
    
    CU_ASSERT(NULL!=strstr(grasp_output_stream_filename_generate_by_pwd("pwd", NULL, settings, segment, output, &tile_description.dimensions, 1,1,1), "/src/test/"));
    CU_ASSERT_STRING_EQUAL("/src/test/",strstr(grasp_output_stream_filename_generate_by_pwd("pwd", NULL, settings, segment, output, &tile_description.dimensions, 1,1,1), "/src/test/"));
    CU_ASSERT_STRING_EQUAL("",grasp_output_stream_filename_generate_by_pwd("a", NULL, settings, segment, output, &tile_description.dimensions, 1,1,1));
    
    CU_ASSERT_STRING_EQUAL("sdata/",grasp_output_stream_filename_generate_by_yml("yml", NULL, settings, segment, output, &tile_description.dimensions, 1,1,1));
    CU_ASSERT_STRING_EQUAL("",grasp_output_stream_filename_generate_by_pwd("a", NULL, settings, segment, output, &tile_description.dimensions, 1,1,1));
    
    CU_ASSERT_STRING_EQUAL("david",grasp_output_stream_filename_generate_by_token("david", NULL, settings, segment, output, &tile_description.dimensions, 1,1,1));
    CU_ASSERT_STRING_EQUAL("david(3)",grasp_output_stream_filename_generate_by_token("david(3)", NULL, settings, segment, output, &tile_description.dimensions, 1,1,1));
    CU_ASSERT_STRING_EQUAL("david(3)",grasp_output_stream_filename_generate_by_token("david", "3", settings, segment, output, &tile_description.dimensions, 1,1,1));    
    CU_ASSERT_STRING_EQUAL("",grasp_output_stream_filename_generate_by_token("", NULL, settings, segment, output, &tile_description.dimensions, 1,1,1));
    
    
    // 5) Clean up 
    grasp_controller_clean_memory(settings, &tile_description, &results, &functions);
    free(segment);
    free(output);
    
}

int main(int argc, char *argv[]) {
    CU_pSuite grasp_Suite = NULL;

    /* Initialize the CUnit test registry */
    if (CUE_SUCCESS != CU_initialize_registry())
        return CU_get_error();

    /* Add a suite to the registry */
    grasp_Suite = CU_add_suite("grasp_test_output", init_suite, clean_suite);
    if (NULL == grasp_Suite) {
        CU_cleanup_registry();
        return CU_get_error();
    }
   if (
       (NULL == CU_add_test(grasp_Suite, "grasp_test_output_stream", grasp_test_output_stream))
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
