/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "grasp_controller.h"
#include "grasp_mpi_engine.h"
#include "yamlsettings/yamlsettings.h"
#include "../input/grasp_input.h"
#include "../output/grasp_output.h"
#include "../settings/grasp_settings.h"
#include <grasp/utils.h>
#include "grasp_controller.h"
#include "mo_grasp_controller.h"
#include "../input/grasp_input_load_functions.h"
#include "../output/grasp_output_load_function.h"
#include "../settings/grasp_settings_data_types.h"
#include "../output/grasp_output_stream.h"
#include "../global/grasp_compilation_information.h"
#include "../global/grasp_runtime_information.h"

grasp_output_stream controller_stream;
grasp_output_stream controller_trackmem_stream;

// Real time in seconds used by retrieval algorithm
static double algorithm_ut;
// User time in seconds used by retrieval algorithm
static double algorithm_ct;
// Segment problem statistics
static int nerror_segment;
static int nerror_pixel;


grasp_output_stream *grasp_controller_get_stream(){
    return &controller_stream;
}

double grasp_controller_get_algorithm_ut(){
    return algorithm_ut;
}

double grasp_controller_get_algorithm_ct(){
    return algorithm_ct;
}

int grasp_controller_get_nerror_segment(){
    return nerror_segment;
}

int grasp_controller_get_nerror_pixel(){
    return nerror_pixel;
}


grasp_settings *grasp_controller_read_settings(int argc, char** argv) {
    yamlsettings_dictionary_t *dictionary;
    grasp_settings *settings;
    
    // Initialize runtime otpions
    grasp_runtime_initialize(argv);

    if (argc < 2) {
        strcpy(grasp_exec_file, argv[0]);
        grasp_compilation_information_print(stdout);
        exit(0);
    }

    // Read config file
    grasp_settings_read(&controller_stream, &dictionary, argc - 1, (const char **)&argv[1]);

    // Initialize runtime information
    grasp_runtime_set(argv[0], dictionary->files[YAMLSETTINGS_FILE_COMMAND_LINE], dictionary->files[YAMLSETTINGS_FILE_MAIN_CONFIGURATION]);

    // Controller will check if it should do something
    grasp_controller_process_options((grasp_settings *) dictionary->settings);

    // Prepare data for return
    settings = (grasp_settings *) dictionary->settings;

    gos_fprintf_flushed(&controller_stream, "Config file read successfully\n");

    yamlsettings_dictionary_destroy(dictionary);

    return settings;
}

int grasp_controller_initialize_inversion(grasp_settings *settings, grasp_tile_description_t *tile_description, grasp_processing_functions_t *functions, grasp_results_t *results) {
    int status;

    status = grasp_input_initialize_tile_description(settings, tile_description);
    assert(status == 0);
    status = grasp_controller_initialize_functions(settings, tile_description, functions);
    assert(status == 0);
    status = grasp_output_initialize_results(settings, tile_description, results);
    assert(status == 0);

    // Initilize times take by the framework
    algorithm_ut=0.0;
    algorithm_ct=0.0;
    nerror_segment=0;
    nerror_pixel=0;
    
    return status;
}

int grasp_controller_get_next_inversion(grasp_tile_dimensions_t *tile_dimensions, int iinversion){
    if(iinversion<0){
        return -1;
    }
    if(iinversion>tile_dimensions->segment_ntimes*tile_dimensions->segment_nrows*tile_dimensions->segment_ncols){
        return -1;
    }
    
    return iinversion-1;
    /*int indexes[3];
    int maximums[3];
    
    maximums[0]=tile_dimensions->segment_ntimes;
    maximums[1]=tile_dimensions->segment_nrows;
    maximums[2]=tile_dimensions->segment_ncols;
             
    if(iinversion > maximums[0]*maximums[1]*maximums[2]){ // If the index is out bound it return -1
        return -1;
    }
    
    indexesofND(iinversion-1, indexes, 3,  maximums);
    
    *itime=indexes[0];
    *irow=indexes[1];
    *icol=indexes[2];
    
    return 0;*/
}

void grasp_controller_invert_tile_sequential(grasp_settings *settings, grasp_tile_description_t *tile_description, grasp_results_t *results, grasp_processing_functions_t *functions) {
    int iinversion;
    grasp_segment_t *segment;
    output_segment_general *output;
    int id_inversion;
    
    segment = (grasp_segment_t *) trackmem_malloc(sizeof (grasp_segment_t)*1);
    assert(segment!=NULL);
    output = (output_segment_general *) trackmem_malloc(sizeof (output_segment_general)*1);
    assert(output!=NULL);       
    
    grasp_init_inversion(&settings->retrieval);

    iinversion = 1;
    // Iterate over segments
    while((id_inversion=grasp_controller_get_next_inversion(&tile_description->dimensions, iinversion))>=0){
        grasp_controller_call_inversion(settings, segment, output, tile_description, results, functions, iinversion, id_inversion);
        iinversion++; 
    }


    grasp_finalize_inversion(&settings->retrieval);
    
    trackmem_free(segment);
    trackmem_free(output);
}
//#define USE_MPI
#ifdef USE_MPI

typedef struct my_master_info_t_ {
    grasp_settings *settings;
    grasp_tile_description_t *tile_description;
    grasp_processing_functions_t *functions;
    grasp_results_t *results;
} my_master_info_t;

typedef struct my_worker_info_t_ {
    grasp_settings *settings;
} my_worker_info_t;

typedef struct my_order_payload_t_ {
    grasp_segment_t input_segment;
    grasp_tile_dimensions_t tile_dimensions;
    int id_inversion;
    int ninversion;
    int npixel;
} my_order_payload_t;

typedef struct my_result_payload_t_ {
    grasp_segment_t input_segment;
    output_segment_general output_segment;
    grasp_tile_dimensions_t tile_dimensions;
    int id_inversion;
    int ninversion;
    int npixel;
    float algorithm_ut;
    float algorithm_ct;
    int status;
} my_result_payload_t;

order_t *grasp_controller_order_callback(int worker_rank, const void *master_info, bool *no_more_order) {
    order_t *order;
    char label[255 + 1];
    int id_inversion;
    my_master_info_t *my_master_info;
    grasp_tile_dimensions_t *tile_dimensions;
    const grasp_settings *settings;
    const grasp_processing_functions_t *functions;
    static int ninversion = 1;
    static bool settings_are_already_set=false;
    
    assert(master_info != NULL);
    assert(no_more_order != NULL);

    my_master_info = (my_master_info_t *) master_info;
    tile_dimensions = &my_master_info->tile_description->dimensions;
    assert(tile_dimensions != NULL);
    settings = my_master_info->settings;
    assert(settings != NULL);
    functions = my_master_info->functions;
    assert(functions != NULL);

    /* obtain the next segment to process */
    id_inversion = grasp_controller_get_next_inversion((grasp_tile_dimensions_t *) tile_dimensions, ninversion);

    if (id_inversion >= 0) { /* there is a next segment */
        my_order_payload_t *my_order_payload;
        int npixel;

        *no_more_order = false;

#ifdef DEBUG
	fprintf(stderr, "%s:%d: next inversion #%d will be performed ( segment id_inversion: %d  )\n", __FILE__, __LINE__, ninversion, id_inversion);
#endif

        my_order_payload = trackmem_malloc(sizeof (my_order_payload_t));
        assert(my_order_payload!=NULL);
        assert(my_order_payload != NULL);

        my_order_payload->id_inversion = id_inversion;
        my_order_payload->ninversion = ninversion;
        memcpy(&my_order_payload->tile_dimensions, tile_dimensions, sizeof(grasp_tile_dimensions_t));

        npixel = grasp_input_extract_segment((grasp_settings *) settings, (grasp_input_driver_t *) &(functions->driver), 
                 functions->ntransformers, (grasp_input_transformer_t *) functions->transformers, &my_order_payload->input_segment, my_master_info->results, tile_dimensions, id_inversion);
        
        // Trick: In master process we need a settings with all the values. For obtaining that, we need to call grasp_complete_input_settings_and_segment_data
        // function in order to fill the missing information. We call this function with a random segment (we are using first segment with data) because we are
        // not interested in segment data, only we need it for calling the function. 
        if(!settings_are_already_set && npixel>0){
            assert(&(my_order_payload->input_segment)!=NULL);
            grasp_prepare_segment_settings((retr_input *) &settings->retrieval, (sensor_data_t *) &(my_order_payload->input_segment.sdata));
            settings_are_already_set=true;
        }
        
#ifdef DEBUG
	fprintf(stderr, "%s:%d: segment ( ninversion: %d id_inversion: %d npixel: %d ) extracted\n", __FILE__, __LINE__, ninversion, id_inversion, npixel);
#endif

        my_order_payload->npixel = npixel;

        if (npixel > 0) {
	  sprintf(label, "order for worker #%d: process segment ( ninversion: %d id_inversion: %d npixel: %d )", worker_rank, ninversion, id_inversion, npixel);
            order = grasp_mpi_engine_new_order(label, sizeof (my_order_payload_t), my_order_payload);
            assert(order != NULL);
        } 
	else {
#ifdef DEBUG_MPI
	  fprintf(stderr, "%s:%d: Master: no valid pixel found in the segment ( ninversion: %d id_inversion: %d ), it will be discarded\n",
		  __FILE__, __LINE__, ninversion, id_inversion);
#endif
	  order = NULL;
        }
	
        trackmem_free(my_order_payload);

	ninversion++;
        return order;
    } 
    else {
#ifdef DEBUG_MPI
      fprintf(stderr, "%s:%d: Master: no more segment to process\n", __FILE__, __LINE__);
#endif
      
      /* this is for checking that workers are correctly stopped at the end of the processing */
#ifdef DEBUG_MPI_WORK_DONE
      grasp_mpi_engine_set_debug_level(2);
#endif
      
      *no_more_order = true;
    }

    return NULL;
}

result_t *grasp_controller_task_callback(const order_t *order, const void *worker_info) {
    result_t *result;
    my_worker_info_t *my_worker_info;
    const grasp_settings *settings;
    my_order_payload_t *my_order_payload = NULL;
    my_result_payload_t *my_result_payload = NULL; /* the size of the result payload depends on the application and may be too large for stack allocation */
    output_segment_general *output_segment = NULL;
    char order_label[255 + 1];
    char result_label[255 + 1];
    benchmark_t benchmark;
    
    assert(order != NULL);
    assert(worker_info != NULL);
    my_worker_info = (my_worker_info_t *) worker_info;
    settings = my_worker_info->settings;
    assert(settings != NULL);

    if (settings->controller.perform_retrieval == true) {
        /* controller.debug.perform_retrieval=true (by default) */
        grasp_segment_t *input_segment;
        const grasp_tile_dimensions_t *tile_dimensions;
        int id_inversion, ninversion;
        int status;

        my_order_payload = grasp_mpi_engine_get_order_payload(order);
        assert(my_order_payload != NULL);
        input_segment = &my_order_payload->input_segment;
        
        id_inversion = my_order_payload->id_inversion;
        ninversion = my_order_payload->ninversion;
        tile_dimensions = &my_order_payload->tile_dimensions;

        benchmark_start(&benchmark);
        output_segment = trackmem_malloc(sizeof(output_segment_general));
        assert(output_segment != NULL);
        status=grasp_controller_processor_unit(settings, input_segment, output_segment, 
					 tile_dimensions, ninversion);
        benchmark_stop(&benchmark); 
          
	my_result_payload = trackmem_malloc(sizeof (my_result_payload_t));
        assert(my_result_payload!=NULL);
        assert(my_result_payload != NULL);
	my_result_payload->id_inversion = id_inversion;
	my_result_payload->ninversion = ninversion;
	my_result_payload->npixel = my_order_payload->npixel;
        my_result_payload->algorithm_ut=benchmark.delta_ut;
        my_result_payload->algorithm_ct=benchmark.delta_ct;
        my_result_payload->status = status;
	memcpy(&my_result_payload->tile_dimensions, &my_order_payload->tile_dimensions, sizeof (grasp_tile_dimensions_t));
        memcpy(&my_result_payload->input_segment, &my_order_payload->input_segment, sizeof (grasp_segment_t));
	memcpy(&my_result_payload->output_segment, output_segment, sizeof (output_segment_general));
	trackmem_free(output_segment);
	
    } else {
        fprintf(stderr, "%s:%d: inversion disabled (controller.debug.perform_retrieval set to false)\n", __FILE__, __LINE__);
    }

    grasp_mpi_engine_get_order_label(order, sizeof (order_label) - 1, order_label);
    snprintf(result_label, sizeof (result_label) - 1, "%s (completed)", order_label);
    result = grasp_mpi_engine_new_result(result_label, sizeof (my_result_payload_t), my_result_payload);
    assert(result != NULL);
    trackmem_free(my_result_payload);
    return result;
}

void grasp_controller_collect_callback(const result_t *result, const void *master_info) {
    my_master_info_t *my_master_info;
    const grasp_tile_description_t *tile_description;
    const grasp_settings *settings;
    grasp_segment_t *input_segment;
    const output_segment_general *output_segment;
    my_result_payload_t *my_result_payload;
    const grasp_results_t *results;
    const grasp_processing_functions_t *functions;
    int npixel;

    assert(master_info != NULL);
    my_master_info = (my_master_info_t *) master_info;
    tile_description = my_master_info->tile_description;
    assert(tile_description != NULL);
    settings = my_master_info->settings;
    assert(settings != NULL);
    results = my_master_info->results;
    assert(results != NULL);
    functions = my_master_info->functions;
    assert(functions != NULL);

    my_result_payload = grasp_mpi_engine_get_result_payload(result);
    assert(my_result_payload != NULL);

    input_segment = &my_result_payload->input_segment;
    output_segment = &my_result_payload->output_segment;
    npixel = my_result_payload->npixel;
    algorithm_ut+=my_result_payload->algorithm_ut;
    algorithm_ct+=my_result_payload->algorithm_ct;

    if (npixel == 0) {
#ifdef DEBUG
      fprintf(stderr, "%s:%d: %s: empty result found (no valid data or interrupted task)\n", __FILE__, __LINE__, __func__);
#endif
      return;
    }

#ifdef DEBUG
    fprintf(stderr, "%s:%d: %s: results from segment ( ninversion: %d id_inversion: %d npixel: %d ) collected\n", 
	    __FILE__, __LINE__, __func__, my_result_payload->ninversion, my_result_payload->id_inversion, npixel);
#endif

    if(my_result_payload->status==-2){
        abort();
    }
    if(my_result_payload->status==-1){
        nerror_segment++;
        nerror_pixel+=my_result_payload->input_segment.sdata.npixels;
        my_result_payload->input_segment.sdata.npixels=0; // Ignoring what we have tried to retrieve because it return error
    }
    
    grasp_controller_post_process_segment((grasp_settings *) settings, (grasp_segment_t *) input_segment, (output_segment_general *) output_segment, 
					  (grasp_tile_description_t *) tile_description, (grasp_results_t *) results, (grasp_processing_functions_t *) functions 
                                          );
}

void grasp_controller_progress_info_callback(int num_orders_performed, const void *master_info) {
    my_master_info_t *my_master_info;
    const grasp_tile_dimensions_t *tile_dimensions;
    double progress_ratio;
    time_t timestamp;
    struct tm *timeinfo;
    char cstr_time[255 + 1];
    int num_orders_expected;
    int ncols, nrows, ntimes;

    assert(num_orders_performed >= 0);
    assert(master_info != NULL);
    my_master_info = (my_master_info_t *) master_info;
    tile_dimensions = &my_master_info->tile_description->dimensions;
    assert(tile_dimensions != NULL);
    ncols = tile_dimensions->segment_ncols;
    nrows = tile_dimensions->segment_nrows;
    ntimes = tile_dimensions->segment_ntimes;

    num_orders_expected = ncols * nrows * ntimes;
    if (num_orders_expected == 0) {
#ifdef DEBUG
      gos_fprintf_flushed(&controller_stream, "%s:%d: ", __FILE__, __LINE__);
#endif
      gos_fprintf_flushed(&controller_stream, "Master: the tile is empty, nothing to do");
      return;
    }

    time(&timestamp);
    timeinfo = localtime(&timestamp);
    asctime_r(timeinfo, cstr_time); /* contains newline */

    progress_ratio = 100. * num_orders_performed / num_orders_expected;

#ifdef DEBUG
    gos_fprintf_flushed(&controller_stream, "%s:%d: ", __FILE__, __LINE__);
#endif
    /* no newline character to add, it is included already in the output of asctime[_r] */
    gos_fprintf_flushed(&controller_stream, "Master: progression status %.2f%% (%d/%d) at %s",
            progress_ratio, num_orders_performed, num_orders_expected, cstr_time);
}

void grasp_controller_invert_tile_mpi(grasp_settings *settings, grasp_tile_description_t *tile_description, grasp_results_t *results, grasp_processing_functions_t *functions) {
    my_master_info_t my_master_info;
    my_worker_info_t my_worker_info;
    
    assert(settings != NULL);
    if (grasp_mpi_engine_is_master()) {
      assert(tile_description != NULL);
      assert(results != NULL);
      assert(functions != NULL);
    }
    else {
      assert(tile_description == NULL);
      assert(results == NULL);
      assert(functions == NULL);
    }

    grasp_mpi_engine_set_appname("grasp");
    my_master_info.tile_description = tile_description;
    my_master_info.settings = settings;
    my_master_info.results = results;
    my_master_info.functions = functions;
    my_worker_info.settings = settings;

    grasp_mpi_engine_set_polling_time(settings->controller.polling_time);
    grasp_mpi_engine_set_maximum_job_time(settings->controller.maximum_job_time);

    /* perform the actual processing */
    grasp_mpi_engine_set_order_callback(grasp_controller_order_callback, sizeof (my_order_payload_t));
    grasp_mpi_engine_set_task_callback(grasp_controller_task_callback, sizeof (my_result_payload_t));
    grasp_mpi_engine_set_collect_callback(grasp_controller_collect_callback);
    grasp_mpi_engine_set_progress_info_callback(grasp_controller_progress_info_callback);

    if (grasp_mpi_engine_is_worker()) {
        grasp_init_inversion(&settings->retrieval);
    }

    grasp_mpi_engine_main_loop(&my_master_info, &my_worker_info);

    if (grasp_mpi_engine_is_worker()) {
        grasp_finalize_inversion(&settings->retrieval);
    }
}
#endif /* USE_MPI */

void grasp_controller_invert_tile(grasp_settings *settings, grasp_tile_description_t *tile_description, grasp_results_t *results, grasp_processing_functions_t *functions) {

#ifdef USE_MPI
    if (grasp_mpi_engine_is_master()) {
        gos_fprintf_flushed(&controller_stream, "The tile is divided in segments with %d rows, %d cols and %d times. %d inversions will be performed (mpi version)\n", tile_description->dimensions.segment_nrows, tile_description->dimensions.segment_ncols, tile_description->dimensions.segment_ntimes, tile_description->dimensions.segment_nrows * tile_description->dimensions.segment_ncols * tile_description->dimensions.segment_ntimes);
    }
    grasp_controller_invert_tile_mpi(settings, tile_description, results, functions);
#else
    gos_fprintf_flushed(&controller_stream, "The tile is divided in segments with %d rows, %d cols and %d times. %d inversions will be performed (sequential version)\n", tile_description->dimensions.segment_nrows, tile_description->dimensions.segment_ncols, tile_description->dimensions.segment_ntimes, tile_description->dimensions.segment_nrows * tile_description->dimensions.segment_ncols * tile_description->dimensions.segment_ntimes);
    grasp_controller_invert_tile_sequential(settings, tile_description, results, functions);
#endif
}

bool grasp_controller_segment_is_invertible(int id_inversion, const grasp_settings *settings, const grasp_tile_dimensions_t *tile_description){
    int min, max;
    
    min=0;
    max=tile_description->segment_ntimes * tile_description->segment_ncols * tile_description->segment_nrows;
    max=max-1; //because last check of function will check with "lower or equal than"
    
    if(settings->controller.segment_range[0]>=0){ // If minimum is defined
        min=settings->controller.segment_range[0];
    }
    
    if(settings->controller.segment_range[1]>=0){ // If maximum is defined
        max=settings->controller.segment_range[1];
    }
    
    if(settings->controller.nsegment_range==1){ // If only a number if defined we set maximum like minimum so only one retrieval will be performed
        max=settings->controller.segment_range[0];
    }
    
    if(id_inversion>=min && id_inversion<=max){
        return true;
    }else{
        return false;
    }
}

int grasp_controller_processor_unit(const grasp_settings *settings, const grasp_segment_t *segment, output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int iinversion) {
    grasp_output_stream grasp_stream;
    benchmark_t benchmark;
    int total_inversions;
    FILE *f;
    int icol, irow,itime;
    int inversion_result;
    int result=0;
    
    grasp_input_position_of_inversion(tile_description,segment->sdata.id_inversion, &icol, &irow,&itime);
    
    total_inversions=tile_description->segment_ncols*tile_description->segment_nrows*tile_description->segment_ntimes;
    
    if (grasp_controller_segment_is_invertible(segment->sdata.id_inversion, settings, tile_description)==true){ //settings->controller.do_retrieval_n < 0 || settings->controller.do_retrieval_n == ninversion) {
        gos_fprintf_flushed(&controller_stream, "Retrieval #%d (%d/%d): %.2f%%: ", 
                segment->sdata.id_inversion,
                iinversion, 
                total_inversions,
                (iinversion*100.0)/(total_inversions));
        if (segment->sdata.npixels < 0) {
            gos_fprintf_flushed(&controller_stream, "there is an error retrieving pixels and this inversion will be skipped\n");
        } else {
            gos_fprintf_flushed(&controller_stream, "%d pixels will be processed\n", segment->sdata.npixels);

            if (segment->sdata.npixels > 0) {
                // Dumping SDATA file
                grasp_output_stream_initialize(settings->input.sdata_stream, &grasp_stream);
                f = grasp_output_stream_open(&grasp_stream, settings, segment, output, tile_description, icol, irow, itime);
                if (grasp_output_stream_writable(&grasp_stream)) {
                    grasp_input_dump_segment(f, segment);
                }
                grasp_output_stream_close(&grasp_stream);
                grasp_output_stream_destroy(&grasp_stream);
                // Dumping imagedat file
                grasp_output_stream_initialize(settings->input.imagedat_stream, &grasp_stream);
                f = grasp_output_stream_open(&grasp_stream, settings, segment, output, tile_description, icol, irow, itime);
                if (grasp_output_stream_writable(&grasp_stream)) {
                    grasp_input_dump_iguess(f, settings, segment);
                }
                grasp_output_stream_close(&grasp_stream);
                grasp_output_stream_destroy(&grasp_stream);
                    
                if (settings->controller.perform_retrieval == true) {
                    // Call inversion
                    benchmark_start(&benchmark);

                    inversion_result=grasp_input_inversion(&settings->retrieval, &segment->sdata, segment->iguess, &segment->edges, output);

                    benchmark_stop(&benchmark);
                    algorithm_ut+=benchmark.delta_ut;
                    algorithm_ct+=benchmark.delta_ct;
                    
                    if(inversion_result==-1){
                        gos_fprintf_flushed(&controller_stream, "Retrieval #%d (%d/%d): Segment error. It is set as empty segment (0 pixels) and it is ignored. Retrieval code keep processing the rest of segments\n", segment->sdata.id_inversion, iinversion,total_inversions);
                        result=-1;
                    }else if(inversion_result==-2){
                        fprintf(stderr,"FATAL ERROR: Retrieval code can not continue and it is stopped.\n");
                        result=-2;
                    }else{
                        gos_fprintf_flushed(&controller_stream, "Retrieval #%d (%d/%d): %.2f%%: %d pixels processed in %f seconds (cpu time: %f). Average per pixel: %f (cpu time: %f)\n", segment->sdata.id_inversion, iinversion,total_inversions, (iinversion*100.0)/(total_inversions), segment->sdata.npixels, benchmark.delta_ct, benchmark.delta_ut, benchmark.delta_ct / segment->sdata.npixels, benchmark.delta_ut / segment->sdata.npixels);
                    }
                }


            }
        }
        if (result!=-2){
            gos_fprintf_flushed(&controller_stream, "Retrieval #%d (%d/%d): %.2f%%: finished\n", segment->sdata.id_inversion, iinversion,total_inversions, (iinversion*100.0)/(total_inversions));
        }
    }
    
    return result;
}

void grasp_controller_post_process_segment(grasp_settings *settings, grasp_segment_t *segment, output_segment_general *output, grasp_tile_description_t *tile_description, grasp_results_t *results, grasp_processing_functions_t *functions) {
    int i;
    grasp_output_stream stream;
    int icol, irow,itime;
    
    grasp_input_position_of_inversion(&tile_description->dimensions,segment->sdata.id_inversion, &icol, &irow, &itime);

    if (grasp_controller_segment_is_invertible(segment->sdata.id_inversion, settings, &(tile_description->dimensions))==true && settings->controller.perform_retrieval==true){
        // Processing output
        // Call output segment drivers
        for (i = 0; i < settings->output.nsegment_output_function; i++) {
            grasp_output_stream_initialize(settings->output.segment_stream[i], &stream);
            functions->segment_output_functions[i].function(&stream, settings, segment, output, tile_description, icol, irow, itime);
            grasp_output_stream_destroy(&stream);
        }

        // Allocate and process output for generating tile results.
        grasp_output_process_output(settings, segment, output, tile_description, results, icol, irow, itime);

        // Call output current functions
        for (i = 0; i < settings->output.ncurrent_output_function; i++) {
            grasp_output_stream_initialize(settings->output.current_stream[i], &stream);
            functions->current_output_functions[i].function(&stream, settings, segment, output, results, tile_description, icol, irow, itime);
            grasp_output_stream_destroy(&stream);
        }
    }

}

void grasp_controller_call_inversion(grasp_settings *settings, grasp_segment_t *segment, output_segment_general *output, grasp_tile_description_t *tile_description, grasp_results_t *results, grasp_processing_functions_t *functions, int iinversion, int id_inversion) {
    int npixel;
    (void)npixel; // This avoid -Wunused-but-set-variable warning
    int status;
    
    // Getting pixels
    npixel = grasp_input_extract_segment(settings, &(functions->driver), functions->ntransformers, functions->transformers, segment, results, &tile_description->dimensions, id_inversion);
    assert(npixel>=0); 
    
    // Calling inversions
    status=grasp_controller_processor_unit(settings, segment, output, &tile_description->dimensions, iinversion); //call works from master and obtain the result              
    if(status==-2){
        abort();
    }
    if(status==-1){
        nerror_segment++;
        nerror_pixel+=segment->sdata.npixels;
        segment->sdata.npixels=0; // Ignoring what we have tried to retrieve because it return error
    }
    grasp_controller_post_process_segment(settings, segment, output, tile_description, results, functions);
}

void grasp_controller_process_options(grasp_settings *settings) {

    if (settings->controller.compilation_information == true) {
        grasp_compilation_information_print(stdout);
    }
    
    grasp_output_stream_initialize(settings->controller.track_mem_stream, &controller_trackmem_stream);    
    trackmem_set_debug_stream(grasp_output_stream_open(&controller_trackmem_stream, settings, NULL, NULL, NULL, -1, -1, -1));
}

void grasp_controller_manage_tile(grasp_settings *settings, grasp_tile_description_t *tile_description, grasp_results_t *results, grasp_processing_functions_t *functions) {
    int i;
    grasp_output_stream stream;
    
    // Manage output
    for (i = 0; i < settings->output.ntile_output_function; i++) {
        grasp_output_stream_initialize(settings->output.tile_stream[i], &stream);
        functions->tile_output_functions[i].function(&stream, settings, tile_description, results);
        grasp_output_stream_destroy(&stream);
    }
    
    grasp_output_stream_initialize(settings->input.print_used_files, &stream);
    grasp_output_stream_open(&stream, settings, NULL, NULL, &tile_description->dimensions, -1, -1, -1);
    if (grasp_output_stream_writable(&stream)) {
        gos_fprintf(&stream, "The following files was used\n");
        for (i = 0; i < tile_description->nused_files; i++) {
            gos_fprintf(&stream," %s\n" ,tile_description->used_files[i]);
        }
    }
    grasp_output_stream_close(&stream);
    grasp_output_stream_destroy(&stream);         
}

void grasp_controller_clean_memory(grasp_settings *settings, grasp_tile_description_t *tile_description, grasp_results_t *results, grasp_processing_functions_t *functions) {
    int i;
    int status;
    (void)status; // This avoid -Wunused-but-set-variable warning
#ifdef USE_MPI
    if (grasp_mpi_engine_is_master()) {
        /* in MPI mode, only the master has to release the tile */
#endif

        status = functions->driver.close();

        assert(status == 0);

        grasp_input_deallocate_expanded_files(tile_description->input_files, tile_description->ninput_files);

        if (tile_description->used_files != NULL) {
            for (i = 0; i < tile_description->nused_files; i++) {
                trackmem_free(tile_description->used_files[i]);
            }
            trackmem_free(tile_description->used_files);
        }

        grasp_output_destroy_result(tile_description, results);

        for (i = 0; i < functions->nsegment_output_functions; i++) {
            status = functions->segment_output_functions[i].close();
            assert(status == 0);
        }
        for (i = 0; i < functions->ntile_output_functions; i++) {
            status = functions->tile_output_functions[i].close();
            assert(status == 0);
        }
        for (i = 0; i < functions->ncurrent_output_functions; i++) {
            status = functions->current_output_functions[i].close();
            assert(status == 0);
        }        

        for (i = 0; i < functions->ntransformers; i++) {
            status = functions->transformers[i].close();
            assert(status == 0);
        }
        trackmem_free(functions->transformers);

        trackmem_free(functions->tile_output_functions);
        trackmem_free(functions->segment_output_functions);
        trackmem_free(functions->current_output_functions);

#ifdef USE_MPI
    } /* if (grasp_mpi_engine_is_master()) */
#endif

    trackmem_free(settings);
    
    grasp_output_stream_close(&controller_trackmem_stream);
    grasp_output_stream_destroy(&controller_trackmem_stream);    
}

int grasp_controller_initialize_functions(grasp_settings *settings, grasp_tile_description_t *tile_description, grasp_processing_functions_t *functions) {
    int i;
    int status=0;

    // Get driver
    functions->driver = grasp_input_driver_get_function(settings->input.driver_name);
    // Initialize 
    status += functions->driver.init(settings, tile_description);

    // Prepare transformers
    functions->ntransformers = settings->input.ntransformers;
    functions->transformers = (grasp_input_transformer_t *) trackmem_malloc(sizeof (grasp_input_transformer_t) * functions->ntransformers);
    assert(functions->transformers!=NULL);
    for (i = 0; i < functions->ntransformers; i++) {
        functions->transformers[i] = grasp_input_transformer_get_function(settings->input.transformers_name[i]);
        status += functions->transformers[i].init(settings, tile_description);
    }

    // For output...
    functions->nsegment_output_functions = settings->output.nsegment_output_function;
    functions->segment_output_functions = (grasp_output_segment_function_t *) trackmem_malloc(sizeof (grasp_output_segment_function_t) * functions->nsegment_output_functions);
    assert(functions->segment_output_functions!=NULL);
    for (i = 0; i < functions->nsegment_output_functions; i++) {
        functions->segment_output_functions[i] = (grasp_output_segment_function_t) grasp_output_get_segment_output_funtion(settings->output.segment_output_function[i]);
        status += functions->segment_output_functions[i].init(settings, tile_description);
    }


    functions->ntile_output_functions = settings->output.ntile_output_function;
    functions->tile_output_functions = (grasp_output_tile_function_t *) trackmem_malloc(sizeof (grasp_output_tile_function_t) * functions->ntile_output_functions);
    assert(functions->tile_output_functions!=NULL);
    for (i = 0; i < functions->ntile_output_functions; i++) {
        functions->tile_output_functions[i] = (grasp_output_tile_function_t) grasp_output_get_tile_output_funtion(settings->output.tile_output_function[i]);
        status += functions->tile_output_functions[i].init(settings, tile_description);
    }

    functions->ncurrent_output_functions = settings->output.ncurrent_output_function;
    functions->current_output_functions = (grasp_output_current_function_t *) trackmem_malloc(sizeof (grasp_output_current_function_t) * functions->ncurrent_output_functions);
    assert(functions->current_output_functions!=NULL);
    for (i = 0; i < functions->ncurrent_output_functions; i++) {
        functions->current_output_functions[i] = (grasp_output_current_function_t) grasp_output_get_current_output_funtion(settings->output.current_output_function[i]);
        status += functions->current_output_functions[i].init(settings, tile_description);
    }

    return status;
}


