/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

/* grasp_mpi_engine.h */

/* Fabrice Ducos <fabrice.ducos@univ-lille1.fr> 2014 */

#ifndef GRASP_MPI_ENGINE_H
#define GRASP_MPI_ENGINE_H

#include <stdbool.h>

typedef struct order_t_ order_t;
typedef struct result_t_ result_t;

/* master_info and worker_info are opaque data that can be interpreted only by the callback functions 
 * If the author of callback functions don't need to pass extra information, he or she can leave these
 * arguments to NULL.
 */
typedef order_t * (*order_callback_t)(int worker_rank, const void *master_info, bool *no_more_order); /* master_info is optional and can be left to NULL 
												       * no_more_order must be set to true by the callback routine when
												       * the production of orders is exhausted (with a NULL return value)
												       */
typedef result_t * (*task_callback_t)(const order_t *order, const void *worker_info); /* worker_info is optional and can be left to NULL */
typedef void (*collect_callback_t)(const result_t *result, const void *master_info);
typedef void (*progress_info_callback_t)(int num_orders_performed, const void *master_info);

extern void grasp_mpi_engine_set_appname(const char *appname);

/* The user can provide callback routines for preparing orders, preparing results, collecting results and ending the processing */
extern void grasp_mpi_engine_set_order_callback(order_callback_t order_callback, int order_payload_max_size);
extern void grasp_mpi_engine_set_task_callback(task_callback_t task_callback, int result_payload_max_size);
extern void grasp_mpi_engine_set_collect_callback(collect_callback_t collect_callback);
extern void grasp_mpi_engine_set_progress_info_callback(progress_info_callback_t progress_info_callback);

extern order_t *grasp_mpi_engine_new_order(const char *label, size_t payload_size, const void *payload); /* to be called by an order callback function (see dummy example in grasp_mpi_engine.c) */
extern result_t *grasp_mpi_engine_new_result(const char *label, size_t payload_size, const void *payload); /* to be called by a result callback function (see dummy example in grasp_mpi_engine.c) */
extern void grasp_mpi_engine_delete_order(order_t *order); /* to be called only by the MPI workflow (implemented by grasp_mpi_engine.c) */
extern void grasp_mpi_engine_delete_result(result_t *result); /* to be called only by the MPI workflow (implemented by grasp_mpi_engine.c) */
extern void grasp_mpi_engine_get_order_label(const order_t *order, size_t label_max_size, char *label); /* label is assumed to be allocated by the caller */
extern void grasp_mpi_engine_get_result_labbel(const result_t *result, size_t label_max_size, char *label); /* label is assumed to be allocated by the caller */

extern void *grasp_mpi_engine_get_order_payload(const order_t *order); /* returns a read-only reference to the payload of an order */
extern void *grasp_mpi_engine_get_result_payload(const result_t *result); /* returns a read-only reference to the payload of a result */

extern void grasp_mpi_engine_set_debug_level(int debug_level);
extern bool grasp_mpi_engine_is_master(void); /* returns true if the current processing is the master */
extern bool grasp_mpi_engine_is_worker(void); /* returns true if the current processing is a worker */
extern int  grasp_mpi_engine_main_loop(const void *master_info, const void *worker_info); /* master_info and worker_info are optional information that can be passed
											       * to the callback functions (master_info will be passed to order_callback
											       * and worker_info to task_callback). They can be left to NULL if
											       * the callback functions don't need extra information.
											       */

extern void grasp_mpi_engine_set_polling_time(int polling_time); /* sets the polling time, in seconds (the master will sleep during polling_time seconds before checking the list of available workers). Default value is 1 second. Setting polling_time to 0 is possible but not recommended (it would produce busy waiting and increase considerably the overhead on the system). A too large value will degrade the performance of the application. */

extern void grasp_mpi_engine_set_maximum_job_time(int maximum_job_time); /* sets the maximum time, in seconds, allowed for a task. If the task takes longer, it will be interrupted (default value is 0, meaning no interruption) */

/* these are not really part of the API but dummy callback functions for quick tests, debugging and example of usage */
extern order_t *grasp_mpi_engine_dummy_order_callback(int worker_rank, const void *master_info, bool *no_more_order); /* called by the master to prepare an order for a worker */
extern result_t *grasp_mpi_engine_dummy_task_callback(const order_t *order, const void *worker_info); /* called by a worker in response to an order */
extern bool grasp_mpi_engine_dummy_termination_condition_callback(int num_results_received, const void *master_info); /* called by the master to stop the processing */

#endif /* GRASP_MPI_ENGINE_H */
