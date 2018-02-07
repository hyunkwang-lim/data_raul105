/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

/* grasp_mpi_engine.c */

/* Fabrice Ducos <fabrice.ducos@univ-lille1.fr> 2014, 2015 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <grasp/utils.h>
#include <unistd.h> /* getpid, getppid, alarm, sleep (for simulating processing) */
#include <assert.h>

/* this is for interruptible_task (and also <unistd.h>, already included) */
#include <signal.h>
#include <setjmp.h>
#include <stdbool.h>

#ifdef USE_MPI
#include <mpi.h>
#endif

#ifdef DEBUG_MPI
static int g_debug_level = 1;
#else
static int g_debug_level = 0;
#endif

#include "grasp_mpi_engine.h"

static char g_appname[255 + 1];
static int g_order_payload_max_size;
static int g_result_payload_max_size;

static order_callback_t prepare_order_callback; /* routine set by the user for preparing orders */
static task_callback_t perform_task_callback; /* routine set by the user for performing a task */
static collect_callback_t collect_result_callback = NULL; /* routine set by the user for collecting results (will be ignored if set to NULL) */
static progress_info_callback_t progress_info_callback = NULL;

/* The labels in order_t and result_t are not necessary for the processing (only the payload is).
 * They contain strings with readable information about the content of an order or
 * a result. They are useful for debugging.
 */

struct order_t_ {
  unsigned char *message; /* this pointer references a sequence of bytes holding the metadata and the payload to be passed through MPI */
  char *label; /* a reference to message[ORDER_OFFSET_LABEL], for ease of use. */
  void *payload; /* a reference to message[ORDER_OFFSET_PAYLOAD], for ease of use. */
  bool *ptr_no_more_job; /* a reference to message[ORDER_OFFSET_NO_MORE_JOB] (its address, not its value) for ease of use.  */
  size_t *ptr_payload_size; /* a reference to message[ORDER_OFFSET_PAYLOAD_SIZE], for ease of use */
  size_t *ptr_message_size; /* a reference to message[ORDER_OFFSET_MESSAGE_SIZE], for ease of use */
};

struct result_t_ {
  unsigned char *message; /* this pointer references a sequence of bytes holding the metadata and the payload to be passed through MPI */
  char *label; /* a reference to message[RESULT_OFFSET_LABEL] for ease of use */
  void *payload; /* a reference to message[RESULT_OFFSET_PAYLOAD] for ease of use */
  size_t *ptr_payload_size; /* a reference to message[RESULT_OFFSET_PAYLOAD_SIZE], for ease of use */
  size_t *ptr_message_size; /* a reference to message[RESULT_OFFSET_MESSAGE_SIZE], for ease of use */
};

/* const variables (const int or const size_t) can't be used for declaring sizes of global arrays in plain C, so one uses enums (#define's would work too) */
enum {
  ORDER_LABEL_MAX_LEN = 255 /* this does not include the ending null character */
};

enum {
  RESULT_LABEL_MAX_LEN = 255 /* this does not include the ending null character */
};

enum {
  ORDER_WRAPPER_SIZE = 1024
};

enum {
  RESULT_WRAPPER_SIZE = 1024
};

enum {
  ORDER_OFFSET_NO_MORE_JOB = 0,
  ORDER_OFFSET_MESSAGE_SIZE = 8,
  ORDER_OFFSET_PAYLOAD_SIZE = 16,
  ORDER_OFFSET_LABEL = 24, /* the (ORDER_LABEL_MAX_LEN + 1) bytes starting from this offset are reserved for the order label */
  ORDER_OFFSET_PAYLOAD = ORDER_WRAPPER_SIZE
};

enum {
  RESULT_OFFSET_MESSAGE_SIZE = 0,
  RESULT_OFFSET_PAYLOAD_SIZE = 8,
  RESULT_OFFSET_LABEL = 16, /* the (RESULT_LABEL_MAX_LEN + 1) bytes starting from this offset are reserved for the result label */
  RESULT_OFFSET_PAYLOAD = RESULT_WRAPPER_SIZE
};

enum {
  WORKER_STATUS_TAG,
  ORDER_TAG,
  RESULT_IS_READY_TAG
};

enum {
  NO_AVAILABLE_WORKER = -1
};

typedef enum worker_status_t_ {
  WORKER_STATUS_BUSY, /* the worker has received an order and is busy on its task */
  WORKER_STATUS_AVAILABLE /* the worker is available for a new task. */
} worker_status_t;

typedef struct process_info_t_ {
#ifdef USE_MPI
  MPI_Comm communicator;
#endif
  char label[255 + 1];
  int pid;
  int ppid;
  int master_rank;
  int my_rank;
  int num_processes;
  int num_workers;
} process_info_t;

#ifdef USE_MPI
static void print_order_aux(const char *file, int line, FILE *stream, const char *print_label, const order_t *order, const process_info_t *process_info) {
  const char *process_info_label;

  assert(stream != NULL);
  assert(order != NULL);

  if (process_info == NULL) {
    process_info_label = "unknown";
  }
  else {
    process_info_label = process_info->label;
  }

  if (print_label == NULL) {
    print_label = "order";
  }

  fprintf(stream, "%s:%d: %s: %s->payload_size: %lu\n", file, line, process_info_label, print_label, (long unsigned) *(order->ptr_payload_size));
  fprintf(stream, "%s:%d: %s: %s->message_size: %lu\n", file, line, process_info_label, print_label, (long unsigned) *(order->ptr_message_size));
  fprintf(stream, "%s:%d: %s: %s->label: %s\n", file, line, process_info_label, print_label, order->label);
  fprintf(stream, "%s:%d: %s: %s->no_more_job: %s\n\n", file, line, process_info_label, print_label, *(order->ptr_no_more_job) ? "true" : "false");
  
}

#define print_order(stream, print_label, order, process_info) print_order_aux(__FILE__, __LINE__, stream, print_label, order, process_info)

static void print_result_aux(const char *file, int line, FILE *stream, const char *print_label, const result_t *result, const process_info_t *process_info) {
  const char *process_info_label;

  assert(stream != NULL);
  assert(result != NULL);

  if (process_info == NULL) {
    process_info_label = "unknown";
  }
  else {
    process_info_label = process_info->label;
  }

  if (print_label == NULL) {
    print_label = "result";
  }

  fprintf(stream, "%s:%d: %s: %s->payload_size: %lu\n", file, line, process_info_label, print_label, (long unsigned) *(result->ptr_payload_size));
  fprintf(stream, "%s:%d: %s: %s->message_size: %lu\n", file, line, process_info_label, print_label, (long unsigned) *(result->ptr_message_size));
  fprintf(stream, "%s:%d: %s: %s->label: %s\n", file, line, process_info_label, print_label, result->label);
}

#define print_result(stream, print_label, result, process_info) print_result_aux(__FILE__, __LINE__, stream, print_label, result, process_info)
#endif /* USE_MPI */

/* this internal routine sets derived fields from an order object to their correct value (derived fields are all the fields but message);
 * to be called each time order is updated
 */
static void set_order_references(order_t *order) {
  assert(order != NULL);
  assert(order->message != NULL);
  
  order->label = (char *) &order->message[ORDER_OFFSET_LABEL];
  order->payload = &order->message[ORDER_OFFSET_PAYLOAD];
  order->ptr_payload_size = (size_t *) &order->message[ORDER_OFFSET_PAYLOAD_SIZE];
  order->ptr_message_size = (size_t *) &order->message[ORDER_OFFSET_MESSAGE_SIZE];
  order->ptr_no_more_job = (bool *) &order->message[ORDER_OFFSET_NO_MORE_JOB];
}

/* this internal routine sets derived fields from a result object to their correct value (derived fields are all the fields but message);
 * to be called each time result is updated
 */
static void set_result_references(result_t *result) {
  assert(result != NULL);
  assert(result->message != NULL);
  
  result->label = (char *) &result->message[RESULT_OFFSET_LABEL];
  result->payload = &result->message[RESULT_OFFSET_PAYLOAD];
  result->ptr_payload_size = (size_t *) &result->message[RESULT_OFFSET_PAYLOAD_SIZE];
  result->ptr_message_size = (size_t *) &result->message[RESULT_OFFSET_MESSAGE_SIZE];
  
}

order_t *grasp_mpi_engine_new_order(const char *label, size_t payload_size, const void *payload) {
  order_t *order;
  size_t message_size = ORDER_WRAPPER_SIZE + payload_size;

  order = calloc(1, sizeof(order_t));
  assert(order != NULL); /* if the memory is exhausted, it's useless to continue */
  
  if (! (payload_size <= g_order_payload_max_size)) {
    fprintf(stderr, "%s:%d: %s: the argument payload_size (%d) is out of range [0..%d]. Note that you can change the upper bound with grasp_mpi_engine_set_order_callback()\n",
	    __FILE__, __LINE__, __func__, (int) payload_size, (int) g_order_payload_max_size);
  }
  assert(payload_size <= g_order_payload_max_size);
  order->message = calloc(1, message_size);
  assert(order->message != NULL);
  set_order_references(order);
  *(order->ptr_no_more_job) = false;
  *(order->ptr_payload_size) = payload_size;
  *(order->ptr_message_size) = message_size;
  
  if (label != NULL) {
    strncpy(order->label, label, ORDER_LABEL_MAX_LEN);
  }
  
  if (payload_size > 0 && payload != NULL) {
    memcpy(order->payload, payload, payload_size);
  }
  
  return order;
}

result_t *grasp_mpi_engine_new_result(const char *label, size_t payload_size, const void *payload) {
  result_t *result;
  size_t message_size = RESULT_WRAPPER_SIZE + payload_size;
  
  result = calloc(1, sizeof(result_t));
  assert(result != NULL);  /* if the memory is exhausted, it's useless to continue */
  
  if (! (payload_size <= g_result_payload_max_size)) {
    fprintf(stderr, "%s:%d: %s: the argument payload_size (%d) is out of range [0..%d]. Note that you can change the upper bound with grasp_mpi_engine_set_task_callback()\n",
	    __FILE__, __LINE__, __func__, (int) payload_size, (int) g_result_payload_max_size);
  }
  
  assert(payload_size <= g_result_payload_max_size);
  result->message = calloc(1, message_size);
  assert(result->message != NULL);
  set_result_references(result);
  *(result->ptr_payload_size) = payload_size;
  *(result->ptr_message_size) = message_size;

  if (label != NULL) {
    strncpy(result->label, label, RESULT_LABEL_MAX_LEN);
  }

  if (payload_size > 0 && payload != NULL) {
    memcpy(result->payload, payload, payload_size);
  }

  return result;
}

void grasp_mpi_engine_delete_order(order_t *order) {
  assert(order != NULL);
  
  trackmem_free(order->message);
  trackmem_free(order);
}

void grasp_mpi_engine_delete_result(result_t *result) {
  assert(result != NULL);

  trackmem_free(result->message);
  trackmem_free(result);
}

void grasp_mpi_engine_get_order_label(const order_t *order, size_t label_max_size, char *label) {
  assert(order != NULL);
  assert(label != NULL);
  strncpy(label, (const char *) &order->message[ORDER_OFFSET_LABEL], label_max_size);
}

void grasp_mpi_engine_get_result_label(const result_t *result, size_t label_max_size, char *label) {
  assert(result != NULL);
  assert(label != NULL);
  strncpy(label, (const char *) &result->message[RESULT_OFFSET_LABEL], label_max_size);
}

#ifdef USE_MPI
static void get_process_info(process_info_t *process_info) {
  MPI_Comm communicator = MPI_COMM_WORLD;
  assert(process_info != NULL);
  
  memset(process_info, 0, sizeof(process_info_t));
  process_info->pid = getpid();
  process_info->ppid = getppid();
  MPI_Comm_size(communicator, &process_info->num_processes);
  MPI_Comm_rank(communicator, &process_info->my_rank);

  process_info->communicator = communicator;
  process_info->num_workers = process_info->num_processes - 1;
  process_info->master_rank = process_info->num_processes - 1;

  if (process_info->my_rank == process_info->master_rank) {
    strcpy(process_info->label, "master");
  }
  else {
    sprintf(process_info->label, "worker #%d", process_info->my_rank);
  }
}

static void print_process_info_aux(const char *file, int line, FILE *stream, const char *print_label, const process_info_t *process_info) {
  assert(stream != NULL);
  assert(process_info != NULL);

  if (print_label == NULL) {
    print_label = "process_info";
  }
  
  fprintf(stream, "%s:%d: %s: %s->label:         %s\n", file, line, process_info->label, print_label, process_info->label);
  fprintf(stream, "%s:%d: %s: %s->pid:           %d\n", file, line, process_info->label, print_label, process_info->pid);
  fprintf(stream, "%s:%d: %s: %s->ppid:          %d\n", file, line, process_info->label, print_label, process_info->ppid); 
  fprintf(stream, "%s:%d: %s: %s->master_rank:   %d\n", file, line, process_info->label, print_label, process_info->master_rank);
  fprintf(stream, "%s:%d: %s: %s->my_rank:       %d\n", file, line, process_info->label, print_label, process_info->my_rank);
  fprintf(stream, "%s:%d: %s: %s->num_processes: %d\n", file, line, process_info->label, print_label, process_info->num_processes);
  fprintf(stream, "%s:%d: %s: %s->num_workers:   %d\n\n", file, line, process_info->label, print_label, process_info->num_workers);
}

#define print_process_info(stream, print_label, process_info) print_process_info_aux(__FILE__, __LINE__, stream, print_label, process_info)

static bool I_am_a_worker(const process_info_t *process_info) {
  assert(process_info != NULL);
  return process_info->my_rank != process_info->master_rank;
}

static bool I_am_the_master(const process_info_t *process_info) {
  assert(process_info != NULL);
  return process_info->my_rank == process_info->master_rank;
}

/********************** Worker's workflow *********************/

static void notify_master_of_availability(const process_info_t *process_info, int worker_status) {
  const char *cstr_status;

  assert(process_info != NULL);
  assert(I_am_a_worker(process_info) == true);
  
  switch (worker_status) {
  case WORKER_STATUS_BUSY:
    cstr_status = "busy";
    break;
  case WORKER_STATUS_AVAILABLE:
    cstr_status = "available";
    break;
  default:
    fprintf(stderr, "%s:%d: %s: unexpected value for worker_status: %d\n", __FILE__, __LINE__, process_info->label, worker_status);
    abort();
  }
  
  if (g_debug_level > 0) {
    fprintf(stderr, "%s:%d: %s: I notify the master that I'm %s\n", __FILE__, __LINE__, process_info->label, cstr_status);
  }
  MPI_Ssend(&worker_status, 1, MPI_INT, process_info->master_rank, WORKER_STATUS_TAG, process_info->communicator);
}

static order_t *wait_for_next_order(const process_info_t *process_info) {
  order_t *order;
  unsigned char *recv_buffer;
  const size_t recv_buffer_size = ORDER_WRAPPER_SIZE + g_order_payload_max_size;
  const char *label;
  void *payload;
  bool no_more_job, *ptr_no_more_job;
  size_t payload_size, *ptr_payload_size;
  size_t message_size, *ptr_message_size;
  
  assert(process_info != NULL);
  assert(I_am_a_worker(process_info) == true);
  assert(recv_buffer_size > 0);
  
  if (g_debug_level > 0) {
    fprintf(stderr, "%s:%d: %s: I'm waiting for a new order\n", __FILE__, __LINE__, process_info->label);
  }
  
  recv_buffer = calloc(1, recv_buffer_size);
  assert(recv_buffer != NULL);
  
  MPI_Recv(recv_buffer, recv_buffer_size, MPI_BYTE, process_info->master_rank, ORDER_TAG, process_info->communicator, MPI_STATUS_IGNORE);
  if (g_debug_level > 0) {
    fprintf(stderr, "%s:%d: %s: I've just received a new order\n", __FILE__, __LINE__, process_info->label);
  }
  
  label = (const char *) &recv_buffer[ORDER_OFFSET_LABEL];
  payload = &recv_buffer[ORDER_OFFSET_PAYLOAD];
  ptr_no_more_job  = (bool *) &recv_buffer[ORDER_OFFSET_NO_MORE_JOB];
  ptr_payload_size = (size_t *) &recv_buffer[ORDER_OFFSET_PAYLOAD_SIZE];
  ptr_message_size = (size_t *) &recv_buffer[ORDER_OFFSET_MESSAGE_SIZE];
  no_more_job  = *ptr_no_more_job;
  payload_size = *ptr_payload_size;
  message_size = *ptr_message_size;
  order = grasp_mpi_engine_new_order(label, payload_size, payload);
  set_order_references(order);
  *(order->ptr_no_more_job) = no_more_job;
  *(order->ptr_payload_size) = payload_size;
  *(order->ptr_message_size) = message_size;

  trackmem_free(recv_buffer);

  if (g_debug_level > 0) {
    print_order(stderr, "order[received from master]", order, process_info);
  }
  
  return order;
}

static int g_maximum_job_time = 0; // 30*60; /* Maximal time allowed for a task before interruption, in seconds (0 means no interruption) */
void grasp_mpi_engine_set_maximum_job_time(int maximum_job_time) {
  if (g_debug_level > 0) {
    fprintf(stderr, "%s:%d: %s(maximum_job_time = %d) called\n", __FILE__, __LINE__, __func__, maximum_job_time);
  }

  if (maximum_job_time < 0) {
    fprintf(stderr, "%s: unexpected value for the maximum job time (%d s). Should be positive or null. Will be set to 0.\n", 
	    g_appname, maximum_job_time);
    maximum_job_time = 0;
  }

  g_maximum_job_time = maximum_job_time;
}

static sigjmp_buf jmp_buffer;

void alarm_handler(int signum) {
  siglongjmp(jmp_buffer, 1 /* return value */);
}

/* a wrapper that prevents a task from lasting too long and blocking a worker, possibly forever.
 * This may happen (and did), when the application algorithm doesn't converge, or converges too slowly.
 * The task will be interrupted after maximum_job_time seconds.
 */
result_t *interruptible_task(const order_t *order, const void *worker_info, int maximum_job_time /* in seconds */) {
  result_t *result;

  assert(order != NULL);

  signal(SIGALRM, alarm_handler); // handler calls siglongjmp
  alarm(maximum_job_time);
  if ( sigsetjmp(jmp_buffer, 0 /* savemask */) == 0 ) {
    /* return from sigsetjmp() */
    result = perform_task_callback(order, worker_info); /* do the work */
    alarm(0); /* cancel alarm */
    return result; /* the task arrived to completion */
  } else {
    /* return from siglongjmp() */
    return NULL; /* it took too long */
  }
}

/* this routine should always return a valid result (even empty) */
static result_t *execute_order(const process_info_t *process_info, const void *worker_info, const order_t *order) {
  result_t *result;

  assert(process_info != NULL);
  assert(order != NULL);
  assert(I_am_a_worker(process_info) == true);

  /* perform the real work here */
  if (g_debug_level > 0) {
    fprintf(stderr, "%s:%d: %s: I'm working on %s\n", __FILE__, __LINE__, process_info->label, order->label);
  }

  assert(perform_task_callback != NULL);

  if (g_maximum_job_time > 0) { /* interruptible tasks are enabled */  
    result = interruptible_task(order, worker_info, g_maximum_job_time); /* the task may not last more the g_maximum_job_time seconds. */
    if (result == NULL) { /* the task has been interrupted*/
      char label[RESULT_LABEL_MAX_LEN + 1];

      fprintf(stderr, "%s:%d: %s: %s took more than %d seconds (%d min and %d sec) to process, there must be a convergence problem. I give up with this task\n", 
	      __FILE__, __LINE__, process_info->label, order->label, g_maximum_job_time, g_maximum_job_time / 60, g_maximum_job_time % 60);

      snprintf(label, RESULT_LABEL_MAX_LEN, "%s could not be completed (more than the time limit: %d seconds, i.e. %d min and %d sec)", 
	       order->label, g_maximum_job_time, g_maximum_job_time / 60, g_maximum_job_time % 60);

      result = grasp_mpi_engine_new_result(label, 0, NULL); /* create an empty result (the routine guarantees to return a valid result, even empty) */
    } /* result == NULL */
  } /* g_maximum_job_time > 0 */
  else { /* interruptions are disabled. This can be a problem when a task never finishes. */
    result = perform_task_callback(order, worker_info); /* the task will be performed directly without possibility of interruption */
  }

  if (g_debug_level > 0) {
    fprintf(stderr, "%s:%d: %s: I've finished my work on %s\n", __FILE__, __LINE__, process_info->label, order->label);
  }

  assert(result != NULL);
  return result;
}

static void send_result_to_master(const process_info_t *process_info, result_t *result) {
  int count;

  assert(process_info != NULL);
  assert(result != NULL);
  assert(I_am_a_worker(process_info) == true);
  
  if (g_debug_level > 0) {
    print_result(stderr, "result[I'm about to send]", result, process_info);
  }

  count = *(result->ptr_message_size);
  assert(count > 0);
  assert(result->message != NULL);
  assert(process_info->communicator == MPI_COMM_WORLD);
  MPI_Send(result->message, count, MPI_BYTE, process_info->master_rank, RESULT_IS_READY_TAG, process_info->communicator);
}

static void worker_main_loop(const process_info_t *process_info, const void *worker_info) {
  bool done = false;
  order_t *order = NULL;
  result_t *result = NULL;

  assert(process_info != NULL);
  assert(I_am_a_worker(process_info) == true);

  while (! done) {
    notify_master_of_availability(process_info, WORKER_STATUS_AVAILABLE);
    order = wait_for_next_order(process_info);
    done = *(order->ptr_no_more_job);
    
    if (! done) {
      notify_master_of_availability(process_info, WORKER_STATUS_BUSY);
      result = execute_order(process_info, worker_info, order);
      send_result_to_master(process_info, result);
      grasp_mpi_engine_delete_result(result);
    }

    grasp_mpi_engine_delete_order(order);
  }

  if (g_debug_level > 0) {
    fprintf(stderr, "%s:%d: %s: I'm stopping now\n", __FILE__, __LINE__, process_info->label);
  }
}

/********************** Master's workflow *********************/
static void send_order_to_worker(const process_info_t *process_info, const order_t *order, int worker_rank) {
  assert(process_info != NULL);
  assert(order != NULL);
  
  assert(I_am_the_master(process_info));
  if (g_debug_level > 0) {
    fprintf(stderr, "%s:%d: %s: I'm sending an order to worker #%d\n", __FILE__, __LINE__, process_info->label, worker_rank);
  }
  assert(0 <= worker_rank && worker_rank < process_info->num_workers);
  if (g_debug_level > 0) {
    char label[255 + 1];
    sprintf(label, "order[sent to worker #%d]", worker_rank);
    print_order(stderr, label, order, process_info);
  }
  
  MPI_Ssend((void *) order->message, *(order->ptr_message_size), MPI_BYTE, worker_rank, ORDER_TAG, process_info->communicator);
  if (g_debug_level > 0) {
    fprintf(stderr, "%s:%d: %s: order has been sent to worker #%d\n", __FILE__, __LINE__, process_info->label, worker_rank);
  }
}

static void print_worker_status_aux(const char *file, int line, FILE *stream, const char *label, 
				    int num_workers, const int *worker_status, const process_info_t *process_info) {
  int worker_rank;
  
  assert(file != NULL);
  assert(stream != NULL);
  assert(num_workers > 0);
  assert(worker_status != NULL);
  assert(process_info != NULL);
  
  fprintf(stream, "%s:%d: %s: ", file, line, process_info->label);
  if (label == NULL) {
    label = "worker_status";
  }

  fprintf(stream, "%s (num_workers: %d): [", label, num_workers);
  for (worker_rank = 0 ; worker_rank < num_workers ; worker_rank++) {
    fprintf(stream, " #%d(%d) ", worker_rank, worker_status[worker_rank]);
  }
  fprintf(stream, "]\n");
}

#define print_worker_status(stream, label, num_workers, worker_status, process_info) \
  print_worker_status_aux(__FILE__, __LINE__, stream, label, num_workers, worker_status, process_info)

/* a listen_* is non blocking. Its purpose is to initiate an asynchronous reception */
static void listen_to_worker_availability(const process_info_t *process_info, MPI_Request *requests_for_available_workers, int *worker_status, int worker_rank) {
  assert(process_info != NULL);
  assert(requests_for_available_workers != NULL);
  assert(worker_status != NULL);
  assert(I_am_the_master(process_info));

  MPI_Irecv(&worker_status[worker_rank], 1, MPI_INT, worker_rank, WORKER_STATUS_TAG, 
	      process_info->communicator, &requests_for_available_workers[worker_rank]);
}

/* a listen_* is non blocking. Its purpose is to initiate an asynchronous reception */
static void listen_to_workers_availability(const process_info_t *process_info, MPI_Request *requests_for_available_workers, int *worker_status) {
  int worker_rank;

  for (worker_rank = 0 ; worker_rank < process_info->num_workers ; worker_rank++) {
    listen_to_worker_availability(process_info, requests_for_available_workers, worker_status, worker_rank);
  }
}

/* a listen_* is non blocking. Its purpose is to initiate an asynchronous reception */
static void listen_to_result_availability(const process_info_t *process_info, MPI_Request *requests_for_available_results, result_t **results, int worker_rank) {
  assert(process_info != NULL);
  assert(requests_for_available_results != NULL);
  assert(I_am_the_master(process_info));
  
  MPI_Irecv(results[worker_rank]->message, *(results[worker_rank]->ptr_message_size), MPI_BYTE, worker_rank, RESULT_IS_READY_TAG, process_info->communicator, &requests_for_available_results[worker_rank]);
}

/* a listen_* is non blocking. Its purpose is to initiate an asynchronous reception */
static void listen_to_results_availability(const process_info_t *process_info, MPI_Request *requests_for_available_results, result_t **results) {
  int worker_rank;
  
  for (worker_rank = 0 ; worker_rank < process_info->num_workers ; worker_rank++) {
    listen_to_result_availability(process_info, requests_for_available_results, results, worker_rank);
  }
}

static int get_first_available_worker(const process_info_t *process_info, const int *worker_status) {
  int worker_rank;

  assert(process_info != NULL);
  assert(worker_status != NULL);
  
  for (worker_rank = 0 ; worker_rank < process_info->num_workers ; worker_rank++) {
    if (worker_status[worker_rank] == WORKER_STATUS_AVAILABLE) {
      return worker_rank;
    }
  }

  return NO_AVAILABLE_WORKER;
}

static int g_polling_time = 1; /* the master process will sleep during g_polling_time seconds before checking for available workers.
				* This is done to avoid spinning (busy waiting) of the application. The default polling time
				* is one second and is generally long enough for reducing considerably the system overhead without degrading
				* the performances of the application.
				*/

void grasp_mpi_engine_set_polling_time(int polling_time) {
  if (g_debug_level > 0) {
    fprintf(stderr, "%s:%d: %s(polling_time = %d) called\n", __FILE__, __LINE__, __func__, polling_time);
  }

  if (polling_time < 0) {
    fprintf(stderr, "%s: unexpected value for the polling time (%d s). Should be positive or null. Will be set to 0.\n", 
	    g_appname, polling_time);
    polling_time = 0;
  }

  g_polling_time = polling_time;
}

/* a wait_* routine is blocking */
static int wait_for_available_worker(const process_info_t *process_info, MPI_Request *requests_for_available_workers, int *worker_status) {
  int index = MPI_UNDEFINED;
  int num_workers; 
  int available_worker_rank = -1;
  int worker_status_has_been_updated = 0;
 
  assert(process_info != NULL);
  assert(requests_for_available_workers != NULL);
  assert(worker_status != NULL);
  assert(I_am_the_master(process_info));
  
  num_workers = process_info->num_workers;
  
  if (g_debug_level > 0) {
    fprintf(stderr, "%s:%d: %s: I'm checking the available workers\n", __FILE__, __LINE__, process_info->label);
  }

  if (get_first_available_worker(process_info, worker_status) == NO_AVAILABLE_WORKER) {
    /* all the workers are busy, wait for one */
    MPI_Waitany(num_workers, requests_for_available_workers, &index, MPI_STATUS_IGNORE);
    worker_status_has_been_updated = 1;
  }
  else {
    /* some workers are available, but update the availability status of all the workers before going on */
    MPI_Testany(num_workers, requests_for_available_workers, &index, &worker_status_has_been_updated, MPI_STATUS_IGNORE);
    
    if (g_polling_time > 0) {
      sleep(g_polling_time); /* allowing the master process to sleep can reduce the load on the system (against a performance cost) */
    }
  }
  
  if (worker_status_has_been_updated) {
    assert(index != MPI_UNDEFINED);
    if (g_debug_level > 1) {
      fprintf(stderr, "%s:%d: %s: the worker #%d has updated its availability status\n", __FILE__, __LINE__, process_info->label, index);
      print_worker_status(stderr, "workers availability status", num_workers, worker_status, process_info);
    }
  
    assert(0 <= index && index < process_info->num_workers);
    listen_to_worker_availability(process_info, requests_for_available_workers, worker_status, index);
    
    if (worker_status[index] == WORKER_STATUS_AVAILABLE) {
      available_worker_rank = index;
      return available_worker_rank;
    }
  }
  
  /* no availability update has occured, simply check the current worker_status table */
  available_worker_rank = get_first_available_worker(process_info, worker_status);
  
  if (g_debug_level > 0) {
    if (available_worker_rank >= 0) {
      fprintf(stderr, "%s:%d: %s: the worker #%d is available\n", __FILE__, __LINE__, process_info->label, available_worker_rank);
    }
    else {
      fprintf(stderr, "%s:%d: %s: all the workers are busy\n", __FILE__, __LINE__, process_info->label);
    }
  }

  return available_worker_rank;
}

/* a peek_* routine is non blocking */
static int peek_for_available_result(const process_info_t *process_info, MPI_Request *requests_for_available_results, result_t **results) {
  int index;
  int flag;
  
  assert(process_info != NULL);
  assert(requests_for_available_results != NULL);
  assert(I_am_the_master(process_info));

  MPI_Testany(process_info->num_workers, requests_for_available_results, &index, &flag, MPI_STATUS_IGNORE);
  if (flag == true) {
    listen_to_result_availability(process_info, requests_for_available_results, results, index);
    return index;
  }
  else {
    return -1;
  }
}

/* the master sends stop orders for asking workers to stop themselves */
static order_t *new_stop_order(int worker_rank) {
  order_t *order;
  char label[ORDER_LABEL_MAX_LEN + 1];

  sprintf(label, "stop order for worker #%d", worker_rank);
  order = grasp_mpi_engine_new_order(label, /* payload_size */ 0, /* payload */ NULL);
  assert(order != NULL);
  order->message[ORDER_OFFSET_NO_MORE_JOB] = true;
  
  set_order_references(order);

  return order;
}

static bool everything_is_done(bool no_more_order, int pending_orders) {
  assert(pending_orders >= 0);
  return no_more_order && (pending_orders == 0);
}

static void master_main_loop(const process_info_t *process_info, const void *master_info) {
  MPI_Request *requests_for_available_workers;
  MPI_Request *requests_for_available_results;
  result_t    **results;
  int *worker_status; /* availability status of workers */
  int num_orders_performed = 0;
  int worker_rank;
  order_t *order;
  int pending_orders = 0; /* orders received but not yet executed */
  bool no_more_order = false;
  bool done;

  assert(process_info != NULL);
  assert(I_am_the_master(process_info));
  assert(prepare_order_callback != NULL);

  if (process_info->num_workers < 1) {
    fprintf(stderr, "%s:%d: no worker found. Please run the code with mpiexec -n <num_processes>, where <num_processes> == <num_workers> + 1\n",
	    __FILE__, __LINE__);
    exit (EXIT_FAILURE);
  }

  requests_for_available_workers = trackmem_malloc(process_info->num_workers * sizeof(MPI_Request));
  assert(requests_for_available_workers!=NULL);
  worker_status = trackmem_malloc(process_info->num_workers * sizeof(*worker_status));
  assert(worker_status!=NULL);
  requests_for_available_results = trackmem_malloc(process_info->num_workers * sizeof(MPI_Request));
  assert(requests_for_available_results!=NULL);
  
  results = trackmem_malloc(process_info->num_workers * sizeof(result_t *));
  assert(results!=NULL);
  for (worker_rank = 0 ; worker_rank < process_info->num_workers ; worker_rank++) {
    result_t *result = grasp_mpi_engine_new_result(/* label */ NULL, g_result_payload_max_size, /* payload */ NULL);
    results[worker_rank] = result;
    worker_status[worker_rank] = WORKER_STATUS_BUSY;
  }

  listen_to_workers_availability(process_info, requests_for_available_workers, worker_status);
  listen_to_results_availability(process_info, requests_for_available_results, results);
  
  /* issue orders and retrieve results */
  done = everything_is_done(no_more_order, pending_orders); /* do not enter the loop if no result is expected (empty data set) */
  while (! done) {
    done = everything_is_done(no_more_order, pending_orders);
    if (g_debug_level > 0) {
        if (no_more_order) {
          fprintf(stderr, "%s:%d: %s: no more order to process, %d pending order(s)\n", __FILE__, __LINE__, process_info->label,
                  pending_orders);
        }
    }
    if (done) break;

    worker_rank = wait_for_available_worker(process_info, requests_for_available_workers, worker_status);
    if (! (0 <= worker_rank && worker_rank < process_info->num_workers)) {
      continue; /* no available worker has been found */
    }
    
    if (worker_status[worker_rank] == WORKER_STATUS_AVAILABLE) {
      order = prepare_order_callback(worker_rank, master_info, &no_more_order);
      if (order != NULL) {
	assert(no_more_order == false);
	send_order_to_worker(process_info, order, worker_rank);
	grasp_mpi_engine_delete_order(order);
	pending_orders++;
	worker_status[worker_rank] = WORKER_STATUS_BUSY;
      }
      /* if no_more_order == true, the production of orders is exhausted and one must exit the loop */
    }

    worker_rank = peek_for_available_result(process_info, requests_for_available_results, results);
    if (worker_rank >= 0) {
      num_orders_performed++;
      pending_orders--;

      if (collect_result_callback != NULL) {
        collect_result_callback(results[worker_rank], master_info);
      }

      if (g_debug_level > 0) {
	char print_label[255 + 1];
	sprintf(print_label, "result[from worker #%d]", worker_rank);
	fprintf(stderr, "%s:%d: %s: I've received results from worker #%d\n", __FILE__, __LINE__, process_info->label,
		worker_rank);
	print_result(stderr, print_label, results[worker_rank], process_info);
      } /* g_debug_level > 0 */

      if (progress_info_callback != NULL) {
	progress_info_callback(num_orders_performed, master_info);
      }
    }
  }
  
  if (g_debug_level > 0) {
    fprintf(stderr, "%s:%d: %s: all the work is done (no new order, no more pending orders). \n", __FILE__, __LINE__, process_info->label);
  }

  /* ask the workers to stop their work */
  for (worker_rank = 0 ; worker_rank < process_info->num_workers ; worker_rank++) {
    order = new_stop_order(worker_rank);
    assert(order != NULL);
    send_order_to_worker(process_info, order, worker_rank);
    if (g_debug_level > 0) {
      fprintf(stderr, "%s:%d: %s: the worker #%d has been told to stop its work\n", __FILE__, __LINE__, process_info->label, worker_rank);
    }
    grasp_mpi_engine_delete_order(order);
  } /* for (worker_rank) */
  
  trackmem_free(requests_for_available_workers);
  trackmem_free(requests_for_available_results);
  trackmem_free(worker_status);
  
  for (worker_rank = 0 ; worker_rank < process_info->num_workers ; worker_rank++) {
    grasp_mpi_engine_delete_result(results[worker_rank]);
  }
  trackmem_free(results);
  
}
#endif /* #ifdef USE_MPI */

/* implementation of the interface */

void grasp_mpi_engine_set_debug_level(int debug_level) {
  g_debug_level = debug_level;
}

void grasp_mpi_engine_set_appname(const char *appname) {
  strncpy(g_appname, appname, sizeof(g_appname) - 1);
}

void *grasp_mpi_engine_get_order_payload(const order_t *order) {
  assert(order != NULL);
  return order->payload;
}

void *grasp_mpi_engine_get_result_payload(const result_t *result) {
  assert(result != NULL);
  return result->payload;
}

order_t *grasp_mpi_engine_dummy_order_callback(int worker_rank, const void *master_info, bool *no_more_order) {
  order_t *order;
  char label[255 + 1];
  const void *payload = NULL;
  size_t payload_size = 0;
  const int NUM_ORDERS_TOTAL = 3; /* number of orders to issue */
  static int num_orders_sent = 0;
  
  if (num_orders_sent < NUM_ORDERS_TOTAL) {
    sprintf(label, "order for worker #%d", worker_rank);
    order = grasp_mpi_engine_new_order(label, payload_size, payload);
    assert(order != NULL);
    num_orders_sent++;
    *no_more_order = false;
    return order;
  }
  else {
    /* no more order will be issued */
    *no_more_order = true;
    return NULL;
  }
}

result_t *grasp_mpi_engine_dummy_task_callback(const order_t *order, const void *worker_info) {
  result_t *result;
  const void *payload = NULL;
  size_t payload_size = 0;
  char order_label[255 + 1];
  char result_label[255 + 1];
  
  assert(order != NULL);
  
  /* SIMULATE THE EXECUTION OF AN ORDER */
  sleep(5);
  
  grasp_mpi_engine_get_order_label(order, sizeof(order_label) - 1, order_label);
  snprintf(result_label, sizeof(result_label) - 1, "%s (completed)", order_label);
  result = grasp_mpi_engine_new_result(result_label, payload_size, payload);
  assert(result != NULL);
  return result;
}

void grasp_mpi_engine_set_order_callback(order_callback_t order_callback, int order_payload_max_size) {
  assert(order_payload_max_size >= 0);
  
  if (order_callback != NULL) {
    g_order_payload_max_size = order_payload_max_size;
    prepare_order_callback = order_callback;
  }
  else {
    g_order_payload_max_size = 0;
    prepare_order_callback = grasp_mpi_engine_dummy_order_callback;
  }
}

void grasp_mpi_engine_set_task_callback(task_callback_t task_callback, int result_payload_max_size) {
  assert(result_payload_max_size >= 0);
  
  if (task_callback != NULL) {
    g_result_payload_max_size = result_payload_max_size;
    perform_task_callback = task_callback;
    
  }
  else {
    g_result_payload_max_size = 0;
    perform_task_callback = grasp_mpi_engine_dummy_task_callback;
  }
}

void grasp_mpi_engine_set_collect_callback(collect_callback_t collect_callback) {
  collect_result_callback = collect_callback; /* will be ignored if set to NULL */
}

void grasp_mpi_engine_set_progress_info_callback(progress_info_callback_t progress_callback) {
  progress_info_callback = progress_callback;
  
}

#ifdef USE_MPI
bool grasp_mpi_engine_is_master(void) {
  process_info_t process_info;

  get_process_info(&process_info);
  return I_am_the_master(&process_info);
}

bool grasp_mpi_engine_is_worker(void) {
  process_info_t process_info;

  get_process_info(&process_info);
  return I_am_a_worker(&process_info);
}

int grasp_mpi_engine_main_loop(const void *master_info, const void *worker_info) {
  process_info_t process_info;
  
  get_process_info(&process_info);
  if (g_debug_level > 0) {
    char label[255 + 1];
    sprintf(label, "process_info[#%d]", process_info.my_rank);
    print_process_info(stderr, label, &process_info);
  }
  
  if (I_am_the_master(&process_info)) {
    master_main_loop(&process_info, master_info);
  }
  else {
    worker_main_loop(&process_info, worker_info);
  }
  
  return EXIT_SUCCESS;
}
#else /* #ifdef USE_MPI is false */
static void grasp_mpi_engine_print_error(void) {
  fprintf(stderr, "%s: error: not an MPI-enabled application (please use BUILD=mpi|mpi-debug|mpi-prod)\n", g_appname);
}

bool grasp_mpi_engine_is_master(void) {
  return false;
}

bool grasp_mpi_engine_is_worker(void) {
  return false;
}

int grasp_mpi_engine_main_loop(const void *master_info, const void *worker_info) {
  grasp_mpi_engine_print_error();

  return EXIT_FAILURE;
}

#endif /* #ifdef USE_MPI */
