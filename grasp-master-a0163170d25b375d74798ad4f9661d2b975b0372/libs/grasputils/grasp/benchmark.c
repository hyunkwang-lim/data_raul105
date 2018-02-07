#include "utils/benchmark.h"

double benchmark_convert_timeval_to_seconds(struct timeval *tv) {
  return tv->tv_sec + 1.0*tv->tv_usec / 1000000.;
}


int benchmark_start(benchmark_t *benchmark){
  struct timeval tv;
  
  gettimeofday(&tv, NULL);
  benchmark->ct0 = benchmark_convert_timeval_to_seconds(&tv);
  benchmark->ut0 = clock();
  benchmark->ct1=benchmark->ct0;
  benchmark->ut1=benchmark->ut0;
  benchmark->delta_ut=0;
  benchmark->delta_ct=0;

  return 0;
}

int benchmark_stop(benchmark_t *benchmark){
  struct timeval tv;

  gettimeofday(&tv, NULL);
  benchmark->ct1 = benchmark_convert_timeval_to_seconds(&tv);
  benchmark->ut1 = clock();

  benchmark->delta_ut = 1.0*(benchmark->ut1 - benchmark->ut0) / CLOCKS_PER_SEC;
  benchmark->delta_ct = 1.0*(benchmark->ct1 - benchmark->ct0);

  return 0;
}

/*
int main(void) {
  // Define benchmark variable
  benchmark_t benchmark;

  // Start benchmark
  benchmark_start(&benchmark);  

  // Do something
  //long_processing();

  // Stop benchmark
  benchmark_stop(&benchmark);  

  // Print result
  benchmark_print(test, &benchmark);

  return 0;
}
*/

