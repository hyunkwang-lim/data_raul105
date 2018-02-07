#include <sys/time.h> /* gettimeofday(): SVr4, 4.3BDS, POSIX.1-2001 */
#include <time.h> /* clock(): ANSI C89, C99, POSIX.1-2001 */
#include <unistd.h>

// Main benchmark structure which store initial, final and delta values.
typedef struct benchmark_t_{
  // Initial timeval
  double ct0;
  // Initial clock
  clock_t ut0;
  // Final timeval
  double ct1;
  // Final clock
  clock_t ut1;
  // Real time in seconds
  double delta_ut;
  // User time in seconds
  double delta_ct;
} benchmark_t;

#define benchmark_print(label,benchmark) printf("%s:%d: %s: clock_time: %f cpu_time: %f\n", __FILE__,__LINE__,label,(benchmark)->delta_ct, (benchmark)->delta_ut);

#ifdef	__cplusplus
extern "C" {
#endif

// Start benchmark variable
int benchmark_start(benchmark_t *benchmark);

// Stop benchmark variable
int benchmark_stop(benchmark_t *benchmark);

#ifdef	__cplusplus
}
#endif
