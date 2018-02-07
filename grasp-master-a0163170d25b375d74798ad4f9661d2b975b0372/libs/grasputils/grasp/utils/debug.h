/**
 * @file:  debug.h
 * @Author: Fabrice Ducos <fabrice.ducos@univ-lille1.fr>
 *
 */

#ifndef DEBUG_H
#define DEBUG_H

#include <stdio.h>
#include <assert.h>

#define DEBUG_STREAM stderr

#ifdef APPNAME
#define DISPLAY_APPNAME_IF_AVAILABLE do { fprintf(DEBUG_STREAM, "%s: ", APPNAME); } while (0)
#else
#define DISPLAY_APPNAME_IF_AVAILABLE
#endif

/* the debugging code will be compiled only if DEBUG is defined 
 * CAUTIOUS: if variables are declared in the wrapped code (especially
 * C++ or C99 style declarations that can be put anywhere), their scope will
 * be reduced to the body of the 'do' statement.
 */
#ifdef DEBUG
#define Debug(code) do { code } while (0);
#else
#define Debug(code)
#endif

/* a way to comment out a part of code (to avoid imbrications of C comments). */
#ifdef DISABLE
#define Disable(code)
#endif

/* measures the time spent in a part of code */
#ifdef BENCH
#include <time.h>
#define Bench(code) \
do { \
    clock_t t0 = clock(); \
    code; \
    clock_t t1 = clock(); \
    double t01 = (t1 - t0)/((double) CLOCKS_PER_SEC); \
    fprintf(DEBUG_STREAM, "%s:%d: %f s (%ld us)\n", __FILE__, __LINE__, t01, t1 - t0); \
} while (0);
#else
#define Bench(code) do { code } while (0)
#endif

/* PCALL displays each call to the current function in DEBUG mode
 * (to be put in the body of the function, usually at the beginning) */
#ifdef DEBUG
#define PCALL \
  do { \
    DISPLAY_APPNAME_IF_AVAILABLE; \
    fprintf(DEBUG_STREAM, "%s:%d: %s called\n", __FILE__, __LINE__, __FUNCTION__); \
    } while (0);
#else
#define PCALL
#endif

/* like PCALL but displays only the first call */
#ifdef DEBUG
#define PCALL_ONCE \
  do { \
    static int warn_trigger = 1;   \
    if (warn_trigger == 1) {	   \
      DISPLAY_APPNAME_IF_AVAILABLE; \
      fprintf(DEBUG_STREAM, "%s:%d: %s called\n", __FILE__, __LINE__, __FUNCTION__); \
      warn_trigger = 0; \
    } \
  } while(0);
#else
#define PCALL_ONCE
#endif

/* displays the call stack in the current scope ;
 * the function names and offset will be available only if the
 * code is compiled with the -rdynamic option. Otherwise only
 * the addresses will be displayed.
 */
#if defined(_GNU_SOURCE) && ! defined(__CYGWIN32__) /* execinfo.h not available on Cygwin */
#include <execinfo.h>
#define PRINT_CALLSTACK do {						\
    void *array[10];							\
    size_t size;							\
    char **strings;							\
    size_t i;								\
    									\
    size = backtrace (array, 10);					\
    strings = backtrace_symbols (array, size);				\
    									\
    for (i = 0; i < size; i++) {					\
      DISPLAY_APPNAME_IF_AVAILABLE;					\
      fprintf (DEBUG_STREAM, "%s\n", strings[i]);				\
    }									\
    free (strings);							\
  } while (0);
#else
#define PRINT_CALLSTACK
#endif /* _GNU_SOURCE && __CYGWIN32__ */

/* displays fatal errors, due to a bug or a deficiency in the software. Displays the call stack
 * and aborts the program. */
#define PRINT_FATAL(...) do {	\
    DISPLAY_APPNAME_IF_AVAILABLE; \
    fprintf(DEBUG_STREAM, "%s:%d: fatal error: ", __FILE__, __LINE__); \
    fprintf(DEBUG_STREAM, __VA_ARGS__); \
    PRINT_CALLSTACK;		  \
    abort();			  \
  } while (0);

/* displays errors due to a bad usage of the software, or a bad user input. Should be dealt with properly. */
#define PRINT_ERROR(...) do {	\
    DISPLAY_APPNAME_IF_AVAILABLE; \
    fprintf(DEBUG_STREAM, __VA_ARGS__); \
  } while (0);

/* displays system errors (a format string with printf-like arguments can be given to specify the source of the error,
 * do not add neither a trailing newline nor columns. Columns and the system error message will be appended to the user message) 
 */
#define PRINT_SYS_ERROR(...) do {	\
    DISPLAY_APPNAME_IF_AVAILABLE; \
    fprintf(DEBUG_STREAM, __VA_ARGS__); \
    fprintf(DEBUG_STREAM, ": "); \
    perror(NULL); \
  } while (0);

/* displays an informational message in debug mode (macro DEBUG defined) */
#ifdef DEBUG
#define PRINT_DEBUG(...) do { \
    DISPLAY_APPNAME_IF_AVAILABLE; \
    fprintf(DEBUG_STREAM, "%s:%d: ", __FILE__, __LINE__);	\
    fprintf(DEBUG_STREAM, __VA_ARGS__);  \
  } while (0);
#else
#define PRINT_DEBUG(...)
#endif

/* displays a user-defined warning, once */
#define WARN_ONCE(...) \
  do { \
    static int warn_trigger = 1;   \
    if (warn_trigger == 1) {	   \
      DISPLAY_APPNAME_IF_AVAILABLE; \
      fprintf(DEBUG_STREAM, __VA_ARGS__);  \
      warn_trigger = 0; \
    } \
  } while(0);

/* to be put in functions still to be implemented (will be displayed on the first call only) */
#define NOT_YET_IMPLEMENTED() \
  do { \
    static int warn_trigger = 1;   \
    if (warn_trigger == 1) {	   \
      DISPLAY_APPNAME_IF_AVAILABLE; \
      fprintf(DEBUG_STREAM, "%s:%d: %s: not yet implemented\n", __FILE__, __LINE__, __FUNCTION__); \
      warn_trigger = 0; \
    } \
  } while(0);

/* to be put in functions with quick and dirty implementations (working in most common cases,
 * but whose performance is not optimal, or whose error handling is often to be reworked or
 * whose behaviour is not guaranteed in all use cases)
 */
#define QUICK_AND_DIRTY() \
  do { \
    static int warn_trigger = 1;   \
    if (warn_trigger == 1) {	   \
      DISPLAY_APPNAME_IF_AVAILABLE; \
      fprintf(DEBUG_STREAM, "%s:%d: %s: quick and dirty (possibly buggy or not optimal)\n", __FILE__, __LINE__, __FUNCTION__); \
      warn_trigger = 0; \
    } \
  } while(0);

#ifdef __GNUC__
#define DEPRECATED(func) func __attribute__ ((deprecated))
#elif defined(_MSC_VER)
#define DEPRECATED(func) __declspec(deprecated) func
#else
#pragma message("WARNING: You need to implement DEPRECATED for this compiler")
#define DEPRECATED(func) func
#endif

//don't use me any more
// DEPRECATED(void OldFunc(int a, float b));

//use me instead
// void NewFunc(int a, double b);


#endif /* DEBUG_H */
