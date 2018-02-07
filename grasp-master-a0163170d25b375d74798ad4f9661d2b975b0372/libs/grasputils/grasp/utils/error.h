/**
 * @file: error.h 
 * @Author: Fabrice Ducos <fabrice.ducos@univ-lille1.fr>
 *
 */


#ifndef ERROR_H_
#define ERROR_H_

#include <errno.h>
#include <stdarg.h>

/* the function to set the prefix of the error messages;
 * usually the prefix will be an application name;
 * e.g., if the prefix is set "program", the error messages will
 * look like this: "program: message" 
 */
extern void set_error_prefix(const char *prefix);

/* the macro to set a new error message. 
 * Granted, macros are usually evil when used with
 * efficiency reasons in mind, but the macro here
 * simply relieves the user of setting himself/herself the
 * first arguments of set_last_error_()
 */
#define SET_LAST_ERROR(sys_errnum, format, ...) do {			\
    set_last_error_( __FILE__, __LINE__, sys_errnum, format, __VA_ARGS__);	\
  } while (0);


extern const char *get_last_error(void);
extern void print_last_error(void);   /* prints the last error without adding a newline character */
extern void println_last_error(void); /* prints the last error and a newline character */
extern void reset_last_error(void);

/* functions below (with a trailing underscore), though part of the interface,
 * should only be called through the SET_LAST_ERROR macro, and never
 * directly.
 */

/* file: the file containing the call (usually set to __FILE__)
   line: the line of the call (usually set to __LINE__)
   sys_errnum: should be set to errno for a system error, and to 0 for a user error
   format: a format string in the printf-style
   other arguments: extra arguments for the format string
 */
extern void set_last_error_(const char *file, int line, int sys_errnum, const char *format, ...);
extern void set_last_error_va_(const char *file, int line, int sys_errnum, const char *format, va_list arg);

#endif /* ERROR_H_ */
