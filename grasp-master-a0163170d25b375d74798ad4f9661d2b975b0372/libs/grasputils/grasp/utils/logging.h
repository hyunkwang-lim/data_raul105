/**
 * @file:  logging.h
 * @Author: Fabrice Ducos <fabrice.ducos@univ-lille1.fr>
 *
 */

#ifndef LOGGING_H_
#define LOGGING_H_

#include <stdio.h> /* FILE * */
#include <stdbool.h> /* bool */

extern void set_error_stream(FILE *stream);
extern void set_debug_stream(FILE *stream);
extern void set_debug_level(int level);
extern FILE *get_error_stream(void);
extern FILE *get_debug_stream(void);
extern int  get_debug_level(void);
extern bool debug_is_enabled(void);
extern void perror2(FILE *stream, const char *s);

#endif /* LOGGING_H_ */
