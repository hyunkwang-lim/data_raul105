#include <stdio.h>

// Use this macro if you want to add constant value like a string (e.g.: "3")
#define CONST2STR(s) STR(s)
#define STR(s) #s

// Return a copy of the data string in lowercase
char *strtolower(const char *data);   

// Transform safely string to intenger. Return 0 if everything is ok or number lower than 0 otherwise
// return -1 if number is outside int range
// return -2 if s is unrecognise string
// return -3 if after transform s to int there are more characters in s
// result is the output value
int safe_atoi (const char *s, int *result);

// Transform safely string to float. Return 0 if everything is ok or a number lower than 0 otherwise
// return -1 if number is outside float range
// return -2 if s is unrecognise string or after transform s to int there are more characters in s
// return -3 if can not perform transformation
// result is the output value
int safe_atofloat (const char *s, float *result);

// Move read point of fp to next new line
void discard_rest_of_line(FILE *fp) ;
