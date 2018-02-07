#include <math.h>
#include <limits.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <stdio.h>

void discard_rest_of_line(FILE *fp) {
  while (fgetc(fp) != '\n') {
    if (feof(fp)) break;
  }
}

char *strtolower(const char *data) {
    char *c;
    int i, max;

    max = strlen(data);
    // I make a copy of input string converting it to lower case
    c = (char *) malloc(sizeof (char)*max + 1);
    for (i = 0; i < max; i++) {
        c[i] = tolower(data[i]);
    }
    c[i] = '\0';

    return c;
}

int safe_atoi (const char *s, int *result){
   long lval;
   char *end;

   lval = strtol (s, &end, 10);
   if (lval > INT_MAX || lval < INT_MIN){
      // The number was not in range int.
      return -1;
   }else if (end == s){
      // Unrecognised string
      return -2;
   }else if (*end){
      // This string has more than a integer number
      return -3;
   }

   *result= (int)lval;
   return 0;
}

int safe_atofloat (const char *s, float *result){
    char *endptr;
    
    float d = strtof(s, &endptr);

    if (s == endptr){
        // Unrecognised string
        return -2;
    }else{
        if(isinf(d)){
            // Is grater than the range
            return -1;
        }else{
            if(isnan(d)){
                return -3;
            }else{
                *result=d;
                return 0;            
            }
        }
    }
}
