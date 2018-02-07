#include <stdlib.h>

int strfortranlen(const char *x, int maxsize) {
    int i;

    for (i = maxsize - 1; i >= 0; i--) {
        if (x[i] != ' ') return i + 1;
    }

    return 0;
}


char *fstr2c(const char *x, int maxsize) {
    int i, length;
    char *result;

    length = strfortranlen(x, maxsize);
    result = (char *) malloc(sizeof (char)*length + 1);

    for (i = 0; i < length; i++) {
        result[i] = x[i];
    }
    result[i] = '\0';

    return result;
}

