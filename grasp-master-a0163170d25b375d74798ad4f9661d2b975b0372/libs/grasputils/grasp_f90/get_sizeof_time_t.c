#include <stdio.h>
#include <time.h>

int main(void) {
  printf("%d\n", (int) (8*sizeof(time_t)));
  return 0;
}
