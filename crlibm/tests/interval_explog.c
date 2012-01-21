/* 
gcc  interval_explog.c ../libcrlibm.a

*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>

#include "../crlibm.h"


typedef struct interval { double INF, SUP;} interval ;

#define LOW(x) x.INF
#define UP(x) x.SUP

interval j_log(interval x);
interval j_exp(interval x);
int i;

int main() {
  interval x, y;
  crlibm_init();
  x.INF = 1.5;
  x.SUP = 1.5;

  //  x=j_log(x);
  //printf(" [ %1.30e  ,\n   %1.30e]", x.INF, x.SUP);
  //exit(0);


  for (i=0; i<100; i++) {
    y=j_log(x);
    x=j_exp(y);
    printHexa("yinf:",y.INF);
    printHexa("ysup:",y.SUP);
    printHexa("xinf:",x.INF);
    printHexa("xsup:",x.SUP);
  }
  printf("y: [%1.30e  ,   %1.30e]     ", y.INF, y.SUP);
  printf("Result: [%1.30e  ,   %1.30e]\n", x.INF, x.SUP);
}
