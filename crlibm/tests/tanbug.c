#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main (int argc, char *argv[]) 
{ 
  double input, result;
  double (*testfun_libm)   () = NULL;

  sscanf(argv[1],"%le", &input);

  testfun_libm=sin;
  printf("Computing sine of %.50e:\n",input);
  fflush(stdout);
  result=testfun_libm(input);
  printf("result is %.50e\n", result); 

  testfun_libm=tan;
  printf("Computing tangent of %.50e:\n",input);
  fflush(stdout);
  result=testfun_libm(input);
  printf("result is %.50e \n", result); 

  return 0;
}

