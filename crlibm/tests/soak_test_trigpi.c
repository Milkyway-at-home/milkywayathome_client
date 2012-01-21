#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "crlibm.h"
#include "crlibm_private.h"
#include "log-interval.h"
#ifdef HAVE_MPFR_H  /* stop here if MPFR not present */
#include "test_common.h"
#include <gmp.h>
#include <mpfr.h>

#ifdef HAVE_MATHLIB_H
#include <MathLib.h>
#endif

/* Stupidely soak-tests a function against mpfr */

/* if set to 1, print out the worst error (wiht values of x)
   if set to 0, don't print out worst error information */
#define WORST_ERROR_REPORT 1 

/* if set to 1, print out detailed errors (vith values of x)
   if set to 0, only count the errors and print out the count */
#define DETAILED_REPORT 1

/* if set to 1, print out errors due to NaNs
   if set to 0, don't print NaNs' errors */
#define PRINT_NAN 0
/* What we are looking for here is misrounding. Therefore these tests
   concern only intervals on which the function does something
   useful. For example for exp, it is enough to soaktest on -1024,1024.

   To achieve this we have several random generator function, tailored
   for each function or group of function.

   The other cases (including all the special cases) are supposedly
   tested exhaustively by the other programs of this directory. */ 


/* Basic-like programming with global variables: */
db_number input, input2, res_crlibm, res_crlibm_low, res_crlibm_up, res_mpfr, res_libultim_low, res_libultim_up, res_libmcr_low, res_libmcr_up, res_libm_low, res_libm_up;
interval input_i;
mpfr_t mp_res, mp_inpt, mp_inpt2, mp_pi, mp_piinpt; 

#define PRINT_INPUT_ERROR\
  printf("  x =%.70e \n         (%08x %08x) \n", input.d, input.i[HI], input.i[LO]);\


double sinpi_rn();
double sinpi_rd();
double sinpi_ru();
double sinpi_rz();

double rand_for_trig_perf();

double (*randfun)       () = NULL;
double (*testfun_crlibm)() = NULL;
int    (*testfun_mpfr)  () = NULL;
mp_rnd_t mpfr_rnd_mode;

/*
 * Give the number of missrounded results 
 */

void test_all() {
  long long int 
    failures_libm=0,
#ifdef HAVE_MATHLIB_H
    failures_libultim=0,
#endif
#ifdef HAVE_LIBMCR_H
    failures_libmcr=0,
#endif
    failures_crlibm=0;
  long long int i;
  double worst_err, global_worst_err=-200;
  db_number global_worst_inpt, global_worst_inpt2;

  i=0; 
  while(1+1==2)
  {
    input.d = randfun();
    res_crlibm.d = testfun_crlibm(input.d);
    mpfr_const_pi (mp_pi,GMP_RNDN);
    mpfr_set_d(mp_inpt, input.d, GMP_RNDN);
    mpfr_mul(mp_piinpt,mp_pi,mp_inpt,GMP_RNDN);
    testfun_mpfr(mp_res, mp_piinpt, mpfr_rnd_mode);
    res_mpfr.d = mpfr_get_d(mp_res, mpfr_rnd_mode);
/*    printHexa("resul crlibm low:",res_crlibm_low.d);
    printHexa("resul crlibm up:",res_crlibm_up.d);
    printHexa("resul mpfr low:",res_mpfr_low.d);
    printHexa("resul mpfr up:",res_mpfr_up.d);*/
    

#if PRINT_NAN
    if(1){
#else
      if( ((res_mpfr.i[HI] & 0x7ff00000) != 0x7ff00000) ) {
#endif
	if( (res_crlibm.i[LO] != res_mpfr.i[LO]) 
	    || (res_crlibm.i[HI] != res_mpfr.i[HI]) )
	    {
#if DETAILED_REPORT
	  printf("*** CRLIBM ERROR ***\n");
	  PRINT_INPUT_ERROR;
	  printf("crlibm gives    %.50e \n         (%08x %08x) \n", 
		 res_crlibm.d,
		 res_crlibm.i[HI],
		 res_crlibm.i[LO]);
	  printf("MPFR gives      %.50e \n         (%08x %08x) \n\n", 
		 res_mpfr.d, 
		 res_mpfr.i[HI], 
		 res_mpfr.i[LO]);
#endif
#if WORST_ERROR_REPORT
	  mpfr_set_d(mp_inpt, res_crlibm.d, GMP_RNDN);  
	  mpfr_sub(mp_inpt, mp_inpt, mp_res, GMP_RNDN);  
	  mpfr_div(mp_inpt, mp_inpt, mp_res, GMP_RNDN);
	  mpfr_abs(mp_inpt, mp_inpt, GMP_RNDN);
	  mpfr_log2(mp_inpt, mp_inpt, GMP_RNDN);
	  worst_err=mpfr_get_d(mp_inpt, GMP_RNDN);

	  if (worst_err>global_worst_err){
	    global_worst_err=worst_err;
	    global_worst_inpt.d  = input.d;
	    global_worst_inpt2.d = input2.d;
	  }
	  printf("Worst crlibm relative error so far : 2^(%f)\n",global_worst_err);
	  printf(" for x =%.50e     (%08x %08x) \n", 
		 global_worst_inpt.d, 
		 global_worst_inpt.i[HI], 
		 global_worst_inpt.i[LO]);
#endif

	  failures_crlibm++;
	}
	
      }
      i++;
      if((i % 10000)==0) 
      {
	printf(" CRLIBM       : %lld failures out of %lld (ratio %e) \n",failures_crlibm, i,
	       ((double)failures_crlibm)/(double)i);
	printf("\n");
	
      }
  }
}




void usage(char *fct_name){
  /* fprintf (stderr, "\n%s: Soak-test for crlibm and other mathematical libraries \n", fct_name); */
  fprintf (stderr, "\nUsage: %s function (RN|RU|RD|RZ) seed \n", fct_name);
  fprintf (stderr, " function      : name of function to test \n");
  fprintf (stderr, " (RN|RU|RD|RZ) : rounding mode, \n");
  fprintf (stderr, " seed          : integer seed for the random number generator \n");
  exit (1);
}



int main (int argc, char *argv[]) 
{ 
  int seed;
  char* function_name;
  char* rounding_mode;


  crlibm_init();
  mpfr_init2(mp_res,  130);
  mpfr_init2(mp_inpt, 53);
  mpfr_init2(mp_inpt2, 53);
  mpfr_init2(mp_pi, 5000);
  mpfr_init2(mp_piinpt, 10000);


  if ((argc != 4)) usage(argv[0]);
  else{
    function_name = argv[1];
    rounding_mode = argv[2];
    sscanf(argv[2],"%d", &seed);
    if ((strcmp(function_name,"sinpi")!=0))
    {
      fprintf (stderr, "\nUnknown function:  %s \n", function_name);
      return 1;
    }
    if ((strcmp(function_name,"sinpi")==0))
    {
      randfun = rand_for_trig_perf;
      testfun_mpfr = mpfr_sin;
    }
    if      (strcmp(rounding_mode,"RU")==0) { mpfr_rnd_mode = GMP_RNDU; testfun_crlibm = sinpi_ru; }
    else if (strcmp(rounding_mode,"RD")==0) { mpfr_rnd_mode = GMP_RNDD; testfun_crlibm = sinpi_rd; }
    else if (strcmp(rounding_mode,"RZ")==0) { mpfr_rnd_mode = GMP_RNDZ; testfun_crlibm = sinpi_rz; }
    else {
      mpfr_rnd_mode = GMP_RNDN;
      testfun_crlibm = sinpi_rn;
      rounding_mode="RN" ;
    }
    printf("Testing %s function with rounding mode : %s \n", function_name, rounding_mode);
    srand(seed);
    test_all();
  }
  return 0;
}


#else
int main (int argc, char *argv[]) 
{ 
  printf("Sorry, I need to be compiled against MPFR\n  (get it from www.mpfr.org, then:   configure --enable-mpfr)\n");
  return 0;
}
#endif
