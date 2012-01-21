#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "crlibm.h"
#include "crlibm_private.h"
#include "interval.h"
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
db_number input, input2, res_crlibm_low, res_crlibm_up, res_mpfr_low, res_mpfr_up, res_libultim_low, res_libultim_up, res_libmcr_low, res_libmcr_up, res_libm_low, res_libm_up;
interval input_i, res_crlibm;
mpfr_t mp_res, mp_inpt, mp_inpt2; 

interval j_exp(interval x);
interval j_log(interval x);
interval j_expm1(interval x);
interval j_log2(interval x);

double (*randfun)       () = NULL;
double (*randfun_perf)       () = NULL;
interval (*testfun_crlibm_interval)() = NULL;
double (*testfun_crlibm_low)() = NULL;
double (*testfun_crlibm_up)() = NULL;
int    (*testfun_mpfr)  () = NULL;


#define PRINT_INPUT_ERROR\
  printf("  x =%.70e \n         (%08x %08x) \n", input.d, input.i[HI], input.i[LO]);\
  printf("  y =%.70e \n         (%08x %08x) \n",input2.d,input2.i[HI],input2.i[LO]);\



void test_filter_cases(){
  interval input,res,waited;
  db_number res_low, res_up, waited_low, waited_up;
  printf("When low>up... ");
  ASSIGN_LOW(input,3);
  ASSIGN_UP(input,2);
  ASSIGN_LOW(waited,0.0/0.0);
  ASSIGN_UP(waited,0.0/0.0);
  res=j_log(input);
  res_low.d=LOW(res);
  res_up.d=UP(res);
  waited_low.d=LOW(waited);
  waited_up.d=UP(waited);
  if (  (waited_low.i[LO]==res_low.i[LO])
     && (waited_low.i[HI]==res_low.i[HI])
     && (waited_up.i[LO]==res_up.i[LO])
     && (waited_up.i[HI]==res_up.i[HI]))
  {
    printf("OK\n");
  }
   else
  {
    printHexa("Waited low:",waited_low.d);
    printHexa("Waited up:",waited_up.d);
    printHexa("Given low:",res_low.d);
    printHexa("Given up:",res_up.d);
  }

  printf("When up<0... ");
  ASSIGN_LOW(input,3);
  ASSIGN_UP(input,-2);
  ASSIGN_LOW(waited,0.0/0.0);
  ASSIGN_UP(waited,0.0/0.0);
  res=j_log(input);
  res_low.d=LOW(res);
  res_up.d=UP(res);
  waited_low.d=LOW(waited);
  waited_up.d=UP(waited);
  if (  (waited_low.i[LO]==res_low.i[LO])
     && (waited_low.i[HI]==res_low.i[HI])
     && (waited_up.i[LO]==res_up.i[LO])
     && (waited_up.i[HI]==res_up.i[HI]))
  {
    printf("OK\n");
  }
   else
  {
    printHexa("Waited low:",waited_low.d);
    printHexa("Waited up:",waited_up.d);
    printHexa("Given low:",res_low.d);
    printHexa("Given up:",res_up.d);
  }

  printf("When low=1, up=1... ");
  ASSIGN_LOW(input,1.0);
  ASSIGN_UP(input,1.0);
  ASSIGN_LOW(waited,0.0);
  ASSIGN_UP(waited,0.0);
  res=j_log(input);
  res_low.d=LOW(res);
  res_up.d=UP(res);
  waited_low.d=LOW(waited);
  waited_up.d=UP(waited);
  if (  (waited_low.i[LO]==res_low.i[LO])
     && (waited_low.i[HI]==res_low.i[HI])
     && (waited_up.i[LO]==res_up.i[LO])
     && (waited_up.i[HI]==res_up.i[HI]))
  {
    printf("OK\n");
  }
   else
  {
    printHexa("Waited low:",waited_low.d);
    printHexa("Waited up:",waited_up.d);
    printHexa("Given low:",res_low.d);
    printHexa("Given up:",res_up.d);
  }
  printf("When low<0 and up>0... ");
  ASSIGN_LOW(input,-100.0);
  ASSIGN_UP(input,2);
  ASSIGN_LOW(waited,-1.0/0.0);
  ASSIGN_UP(waited,-1.0/0.0);
  res=j_log(input);
  res_low.d=LOW(res);
  res_up.d=UP(res);
  waited_low.d=LOW(waited);
  waited_up.d=UP(waited);
  if (  (waited_low.i[LO]==res_low.i[LO])
     && (waited_low.i[HI]==res_low.i[HI]))
   {
    printf("OK\n");
  }
   else
  {
    printHexa("Waited low:",waited_low.d);
    printHexa("Waited up:",waited_up.d);
    printHexa("Given low:",res_low.d);
    printHexa("Given up:",res_up.d);
  }

}

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
    input2.d = randfun();
    if (input.d>input2.d)
    {
      double temp=input.d;
      input.d=input2.d;
      input2.d=temp;
    }
/*    db_number ia,ib;
    ia.i[HI]=0x31100afb;
    ia.i[LO]=0x198a95fe;
    ib.i[HI]=0x3d42897b;
    ib.i[LO]=0x84591a4e;
    input.d=ia.d; input2.d=ib.d;*/
    ASSIGN_LOW(input_i,input.d);
    ASSIGN_UP(input_i,input2.d);
   

    res_crlibm = testfun_crlibm_interval(input_i);
    res_crlibm_low.d=LOW(res_crlibm);
    res_crlibm_up.d=UP(res_crlibm);
    mpfr_set_d(mp_inpt, input.d, GMP_RNDN);
    testfun_mpfr(mp_res, mp_inpt, GMP_RNDD);
    res_mpfr_low.d = mpfr_get_d(mp_res, GMP_RNDD);
    mpfr_set_d(mp_inpt, input2.d, GMP_RNDN);
    testfun_mpfr(mp_res, mp_inpt, GMP_RNDU);
    res_mpfr_up.d = mpfr_get_d(mp_res, GMP_RNDU);
/*    printHexa("resul crlibm low:",res_crlibm_low.d);
    printHexa("resul crlibm up:",res_crlibm_up.d);
    printHexa("resul mpfr low:",res_mpfr_low.d);
    printHexa("resul mpfr up:",res_mpfr_up.d);*/
    

#if PRINT_NAN
    if(1){
#else
      if( ((res_mpfr_low.i[HI] & 0x7ff00000) != 0x7ff00000)
       && ((res_mpfr_up.i[HI] & 0x7ff00000) != 0x7ff00000) ) {
#endif
	if( (res_crlibm_low.i[LO] != res_mpfr_low.i[LO]) 
	    || (res_crlibm_low.i[HI] != res_mpfr_low.i[HI])
	    || (res_crlibm_up.i[LO] != res_mpfr_up.i[LO])
	    || (res_crlibm_up.i[HI] != res_mpfr_up.i[HI]) )
	    {
#if DETAILED_REPORT
	  printf("*** CRLIBM ERROR ***\n");
	  PRINT_INPUT_ERROR;
	  printf("crlibm gives    [%.50e,%.50e] \n         [(%08x %08x),(%08x %08x)] \n", 
		 res_crlibm_low.d, res_crlibm_up.d,
		 res_crlibm_low.i[HI], res_crlibm_low.i[LO],
		 res_crlibm_up.i[HI], res_crlibm_up.i[LO]);
	  printf("MPFR gives      [%.50e,%.50e] \n         [(%08x %08x),(%08x %08x)] \n\n", 
		 res_mpfr_low.d, 
		 res_mpfr_up.d,
		 res_mpfr_low.i[HI], 
		 res_mpfr_low.i[LO],
		 res_mpfr_up.i[HI],
		 res_mpfr_up.i[LO]
		 );
#endif
#if WORST_ERROR_REPORT
	  mpfr_set_d(mp_inpt, res_crlibm_low.d, GMP_RNDN);  
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
	  mpfr_set_d(mp_inpt, res_crlibm_up.d, GMP_RNDN);
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
  fprintf (stderr, "\nUsage: %s function seed \n", fct_name);
  fprintf (stderr, " function      : name of function to test \n");
  fprintf (stderr, " seed          : integer seed for the random number generator \n");
  exit (1);
}

int main (int argc, char *argv[]) 
{ 
  char* function_name;
  int seed;
  if ((argc != 3)) usage(argv[0]);
  else{
    function_name = argv[1];
    sscanf(argv[2],"%d", &seed);
    if ((strcmp(function_name,"exp")!=0) && (strcmp(function_name,"j_exp")!=0) && (strcmp(function_name,"log")!=0) && (strcmp(function_name,"j_log")!=0) && (strcmp(function_name,"expm1")!=0) && (strcmp(function_name,"j_expm1")!=0) && (strcmp(function_name,"log2")!=0) && (strcmp(function_name,"j_log2")!=0))
    {
      fprintf (stderr, "\nUnknown function:  %s \n", function_name);
      return 1;
    }
    if ((strcmp(function_name,"log")==0) || (strcmp(function_name,"j_log")==0))
    {
      randfun = rand_for_log;
      testfun_crlibm_interval = j_log;
      testfun_crlibm_low = log_rd;
      testfun_crlibm_up = log_ru;
      testfun_mpfr = mpfr_log;
    }
    if ((strcmp(function_name,"log2")==0) || (strcmp(function_name,"j_log2")==0))
    {
      randfun = rand_for_log;
      testfun_crlibm_interval = j_log2;
      testfun_crlibm_low = log2_rd;
      testfun_crlibm_up = log2_ru;
      testfun_mpfr = mpfr_log2;
    }
    if ((strcmp(function_name,"exp")==0) || (strcmp(function_name,"j_exp")==0))
    {
      randfun = rand_generic;
      testfun_crlibm_interval = j_exp;
      testfun_crlibm_low = exp_rd;
      testfun_crlibm_up = exp_ru;
      testfun_mpfr = mpfr_exp;
    }
    if ((strcmp(function_name,"expm1")==0) || (strcmp(function_name,"j_expm1")==0))
    {
      randfun = rand_for_expm1_soaktest;
      testfun_crlibm_interval = j_expm1;
      testfun_crlibm_low = expm1_rd;
      testfun_crlibm_up = expm1_ru;
      testfun_mpfr = mpfr_expm1;
    }

    crlibm_init();


    mpfr_init2(mp_res,  200);
    mpfr_init2(mp_inpt, 53);
    mpfr_init2(mp_inpt2, 53);

  printf("Testing %s function \n",function_name);


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
