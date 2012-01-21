/* 
Beware to compile without optimizations

How to link with fi_lib:

gcc  -O0 -std=c99    -o crlibm_testperf_interval  test_perf-interval.o test_common.o ../fi_lib.a ../libcrlibm.a -lm -z muldefs

However there is a modification to make to rand_for_log, see in this function.

*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>

#include "crlibm.h"
#include "crlibm_private.h"
#include "test_common.h"
#include "interval.h"

#include "scs_lib/tests/tbx_timing.h"

#ifdef HAVE_MATHLIB_H
#include <MathLib.h>
#endif

#ifdef HAVE_LIBMCR_H
#include <libmcr.h>
#endif

#ifdef HAVE_MPFR_H
# include <gmp.h>
# include <mpfr.h>
# ifdef MPFR_VERSION
#  if MPFR_VERSION < MPFR_VERSION_NUM(2,2,0)
#   define mpfr_subnormalize(mp_res, inexact, mpfr_rnd_mode) 
#  endif
# else
#  define mpfr_get_version() "<2.1.0"
# endif
#endif


#define N1 20

#define TIMING_ITER 200

#define DETAILED_REPORT 0

/* If set, the behaviour of the function with respect to cache memory
   will be tested*/
#define TEST_CACHE 0

/*
 * Rounding mode to test
 */




#if EVAL_PERF==1  
/* counter of calls to the second step */
extern int crlibm_second_step_taken; 
int crlibm_first_step_taken;
#endif

#ifdef HAVE_MPFR_H  
/* The rounding mode for mpfr function */
mp_rnd_t mpfr_rnd_mode;
#endif

/* Global variable for time stamping */
static unsigned long long tbx_time;

/* Unused random number generator*/
double (*randfun_soaktest) () = NULL;
/* The random number generator*/
double (*randfun)       () = NULL;
/* The function we test */
interval (*testfun_crlibm_interval)() = NULL;
double (*testfun_crlibm_low)() = NULL;
double (*testfun_crlibm_up)() = NULL;
/* The function we trust */
int    (*testfun_mpfr)  () = NULL;
/* The function to show off against for accuracy  */
double (*testfun_libm)  () = NULL;
/* The function to show off against for performance */
double (*testfun_libultim)   () = NULL;
/* The last competitor  */
double (*testfun_libmcr)   () = NULL;


/* TESTSIZE doubles should be enough to flush the cache */
#define TESTSIZE 200000
 
#if TEST_CACHE
static double inputs[TESTSIZE];
static double inputs2[TESTSIZE];
#endif /* TEST_CACHE */


interval j_log(interval x);
interval j_exp(interval x);
interval j_expm1(interval x);
interval j_log2(interval x);
/* indicate the number of argument taken by the function */
static int nbarg;          






static void usage(char *fct_name){
  /*fprintf (stderr, "\n%s: Performance test for crlibm and other mathematical libraries\n", fct_name);*/
  fprintf (stderr, "\nUsage: %s function iterations \n", fct_name);
  fprintf (stderr, " function      : name of function to test \n");
  fprintf (stderr, " iterations    : number of iterations, also seed for the random number generator \n");
  exit (EXIT_FAILURE);
}


#if TEST_CACHE

static void fill_and_flush(int seed) {
  int i;
  srandom(seed);
  for (i=0; i<TESTSIZE; i++){
    inputs[i] = randfun();
    inputs2[i] = randfun();
  }
}

static void test_with_cache(const char *name, double (*testfun_low)(), double (*testfun_up)(), int n){
  int i, j, k;
  double i1, i2, rd;
  tbx_tick_t   t1, t2; 
  unsigned long long dt, min_dtsum, dtsum;

  if(testfun_low!=NULL && testfun_up!=NULL) { /* test if some functions are missing */
    printf("\n%s\n",name);
    for(i=1; i<=10000; i*=10) { /* i=1,10,100...*/
      min_dtsum=1<<30; 
      for(k=0;k<10;k++) { /* do the whole test 10 times and take the min */
	fill_and_flush(n);
	dtsum=0;
	for (j=0; j<i; j++) {
	  i1 = inputs[i];
	  i2 = inputs2[i];
	  if (i1>i2)
	  {
	    double temp=i1;
	    i1=i2;
	    i2=temp;
	  }

	  if (nbarg==1){
	    TBX_GET_TICK(t1);
	    rd = testfun(i1);
	    rd = testfun(i2);
	    TBX_GET_TICK(t2);
	  }else{
	    TBX_GET_TICK(t1);
	    rd = testfun(i1, i2);
	    TBX_GET_TICK(t2);
	  }

	  dt = TBX_TICK_RAW_DIFF(t1, t2)-tbx_time; 
	  dtsum += dt;
	}
	if (dtsum < min_dtsum) min_dtsum=dtsum; 
      }
      printf("  %d loops: \t avg time = %f ticks\n",i, ((double)min_dtsum)/i);
    }
  }
}
#endif /* TEST_CACHE */

static void test_without_cache(const char *name, 
			double i1,
			double i2,
			unsigned long long *lib_dtmin,
			unsigned long long *lib_dtmax,
			unsigned long long *lib_dtsum,
			int func_type){
  double result;
  unsigned long long dt, dtmin;
  tbx_tick_t   t1, t2; 
  int j;
#ifdef TIMING_USES_GETTIMEOFDAY
  int k;
#endif 
#ifdef HAVE_MPFR_H  
  mpfr_t mp_res, mp_inpt;
  mpfr_t mp_inpt2; /* For the pow function */
  int inexact;
  
  if (i1>i2)
  {
    double temp=i1;
    i1=i2;
    i2=temp;
  }


  /*  db_number ia,ib;
  ia.i[HI]=0x31100afb;
  ia.i[LO]=0x198a95fe;
  ib.i[HI]=0x3d42897b;
  ib.i[LO]=0x84591a4e;
  i1=ia.d; i2=ib.d;*/

  mpfr_init2(mp_res,  53);
  mpfr_init2(mp_inpt, 53);
  mpfr_init2(mp_inpt2, 53);
#endif

  if(testfun_crlibm_low!=NULL && testfun_crlibm_up!=NULL) { /* test if some functions are missing */
    dtmin=1<<30;
    /* take the min of N1 consecutive calls */
    for(j=0; j<N1; j++) {

      if (func_type == 0){              /* func_type = normal function */
	if (nbarg==1){
	  TBX_GET_TICK(t1);
#ifdef TIMING_USES_GETTIMEOFDAY /* use inaccurate timer, do many loops */
	  for(k=0; k<TIMING_ITER;k++)
#endif
	    result = testfun_crlibm_low(i1);
	    result = testfun_crlibm_up(i2);
	  TBX_GET_TICK(t2);
	}else{
	  TBX_GET_TICK(t1);
#ifdef TIMING_USES_GETTIMEOFDAY /* use inaccurate timer, do many loops */
	  for(k=0; k<TIMING_ITER;k++)
#endif
//	    result = testfun(i1,i2);
	  TBX_GET_TICK(t2);	  
	}
      }else{                             /* func_type = MPFR function  */
#ifdef   HAVE_MPFR_H 
	if (nbarg==1){
	  TBX_GET_TICK(t1);
#ifdef TIMING_USES_GETTIMEOFDAY /* use inaccurate timer, do many loops */
	  for(k=0; k<TIMING_ITER;k++){
#endif    
	    mpfr_set_d(mp_inpt, i1, GMP_RNDN);
//	    inexact = testfun(mp_res, mp_inpt, mpfr_rnd_mode);
            mpfr_subnormalize (mp_res, inexact, mpfr_rnd_mode);
	    result = mpfr_get_d(mp_res, GMP_RNDN);
#ifdef TIMING_USES_GETTIMEOFDAY /* use inaccurate timer, do many loops */
	  }
#endif
	  TBX_GET_TICK(t2);
	}else{
	  TBX_GET_TICK(t1);
#ifdef TIMING_USES_GETTIMEOFDAY /* use inaccurate timer, do many loops */
	  for(k=0; k<TIMING_ITER;k++){
#endif    
	    mpfr_set_d(mp_inpt, i1, GMP_RNDN);
	    mpfr_set_d(mp_inpt2, i2, GMP_RNDN);
//	    inexact = testfun(mp_res, mp_inpt, mp_inpt2, mpfr_rnd_mode);
            mpfr_subnormalize (mp_res, inexact, mpfr_rnd_mode);
	    result = mpfr_get_d(mp_res, GMP_RNDN);
#ifdef TIMING_USES_GETTIMEOFDAY /* use inaccurate timer, do many loops */
	  }
#endif 
	  TBX_GET_TICK(t2);
	}
#endif /*HAVE_MPFR_H*/
      }

      dt = TBX_TICK_RAW_DIFF(t1, t2)-tbx_time;
      if (dt<dtmin)  dtmin=dt;
    }
    *lib_dtsum+=dtmin;
    if (dtmin<*lib_dtmin)  *lib_dtmin=dtmin;
    if (dtmin>*lib_dtmax)  *lib_dtmax=dtmin;
#if      DETAILED_REPORT
    printf("\n input=%1.15e\tT%s=%lld", i1, name, dtmin);
#endif /*DETAILED_REPORT*/
  }

  /* release memory */
#ifdef   HAVE_MPFR_H
  mpfr_clear(mp_inpt2);
  mpfr_clear(mp_inpt);
  mpfr_clear(mp_res);
#endif /*HAVE_MPFR_H*/

}

static void test_without_cache_inter(const char *name, 
			double i1,
			double i2,
			unsigned long long *lib_dtmin,
			unsigned long long *lib_dtmax,
			unsigned long long *lib_dtsum,
			int func_type){
  interval input,result;
  ASSIGN_LOW(input,i1);
  ASSIGN_UP(input,i2);
  unsigned long long dt, dtmin;
  tbx_tick_t   t1, t2; 
  int j;
#ifdef TIMING_USES_GETTIMEOFDAY
  int k;
#endif 
#ifdef HAVE_MPFR_H  
  mpfr_t mp_res, mp_inpt;
  mpfr_t mp_inpt2; /* For the pow function */
  int inexact;

  mpfr_init2(mp_res,  53);
  mpfr_init2(mp_inpt, 53);
  mpfr_init2(mp_inpt2, 53);
#endif

  if(testfun_crlibm_interval!=NULL) { /* test if some functions are missing */
    dtmin=1<<30;
    /* take the min of N1 consecutive calls */
    for(j=0; j<N1; j++) {

      if (func_type == 0){              /* func_type = normal function */
	if (nbarg==1){
	  TBX_GET_TICK(t1);
#ifdef TIMING_USES_GETTIMEOFDAY /* use inaccurate timer, do many loops */
	  for(k=0; k<TIMING_ITER;k++)
#endif
	    result = testfun_crlibm_interval(input);
/*	  printHexa("low res:",LOW(result));
          printHexa("up res:",UP(result));*/
	  TBX_GET_TICK(t2);
	}
      }

      dt = TBX_TICK_RAW_DIFF(t1, t2)-tbx_time;
      if (dt<dtmin)  dtmin=dt;
    }
    *lib_dtsum+=dtmin;
    if (dtmin<*lib_dtmin)  *lib_dtmin=dtmin;
    if (dtmin>*lib_dtmax)  *lib_dtmax=dtmin;
#if      DETAILED_REPORT
    printf("\n input=%1.15e\tT%s=%lld", i1, name, dtmin);
#endif /*DETAILED_REPORT*/
  }

  /* release memory */
#ifdef   HAVE_MPFR_H
  mpfr_clear(mp_inpt2);
  mpfr_clear(mp_inpt);
  mpfr_clear(mp_res);
#endif /*HAVE_MPFR_H*/

}




static void test_worst_case(double (*testfun)(), 
		     double i1, 
		     double i2, 
		     unsigned long long *lib_dtwc, 
		     int func_type){


  double res;
  tbx_tick_t   t1, t2; 
  unsigned long long dtmin, dt;
  int j;
#ifdef TIMING_USES_GETTIMEOFDAY
  int k;
#endif 
#ifdef HAVE_MPFR_H  
  mpfr_t mp_res, mp_inpt;
  mpfr_t mp_inpt2; /* For the pow function */
  int inexact;

  mpfr_init2(mp_res,  53);
  mpfr_init2(mp_inpt, 53);
  mpfr_init2(mp_inpt2, 53);
#endif

  if(testfun!=NULL) { /* test if some functions are missing  */
    dtmin=1<<30;
    for(j=0; j<N1; j++) {
      if (func_type == 0){    /* func_type = normal function */
	if (nbarg==1){
	  TBX_GET_TICK(t1);
#ifdef TIMING_USES_GETTIMEOFDAY
	  for(k=0; k<TIMING_ITER;k++)
#endif
	    res = testfun(i1);
	  TBX_GET_TICK(t2);
	}else{
	  TBX_GET_TICK(t1);
#ifdef TIMING_USES_GETTIMEOFDAY
	  for(k=0; k<TIMING_ITER;k++)
#endif
	    res = testfun(i1,i2);
	  TBX_GET_TICK(t2);
	}
      }else{                   /* func_type = MPFR function  */
#ifdef HAVE_MPFR_H 
	if (nbarg==1){
	  TBX_GET_TICK(t1);
#ifdef TIMING_USES_GETTIMEOFDAY
	  for(k=0; k<TIMING_ITER;k++){
#endif
	    mpfr_set_d (mp_inpt, i1, GMP_RNDN);
	    inexact = testfun (mp_res, mp_inpt, mpfr_rnd_mode);
            mpfr_subnormalize (mp_res, inexact, mpfr_rnd_mode);
	    res = mpfr_get_d (mp_res, GMP_RNDN);
#ifdef TIMING_USES_GETTIMEOFDAY
	  }
#endif
	  TBX_GET_TICK(t2);
	}else{
	  TBX_GET_TICK(t1);
#ifdef TIMING_USES_GETTIMEOFDAY
	  for(k=0; k<TIMING_ITER;k++){
#endif
	    mpfr_set_d (mp_inpt, i1, GMP_RNDN);
	    mpfr_set_d (mp_inpt2, i2, GMP_RNDN);
	    inexact = testfun (mp_res, mp_inpt, mp_inpt2, mpfr_rnd_mode);
            mpfr_subnormalize (mp_res, inexact, mpfr_rnd_mode);
	    res = mpfr_get_d (mp_res, GMP_RNDN);
#ifdef TIMING_USES_GETTIMEOFDAY
	  }
#endif
	  TBX_GET_TICK(t2);
	}
#endif /*HAVE_MPFR_H*/
      }
      dt = TBX_TICK_RAW_DIFF(t1, t2)-tbx_time; 
      if (dt<dtmin)  dtmin=dt;
    }
    *lib_dtwc = dtmin;
  }

  /* release memory */
#ifdef   HAVE_MPFR_H
  mpfr_clear(mp_inpt2);
  mpfr_clear(mp_inpt);
  mpfr_clear(mp_res);
#endif /*HAVE_MPFR_H*/
}



static void normal_output(const char *name,
		   double (*testfun)(),
		   unsigned long long lib_dtmin,
		   unsigned long long lib_dtmax,
		   unsigned long long lib_dtsum,
		   unsigned long long lib_dtwc,
		   int n){
  if(testfun!=NULL) { /* some functions are missing in libultim (cosh, ...  */
    printf("\n%s\nTmin = %lld ticks,\t Tmax = %lld ticks\t avg = %f\tT worst case = %lld\n",
	   name, lib_dtmin, lib_dtmax, ((double)lib_dtsum) / ((double) n), lib_dtwc);
  }
}

static void latex_output(const char *name,
		  double (*testfun)(),
		  unsigned long long lib_dtmin,
		  unsigned long long lib_dtmax,
		  unsigned long long lib_dtsum,
		  unsigned long long lib_dtwc,
		  int n){
    if(testfun!=NULL) { /* some functions are missing in libultim (cosh, ...  */
      if (lib_dtwc > lib_dtmax) lib_dtmax=lib_dtwc;
      printf(" %s  \t& %lld    \t& %10.0f   \t& %lld      \\\\ \n \\hline\n",  
	     name, lib_dtmin, ((double)lib_dtsum) / ((double) n), lib_dtmax);
    }
}


int main (int argc, char *argv[]){ 
  nbarg=1;
  int i, j, n;
  double i1, i2;
  char* rounding_mode;
  char* function_name;
  double worstcase;
  tbx_tick_t   t1, t2; 
  unsigned long long 
    dt,
    libm_dtmin, libm_dtmax, libm_dtsum, libm_dtwc,
    crlibm_dtmin, crlibm_dtmax, crlibm_dtsum, crlibm_dtwc,
    crlibm_inter_dtmin, crlibm_inter_dtmax, crlibm_inter_dtsum, crlibm_inter_dtwc,
    mpfr_dtmin, mpfr_dtmax, mpfr_dtsum, mpfr_dtwc,
    libultim_dtmin, libultim_dtmax, libultim_dtsum, libultim_dtwc,
    libmcr_dtmin, libmcr_dtmax, libmcr_dtsum, libmcr_dtwc;
/*
#ifdef   HAVE_MATHLIB_H
  short Original_Mode;
#endif*/


/*  if ((argc != 3)) usage(argv[0]);*/
  {
    function_name = argv[1];
    sscanf(argv[2],"%d", &n);
    if ((strcmp(function_name,"exp")!=0) && (strcmp(function_name,"j_exp")!=0)  && (strcmp(function_name,"log")!=0) && (strcmp(function_name,"j_log")!=0) && (strcmp(function_name,"expm1")!=0) && (strcmp(function_name,"j_expm1")!=0) && (strcmp(function_name,"log2")!=0) && (strcmp(function_name,"j_log2")!=0))
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
    }
    if ((strcmp(function_name,"log2")==0) || (strcmp(function_name,"j_log2")==0))
    {
      randfun = rand_for_log;
      testfun_crlibm_interval = j_log2;
      testfun_crlibm_low = log2_rd;
      testfun_crlibm_up = log2_ru;
    }
    if ((strcmp(function_name,"exp")==0) || (strcmp(function_name,"j_exp")==0))
    {
      randfun = rand_for_exp_perf;
      testfun_crlibm_interval = j_exp;
      testfun_crlibm_low = exp_rd;
      testfun_crlibm_up = exp_ru;
    }
    if ((strcmp(function_name,"expm1")==0) || (strcmp(function_name,"j_expm1")==0))
    {
      randfun = rand_for_expm1_soaktest;
      testfun_crlibm_interval = j_expm1;
      testfun_crlibm_low = expm1_rd;
      testfun_crlibm_up = expm1_ru;
    }

  crlibm_init();
  

  crlibm_dtmin=1<<30; crlibm_dtmax=0; crlibm_dtsum=0; crlibm_dtwc=0;
  crlibm_inter_dtmin=1<<30; crlibm_inter_dtmax=0; crlibm_inter_dtsum=0; crlibm_inter_dtwc=0;
  libm_dtmin=1<<30;   libm_dtmax=0;   libm_dtsum=0; libm_dtwc=0;
  libultim_dtmin=1<<30;    libultim_dtmax=0;    libultim_dtsum=0; libultim_dtwc=0;
  libmcr_dtmin=1<<30;    libmcr_dtmax=0;    libmcr_dtsum=0; libmcr_dtwc=0;
  mpfr_dtmin=1<<30;   mpfr_dtmax=0;   mpfr_dtsum=0; mpfr_dtwc=0;


  /* take the min of N1 consecutive calls */
  tbx_time=1<<30;
  for(j=0; j<1000*N1; j++) {
    TBX_GET_TICK(t1);
    TBX_GET_TICK(t2);
    dt = TBX_TICK_RAW_DIFF(t1, t2);
    if(dt<tbx_time) tbx_time = dt;
  }
  printf("tbx_time=%llu\n", tbx_time);

  /************  TESTS WITHOUT CACHES  *******************/
  srandom(n);
#if EVAL_PERF==1  
  crlibm_second_step_taken=0; 
#endif

  /* take the min of N1 identical calls to leverage interruptions */
  /* As a consequence, the cache impact of these calls disappear...*/
  for(i=0; i< n; i++){ 
    i1 = randfun();
    i2 = randfun();
    if(!(i1<=i2))
    {
      double temp=i1;
      i1=i2;
      i2=temp;
    }
    /*    db_number ia,ib;
    ia.i[HI]=0x31100afb;
    ia.i[LO]=0x198a95fe;
    ib.i[HI]=0x3d42897b;
    ib.i[LO]=0x84591a4e;
    i1=ia.d; i2=ib.d;*/
    test_without_cache("crlibm", i1, i2, &crlibm_dtmin, &crlibm_dtmax, &crlibm_dtsum, 0);
    test_without_cache_inter("crlibm interval", i1, i2, &crlibm_inter_dtmin, &crlibm_inter_dtmax, &crlibm_inter_dtsum, 0);
    
  } 

#if EVAL_PERF==1  
#ifdef TIMING_USES_GETTIMEOFDAY /* use inaccurate timer, do many loops */
	 printf("\nCRLIBM : Second step taken %d times out of %d\n",
		crlibm_second_step_taken/(N1 * TIMING_ITER), n );
#else
	 printf("\nCRLIBM : Second step taken %d times out of %d\n",
		crlibm_second_step_taken/N1, n );
#endif

#endif


  /************  WORST CASE TESTS   *********************/
  /* worst case test */
  i1 = worstcase;
  i2 = 0 ; /* TODO when we have worst cases for power...*/

//  test_worst_case(log_rd(), i1, i2, &libm_dtwc, 0);
  test_worst_case(log_rd, i1, i2, &crlibm_dtwc, 0);
/*#ifdef   HAVE_MPFR_H
  test_worst_case((double(*)())testfun_mpfr, i1, i2, &mpfr_dtwc, 1);
#endif*/ /*HAVE_MPFR_H*/
/*#ifdef   HAVE_MATHLIB_H
  Original_Mode = Init_Lib(); 
  test_worst_case(testfun_libultim, i1, i2, &libultim_dtwc, 0);
  Exit_Lib(Original_Mode);
#endif*/ /*HAVE_MATHLIB_H*/
/*#ifdef   HAVE_LIBMCR_H
  test_worst_case(testfun_libmcr, i1, i2, &libmcr_dtwc, 0);
#endif*/ /*HAVE_LIBMCR_H*/

    /*************Normal output*************************/
//  normal_output("LIBM", testfun_libm, libm_dtmin, libm_dtmax, libm_dtsum, libm_dtwc, n);
  normal_output("CRLIBM", testfun_crlibm_low, crlibm_dtmin, crlibm_dtmax, crlibm_dtsum, crlibm_dtwc, n);
  normal_output("CRLIBM interval", testfun_crlibm_interval, crlibm_inter_dtmin, crlibm_inter_dtmax, crlibm_inter_dtsum, crlibm_inter_dtwc, n);  


  /******************* Latex output ****************/
  printf("\\multicolumn{4}{|c|}{Processor / system / compiler}   \\\\ \n \\hline");
  printf("\n                             & min time \t & avg time \t& max time \t  \\\\ \n \\hline\n");
  latex_output("\\texttt{crlibm}        ", testfun_crlibm_low, crlibm_dtmin, crlibm_dtmax, crlibm_dtsum, crlibm_dtwc, n);
  latex_output("\\texttt{crlibm inter}        ", testfun_crlibm_interval, crlibm_inter_dtmin, crlibm_inter_dtmax, crlibm_inter_dtsum, crlibm_inter_dtwc, n);  
  }
  return 0;
}


