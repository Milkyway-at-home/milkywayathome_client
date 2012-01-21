#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "crlibm.h"
#include "crlibm_private.h"
#include "test_common.h"

#ifdef HAVE_MPFR_H
#include <gmp.h>
#include <mpfr.h>
#endif

#ifdef HAVE_MATHLIB_H
#include <MathLib.h>
#endif

#ifdef HAVE_LIBMCR_H
#include <libmcr.h>
#endif

#if defined(_MSC_VER) || defined(__MINGW32__)
  #define random rand
#endif



#define RN 1
#define RU 2
#define RD 3
#define RZ 4  

/* A variable equal to zero, stored here so that the compiler doesn't
   know its value in the other functions, which allows to prevent some
   optimizations  */


double zero ;

/* Here come the various random number generators. They all use the
   rand() function.  

   We may have two rand functions for each function under
   test. The first is for the soaktest, the second for testing the
   performance under what is supposed the main domain of use the function. 

   Typical examples:

   log has identical functions for soaktest and performance: random
   positive numbers. This means that negative numbers are not tested
   by soaktest, though.

   sin soaktests on all the floats, but tests for perf on a small
   interval around zero, shamelessely chosen as the one on which crlibm
   is the fastest.
*/


/**/


/* Return 'sizeof(int)' random bits    */
int rand_int(){
  int val;
  int i;
  val = (random() & 0x000000ff);
  for(i=0; i<(sizeof(int)); i++){
    val = val << 8;
    val += (random() & 0x000000ff ); /* we keep only 8 bits */
  }
  return val;
}




/* Return a completely random double  */

double rand_generic(){
  db_number result;
  
  result.i[LO]=rand_int();
  result.i[HI]=rand_int();
 
  return result.d;
}


/* Return a random double between 0 and 1, with a normal law on the
   exponent  */

double rand_double(){
  db_number result;
  int e;
  /*first the low bits of the mantissa*/
  result.i[LO]=rand_int();
  /* then the high bits of the mantissa, and the sign bit */
  result.i[HI]=  rand_int() & 0x000fffff;
  /* Now set the exponent (negative value) */
  e = rand() & 0x000003ff; 
  if (e>0) e-=1;
  result.i[HI] += e<<20;
  return (result.d);
}


/* Return a random double between 1 and 2, with a normal law on the
   mantissa and a constant exponent  */

double rand_double_normal(){
  db_number result;
  int e;
  /*first the low bits of the mantissa*/
  result.i[LO]=rand_int();
  /* then the high bits of the mantissa, and the sign bit */
  result.i[HI]=  rand_int() & 0x000fffff;
  /* Now set the exponent */
  e = 1023; 
  result.i[HI] += e<<20;
  return (result.d);
}



/* For exp we will test perf on numbers with a random sign, a random mantissa, and
   a random exponent between -9 and 9. And we soaktest on all the doubles */

#define rand_for_exp_soaktest rand_generic

double rand_for_exp_perf(){
  db_number result;
  int e;

  /*first the low bits of the mantissa*/
  result.i[LO]=rand_int();
  /* then the high bits of the mantissa, and the sign bit */
  result.i[HI]=  rand_int() & 0x800fffff;
  /* Now set the exponent between -9 and 9, enough to cover the useful range  */
  e =  (int) ( (rand_double_normal()-1) * 18 );
  result.i[HI] += (1023 + e -9)<<20;
  return result.d;
}


/* a number in the range which never produces over/underflow for the
   exp function. I don't trust the randomness of this function */
double rand_for_exp_normal(){
  return((750+710)*(rand_double_normal()-1)-750);
}


/* Return a number between -2^10 and +2^10 leaving out exponents <= -60
   Attention: no normal law
*/
double rand_for_expm1_soaktest() {
  db_number resdb;
  int expo;

  resdb.i[HI] = rand_int() & 0x000fffff;
  resdb.i[LO] = rand_int();

  while ((expo = (rand_int() & 0x0000007f)) > 70);
  expo -= 60;

  resdb.i[HI] |= (expo + 1023) << 20;

  resdb.i[HI] |= (rand_int() & 0x80000000);

  return resdb.d;
}


#define rand_for_expm1_testperf rand_for_expm1_soaktest

#define rand_for_csh_perf rand_for_exp_perf
/* I wish we could soaktest using rand_generic, but current MPFR is
   very slow for small and large arguments (up to a few minutes per
   call !). To check regularly, this is bound to improve. */

#define rand_for_csh_soaktest  rand_for_exp_perf

/* For log we only test the positive numbers*/
double rand_for_log(){
  db_number result;

  /*first the low bits of the mantissa*/
  result.i[LO]=rand_int();
  /* then the high bits of the mantissa, and the sign bit */
  result.i[HI]=  rand_int() & 0x7fffffff;
  /* printf("x = %1.5e\n", result.d);*/
#if 0 /* to avoid NaNs and subnormals which kill fi_lib */ 
  if(0x7ff==result.i[HI]>>20  ||  0==result.i[HI]>>20) return rand_for_log();
#endif
  return result.d;
}

/* Produce a domain ]-1;+infty[ 

   Attention: we have no normal law on the FP results.

*/
double rand_for_log1p() {
  db_number result;
  /* int e; */
  
  result.d = rand_for_log();

  /* If you find the MPFR errors in directed rounding modes for small values
     annoying replace the previous line by the following ones: */

  /*
  result.i[LO] = rand_int();
  result.i[HI] = 0x000fffff & rand_int();
  
  e = 0x3c9 + (0x3ff & rand_int());

  result.i[HI] |= e << 20;
  */

  if (result.d < 1.0) {
    if (rand_int() && 0x1) 
      result.i[HI] = 0x80000000 | result.i[HI];
  }

  return result.d;
}


/* For trigonometric functions it is difficult to tell what the test
   function should be. We chose an exponent between -20 and 40... */ 

double rand_for_trig_perf(){
  db_number result;
  int e;
  /*first the low bits of the mantissa*/
  result.i[LO]=rand_int();
  /* then the high bits of the mantissa, and the sign bit */
  result.i[HI]=  rand_int() & 0x800fffff;
  /* Now set the exponent  between -20 and 40 */
    e =  (int) ( (rand_double_normal()-1) * 60 ); 
   result.i[HI] += (1023 + e -20)<<20;
  return result.d;

}


/* For trigpi functions there is nobody to compare to, so who cares... */ 
#define rand_for_trigpi_soaktest rand_generic


#define rand_for_trigpi_perf rand_generic
/* Actually we might exhibit better results if we spent less time near
   zero */ 



#if 0
#define rand_for_trig_soaktest rand_generic
#else
#define rand_for_trig_soaktest rand_for_trig_perf
#endif


#define PIH 3.14159265358979311599796346854418516159057617187500e+00


double rand_for_atan_perf(){
  db_number result;
  int e;
  /*first the low bits of the mantissa*/
  result.i[LO]=rand_int();
  /* then the high bits of the mantissa, and the sign bit */
  result.i[HI]=  rand_int() & 0x800fffff;
  /* Now set the exponent between -20 and 50, enough to cover the useful range  */
  e =  (int) ( (rand_double_normal()-1) * 70 );
  result.i[HI] += (1023 + e -20)<<20;
  return result.d;

}

double rand_for_atan_soaktest(){
  db_number result;
  int e;

  /*first the low bits of the mantissa*/
  result.i[LO]=rand_int();
  /* then the high bits of the mantissa, and the sign bit */
  result.i[HI]=  rand_int() & 0x800fffff;
  /* Now set the exponent between -20 and 50, enough to cover the useful range  */
  e =  (int) ( (rand_double_normal()-1) * 70 );
  result.i[HI] += (1023 + e -20)<<20;
  return result.d;

}





/* For power, we must generate two values x and y
   The function returns x and sets y by a sideeffect
   We must cast the function correctly.

   We generate x and y such that over- or
   underflow is unlikely
   The value x is always positive.

   We will have no idea on the randomness

   We generate first x and z on the whole 
   positive double range. 
   We compute than y1 = log_rn(z) /arith log_rn(x)

   We generate than eps in the range [0;2^(-32)]
   If z > 1 then we return y2 = y1 -arith y1 *arith eps
            else we return y2 = y1 +arith y1 *arith eps

*/
double rand_for_pow_perf(double *yr){
  double x,z,y1,y2,logz,logx,delta;
  db_number tempdb;

  tempdb.i[LO] = rand_int();
  tempdb.i[HI] = rand_int() & 0x7fffffff;
  x = tempdb.d;

  tempdb.i[LO] = rand_int();
  tempdb.i[HI] = rand_int() & 0x7fffffff;
  z = tempdb.d;

  logx = log_rn(x);
  logz = log_rn(z);

  y1 = logx / logz;
 
  tempdb.i[LO] = rand_int();
  tempdb.i[HI] = (rand_int() & 0x000fffff) | ((1023 - (rand_int() & 0x1f)) << 20);
 
  delta = y1 * tempdb.d;

  if (z > 1.0) y2 = y1 - delta; else y2 = y1 + delta;
 
  *yr = y2;
  return x;
}

#define rand_for_pow_soaktest rand_for_pow_perf



/* This #if selects among two possible policies for perf test. 
   Default is 0,  the less favorable to crlibm. */
#if 0
/* Produces x with |x| in [2^(-31);1] and a random sign This passes
 through all the interesting parts of any implementation (for x
 smaller than 2^-31 one returns x) but most tested numbers will be
 very small, which is does probably not correspond to real use
 cases. */
double rand_for_asin_testperf(){
  db_number result;
  int e;
  /*first the low bits of the mantissa*/
  result.i[LO]=rand_int();
  /* then the high bits of the mantissa */
  result.i[HI]=  rand_int() & 0x000fffff;
  /* Now set the exponent */
  e = (rand_int() & 0x0000001f) + 0x000003e0;
  if (e == 0x000003ff) {
    result.i[HI] = 0;
    result.i[LO] = 0;
  }
  /* Now set the sign */
  e += (rand_int() & 0x00000800);
  result.i[HI] += e<<20;
  return (result.d);
}

/* Produces x with |x| in [2^(-127);1] and a random sign. Same comment
   as for asin */
double rand_for_acos_testperf(){
  db_number result;
  int e;
  /*first the low bits of the mantissa*/
  result.i[LO]=rand_int();
  /* then the high bits of the mantissa */
  result.i[HI]=  rand_int() & 0x000fffff;
  /* Now set the exponent */
  e = (rand_int() & 0x0000007f) + 0x00000380;
  if (e == 0x000003ff) {
    result.i[HI] = 0;
    result.i[LO] = 0;
  }
  /* Now set the sign */
  e += (rand_int() & 0x00000800);
  result.i[HI] += e<<20;
  return (result.d);
}
#else
/* normal number between -1 and 1 . This can be discussed */
double rand_for_asin_testperf(){
  return (rand_double_normal()-1.5) * 2.0;
}
#define rand_for_acos_testperf rand_for_asin_testperf
#endif



/* Produces x with |x| in [0;1] */
double rand_for_asin_soaktest(){
  db_number result;
  int e;
  /*first the low bits of the mantissa*/
  result.i[LO]=rand_int();
  /* then the high bits of the mantissa */
  result.i[HI]=  rand_int() & 0x000fffff;
  /* Now set the exponent */
  e = (rand_int() & 0x000003ff);
  if (e == 0x000003ff) {
    result.i[HI] = 0;
    result.i[LO] = 0;
  }
  /* Now set the sign */
  e += (rand_int() & 0x00000800);
  result.i[HI] += e<<20;
  return (result.d);
}




/*******************************************************************/
/*             A few tinkered functions for comparisons            */


/* These ones are mostly inaccurate, for perf comparison only  */

double tinkered_sinpi (double x) {
  return sin(PIH*x);
}
double tinkered_cospi (double x) {
  return cos(PIH*x);
}
double tinkered_tanpi (double x) {
  return tan(PIH*x);
}
double tinkered_atanpi (double x) {
  return atan(x)/PIH;
}
double tinkered_asinpi (double x) {
  return asin(x)/PIH;
}
double tinkered_acospi (double x) {
  return acos(x)/PIH;
}


#ifdef HAVE_MPFR_H
/* we might even attempt to prove some day that these are correct for
   RN mode. We compute pi*x on enough bits to avoid any
   cancellation in the argument reduction. 

   Anyway someday the mpfr team will write real ones.

   These functions are probably wrong sometimes in directed rounding
   modes, although this has improved thanks to Rob Clark and Guillaume
   Hanrot.  

   Useful mostly for soaktesting the RN trigpi functions.

*/ 

int tinkered_mpfr_sinpi (mpfr_t mpr, mpfr_t mpx, mp_rnd_t rnd) {
  mpfr_t pi, pix;
  mpfr_init2(pi,  2000);
  mpfr_init2(pix, 2060);

  /* Test for exact cases */
  mpfr_mul_2exp(mpr, mpx, 1, GMP_RNDN); 
  if(mpfr_integer_p(mpr)) { /* Exact cases */ 
    mpfr_div_2exp(mpr, mpx, 1, GMP_RNDN); 
    mpfr_frac(mpr, mpr, GMP_RNDN);
    double d = mpfr_get_d(mpr, GMP_RNDN); /* -1/4, 0, 1/2, 3/4 */
    if (d == 0.0)       { mpfr_set_si(mpr, 0, rnd); }
    else if (d == 0.25) { mpfr_set_si(mpr, 1, rnd); }
    else if (d == 0.5)  { mpfr_set_si(mpr, 0, rnd); }
    else                { mpfr_set_si(mpr, -1, rnd); } 
    return 0; 
  }
  
  mpfr_const_pi(pi,  GMP_RNDN);
  mpfr_mul(pix, pi, mpx, GMP_RNDN);
  mpfr_sin(mpr, pix, rnd);
  mpfr_clear(pi);
  mpfr_clear(pix);
  return 0;
}

int tinkered_mpfr_cospi (mpfr_t mpr, mpfr_t mpx, mp_rnd_t rnd) {
  mpfr_t pi, pix;
  mpfr_init2(pi,  2000);
  mpfr_init2(pix, 2060);

  /* Test for exact cases */
  mpfr_mul_2exp(mpr, mpx, 1, GMP_RNDN); 
  if(mpfr_integer_p(mpr)) { /* Exact cases */ 
    mpfr_div_2exp(mpr, mpx, 1, GMP_RNDN); 
    mpfr_frac(mpr, mpr, GMP_RNDN);
    double d = mpfr_get_d(mpr, GMP_RNDN); /* -1/4, 0, 1/2, 3/4 */
    if (d == 0.0)       { mpfr_set_si(mpr, 1, rnd); }
    else if (d == 0.25) { mpfr_set_si(mpr, 0, rnd); }
    else if (d == 0.5)  { mpfr_set_si(mpr, -1, rnd); }
    else                { mpfr_set_si(mpr, 0, rnd); } 
    return 0; 
  }

  mpfr_const_pi(pi,  GMP_RNDN);
  mpfr_mul(pix, pi, mpx, GMP_RNDN);
  mpfr_cos(mpr, pix, rnd);
  mpfr_clear(pi);
  mpfr_clear(pix);
  return 0;
}

int tinkered_mpfr_tanpi (mpfr_t mpr, mpfr_t mpx, mp_rnd_t rnd) {
  mpfr_t pi, pix;
  mpfr_init2(pi,  2000);
  mpfr_init2(pix, 2060);

  /* Test for exact cases */
  mpfr_mul_2exp(mpr, mpx, 2, GMP_RNDN); 
  if(mpfr_integer_p(mpr)) { /* Exact cases */ 
    mpfr_frac(mpr, mpx, GMP_RNDN); 
    double d = mpfr_get_d(mpr, GMP_RNDN); /* -1/4, 0, 1/4, 1/2 */
    if (d == 0.0)       { mpfr_set_ui(mpr, 0, rnd); }
    else if (d == 0.25) { mpfr_set_ui(mpr, 1, rnd); }
    else if (d == 0.5)  { mpfr_set_nan(mpr); }
    else                { mpfr_set_si(mpr, -1, rnd); } 
    return 0; 
  }
  
  /* Otherwise */
  mpfr_const_pi(pi,  GMP_RNDN);
  mpfr_mul(pix, pi, mpx, GMP_RNDN);
  mpfr_tan(mpr, pix, rnd);
  mpfr_clear(pi);
  mpfr_clear(pix);
  return 0;
}

int tinkered_mpfr_atanpi (mpfr_t mpr, mpfr_t mpx, mp_rnd_t rnd) {
  mpfr_t pi, at;
  mpfr_init2(pi,  2000);
  mpfr_init2(at, 2060);

  mpfr_const_pi(pi,  GMP_RNDN);
  mpfr_atan(at, mpx, GMP_RNDN);
  mpfr_div(mpr, at, pi,  rnd);
  mpfr_clear(pi);
  mpfr_clear(at);
  return 0;
}

int tinkered_mpfr_acospi (mpfr_t mpr, mpfr_t mpx, mp_rnd_t rnd) {
  mpfr_t pi, pix;
  mpfr_init2(pi,  2000);
  mpfr_init2(pix, 2000);

  mpfr_const_pi(pi,  GMP_RNDN);
  mpfr_acos(pix, mpx, GMP_RNDN);
  mpfr_div(mpr, pix, pi, rnd);
  mpfr_clear(pi);
  mpfr_clear(pix);
  return 0;
}

int tinkered_mpfr_asinpi (mpfr_t mpr, mpfr_t mpx, mp_rnd_t rnd) {
  mpfr_t pi, pix;
  mpfr_init2(pi,  2000);
  mpfr_init2(pix, 2000);

  mpfr_const_pi(pi,  GMP_RNDN);
  mpfr_asin(pix, mpx, GMP_RNDN);
  mpfr_div(mpr, pix, pi, rnd);
  mpfr_clear(pi);
  mpfr_clear(pix);
  return 0;
}

#endif



void test_rand()  {
  int i;
  double min=1e300, max=0.0;
  db_number input;
  for(i=0; i< 1000; i++){ 
    input.d = rand_for_exp_perf();
    if (input.d<min) min=input.d;
    if (input.d>max) max=input.d;
    printf("%1.5ex \t%.8x %.8x\t%1.5e \t%1.5e\n", input.d, input.i[HI], input.i[LO],min,max );
  }
}





/* general init function */


void test_init(/* pointers to returned value */
	       double (**randfun_perf)(), 
	       double (**randfun_soaktest)(), 
	       double (**testfun_crlibm)(), 
	       int    (**testfun_mpfr)  (),
	       double (**testfun_libultim)   (),
	       double (**testfun_libmcr)  (),
	       double (**testfun_libm)  (),
	       double* worst_case,
	       /* arguments */
	       char *func_name,
	       char *rnd_mode)  {

  int crlibm_rnd_mode;

  /* We have added the rounding mode designation used in libmcr's test files */
  if      ((strcmp(rnd_mode,"RU")==0) || (strcmp(rnd_mode,"P")==0)) crlibm_rnd_mode = RU;
  else if ((strcmp(rnd_mode,"RD")==0) || (strcmp(rnd_mode,"M")==0)) crlibm_rnd_mode = RD;
  else if ((strcmp(rnd_mode,"RZ")==0) || (strcmp(rnd_mode,"Z")==0)) crlibm_rnd_mode = RZ;
  else if ((strcmp(rnd_mode,"RN")==0) || (strcmp(rnd_mode,"N")==0)) crlibm_rnd_mode = RN;
  else {
    fprintf(stderr, "Unknown rounding mode: %s, exiting\n", rnd_mode);
    exit(EXIT_FAILURE);
  }


  *randfun_perf     = rand_generic; /* the default random function */
  *randfun_soaktest = rand_generic; /* the default random function */
  *testfun_mpfr     = NULL;
  *testfun_libm     = NULL;
  *worst_case=0.;

  if (strcmp (func_name, "exp") == 0)
    {
      *randfun_perf     = rand_for_exp_perf;
      *randfun_soaktest = rand_for_exp_soaktest;
      *worst_case= .75417527749959590085206221024712557043923055744016892276704311370849609375e-9;
      *testfun_libm   = exp;
      switch(crlibm_rnd_mode){
      case RU:
	*testfun_crlibm = exp_ru;	break;
      case RD:
	*testfun_crlibm = exp_rd;	break;
      case RZ:
	*testfun_crlibm = exp_rz;	break;
      default:
	*testfun_crlibm = exp_rn;
      }
#ifdef HAVE_MATHLIB_H
      *testfun_libultim    = uexp;
#endif
#ifdef HAVE_LIBMCR_H
      *testfun_libmcr    = __libmcr_exp;
#endif
#ifdef HAVE_MPFR_H
      *testfun_mpfr   = mpfr_exp;
#endif
    }


  else  if (strcmp (func_name, "log") == 0)
    {
      *randfun_perf     = rand_for_log;
      *randfun_soaktest = rand_for_log;
      *worst_case=0.4009793462309855760053830468258630076242931610568335144339734234840014178511334897967240437927437320e-115;
      *testfun_libm   = log;
      switch(crlibm_rnd_mode){
      case RU:
	*testfun_crlibm = log_ru;	break;
      case RD:
	*testfun_crlibm = log_rd;	break;
      case RZ:
	*testfun_crlibm = log_rz;	break;
      default:
	*testfun_crlibm = log_rn;
      }
#ifdef HAVE_MATHLIB_H
      *testfun_libultim    = ulog;
#endif
#ifdef HAVE_LIBMCR_H
      *testfun_libmcr    = __libmcr_log;
#endif
#ifdef HAVE_MPFR_H
      *testfun_mpfr   = mpfr_log;
#endif
    }

  else  if (strcmp (func_name, "log2") == 0)
    {
      *randfun_perf     = rand_for_log;
      *randfun_soaktest = rand_for_log;
      *worst_case= 1.706724408218747379706314859504e+00; 
      *testfun_libm   = log2; 
      switch(crlibm_rnd_mode){
      case RU:
	*testfun_crlibm = log2_ru;	break;
      case RD:
	*testfun_crlibm = log2_rd;	break;
      case RZ:
	*testfun_crlibm = log2_rz;	break;  
      default:
	*testfun_crlibm = log2_rn;
      }
#ifdef HAVE_MATHLIB_H
      *testfun_libultim    = ulog2;
#endif
#ifdef HAVE_LIBMCR_H
      *testfun_libmcr    = NULL;
#endif
#ifdef HAVE_MPFR_H
      *testfun_mpfr   = mpfr_log2;
#endif
    }

  else  if (strcmp (func_name, "log10") == 0)
    {
      *randfun_perf     = rand_for_log;
      *randfun_soaktest = rand_for_log;
      *worst_case = 2.60575359533670695497442621444820894404798523211076e+129;
      *testfun_libm   = log10; 
      switch(crlibm_rnd_mode){
      case RU:
	*testfun_crlibm = log10_ru;	break;
      case RD:
	*testfun_crlibm = log10_rd;	break;
      case RZ:
	*testfun_crlibm = log10_rz;	break;  
      default:
	*testfun_crlibm = log10_rn;
      }

      /*
#ifdef HAVE_MATHLIB_H
      *testfun_libultim    = ulog10;
#endif

      */

#ifdef HAVE_LIBMCR_H
      *testfun_libmcr    = NULL;
#endif
#ifdef HAVE_MPFR_H
      *testfun_mpfr   = mpfr_log10;
#endif
    }

  else  if (strcmp (func_name, "expm1") == 0)
    {
      *randfun_perf     = rand_for_expm1_testperf;
      *randfun_soaktest = rand_for_expm1_soaktest;
      *worst_case= 1.312678236466234442186232916905e-07; /* worst case for RN only, there are much worse cases for RZ */
      *testfun_libm   = expm1;
      switch(crlibm_rnd_mode){
      case RU:
	*testfun_crlibm = expm1_ru;	break;
      case RD:
	*testfun_crlibm = expm1_rd;	break;
      case RZ:
	*testfun_crlibm = expm1_rz;	break;
      default:
	*testfun_crlibm = expm1_rn;
      }
#ifdef HAVE_MATHLIB_H
      *testfun_libultim    = NULL; /* TODO */
#endif
#ifdef HAVE_LIBMCR_H
      *testfun_libmcr    = NULL; /* TODO */
#endif
#ifdef HAVE_MPFR_H
      *testfun_mpfr   = mpfr_expm1;
#endif
    }

  else if (strcmp (func_name, "log1p") == 0)
    {
      *randfun_perf     = rand_for_log1p;
      *randfun_soaktest = rand_for_log1p;
      *worst_case= 1.332267629550187256862679085950e-15;
      *testfun_libm   = log1p;
      switch(crlibm_rnd_mode){
      case RU:
	*testfun_crlibm = log1p_ru;	break;
      case RD:
	*testfun_crlibm = log1p_rd;	break;
      case RZ:
	*testfun_crlibm = log1p_rz;	break;
      default:
	*testfun_crlibm = log1p_rn;
      }
#ifdef HAVE_MATHLIB_H
      *testfun_libultim    = NULL;
#endif
#ifdef HAVE_LIBMCR_H
      *testfun_libmcr    = NULL;
#endif
#ifdef HAVE_MPFR_H
      *testfun_mpfr   = mpfr_log1p;
#endif
    }

  else if (strcmp (func_name, "sin") == 0)
    {
      *randfun_perf     = rand_for_trig_perf;
      *randfun_soaktest = rand_for_trig_soaktest;
      *worst_case= 9.24898516520941904595076721307123079895973205566406e-01;
      *testfun_libm   = sin;
      switch(crlibm_rnd_mode){
      case RU:
	*testfun_crlibm = sin_ru;	break;
      case RD:
	*testfun_crlibm = sin_rd;	break;
      case RZ:
	*testfun_crlibm = sin_rz;	break;
      default:
	*testfun_crlibm = sin_rn;
      }
#ifdef HAVE_MATHLIB_H
      *testfun_libultim    = usin;
#endif
#ifdef HAVE_LIBMCR_H
      *testfun_libmcr    = __libmcr_sin;
#endif
#ifdef HAVE_MPFR_H
      *testfun_mpfr   = mpfr_sin;
#endif
    }

  else  if (strcmp (func_name, "cos") == 0)
    {
      *randfun_perf     = rand_for_trig_perf;
      *randfun_soaktest = rand_for_trig_soaktest;
      *worst_case=  8.87406081479789610177988379291491582989692687988281e-01;
      *testfun_libm   = cos;
      switch(crlibm_rnd_mode){
      case RU:
	*testfun_crlibm = cos_ru;	break;
      case RD:
	*testfun_crlibm = cos_rd;	break;
      case RZ:
	*testfun_crlibm = cos_rz;	break;
      default:
	*testfun_crlibm = cos_rn;
      }
#ifdef HAVE_MATHLIB_H
      *testfun_libultim    = ucos;
#endif
#ifdef HAVE_LIBMCR_H
      *testfun_libmcr    = __libmcr_cos;
#endif
#ifdef HAVE_MPFR_H
      *testfun_mpfr   = mpfr_cos;
#endif
    }

  else  if (strcmp (func_name, "tan") == 0)
    {
      *randfun_perf     = rand_for_trig_perf;
      *randfun_soaktest = rand_for_trig_soaktest;
      *worst_case= 1.18008664944477814628953638020902872085571289062500e-01;
      *testfun_libm   = tan; 
      switch(crlibm_rnd_mode){
      case RU:
	*testfun_crlibm = tan_ru;	break;
      case RD:
	*testfun_crlibm = tan_rd;	break;
      case RZ:
	*testfun_crlibm = tan_rz;	break;
      default:
	*testfun_crlibm = tan_rn;
      }
#ifdef HAVE_MATHLIB_H
      *testfun_libultim    = utan;
#endif
#ifdef HAVE_LIBMCR_H
      *testfun_libmcr    = __libmcr_tan;
#endif
#ifdef HAVE_MPFR_H
      *testfun_mpfr   = mpfr_tan;
#endif
    }

#if 0 /* No cotan in the standard math.h ? */
  else  if (strcmp (func_name, "cotan") == 0)
    {
      *randfun_perf     = rand_for_trig_perf;
      *randfun_soaktest = rand_for_trig_soaktest;
      *worst_case= 1.18008664944477814628953638020902872085571289062500e-01;
      *testfun_libm   = cotan; 
      switch(crlibm_rnd_mode){
      case RU:
	*testfun_crlibm = cotan_ru;	break;
      case RD:
	*testfun_crlibm = cotan_rd;	break;
      case RZ:
	*testfun_crlibm = cotan_rz;	break;
      default:
	*testfun_crlibm = cotan_rn;
      }
#ifdef HAVE_MATHLIB_H
      *testfun_libultim    = ucotan;
#endif
#ifdef HAVE_MPFR_H
      *testfun_mpfr   = mpfr_cotan;
#endif
    }
#endif /* no cotan */

  else if (strcmp (func_name, "atan") == 0)
    {
      *randfun_perf     = rand_for_atan_perf;
      *randfun_soaktest = rand_for_atan_soaktest;
      *worst_case= 9.54714164331460501955461950274184346199035644531250e-02; 
      *testfun_libm   = atan;
      switch(crlibm_rnd_mode){
      case RU:
	*testfun_crlibm = atan_ru;	break;
      case RD:
	*testfun_crlibm = atan_rd;	break;
      case RZ:
	*testfun_crlibm = atan_rz;	break;
      default:
        *testfun_crlibm = atan_rn ;
      }
#ifdef HAVE_MATHLIB_H
      *testfun_libultim    = uatan;
#endif
#ifdef HAVE_LIBMCR_H
      *testfun_libmcr    = __libmcr_atan;
#endif
#ifdef HAVE_MPFR_H
      *testfun_mpfr   = mpfr_atan;
#endif
    }

  else if (strcmp (func_name, "atanpi") == 0)
    {
      *randfun_perf     = rand_for_atan_perf;
      *randfun_soaktest = rand_for_atan_soaktest;
      *worst_case= 0 /* TODO */; 
      *testfun_libm   = tinkered_atanpi;
      switch(crlibm_rnd_mode){
      case RU:
	*testfun_crlibm = atanpi_ru;	break;
      case RD:
	*testfun_crlibm = atanpi_rd;	break;
      case RZ:
	*testfun_crlibm = atanpi_rz;	break;
      default:
        *testfun_crlibm = atanpi_rn ;
      }
#ifdef HAVE_MPFR_H
      *testfun_mpfr   = tinkered_mpfr_atanpi;
#endif
    }



  else if (strcmp (func_name, "cosh") == 0)
    {
      *randfun_perf     = rand_for_csh_perf;
      *randfun_soaktest = rand_for_csh_soaktest;
      *worst_case= 3.76323248339103422210882854415103793144226074218750e+00;
      *testfun_libm   = cosh;
      switch(crlibm_rnd_mode){
      case RU:
	*testfun_crlibm = cosh_ru;	break;
      case RD:
	*testfun_crlibm = cosh_rd;	break;
      case RZ:
	*testfun_crlibm = cosh_rz;	break;
      default:
	*testfun_crlibm = cosh_rn;
      }
#ifdef HAVE_MATHLIB_H
      /* No hyperbolic function in Ziv library */ 
      *testfun_libultim    = NULL;
#endif
#ifdef HAVE_MPFR_H
      *testfun_mpfr   = mpfr_cosh;
#endif
    }

  else if (strcmp (func_name, "sinh") == 0)
    {
      *randfun_perf     = rand_for_csh_perf;
      *randfun_soaktest = rand_for_csh_soaktest;
      *worst_case= 5.81191276791475441854117889306508004665374755859375;
      *testfun_libm   = sinh;
      switch(crlibm_rnd_mode){
      case RU:
	*testfun_crlibm = sinh_ru;	break;
      case RD:
	*testfun_crlibm = sinh_rd;	break;
      case RZ:
	*testfun_crlibm = sinh_rz; 	break;
      default:
	*testfun_crlibm = sinh_rn;
      }
#ifdef HAVE_MATHLIB_H
      /* No hyperbolic function in Ziv library */ 
      *testfun_libultim    = NULL;
#endif
#ifdef HAVE_MPFR_H
      *testfun_mpfr   = mpfr_sinh;
#endif
    } 
  else if (strcmp (func_name, "asin") == 0)
    {
      *randfun_perf     = rand_for_asin_testperf;
      *randfun_soaktest = rand_for_asin_soaktest;
      *worst_case = 5.6018392194160160357796485186554491519927978515625000000000000000000000e-01;  
                                /* Provokes high timing on Ziv's lib; might not be the real worst case */
      *testfun_libm   = asin;
      switch(crlibm_rnd_mode){
      case RU:
	*testfun_crlibm = asin_ru;	break;
      case RD:
	*testfun_crlibm = asin_rd;	break;
      case RZ:
	*testfun_crlibm = asin_rz; 	break;
      default:
	*testfun_crlibm = asin_rn;
      }
#ifdef HAVE_MATHLIB_H
      *testfun_libultim  = uasin;  
#endif
#ifdef HAVE_LIBMCR_H
      *testfun_libmcr    = NULL;   /* TODO */
#endif
#ifdef HAVE_MPFR_H
      *testfun_mpfr   = mpfr_asin;
#endif
    }
  else if (strcmp (func_name, "acos") == 0)
    {
      *randfun_perf     = rand_for_acos_testperf;
      *randfun_soaktest = rand_for_asin_soaktest;
      *worst_case = 0.9999999999968665065352979581803 ;
      *testfun_libm   = acos;
      switch(crlibm_rnd_mode){
      case RU:
	*testfun_crlibm = acos_ru;	break;
      case RD:
	*testfun_crlibm = acos_rd;	break;
      case RZ:
	*testfun_crlibm = acos_rz; 	break;
      default:
	*testfun_crlibm = acos_rn;
      }
#ifdef HAVE_MATHLIB_H
      *testfun_libultim  = uacos;  
#endif
#ifdef HAVE_LIBMCR_H
      *testfun_libmcr    = NULL;   /* TODO */
#endif
#ifdef HAVE_MPFR_H
      *testfun_mpfr   = mpfr_acos;
#endif
    }

  else if (strcmp (func_name, "acospi") == 0)
    {
      *randfun_perf     = rand_for_asin_testperf;
      *randfun_soaktest = rand_for_asin_soaktest;
      *worst_case = 0.5; /* TODO */
      *testfun_libm   = tinkered_acospi;
      switch(crlibm_rnd_mode){
      case RU:
	*testfun_crlibm = acospi_ru;	break;
      case RD:
	*testfun_crlibm = acospi_rd;	break;
      case RZ:
	*testfun_crlibm = acospi_rz; 	break;
      default:
	*testfun_crlibm = acospi_rn;
      }
#ifdef HAVE_MATHLIB_H
      *testfun_libultim  = NULL;  
#endif
#ifdef HAVE_LIBMCR_H
      *testfun_libmcr    = NULL;   /* TODO */
#endif
#ifdef HAVE_MPFR_H
      *testfun_mpfr   = tinkered_mpfr_acospi;
#endif
    }

  else if (strcmp (func_name, "asinpi") == 0)
    {
      *randfun_perf     = rand_for_asin_soaktest;
      *randfun_soaktest = rand_for_asin_soaktest;
      *worst_case = 0.5; /* TODO */
      *testfun_libm   = tinkered_asinpi;
      switch(crlibm_rnd_mode){
      case RU:
	*testfun_crlibm = asinpi_ru;	break;
      case RD:
	*testfun_crlibm = asinpi_rd;	break;
      case RZ:
	*testfun_crlibm = asinpi_rz; 	break;
      default:
	*testfun_crlibm = asinpi_rn;
      }
#ifdef HAVE_MATHLIB_H
      *testfun_libultim  = NULL;  
#endif
#ifdef HAVE_LIBMCR_H
      *testfun_libmcr    = NULL;   /* TODO */
#endif
#ifdef HAVE_MPFR_H
      *testfun_mpfr   = tinkered_mpfr_asinpi;
#endif
    }


  else if (strcmp (func_name, "sinpi") == 0)
    {
      *randfun_perf     = rand_for_exp_perf;
      *randfun_soaktest = rand_for_trig_soaktest;
      *worst_case= 8.660177969107699326025337835166e-16;
      /* there is a worst one for directed rounding modes actually*/
      *testfun_libm   = tinkered_sinpi;
#ifdef HAVE_MPFR_H
      *testfun_mpfr   = tinkered_mpfr_sinpi;
#endif
      switch(crlibm_rnd_mode){
      case RU:
	*testfun_crlibm = sinpi_ru;	break;
      case RD:
	*testfun_crlibm = sinpi_rd;	break;
      case RZ:
	*testfun_crlibm = sinpi_rz;	break;
      default:
	*testfun_crlibm = sinpi_rn;
      }
    }



  else if (strcmp (func_name, "cospi") == 0)
    {
      *randfun_perf     = rand_for_exp_perf;
      *randfun_soaktest = rand_for_trig_soaktest;
      *worst_case= 4.496656279439104163951187380192e-08;
      *testfun_libm   = tinkered_cospi;
#ifdef HAVE_MPFR_H
      *testfun_mpfr   = tinkered_mpfr_cospi;
#endif
      switch(crlibm_rnd_mode){
      case RU:
	*testfun_crlibm = cospi_ru;	break;
      case RD:
	*testfun_crlibm = cospi_rd;	break;
      case RZ:
	*testfun_crlibm = cospi_rz;	break;
      default:
	*testfun_crlibm = cospi_rn;
      }
    }



  else if (strcmp (func_name, "tanpi") == 0)
    {
      *randfun_perf     = rand_for_exp_perf;
      *randfun_soaktest = rand_for_trig_soaktest;
      *worst_case= 2.826223498647412316943112325918e-06;
      *testfun_libm   = tinkered_tanpi;
#ifdef HAVE_MPFR_H
      *testfun_mpfr   = tinkered_mpfr_tanpi;
#endif
      switch(crlibm_rnd_mode){
      case RU:
	*testfun_crlibm = tanpi_ru;	break;
      case RD:
	*testfun_crlibm = tanpi_rd;	break;
      case RZ:
	*testfun_crlibm = tanpi_rz;	break;
      default:
	*testfun_crlibm = tanpi_rn;
      }
    }




  else if (strcmp (func_name, "pow") == 0)
    {
      *randfun_perf     = (double (*)()) rand_for_pow_perf;
      *randfun_soaktest = (double (*)()) rand_for_pow_soaktest;
      *worst_case = 2.0;
      *testfun_libm   = pow;
      switch(crlibm_rnd_mode){
      case RU:
	*testfun_crlibm = NULL;	break;
      case RD:
	*testfun_crlibm = NULL;	break;
      case RZ:
	*testfun_crlibm = NULL;	break;
      default:
	/*	*testfun_crlibm = pow_rn;*/
	*testfun_crlibm = pow_rn;

      }
#ifdef HAVE_MATHLIB_H
      *testfun_libultim = upow;
#endif
#ifdef HAVE_LIBMCR_H
      *testfun_libmcr    = __libmcr_pow;
#endif
#ifdef HAVE_MPFR_H
      *testfun_mpfr     = mpfr_pow; 
#endif
    }



  else
    {
      fprintf (stderr, "Unknown function: %s\n", func_name);
      exit (1);
    }
}

