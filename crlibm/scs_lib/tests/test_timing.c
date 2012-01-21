#include "scs.h"
#include "scs_private.h"

/* Compile only if gmp is present */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include "tbx_timing.h"
#ifdef HAVE_GMP_H
#include <gmp.h>
#endif
#ifdef HAVE_MPFR_H
#include <mpfr.h>
#endif


#define LOOPS 1000



/*
 * Used to mesure time taken by the instruction "name"
 */
#define TST_FCT(char, name)\
         /* one untimed loop to load cache */ \
         for(i=0; i< LOOPS-1; i++){ \
           name;\
           } \
         TBX_GET_TICK(t1);\
         for(i=0; i< LOOPS-1; i++){ \
           name;\
           } \
         TBX_GET_TICK(t2);\
	 deltat = TBX_TICK_RAW_DIFF(t1, t2); \
         printf("%28s :  %lld ticks,\t ratio to  FP +:  %f\n", char, deltat, (double) deltat/tadd);


/* Similar to the previous, computes the time for one FP add so that
   we have something to normalize against */
#define COMPUTE_TADD()\
         /* one untimed loop to load cache */ \
         for(i=0; i< LOOPS-1; i++){ \
           d3 = double_table[i]*double_table[i+1];\
           } \
         TBX_GET_TICK(t1);\
         for(i=0; i< LOOPS-1; i++){ \
           d3 = double_table[i]*double_table[i+1];\
           } \
         TBX_GET_TICK(t2);\
	 tadd = TBX_TICK_RAW_DIFF(t1, t2);


/*
 * Fct de test . 
 */
int main(){
  scs_t n1I;
  volatile double d3;
  volatile int int_r;
  unsigned long long deltat,tadd;
  tbx_tick_t   t1, t2; 
  int i;
 
  /* table storing random number */
  scs_t  scs_table[LOOPS];
#ifdef HAVE_GMP_H
  mpf_t  mpf_table[LOOPS];
  mpf_t a;
#endif
#ifdef HAVE_MPFR_H
  mpfr_t  mpfr_table[LOOPS];
  mpfr_t  mpfr_a;
#endif

  double double_table[LOOPS];
  int int_table[LOOPS];



  printf("Random generation  ... ");
  srand(42);  
  for(i=0; i<LOOPS; i++){ 
    scs_rand(scs_table[i], 7);  
#if 1
    scs_table[i]->sign = 1; /* only positive numbers */
#endif
    scs_get_d(&double_table[i], scs_table[i]);
    int_table[i] = double_table[i];
#ifdef HAVE_GMP_H
    mpf_init2(mpf_table[i], (SCS_NB_BITS*SCS_NB_WORDS) );
    scs_get_mpf(scs_table[i], mpf_table[i]);
#endif
#ifdef HAVE_MPFR_H
    mpfr_init2(mpfr_table[i], (SCS_NB_BITS*SCS_NB_WORDS) );
    scs_get_mpfr(scs_table[i], mpfr_table[i]);
#endif
  }
  printf(" done \n\n");


#ifdef HAVE_GMP_H
  mpf_init2(a, (SCS_NB_BITS*SCS_NB_WORDS));
  mpf_set(a, mpf_table[1]);
#endif
#ifdef HAVE_MPFR_H
  mpfr_init2(mpfr_a, (SCS_NB_BITS*SCS_NB_WORDS));
  mpfr_set(mpfr_a, mpfr_table[1], GMP_RNDN);
#endif

  printf("These first timings don't mean much\n");
  COMPUTE_TADD()
  TST_FCT("int a + b ", int_r = int_table[i]+int_table[i+1])
  TST_FCT("int a * b ", d3 = int_table[i]*int_table[i+1])
  TST_FCT("double a + b ", d3 = double_table[i]+double_table[i+1])
  TST_FCT("double a * b ", d3 = double_table[i]*double_table[i+1])
  TST_FCT("double a / b ", d3 = double_table[i]/double_table[i+1])
  printf("\n"); 

  printf("Here come the meaningful timings\n");
  /* scs library test */
  TST_FCT("conversion scs=>doubles ", scs_get_d(&double_table[i], scs_table[i]))
  TST_FCT("conversion doubles=>scs ", scs_set_d(n1I, double_table[i]))
  TST_FCT("scs_add ", scs_add(n1I, scs_table[i], scs_table[i+1]))
  TST_FCT("scs_sub ", scs_sub(n1I, scs_table[i], scs_table[i+1]))
  TST_FCT("scs_add_no_renorm ",scs_add_no_renorm(n1I, scs_table[i], scs_table[i+1]))
  TST_FCT("scs_mul ", scs_mul(n1I, scs_table[i], scs_table[i+1]))
  TST_FCT("scs_mul_ui ", scs_mul_ui(scs_table[i], 31242436))
  TST_FCT("scs_square  ", scs_square(n1I, scs_table[i]))
    //  TST_FCT("scs_fma ", scs_fma(n1I, scs_table[i], scs_table[i], scs_table[i+1]))
  TST_FCT("add + mul scs ", scs_mul(n1I, scs_table[i], scs_table[i+1]); scs_add(n1I, n1I, scs_table[i])) 

  TST_FCT("renormalization scs ", scs_renorm(scs_table[i]))
  TST_FCT("scs_div ", scs_div(n1I, scs_table[i], scs_table[i+1]))
  printf("\n");  

#ifdef HAVE_GMP_H
  /* mpf (gmp) library test */
  TST_FCT("Conversion mpf=>double", double_table[i] = mpf_get_d(mpf_table[i]))
  TST_FCT("Conversion double=>mpf ", mpf_set_d(a, double_table[i]))
  TST_FCT("Addition mpf ", mpf_add(a, mpf_table[i], mpf_table[i+1]))
  TST_FCT("Multiplication mpf ", mpf_mul(a, mpf_table[i], mpf_table[i+1]))
  TST_FCT("Multiplication_with_int mpf ", mpf_mul_ui(a, mpf_table[i], 3254353))
  TST_FCT("Division mpf ", mpf_div(a, mpf_table[i], mpf_table[i+1]))
  printf("\n");  
#endif
 

#ifdef HAVE_MPFR_H
  /* mpf (gmp) library test */
  TST_FCT("Conversion mpfr=>double", double_table[i] = mpfr_get_d(mpfr_table[i],GMP_RNDN))
  TST_FCT("Conversion double=>mpfr ", mpfr_set_d(mpfr_a, double_table[i],GMP_RNDN))
  TST_FCT("Addition mpfr ", mpfr_add(mpfr_a, mpfr_table[i], mpfr_table[i+1],GMP_RNDN))
  TST_FCT("Multiplication mpfr ", mpfr_mul(mpfr_a, mpfr_table[i], mpfr_table[i+1],GMP_RNDN))
  TST_FCT("Multiplication_with_int mpfr ", mpfr_mul_ui(mpfr_a, mpfr_table[i], 3254353,GMP_RNDN))
  TST_FCT("Division mpfr ", mpfr_div(mpfr_a, mpfr_table[i], mpfr_table[i+1],GMP_RNDN))
  printf("\n");  
#endif
 



  return 0;
}

