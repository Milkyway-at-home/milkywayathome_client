/*
 *this function computes log, correctly rounded,
 using experimental techniques based on  double-extended arithmetic

 THIS IS EXPERIMENTAL SOFTWARE

In particular it changes rounding modes all the time without warning
nor restoring.

 *
 * Author :  Florent de Dinechin
 * Florent.de.Dinechin at ens-lyon.fr
 *

 To have it replace the crlibm log, do:

on pentium,
 gcc -DHAVE_CONFIG_H -I.  -fPIC  -O2 -c log-de2.c;   mv log-de2.o log_fast.o; make

on itanium,
icc  -I/users/fdedinex/local/IA64/include -mcpu=itanium2\
 -Qoption,cpp,--extended_float_types -IPF_fp_speculationsafe -c log-de2.c;\
 mv log-de2.o log_fast.o; make


*/


#include <stdio.h>
#include <stdlib.h>
#include "crlibm.h"
#include "crlibm_private.h"
#include "double-extended.h"
#include "log-de2.h"


double log_rn(double x) {
  double wi;
#if defined(CRLIBM_TYPECPU_X86) || defined(CRLIBM_TYPECPU_AMD64)
  db_number y;
  int E, i, index, roundtestmask;
  long double r, logirh, logirl, z, z2, z4, t,evenp,oddp, th, tl, eh,el,p1,p2,p3;
#else
  long int E, i;
  unsigned long int index, roundtestmask;
  __fpreg xe, ye, r, logirh, logirl, z, z2, z4, t,evenp,oddp, th, tl, eh,el, c1,c2,c3,c4,c5,c6,p1,p2,p3;
#endif


   E=0;


#if defined(CRLIBM_TYPECPU_X86) || defined(CRLIBM_TYPECPU_AMD64)
   /* Pentium has few registers, so load coefficients directly from memory */
#define c6 c[0]
#define c5 c[1]
#define c4 c[2]
#define c3 c[3]
#define c2 c[4]
#define c1 c[5]
  /* Special case hadling on a x86  */
   y.d=x;

   /* Filter special cases */
   if (y.i[HI] < 0x00100000){        /* x < 2^(-1022)    */

     if (((y.i[HI] & 0x7fffffff)|y.i[LO])==0){
       return -1.0/0.0;
     }                    		   /* log(+/-0) = -Inf */
     if (y.i[HI] < 0){
       return (x-x)/0;                      /* log(-x) = Nan    */
     }
     /* Subnormal number */
     E = -64; 		
     y.d *= two64; 	  /* make x a normal number    */
   }

   if (y.i[HI] >= 0x7ff00000){
     return  x+x;				 /* Inf or Nan       */
   }

   DOUBLE_EXTENDED_MODE;  /* This one should be overlapped with integer computation */

#define X_NEAR_1 (y.i[HI]>MINYFAST) && (y.i[HI]<MAXYFAST)


#else
   /* prefetch coefficients */
   c6=c[0]; c5=c[1]; c4=c[2]; c3=c[3]; c2=c[4]; c1=c[5];
   /* Special case handling on the Itanium */
   xe=x;
   i =  _Asm_getf(2/*_FR_D*/, xe);

   /* Filter special cases */
   //if (__builtin_expect( (i<0x0010000000000000ULL), (1+1==3))){        /* x < 2^(-1022)    */
   if (i<0x0010000000000000LL){        /* x < 2^(-1022)    */
     if ((i & 0x7fffffffffffffffULL)==0){
       return -1.0/0.0;
     }                    		   /* log(+/-0) = -Inf */
     if (i<0){
       return (x-x)/0;                      /* log(-x) = Nan    */
     }
     /* Subnormal number */
     xe *= two64; 	  /* make x a normal number    */
     E = -64;
     i =  _Asm_getf(2/*_FR_D*/, xe); /* and update i */
   }

   //if (__builtin_expect(  (i >= 0x7ff0000000000000ULL), (1+1==3) )){
   if (i >= 0x7ff0000000000000ULL){
     return  x+x;				 /* Inf or Nan       */
   }
#define X_NEAR_1 __builtin_expect((i>(((uint64_t) MINYFAST)<<32)) && (i<(((uint64_t) MAXYFAST)<<32)), (1+1==3))


#endif


   if(X_NEAR_1) {
     roundtestmask=0x7fc;
     z = x - 1 ; /* Sterbenz exact */
     z2 = z*z;
     evenp = c6;                  /* c6 */
     oddp  = c5;                  /* c5 */
     evenp = c4  +  z2 * evenp;   /* c4 */
     oddp  = c3  +  z2 * oddp;    /* c3 */
     evenp = c2  +  z2 * evenp;   /* c2 */

     t  = c1 + (z * evenp + z2 * oddp);
     t = z*t;

     //printf("z= %1.20e,  t=%1.20e  \n  ", (double)z, (double)t);

   }

   else {

#if defined(CRLIBM_TYPECPU_X86) || defined(CRLIBM_TYPECPU_AMD64)
   /* Extract exponent and mantissa */
     E += (y.i[HI]>>20)-1023;             /* extract the exponent */
     index = (y.i[HI] & 0x000fffff);
     y.i[HI] =  index | 0x3ff00000;	/* do exponent = 0 */
     index = index  >> (20-L);
     /* now y.d holds 1+f, and E is the exponent */

     logirh = argredtable[index].h;
     r = (long double) (argredtable[index].r); /* approx to 1/y.d */
     z = r*(long double)y.d - 1. ; /* even without an FMA, all exact */

     if(E>12 || E<-12) { /* faster and less accurate polynomial evaluation */
       roundtestmask=0x7fe;
       p1 = logirh + z*c1;   p2 = c2 + z*c3;   p3 = c4 + z*c5;   z2 = z*z;
       p1 = p1 + z2*p2;      p3 = p3 + z2*c6;  z4=z2*z2;
       t =  p1 + z4*p3;
       t = t + E*ln2h;
     }
     else {
       roundtestmask=0x7f0;
       p1 = c5 + z*c6;         z2 = z*z;
       p2 = c3 + z*c4;         p3 = c1+z*c2;

       t = p2 + z2*p1;
       t = p3 + z2*t;
       t = logirh + z*t;
       t = t + E*ln2h;
     }

#if 0
     if(E>12 || E<-12)
       roundtestmask=0x7fe;
     else
       roundtestmask=0x7f0;

     p1 = c5 + z*c6;         z2 = z*z;
     p2 = c3 + z*c4;         p3 = c1+z*c2;

     t = p2 + z2*p1;
     t = p3 + z2*t;
     t = logirh + z*t;
#endif

#else /* Itanium here*/
   /* Extract exponent and mantissa */
     E += (i>>52)-1023;
     //printf("\nE    = %llx\n", E);
     i = i & 0x000fffffffffffffULL;  /* keep only mantissa */
     ye = _Asm_setf(2/*_FR_D*/, i | 0x3ff0000000000000ULL ); /* exponent = 0*/
     index = i >> (52-L);
     //printf("\nindex= %lld\n", index);
     //printf("\n  ye = %1.20Le\n",(long double)ye);
     /* now ye holds 1+f, and E is the exponent */

     logirh = argredtable[index].h;

     _Asm_frcpa(&r, 1.0L, ye, 1/*_SF1*/);
     z = r*ye - 1. ; /* even without an FMA, all exact */

     if(E>12 || E<-12) { /* faster and less accurate polynomial evaluation */
       roundtestmask=0x7fe;
       p1 = logirh + z*c1;   p2 = c2 + z*c3;   p3 = c4 + z*c5;   z2 = z*z;
       p1 = p1 + z2*p2;      p3 = p3 + z2*c6;  z4=z2*z2;
       t =  p1 + z4*p3;
       t = t + E*ln2h;
     }
     else {
       roundtestmask=0x7f0;
       p1 = c5 + z*c6;         z2 = z*z;
       p2 = c3 + z*c4;         p3 = c1+z*c2;

       t = p2 + z2*p1;
       t = p3 + z2*t;
       t = logirh + z*t;
       t = t + E*ln2h;
     }


#endif

     //printf("  x=%1.20Le\n  r=%1.20Le\n z=%1.20Le\n  logirh=%1.20Le\n  ",(long double)xe, (long double)r,(long double)z, (long double)logirh);
     /* Polynomial evaluation, unrolled to go through Gappa */

     //printf("t=%1.20e  \n  ", (double)t);


   }





#if 0 /* to time the first step only */
   BACK_TO_DOUBLE_MODE; return (double)t;
#endif


   /* To test the second step only, comment out the following line */
   DE_TEST_AND_RETURN_RN(t, roundtestmask);


   /* Accurate phase */
#if EVAL_PERF
  crlibm_second_step_taken++;
#endif

   t = c13h;
   t = c12h + z*t;
   t = c11h + z*t;
   t = c10h + z*t;
   t = c9h  + z*t;
   t = c8h  + z*t;

   //printf("\n t = %1.20Le\n", (long double)t);

   Mul12_ext(&th, &tl,   z,   t);
   Add22_ext(&th, &tl,   c7h,c7l,  th,tl);
   FMA22_ext(&th, &tl,    z,0, th,tl,  c6h,c6l);
   FMA22_ext(&th, &tl,    z,0, th,tl,  c5h,c5l);
   FMA22_ext(&th, &tl,    z,0, th,tl,  c4h,c4l);
   FMA22_ext(&th, &tl,    z,0, th,tl,  c3h,c3l);
   FMA22_ext(&th, &tl,    z,0, th,tl,  c2h,c2l);
   FMA22_ext(&th, &tl,    z,0, th,tl,  c1h,c1l);

   if((X_NEAR_1)) {
     Mul22_ext(&th, &tl,  z,0,  th,tl);
   }
   else{
   FMA22_ext(&th, &tl,    z,0, th,tl,  logirh, argredtable[index].l);

     /* the following is not a FMA22 (Eln2 > th+tl) */
     Mul22_ext(&eh, &el,   ln2h,ln2l, E, 0);
     Add22_ext(&th, &tl,   eh,el,  th,tl);
   }
   BACK_TO_DOUBLE_MODE;
   return (double) (th+tl); /* The exact sum of these double-extended is rounded to the nearest */
}


double log_ru(double x) { return x;};
double log_rd(double x) { return x;};
double log_rz(double x) { return x;};
