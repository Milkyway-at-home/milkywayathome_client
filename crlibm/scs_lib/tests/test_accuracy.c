/** Test of SCS library against MPFR
@file tests/test_accuracy.c

@author David Defour David.Defour@ens-lyon.fr
@author Florent de Dinechin Florent.de.Dinechin@ens-lyon.fr 

This file is part of the SCS library.
*/

/*
Copyright (C) 2002  David Defour and Florent de Dinechin

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

 */
#include "scs.h"
#include "scs_private.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


/* Compile only if mpfr is present */
#ifdef HAVE_MPFR_H

#include <mpfr.h>


void (* test_scs_fct) () = NULL;
int (* test_mpfr_fct) () = NULL;
mp_rnd_t ROUNDING = GMP_RNDN;


/* Function defined in scs2MPF.c */
//void scs2MPF(scs_t, mpf_t);


void usage(char *namefct){
    printf("Usage :\n \t %s n [k] \n\n", namefct);
    printf("Where k is the number of random tests to generate (default 1000) \n");
    printf("and n is the function to test \n");
    printf(" 1 : scs_add \n");
    printf(" 2 : scs_sub \n");
    printf(" 3 : scs_mul \n");
    printf(" 4 : scs_div \n");
    printf(" 5 : scs_get_d \n");
    printf(" 6 : scs_get_d_minf \n");
    printf(" 7 : scs_get_d_pinf \n");
    printf(" 8 : scs_get_d_zero \n");
    printf(" 9 : scs_set \n");
    printf(" 10: scs_square \n");

    printf("\n");
    exit (1);
}


void test_square(int bcl){
  scs_t scx, scr, scx_maxerr, scr_maxerr;
  mpfr_t mpx, mpr,mpscr,mperr;
  int j;
  double err, maxerr, avgerr;

  mpfr_init2(mpx, SCS_NB_WORDS*SCS_NB_BITS*2);
  mpfr_init2(mpr, SCS_NB_WORDS*SCS_NB_BITS*2);
  mpfr_init2(mpscr, SCS_NB_WORDS*SCS_NB_BITS*2);
  mpfr_init2(mperr, SCS_NB_WORDS*SCS_NB_BITS*2);
  maxerr = 0; avgerr=0;

  for(j=0; j<bcl; j++){
    scs_rand(scx, 20);
    #if 1 
    /* set worst case */
    scx -> h_word[0] = 1; 
    #endif
    scs_get_mpfr(scx, mpx);  

    scs_square(scr, scx);
    mpfr_mul(mpr, mpx, mpx, GMP_RNDN);
    scs_get_mpfr(scr, mpscr);  

    mpfr_sub(mperr, mpscr, mpr, GMP_RNDN); 
    mpfr_div(mperr, mperr, mpr, GMP_RNDN); 
    mpfr_abs(mperr, mperr, GMP_RNDN);

    err = mpfr_get_d(mperr, GMP_RNDN);
    if (err > maxerr){
      maxerr = err;
      scs_set(scx_maxerr, scx);
      scs_set(scr_maxerr, scr);
    }
    avgerr +=  err;
  }
  avgerr = avgerr / bcl;
  printf("Average error : %e nb bit : %d \n", 
	 avgerr, -(int)(log(avgerr)/log(2.)));
  printf("Max     error : %e nb bit : %d \n", 
	 maxerr,  -(int)(log(maxerr)/log(2.)));

  printf("Argument giving max error : \n"); scs_get_std(scx_maxerr);
  printf("Result with max error: \n"); scs_get_std(scr_maxerr);
 
  mpfr_clear(mpx); mpfr_clear(mpr); mpfr_clear(mperr); 
  mpfr_clear(mpscr);
}




void test_one_arg(int bcl){
  scs_t sc1, scex;
  mpfr_t mpf1, mpfex;
  double scsd, mpfrd, scsdmax;
  int j, nberr;


  mpfr_init2(mpf1, SCS_NB_WORDS*SCS_NB_BITS*2);
  mpfr_init2(mpfex, SCS_NB_WORDS*SCS_NB_BITS*2);
  scsd = 0; mpfrd = 0; nberr = 0; scsdmax = 0;
  
  for(j=0; j<bcl; j++){
    scs_rand(sc1, 1500); /* test some special cases */
    test_scs_fct(&scsd, sc1);

    scs_get_mpfr(sc1, mpf1);  

    mpfrd = mpfr_get_d(mpf1, ROUNDING);

    if (mpfrd != scsd){
      scs_set(scex, sc1);
      nberr ++;}
  }
  printf("Number of misrounds: %d\n", nberr);
  if (nberr){
    printf("Original number :\n"); 
    scs_get_std(scex);
    scs_get_mpfr(scex, mpfex); 
    /* commented out because this function has disappeared from recent MPFR versions  */
/*     mpfr_out_str(stdout, 2, 150, mpfex, ROUNDING);  printf("\n"); */
 
    printf("SCS rounding : \n");
    test_scs_fct(&scsd, scex);
    printf("%.40e\n",scsd);
    printf("In binary : \n");
    mpfr_set_d(mpf1, scsd, ROUNDING); 
/*     mpfr_out_str(stdout, 2, 150, mpf1, ROUNDING);  printf("\n");  */
 
    printf("MPFR rounding : \n");
    mpfrd = mpfr_get_d(mpfex, ROUNDING);
    printf("%.40e\n",mpfrd);
    printf("In binary : \n");
    mpfr_set_d(mpf1, mpfrd, ROUNDING); 
/*     mpfr_out_str(stdout, 2, 150, mpf1, ROUNDING);  printf("\n");  */
  }

  mpfr_clear(mpf1); mpfr_clear(mpfex); 
}


void test_two_arg(int bcl){
  scs_t sc1, sc2, sc3, scm1, scm2, scm3;
  mpfr_t mp1, mp2, mp3, mp4, mp5;
  double d1, d2, max;
  int j;

  mpfr_init2(mp1, SCS_NB_WORDS*SCS_NB_BITS*2); 
  mpfr_init2(mp2, SCS_NB_WORDS*SCS_NB_BITS*2);
  mpfr_init2(mp3, SCS_NB_WORDS*SCS_NB_BITS*2);
  mpfr_init2(mp4, SCS_NB_WORDS*SCS_NB_BITS*2);
  mpfr_init2(mp5, SCS_NB_WORDS*SCS_NB_BITS*2); 
  d1 = 0; d2 = 0;  max = 0;
 
  
  for(j=0; j<bcl; j++){
    scs_rand(sc1, 20); scs_rand(sc2, 20); 
    /* You get most worst cases by imposing the following: */
    #if 1
    sc1 -> h_word[0] = 1; 
    sc2 -> h_word[0] = 1; 
    sc1 -> h_word[1] = 1; 
    sc2 -> h_word[1] = 1; 
    #endif
    test_scs_fct(sc3, sc1, sc2);

    scs_get_mpfr(sc1, mp1);  scs_get_mpfr(sc2, mp2);  scs_get_mpfr(sc3, mp3);
 
    test_mpfr_fct(mp4, mp1, mp2, GMP_RNDN);
 
    /* Error Computation */   
    mpfr_sub(mp5, mp4, mp3, GMP_RNDN); mpfr_div(mp5, mp5, mp4, GMP_RNDN); mpfr_abs(mp5, mp5, GMP_RNDN);
    
    d2 = mpfr_get_d(mp5, GMP_RNDN);
    if (d2 > max){
      max = d2;
      scs_set(scm1, sc1);
      scs_set(scm2, sc2);
      scs_set(scm3, sc3); }

    d1 +=  d2;
  }
  
  printf("Average error : %e nb bit : %d \n", d1/bcl, -(int)(log(d1/bcl)/log(2.)));
  printf("Max     error : %e nb bit : %d \n", max,    -(int)(log(max)/log(2.)));


  printf("First  Number : \n"); scs_get_std(scm1);
  printf("Second Number : \n"); scs_get_std(scm2);
  printf("Result Number : \n"); scs_get_std(scm3);

  mpfr_clear(mp1); mpfr_clear(mp2); mpfr_clear(mp3); 
  mpfr_clear(mp4); mpfr_clear(mp5);

}


void test_scs_set_d(int loops) {
  db_number d;
  double d1, d2;
  scs_t s;
  int i,j, errors, ep;
  double errorcases1[20];
  double errorcases2[20];

  errors=0; ep=0;
  
  for (i=0; i<loops; i++) {
    d.l = (rand() & 0x000000ff);
    for(j=0; j<(sizeof(db_number)); j++){
      d.l = d.l << 8;
      d.l += (rand() & 0x000000ff );
    }
    d1 = d.d;
    scs_set_d(s,d1);
    scs_get_d(&d2,s);
    if((d1!=d2) && (!(!(d1==d1)))) { 
      /* second test to filter out NaNs, for which the first is always
	 true */
      errors++;
      if(ep<20) {
	errorcases1[ep]=d1;
	errorcases2[ep]=d2;
	ep++;
      }
    }
  }
  printf("Number of errors: %d\n", errors);
  if(ep!=0) {
    for(i=0; i<ep; i++)
      printf("example:\t%e\t%e\n", errorcases1[i], errorcases2[i]);
  }
}


/*
 * Fct de test . 
 */
int main(int argc, char *argv[]){
  int fct, bcl;
 
 
  if (argc < 2)    usage(argv[0]);

  fct  = atoi(argv[1]);
  bcl  = (argc >= 3)? atoi(argv[2]) : 1000;


  switch(fct){
  case 1:
    test_scs_fct = scs_add;   test_mpfr_fct = mpfr_add;
    test_two_arg(bcl);
    break;
  case 2:
    test_scs_fct = scs_sub;   test_mpfr_fct = mpfr_sub;
    test_two_arg(bcl);
    break;
  case 3:
    test_scs_fct = scs_mul;   test_mpfr_fct = mpfr_mul;    
    test_two_arg(bcl);
    break;
  case 4:
    test_scs_fct = scs_div;   test_mpfr_fct = mpfr_div;    
    test_two_arg(bcl);
   break;
  case 5:
    test_scs_fct = scs_get_d;      ROUNDING = GMP_RNDN;
    test_one_arg(bcl);
    break;
  case 6:
    test_scs_fct = scs_get_d_minf; ROUNDING = GMP_RNDD;
    test_one_arg(bcl);
    break;
  case 7:
    test_scs_fct = scs_get_d_pinf; ROUNDING = GMP_RNDU;
    test_one_arg(bcl);
    break;
  case 8:
    test_scs_fct = scs_get_d_zero; ROUNDING = GMP_RNDZ;
    test_one_arg(bcl);
    break;
  case 9:
    test_scs_set_d(bcl);
    break;
  case 10:
    test_square(bcl);
    break;
  default :
    fprintf(stderr,"Unknown Function \n");
    usage(argv[0]);
  }

  return 0;
}
#else
int main(){
  fprintf(stderr,"No MPFR detected on your system.\n Please install MPFR, and compile scslib with MPFR support\n (see ./configure --help) \n");
  return 0;
}
#endif /*HAVE_LIBMPFR*/
