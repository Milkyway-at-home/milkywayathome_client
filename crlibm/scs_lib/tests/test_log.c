/*
 * Function to compute the logarithm with fully exact rounding
 *
 * Author : Defour David  (David.Defour@ens-lyon.fr)
 *
 * Date of creation : 13/05/2002   
 * Last Modified    : 17/05/2002
 */
#include "scs.h"
#include "scs_private.h"
#include "log.h"
#include <stdio.h>
#ifdef HAVE_GMP_H
#include <gmp.h>
#include "tbx_timing.h"




#define LOOPS 10000

mpf_t mpf_sc_ln2_ptr;
mpf_t mpf_constant_poly_ptr[20];
mpf_t mpf_table_ti_ptr[13];
mpf_t mpf_table_inv_wi_ptr[13];


/*
 * 6) SCS CONSTANTS
 */
#ifdef SCS_TYPECPU_SPARC
static const scs
/* 0   */
   scs_zer ={{0x00000000, 0x00000000, 0x00000000, 0x00000000},
             {{0, 0}},  0,   1 },
/* 1/2 */
   scs_half={{0x02000000, 0x00000000, 0x00000000, 0x00000000},
             DB_ONE, -1,   1 },
/*  1  */  
   scs_one ={{0x00000001, 0x00000000, 0x00000000, 0x00000000},
             DB_ONE,  0,   1 },
/*  2  */
   scs_two ={{0x00000002, 0x00000000, 0x00000000, 0x00000000},
             DB_ONE,  0,   1 };
#else
static const scs
/* 0   */
   scs_zer ={{0x00000000, 0x00000000, 0x00000000, 0x00000000,
             0x00000000, 0x00000000, 0x00000000, 0x00000000},
             {{0, 0}},  0,   1 },
/* 1/2 */
   scs_half={{0x20000000, 0x00000000, 0x00000000, 0x00000000,
             0x00000000, 0x00000000, 0x00000000, 0x00000000},
             DB_ONE, -1,   1 },
/*  1  */  
   scs_one ={{0x00000001, 0x00000000, 0x00000000, 0x00000000,
             0x00000000, 0x00000000, 0x00000000, 0x00000000},
             DB_ONE,  0,   1 },
/*  2  */
   scs_two ={{0x00000002, 0x00000000, 0x00000000, 0x00000000,
             0x00000000, 0x00000000, 0x00000000, 0x00000000},
             DB_ONE,  0,   1 };
#endif

#define SCS_ZERO (scs_ptr)(&scs_zer)
#define SCS_HALF (scs_ptr)(&scs_half)
#define SCS_ONE  (scs_ptr)(&scs_one)
#define SCS_TWO  (scs_ptr)(&scs_two)

double sn_log(double); 
double mpf_sn_log(double);


/*
 *  1) First reduction: exponent extraction      
 *         E  
 *  x = 2^   .(1+f)    with  0 <= f < 1
 *
 *  log(x) = E.log(2) + log(1+f) where:
 *     - log(2)   is tabulated
 *     - log(1+f) need to be evalute 
 *  
 *
 *  2) Avoiding accuracy problem when E=-1 by testing
 *   
 *    if (1+f >= sqrt(2)) then 
 *        1+f = (1+f)/2;  E = E+1; 
 *    and,
 *        log(x) = (E+1).log(2) + log((1+f)/2)
 *
 *    so now:      sqrt(2)/2 <= (1+f) < sqrt(2)
 *
 *
 *  3) Second reduction: tabular reduction
 *                   -4  
 *   wi = 1 + i. 2^
 *                                   1  
 *   log(1+f) = log(wi) + log ( 1 + --- . (1 + f - wi) ) 
 *                                   wi 
 *
 *   then |(1+f-wi)/wi| <= 2^-5 if we use rounded to nearest.
 *
 *  4) Computation:
 *   a) Table lookup of: 
 *        - ti    = log(wi)
 *        - inv_wi = 1/(wi)
 *   b) Polynomial evaluation of:
 *        - P(R) ~ log(1 + R), where R = (1+f-wi) * inv_wi 
 *
 *                 -5 
 *   with  |R| < 2^
 *
 *
 *  5) Reconstruction:
 *   log(x) = E.log(2) + t_i + P(R)
 *
 */
#define SQRT_2 1.4142135623730950489e0 



/*************************************************************
 *************************************************************
 *               ROUNDED  TO NEAREST
 *************************************************************
 *************************************************************/
double sn_log(double x){ 
  scs_t R, sc_ln2_times_E, res1;
  db_number nb, nb2, wi, resd;
  int ti, i, E=0;

  nb.d = x;
  /* Filter cases */
  if (nb.i[HI_ENDIAN] < 0x00100000){        /* x < 2^(-1022)    */
    if (((nb.i[HI_ENDIAN] & 0x7fffffff)|nb.i[LO_ENDIAN])==0)
      return 1.0/0.0;                       /* log(+/-0) = -Inf */
    if (nb.i[HI_ENDIAN] < 0) 
      return (x-x)/0;                       /* log(-x) = Nan    */

    /* Subnormal number */
    E    -= (SCS_NB_BITS*2); /* keep in mind that x is a subnormal number */ 
    nb.d *=SCS_RADIX_TWO_DOUBLE;  /* make x as normal number     */         
    /* We may just want add 2 to the scs number.index */
    /* may be .... we will see */
  }
  if (nb.i[HI_ENDIAN] >= 0x7ff00000)
    return x+x;                             /* Inf or Nan       */

  /* find n, nb.d such that sqrt(2)/2 < nb.d < sqrt(2) */
  E += (nb.i[HI_ENDIAN]>>20)-1023;
  nb.i[HI_ENDIAN] =  (nb.i[HI_ENDIAN] & 0x000fffff) | 0x3ff00000;
  if (nb.d > SQRT_2){
    nb.d *= 0.5;
    E++;
  }

  /* to normalize nb.d and round to nearest      */
  /* + (1-trunc(sqrt(2.)/2 * 2^(4))*2^(-4) )+2.^(-(4+1))*/ 
  nb2.d = nb.d + norm_number.d; 
  i = (nb2.i[HI_ENDIAN] & 0x000fffff);
  i = i >> 16; /* 0<= i <=11 */
  

  wi.d = (11+i)*(double)0.6250e-1;

  /* (1+f-w_i) */
  nb.d -= wi.d; 

  /* Table reduction */
  ti = i; 

  /* R = (1+f-w_i)/w_i */
  scs_set_d(R, nb.d);
  scs_mul(R, R, table_inv_wi_ptr[i]);

  /* sc_ln2_times_E = E*log(2)  */
  scs_set(sc_ln2_times_E, sc_ln2_ptr);

  if (E >= 0){
    scs_mul_ui(sc_ln2_times_E, E);
  }else{
    scs_mul_ui(sc_ln2_times_E, -E);
    sc_ln2_times_E->sign = -1;
  }

  /*
   * Polynomial evaluation of log(1 + R) with an error less than 2^(-130)
   */
  scs_mul(res1, constant_poly_ptr[0], R);
  for(i=1; i<20; i++){
    scs_add(res1, constant_poly_ptr[i], res1);
    scs_mul(res1, res1, R);
  }
  scs_add(res1, res1, table_ti_ptr[ti]);
  scs_add(res1, res1, sc_ln2_times_E);  
  

  scs_get_d(&resd.d, res1);  

  return resd.d;
}



/*************************************************************
 *************************************************************
 *               ROUNDED  TO NEAREST
 *************************************************************
 *************************************************************/
double mpf_sn_log(double x){ 
    mpf_t R, sc_ln2_times_E, res1;
    db_number nb, nb2, wi, resd;
    int ti, i, E=0;

    /* Memory allocation */
    mpf_init2(R, 210);
    mpf_init2(sc_ln2_times_E, 210);
    mpf_init2(res1, 210);

  nb.d = x;
  /* Filter cases */
  if (nb.i[HI_ENDIAN] < 0x00100000){        /* x < 2^(-1022)    */
    if (((nb.i[HI_ENDIAN] & 0x7fffffff)|nb.i[LO_ENDIAN])==0)
      return 1.0/0.0;                       /* log(+/-0) = -Inf */
    if (nb.i[HI_ENDIAN] < 0) 
      return (x-x)/0;                       /* log(-x) = Nan    */

    /* Subnormal number */
    E    -= (SCS_NB_BITS*2); /* keep in mind that x is a subnormal number */ 
    nb.d *=SCS_RADIX_TWO_DOUBLE;  /* make x as normal number     */         
    /* We may just want add 2 to the scs number.index */
    /* may be .... we will see */
  }
  if (nb.i[HI_ENDIAN] >= 0x7ff00000)
    return x+x;                             /* Inf or Nan       */


  /* find n, nb.d such that sqrt(2)/2 < nb.d < sqrt(2) */
  E += (nb.i[HI_ENDIAN]>>20)-1023;
  nb.i[HI_ENDIAN] =  (nb.i[HI_ENDIAN] & 0x000fffff) | 0x3ff00000;
  if (nb.d > SQRT_2){
    nb.d *= 0.5;
    E++;
  }

  /* to normalize nb.d and round to nearest      */
  /* + (1-trunc(sqrt(2.)/2 * 2^(4))*2^(-4) )+2.^(-(4+1))*/ 
  nb2.d = nb.d + norm_number.d; 
  i = (nb2.i[HI_ENDIAN] & 0x000fffff);
  i = i >> 16; /* 0<= i <=11 */
  

  wi.d = (11+i)*(double)0.6250e-1;

  /* (1+f-w_i) */
  nb.d -= wi.d; 


  /* Table reduction */
  ti = i; 
 
  /* R = (1+f-w_i)/w_i */
  mpf_set_d(R, nb.d);
  mpf_mul(R, R, mpf_table_inv_wi_ptr[i]);
 

  /* sc_ln2_times_E = E*log(2)  */
  mpf_set(sc_ln2_times_E, mpf_sc_ln2_ptr);


  if (E >= 0){
    mpf_mul_ui(sc_ln2_times_E, sc_ln2_times_E, E);
  }else{
    mpf_mul_ui(sc_ln2_times_E, sc_ln2_times_E, -E);
    mpf_neg(sc_ln2_times_E, sc_ln2_times_E);
  }


  /*
   * Polynomial evaluation of log(1 + R) with an error less than 2^(-130)
   */
  mpf_mul(res1, mpf_constant_poly_ptr[0], R);
  for(i=1; i<20; i++){
    mpf_add(res1, mpf_constant_poly_ptr[i], res1);
    mpf_mul(res1, res1, R);
  }
  mpf_add(res1, res1, mpf_table_ti_ptr[ti]);
  mpf_add(res1, res1, sc_ln2_times_E);  

  resd.d = mpf_get_d(res1);  

  /* Free Memory */
  mpf_clear(R);  mpf_clear(sc_ln2_times_E);  mpf_clear(res1);

  return resd.d;
}



void free_mpf_cst(){
    int i;

    for (i=0; i<13; i++) mpf_clear(mpf_table_ti_ptr[i]);
    for (i=0; i<13; i++) mpf_clear(mpf_table_inv_wi_ptr[i]);
    for (i=0; i<20; i++) mpf_clear(mpf_constant_poly_ptr[i]);
    mpf_clear(mpf_sc_ln2_ptr);
}

void init_mpf_cst(){

    /*
     * mpf constant initialisation 
     */
 mpf_init_set_str(mpf_sc_ln2_ptr,".6931471805599453094172321214581765680755001343602552541206800094933936",10);
 mpf_init_set_str(mpf_table_ti_ptr[0],"-.3746934494414106936069849078675769724802936835036038412641523288430009",10);
 mpf_init_set_str(mpf_table_ti_ptr[1],"-.2876820724517809274392190059938274315035097108977610565066656853492930",10);
 mpf_init_set_str(mpf_table_ti_ptr[2],"-.2076393647782445016154410442673876674967325926808139000636745273101098",10);
 mpf_init_set_str(mpf_table_ti_ptr[3],"-.1335313926245226231463436209313499745894156734989045739026498785426010",10);
 mpf_init_set_str(mpf_table_ti_ptr[4],"-.06453852113757117167292391568399292812890862534975384283537781286190121",10);
 mpf_init_set_str(mpf_table_ti_ptr[5],"0.0",10);
 mpf_init_set_str(mpf_table_ti_ptr[6],".06062462181643484258060613204042026328620247514472377081451769990871809",10);
 mpf_init_set_str(mpf_table_ti_ptr[7],".1177830356563834545387941094705217050684807125647331411073486387948077",10);
 mpf_init_set_str(mpf_table_ti_ptr[8],".1718502569266592223400989460551472649353787238581078020552401984357182",10);
 mpf_init_set_str(mpf_table_ti_ptr[9],".2231435513142097557662950903098345033746010855480072136712878724873917",10);
 mpf_init_set_str(mpf_table_ti_ptr[10],".2719337154836417588316694945329991619825747499635896237113644456014997",10);
 mpf_init_set_str(mpf_table_ti_ptr[11],".3184537311185346158102472135905995955952064508566514128565276806503928",10);
 mpf_init_set_str(mpf_table_ti_ptr[12],".3629054936893684531378243459774898461403797773994147255159153395094188",10);
  
 mpf_init_set_str(mpf_table_inv_wi_ptr[0],"1.454545454545454545454545454545454545454545454545454545454545454545455",10);
 mpf_init_set_str(mpf_table_inv_wi_ptr[1],"1.333333333333333333333333333333333333333333333333333333333333333333333",10);
 mpf_init_set_str(mpf_table_inv_wi_ptr[2],"1.230769230769230769230769230769230769230769230769230769230769230769231",10);
 mpf_init_set_str(mpf_table_inv_wi_ptr[3],"1.142857142857142857142857142857142857142857142857142857142857142857143",10);
 mpf_init_set_str(mpf_table_inv_wi_ptr[4],"1.066666666666666666666666666666666666666666666666666666666666666666667",10);
 mpf_init_set_str(mpf_table_inv_wi_ptr[5],"1.0",10);
 mpf_init_set_str(mpf_table_inv_wi_ptr[6],".9411764705882352941176470588235294117647058823529411764705882352941176",10);
 mpf_init_set_str(mpf_table_inv_wi_ptr[7],".8888888888888888888888888888888888888888888888888888888888888888888889",10);
 mpf_init_set_str(mpf_table_inv_wi_ptr[8],".8421052631578947368421052631578947368421052631578947368421052631578947",10);
 mpf_init_set_str(mpf_table_inv_wi_ptr[9],".8000000000000000000000000000000000000000000000000000000000000000000000",10);
 mpf_init_set_str(mpf_table_inv_wi_ptr[10],".7619047619047619047619047619047619047619047619047619047619047619047619",10);
 mpf_init_set_str(mpf_table_inv_wi_ptr[11],".7272727272727272727272727272727272727272727272727272727272727272727273",10);
 mpf_init_set_str(mpf_table_inv_wi_ptr[12],".6956521739130434782608695652173913043478260869565217391304347826086957",10);
  


 mpf_init_set_str(mpf_constant_poly_ptr[0],"-.5023367294567568078075892234129085015737046673117638148835793879900684e-1",10);
 mpf_init_set_str(mpf_constant_poly_ptr[1],".5286469364636800162451931056605008699849205395651010015123349866702438e-1",10);
 mpf_init_set_str(mpf_constant_poly_ptr[2],"-.5555504165240301671600703643945025506173351812330050928639639013749408e-1",10);
 mpf_init_set_str(mpf_constant_poly_ptr[3],".5882304523791230393387514237891031448778320732847321758009922427128675e-1",10);
 mpf_init_set_str(mpf_constant_poly_ptr[4],"-.6250000063225403563289268352338121550211573212295444469323271379834995e-1",10);
 mpf_init_set_str(mpf_constant_poly_ptr[5],".6666666722317130353982967043683210254854513387493745794307824211890102e-1",10);
 mpf_init_set_str(mpf_constant_poly_ptr[6],"-.7142857142809460344869308613352646790467720996336436837300223496769771e-1",10);
 mpf_init_set_str(mpf_constant_poly_ptr[7],".7692307692269045639218260467556323317594125412187843530035925183575360e-1",10);
 mpf_init_set_str(mpf_constant_poly_ptr[8],"-.8333333333333356037828752693534976595948243524382699349340678490229276e-1",10);
 mpf_init_set_str(mpf_constant_poly_ptr[9],".9090909090909107517931824416149710386096279143853141933153563186364239e-1",10);
 mpf_init_set_str(mpf_constant_poly_ptr[10],"-.9999999999999999993224242653285578758483071952239027712038970659952908e-1",10);
 mpf_init_set_str(mpf_constant_poly_ptr[11],".1111111111111111110676605039270706067112823008194317757408530067415460",10);
 mpf_init_set_str(mpf_constant_poly_ptr[12],"-.1250000000000000000000121547531034901577892307419239533855249184555657",10);
 mpf_init_set_str(mpf_constant_poly_ptr[13],".1428571428571428571428636714934431388385979630564777290210429555695253",10);
 mpf_init_set_str(mpf_constant_poly_ptr[14],"-.1666666666666666666666666654681759449880186388123105519286024077940397",10);
 mpf_init_set_str(mpf_constant_poly_ptr[15],".1999999999999999999999999995018689489549856935258632875129946570831422",10);
 mpf_init_set_str(mpf_constant_poly_ptr[16],"-.2500000000000000000000000000000541884832055314060603150374550752650562",10);
 mpf_init_set_str(mpf_constant_poly_ptr[17],".3333333333333333333333333333333480752710184240674027421784076979756719",10);
 mpf_init_set_str(mpf_constant_poly_ptr[18],"-.4999999999999999999999999999999999992783495160517279217122998393096071",10);
 mpf_init_set_str(mpf_constant_poly_ptr[19],".9999999999999999999999999999999999999280145118565650165317802247440368",10);

}


/*
 * Used to mesure time taken by the instruction "name"
 */
#define TST_FCT(char, name)\
         TBX_GET_TICK(t1);\
         for(i=0; i< LOOPS-1; i++){ \
           name;\
           } \
         TBX_GET_TICK(t2);\
	 deltat = TBX_TICK_RAW_DIFF(t1, t2); \
	 printf("%28s :  %lld ticks \n", char, deltat);



int main(){
    /* table storing random number */
    double  tableau[LOOPS];
    double res[LOOPS];
    mpf_t  TAB_MPF[LOOPS];
    int i;
    tbx_tick_t   t1, t2; 
    unsigned long long deltat;

    printf("mpf constant initialisation ... \n");
    init_mpf_cst();
    printf("Finished ... \n");

    srand(42);  
    for(i=0; i<LOOPS; i++) tableau[i]=rand();
    printf("End of random number generation \n");

    for(i=0; i< LOOPS-1; i++){ 
      res[i]=sn_log(tableau[i]);
    } 
    for(i=0; i< LOOPS-1; i++){ 
      res[i]=mpf_sn_log(tableau[i]);
    } 


    TBX_GET_TICK(t1);
    for(i=0; i< LOOPS-1; i++){ 
      res[i]=sn_log(tableau[i]);
    } 
    TBX_GET_TICK(t2);
    deltat = TBX_TICK_RAW_DIFF(t1, t2); 
    printf("%28s :  %lld ticks \n", "scs_log ", deltat);
    printf("\n");




    TBX_GET_TICK(t1);
    for(i=0; i< LOOPS-1; i++){ 
      res[i]=mpf_sn_log(tableau[i]);
    } 
    TBX_GET_TICK(t2);
    deltat = TBX_TICK_RAW_DIFF(t1, t2); 
    printf("%28s :  %lld ticks \n", "mpf_log " , deltat);


    free_mpf_cst();
}

#else
main(){
    fprintf(stderr,"No GMP detected on your system\n");
}
#endif /*HAVE_GMP_H*/
