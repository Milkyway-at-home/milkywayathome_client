/*
 * Correctly rounded base 2 logarithm
 *
 * Author : David Defour
 *
 * This file is part of the crlibm library developed by the Arenaire
 * project at Ecole Normale Superieure de Lyon
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
*/

#include "log2_accurate.h"

/*
 *  1) First reduction: exponent extraction
 *         E
 *  x = 2^   .(1+f)    with  0 <= f < 1
 *
 *  log2(x) = E + log2(1+f) where:
 *     - log2(1+f) need to be evalute
 *
 *
 *  2) Avoiding accuracy problem when E=-1 by testing
 *
 *    if (1+f >= sqrt(2)) then
 *        1+f = (1+f)/2;  E = E+1;
 *    and,
 *        log2(x) = (E+1) + log2((1+f)/2)
 *
 *    so now:      sqrt(2)/2 <= (1+f) < sqrt(2)
 *
 *
 *  3) Second reduction: tabular reduction
 *                   -4
 *   wi = 1 + i. 2^
 *                                      1
 *   log2(1+f) = log2(wi) + log2 ( 1 + --- . (1 + f - wi) )
 *                                     wi
 *
 *   then |(1+f-wi)/wi| <= 2^-5 if we use rounded to nearest.
 *
 *  4) Computation:
 *   a) Table lookup of:
 *        - ti    = log2(wi)
 *        - inv_wi = 1/(wi)
 *   b) Polynomial evaluation of:
 *        - P(R) ~ log2(1 + R), where R = (1+f-wi) * inv_wi
 *
 *                 -5
 *   with  |R| < 2^
 *
 *
 *  5) Reconstruction:
 *   log2(x) = E + t_i + P(R)
 *
 *
 *   Note 1:
 *	To guarantee log2(2^n)=n, where 2^n is normal, the rounding
 *	mode must set to Round-to-Nearest.
 *
 *   Special cases:
 *	log2(x) is NaN with signal if x < 0;
 *	log2(+INF) is +INF with no signal; log2(0) is -INF with signal;
 *	log2(NaN) is that NaN with no signal;
 *	log2(2^N) = N
 *
 */
#define SQRT_2 1.4142135623730950489e0


/*************************************************************
 *************************************************************
 *               ROUNDED  TO NEAREST
 *************************************************************
 *************************************************************/
double log2_rn(double x) {
  scs_t R, res1, sc_exp;
  scs_ptr inv_wi, ti;

  db_number nb, nb2, wi, resd;
  int i, E=0;

  nb.d = x;
  /* Filter cases */
  if (nb.i[HI] < 0x00100000){        /* x < 2^(-1022)    */
    if (((nb.i[HI] & 0x7fffffff)|nb.i[LO])==0)
      return 1.0/0.0;                       /* log(+/-0) = -Inf */
    if (nb.i[HI] < 0)
      return (x-x)/0;                       /* log(-x) = Nan    */

    /* Subnormal number */
    E    -= (SCS_NB_BITS*2); /* keep in mind that x is a subnormal number */
    nb.d *=SCS_RADIX_TWO_DOUBLE;  /* make x as normal number     */
    /* We may just want add 2 to the scs number.index */
    /* may be .... we will see */
  }
  if (nb.i[HI] >= 0x7ff00000)
    return x+x;                             /* Inf or Nan       */

  /* find n, nb.d such that sqrt(2)/2 < nb.d < sqrt(2) */
  E += (nb.i[HI]>>20)-1023;
  nb.i[HI] =  (nb.i[HI] & 0x000fffff) | 0x3ff00000;
  if (nb.d > SQRT_2){
    nb.d *= 0.5;
    E++;
  }

  scs_set_si(sc_exp, E);

  /* to normalize nb.d and round to nearest      */
  /* +((2^4 - trunc(sqrt(2)/2) *2^4 )*2 + 1)/2^5 */
  nb2.d = nb.d + norm_number.d;
  i = (nb2.i[HI] & 0x000fffff);
  i = i >> 16; /* 0<= i <=11 */

  wi.d = (11+i)*(double)0.6250e-1;

  /* (1+f-w_i) */
  nb.d -= wi.d;

  /* Table reduction */
  ti     = table_ti_ptr[i];
  inv_wi = table_inv_wi_ptr[i];


  /* R = (1+f-w_i)/w_i */
  scs_set_d(R, nb.d);
  scs_mul(R, R, inv_wi);


  /*
   * Polynomial evaluation of log2(1 + R) with an error less than 2^(-130)
   */
  scs_mul(res1, constant_poly_ptr[0], R);
  for(i=1; i<20; i++){
    scs_add(res1, constant_poly_ptr[i], res1);
    scs_mul(res1, res1, R);
  }
  scs_add(res1, res1, ti);
  scs_add(res1, res1, sc_exp);

  scs_get_d(&resd.d, res1);

  return resd.d;
}





/*************************************************************
 *************************************************************
 *               ROUNDED  TOWARD  -INFINITY
 *************************************************************
 *************************************************************/
double log2_rd(double x) {
  scs_t R, res1, sc_exp;
  scs_ptr inv_wi, ti;

  db_number nb, nb2, wi, resd;
  int i, E=0;

  nb.d = x;
  /* Filter cases */
  if (nb.i[HI] < 0x00100000){        /* x < 2^(-1022)    */
    if (((nb.i[HI] & 0x7fffffff)|nb.i[LO])==0)
      return 1.0/0.0;                       /* log(+/-0) = -Inf */
    if (nb.i[HI] < 0)
      return (x-x)/0;                       /* log(-x) = Nan    */

    /* Subnormal number */
    E    -= (SCS_NB_BITS*2); /* keep in mind that x is a subnormal number */
    nb.d *=SCS_RADIX_TWO_DOUBLE;  /* make x as normal number     */
    /* We may just want add 2 to the scs number.index */
    /* may be .... we will see */
  }
  if (nb.i[HI] >= 0x7ff00000)
    return x+x;                             /* Inf or Nan       */

  /* find n, nb.d such that sqrt(2)/2 < nb.d < sqrt(2) */
  E += (nb.i[HI]>>20)-1023;
  nb.i[HI] =  (nb.i[HI] & 0x000fffff) | 0x3ff00000;
  if (nb.d > SQRT_2){
    nb.d *= 0.5;
    E++;
  }

  scs_set_si(sc_exp, E);

  /* to normalize nb.d and round to nearest      */
  /* +((2^4 - trunc(sqrt(2)/2) *2^4 )*2 + 1)/2^5 */
  nb2.d = nb.d + norm_number.d;
  i = (nb2.i[HI] & 0x000fffff);
  i = i >> 16; /* 0<= i <=11 */

  wi.d = (11+i)*(double)0.6250e-1;

  /* (1+f-w_i) */
  nb.d -= wi.d;

  /* Table reduction */
  ti     = table_ti_ptr[i];
  inv_wi = table_inv_wi_ptr[i];


  /* R = (1+f-w_i)/w_i */
  scs_set_d(R, nb.d);
  scs_mul(R, R, inv_wi);


  /*
   * Polynomial evaluation of log2(1 + R) with an error less than 2^(-130)
   */
  scs_mul(res1, constant_poly_ptr[0], R);
  for(i=1; i<20; i++){
    scs_add(res1, constant_poly_ptr[i], res1);
    scs_mul(res1, res1, R);
  }
  scs_add(res1, res1, ti);
  scs_add(res1, res1, sc_exp);

  scs_get_d_minf(&resd.d, res1);

  return resd.d;
}






/*************************************************************
 *************************************************************
 *               ROUNDED  TOWARD  +INFINITY
 *************************************************************
 *************************************************************/
double log2_ru(double x) {
  scs_t R, res1, sc_exp;
  scs_ptr inv_wi, ti;

  db_number nb, nb2, wi, resd;
  int i, E=0;

  nb.d = x;
  /* Filter cases */
  if (nb.i[HI] < 0x00100000){        /* x < 2^(-1022)    */
    if (((nb.i[HI] & 0x7fffffff)|nb.i[LO])==0)
      return 1.0/0.0;                       /* log(+/-0) = -Inf */
    if (nb.i[HI] < 0)
      return (x-x)/0;                       /* log(-x) = Nan    */

    /* Subnormal number */
    E    -= (SCS_NB_BITS*2); /* keep in mind that x is a subnormal number */
    nb.d *=SCS_RADIX_TWO_DOUBLE;  /* make x as normal number     */
    /* We may just want add 2 to the scs number.index */
    /* may be .... we will see */
  }
  if (nb.i[HI] >= 0x7ff00000)
    return x+x;                             /* Inf or Nan       */

  /* find n, nb.d such that sqrt(2)/2 < nb.d < sqrt(2) */
  E += (nb.i[HI]>>20)-1023;
  nb.i[HI] =  (nb.i[HI] & 0x000fffff) | 0x3ff00000;
  if (nb.d > SQRT_2){
    nb.d *= 0.5;
    E++;
  }

  scs_set_si(sc_exp, E);

  /* to normalize nb.d and round to nearest      */
  /* +((2^4 - trunc(sqrt(2)/2) *2^4 )*2 + 1)/2^5 */
  nb2.d = nb.d + norm_number.d;
  i = (nb2.i[HI] & 0x000fffff);
  i = i >> 16; /* 0<= i <=11 */

  wi.d = (11+i)*(double)0.6250e-1;

  /* (1+f-w_i) */
  nb.d -= wi.d;

  /* Table reduction */
  ti     = table_ti_ptr[i];
  inv_wi = table_inv_wi_ptr[i];


  /* R = (1+f-w_i)/w_i */
  scs_set_d(R, nb.d);
  scs_mul(R, R, inv_wi);


  /*
   * Polynomial evaluation of log2(1 + R) with an error less than 2^(-130)
   */
  scs_mul(res1, constant_poly_ptr[0], R);
  for(i=1; i<20; i++){
    scs_add(res1, constant_poly_ptr[i], res1);
    scs_mul(res1, res1, R);
  }
  scs_add(res1, res1, ti);
  scs_add(res1, res1, sc_exp);

  scs_get_d_pinf(&resd.d, res1);

  return resd.d;
}

/*************************************************************
 *************************************************************
 *               ROUNDED  TOWARDS ZERO
 *************************************************************
 *************************************************************/
double log2_rz(double x) {
  if (x > 1)
    return log2_rd(x);
  else
    return log2_ru(x);
}
