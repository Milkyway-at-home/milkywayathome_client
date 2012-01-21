/*
 * Correctly rounded logarithm
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
#include <stdio.h>
#include "log_accurate.h"

/*
 *  1) First reduction: exponent extraction
 *         E
 *  x = 2^   .(1+f)    with  0 <= f < 1
 *
 *  log(x) = E.log(2) + log(1+f) where:
 *     - log(2)   is tabulated
 *     - log(1+f) need to be evaluated
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




void scs_log(scs_ptr res, db_number y, int E){
  scs_t R, sc_ln2_times_E, res1, addi;
  scs_ptr ti, inv_wi;
  db_number z, wi;
  int i;


#if EVAL_PERF
  crlibm_second_step_taken++;
#endif


  /* to normalize y.d and round to nearest      */
  /* + (1-trunc(sqrt(2.)/2 * 2^(4))*2^(-4) )+2.^(-(4+1))*/
  z.d = y.d + norm_number.d;
  i = (z.i[HI] & 0x000fffff);
  i = i >> 16; /* 0<= i <=11 */


  wi.d = ((double)(11+i))*0.0625;

  /* (1+f-w_i) */
  y.d -= wi.d;

  /* Table reduction */
  ti     = table_ti_ptr[i];
  inv_wi = table_inv_wi_ptr[i];

  /* R = (1+f-w_i)/w_i */
  scs_set_d(R, y.d);
  scs_mul(R, R, inv_wi);


  /*
   * Polynomial evaluation of log(1 + R) with an error less than 2^(-130)
   */

  scs_mul(res1, constant_poly_ptr[0], R);
  for(i=1; i<20; i++){
    scs_add(addi, constant_poly_ptr[i], res1);
    scs_mul(res1, addi, R);
  }

  if(E==0){
    scs_add(res, res1, ti);
  }else{
    /* sc_ln2_times_E = E*log(2)  */
    scs_set(sc_ln2_times_E, sc_ln2_ptr);

    if (E >= 0){
      scs_mul_ui(sc_ln2_times_E, (unsigned int) E);
    }else{
      scs_mul_ui(sc_ln2_times_E, (unsigned int) -E);
      sc_ln2_times_E->sign = -1;
    }
    scs_add(addi, res1, ti);
    scs_add(res, addi, sc_ln2_times_E);
  }
}

