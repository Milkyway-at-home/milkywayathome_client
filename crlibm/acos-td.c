/*
 * Correctly rounded arccosine
 *
 * Author : Christoph Lauter (ENS Lyon)
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
#include <stdlib.h>
#include "crlibm.h"
#include "crlibm_private.h"
#include "triple-double.h"
#include "acos-td.h"

#define AVOID_FMA 1

void acos_accurate_lower(double *acosh, double *acosm, double *acosl, double x, double xSqh, double xSql, double sign) {
  double highPoly;
  double t1h, t1l, t2h, t2l, t3h, t3l, t4h, t4l, t5h, t5l, t6h, t6l, t7h, t7l;
  double tt1h, tt1l;
  double t8h, t8m, t8l, t9h, t9m, t9l, t10h, t10m, t10l, t11h, t11m, t11l, t12h, t12m, t12l;
  double tt8h, tt8m, tt8l, tt9h, tt9m, tt9l, tt10h, tt10m, tt10l, tt11h, tt11m, tt11l, tt12h, tt12m, tt12l;
  double xCubeh, xCubem, xCubel, tt13h, tt13m, tt13l, t13h, t13m, t13l, polyh, polym, polyl;
  double tt11hover, tt11mover, tt11lover;
  double zw1h, zw1m, zw1l, acoshover, acosmover, acoslover;

#if EVAL_PERF
  crlibm_second_step_taken++;
#endif

  /* Evaluate the polynomial of degree 37
     Its coefficients start at tbl[0]

     p(x) = x + x * x^2 * (c3 + x^2 * (c5 + ...

     We receive x^2 as xSqh + xSql = x * x (exactly)
     in argument

     |x| <= 0.185 = 2^(-2.43)

     Compute monomials 27 to 37 in double precision
     monomials 13 to 25 in double-double and
     1 to 11 in triple-double precision in a
     modified Horner form

  */

  /* Double computations */

#if defined(PROCESSOR_HAS_FMA) && !defined(AVOID_FMA)
  highPoly = FMA(FMA(FMA(FMA(FMA(tbl[33],xSqh,tbl[32]),xSqh,tbl[31]),xSqh,tbl[30]),xSqh,tbl[29]),xSqh,tbl[28]);
#else
  highPoly = tbl[28] + xSqh * (tbl[29] + xSqh * (tbl[30] + xSqh * (tbl[31] + xSqh * (tbl[32] + xSqh * tbl[33]))));
#endif

  /* Double-double computations */

  Mul12(&tt1h,&tt1l,xSqh,highPoly);
  Add22(&t1h,&t1l,tbl[27],0,tt1h,tt1l);

  MulAdd22(&t2h,&t2l,tbl[25],tbl[26],xSqh,xSql,t1h,t1l);
  MulAdd22(&t3h,&t3l,tbl[23],tbl[24],xSqh,xSql,t2h,t2l);
  MulAdd22(&t4h,&t4l,tbl[21],tbl[22],xSqh,xSql,t3h,t3l);
  MulAdd22(&t5h,&t5l,tbl[19],tbl[20],xSqh,xSql,t4h,t4l);
  MulAdd22(&t6h,&t6l,tbl[17],tbl[18],xSqh,xSql,t5h,t5l);
  MulAdd22(&t7h,&t7l,tbl[15],tbl[16],xSqh,xSql,t6h,t6l);

  /* Triple-double computations */

  Mul23(&tt8h,&tt8m,&tt8l,xSqh,xSql,t7h,t7l);                             /* 149 - 48/53 */
  Add33(&t8h,&t8m,&t8l,tbl[12],tbl[13],tbl[14],tt8h,tt8m,tt8l);           /* 145 - 43/53 */
  Mul233(&tt9h,&tt9m,&tt9l,xSqh,xSql,t8h,t8m,t8l);                        /* 139 - 39/53 */
  Add33(&t9h,&t9m,&t9l,tbl[9],tbl[10],tbl[11],tt9h,tt9m,tt9l);            /* 136 - 34/53 */
  Mul233(&tt10h,&tt10m,&tt10l,xSqh,xSql,t9h,t9m,t9l);                     /* 130 - 30/53 */
  Add33(&t10h,&t10m,&t10l,tbl[6],tbl[7],tbl[8],tt10h,tt10m,tt10l);        /* 127 - 25/53 */
  Mul233(&tt11hover,&tt11mover,&tt11lover,xSqh,xSql,t10h,t10m,t10l);      /* 121 - 21/53 */

  Renormalize3(&tt11h,&tt11m,&tt11l,tt11hover,tt11mover,tt11lover);       /* infty - 52/53 */

  Add33(&t11h,&t11m,&t11l,tbl[3],tbl[4],tbl[5],tt11h,tt11m,tt11l);        /* 149 - 47/53 */
  Mul233(&tt12h,&tt12m,&tt12l,xSqh,xSql,t11h,t11m,t11l);                  /* 143 - 43/53 */
  Add33(&t12h,&t12m,&t12l,tbl[0],tbl[1],tbl[2],tt12h,tt12m,tt12l);        /* 140 - 38/53 */

  Mul123(&xCubeh,&xCubem,&xCubel,x,xSqh,xSql);                            /* 154 - 47/53 */

  Mul33(&tt13h,&tt13m,&tt13l,xCubeh,xCubem,xCubel,t12h,t12m,t12l);        /* 136 - 34/53 */
  Add133(&t13h,&t13m,&t13l,x,tt13h,tt13m,tt13l);                          /* 138 - 32/53 */

  Renormalize3(&polyh,&polym,&polyl,t13h,t13m,t13l);                      /* infty - 52/53 */

  /* Reconstruction:

     - Multiply by the inverted sign
     - Add Pi/2 in triple-double
     - Renormalize

  */

  zw1h = -sign * polyh;
  zw1m = -sign * polym;
  zw1l = -sign * polyl;

  Add33(&acoshover,&acosmover,&acoslover,PIHALFH,PIHALFM,PIHALFL,zw1h,zw1m,zw1l);

  Renormalize3(acosh,acosm,acosl,acoshover,acosmover,acoslover);

}



void  acos_accurate_middle(double *acosh, double *acosm, double *acosl, double z, int i, double sign) {
  double highPoly;
  double t1h, t1l, t2h, t2l, t3h, t3l, t4h, t4l, t5h, t5l, t6h, t6l, t7h, t7l, t8h, t8l, t9h, t9l;
  double t10h, t10m, t10l, t11h, t11m, t11l, t12h, t12m, t12l, t13h, t13m, t13l, t14h, t14m, t14l;
  double t15h, t15m, t15l, t16h, t16m, t16l;
  double tt1h, tt1l;
  double tt10h, tt10m, tt10l, tt11h, tt11m, tt11l, tt12h, tt12m, tt12l;
  double tt13h, tt13m, tt13l, tt14h, tt14m, tt14l, tt15h, tt15m, tt15l, tt16h, tt16m, tt16l;
  double polyh, polym, polyl, tt13hover, tt13mover, tt13lover;
  double zw1h, zw1m, zw1l, acoshover, acosmover, acoslover;

#if EVAL_PERF
  crlibm_second_step_taken++;
#endif

  /* Evaluate the polynomial of degree 35
     Its coefficients start at tbl[i+1]
     Evaluate degrees 35 to 20 in double precision,
     degrees 20 to 7 in double-double precision and
     finally degrees 6 to 1 in triple-double.
     The constant coefficient is a double-double, the
     computations are nevertheless in triple-double
  */

  /* Double computations */

#if defined(PROCESSOR_HAS_FMA) && !defined(AVOID_FMA)
  highPoly = FMA(FMA(FMA(FMA(FMA(FMA(FMA(FMA(FMA(FMA(FMA(FMA(FMA(FMA(FMA(FMA(FMA(FMA(FMA(
             tbl[i+58] ,z,tbl[i+57]),z,tbl[i+56]),z,tbl[i+55]),z,tbl[i+54]),z,
             tbl[i+53]),z,tbl[i+52]),z,tbl[i+51]),z,tbl[i+50]),z,tbl[i+49]),z,
             tbl[i+48]),z,tbl[i+47]),z,tbl[i+46]),z,tbl[i+45]),z,tbl[i+44]),z,
	     tbl[i+43]),z,tbl[i+42]),z,tbl[i+41]),z,tbl[i+40]),z,tbl[i+39]);

#else
  highPoly = tbl[i+39] + z * (tbl[i+40] + z * (tbl[i+41] + z * (tbl[i+42] + z * (
             tbl[i+43] + z * (tbl[i+44] + z * (tbl[i+45] + z * (tbl[i+46] + z * (
             tbl[i+47] + z * (tbl[i+48] + z * (tbl[i+49] + z * (tbl[i+50] + z * (
             tbl[i+51] + z * (tbl[i+52] + z * (tbl[i+53] + z * (tbl[i+54] + z * (
             tbl[i+55] + z * (tbl[i+56] + z * (tbl[i+57] + z * tbl[i+58]))))))))))))))))));
#endif


  /* Double-double computations */

  Mul12(&tt1h,&tt1l,z,highPoly);
  Add22(&t1h,&t1l,tbl[i+37],tbl[i+38],tt1h,tt1l);

  MulAdd212(&t2h,&t2l,tbl[i+35],tbl[i+36],z,t1h,t1l);
  MulAdd212(&t3h,&t3l,tbl[i+33],tbl[i+34],z,t2h,t2l);
  MulAdd212(&t4h,&t4l,tbl[i+31],tbl[i+32],z,t3h,t3l);
  MulAdd212(&t5h,&t5l,tbl[i+29],tbl[i+30],z,t4h,t4l);
  MulAdd212(&t6h,&t6l,tbl[i+27],tbl[i+28],z,t5h,t5l);
  MulAdd212(&t7h,&t7l,tbl[i+25],tbl[i+26],z,t6h,t6l);
  MulAdd212(&t8h,&t8l,tbl[i+23],tbl[i+24],z,t7h,t7l);
  MulAdd212(&t9h,&t9l,tbl[i+21],tbl[i+22],z,t8h,t8l);

  /* Triple-double computations */

  Mul123(&tt10h,&tt10m,&tt10l,z,t9h,t9l);                                          /* 154 - 47/53 */
  Add33(&t10h,&t10m,&t10l,tbl[i+18],tbl[i+19],tbl[i+20],tt10h,tt10m,tt10l);        /* 144 - 42/53 */
  Mul133(&tt11h,&tt11m,&tt11l,z,t10h,t10m,t10l);                                   /* 142 - 38/53 */
  Add33(&t11h,&t11m,&t11l,tbl[i+15],tbl[i+16],tbl[i+17],tt11h,tt11m,tt11l);        /* 136 - 33/53 */
  Mul133(&tt12h,&tt12m,&tt12l,z,t11h,t11m,t11l);                                   /* 133 - 28/53 */
  Add33(&t12h,&t12m,&t12l,tbl[i+12],tbl[i+13],tbl[i+14],tt12h,tt12m,tt12l);        /* 125 - 23/53 */
  Mul133(&tt13hover,&tt13mover,&tt13lover,z,t12h,t12m,t12l);                       /* 123 - 18/53 */

  Renormalize3(&tt13h,&tt13m,&tt13l,tt13hover,tt13mover,tt13lover);                /* infty - 52/53 */

  Add33(&t13h,&t13m,&t13l,tbl[i+9],tbl[i+10],tbl[i+11],tt13h,tt13m,tt13l);         /* 149 - 47/53 */
  Mul133(&tt14h,&tt14m,&tt14l,z,t13h,t13m,t13l);                                   /* 147 - 42/53 */
  Add33(&t14h,&t14m,&t14l,tbl[i+6],tbl[i+7],tbl[i+8],tt14h,tt14m,tt14l);           /* 139 - 37/53 */
  Mul133(&tt15h,&tt15m,&tt15l,z,t14h,t14m,t14l);                                   /* 137 - 32/53 */
  Add33(&t15h,&t15m,&t15l,tbl[i+3],tbl[i+4],tbl[i+5],tt15h,tt15m,tt15l);           /* 129 - 28/53 */
  Mul133(&tt16h,&tt16m,&tt16l,z,t15h,t15m,t15l);                                   /* 128 - 23/53 */
  Add233(&t16h,&t16m,&t16l,tbl[i+1],tbl[i+2],tt16h,tt16m,tt16l);                   /* 126 - 19/53 */

  Renormalize3(&polyh,&polym,&polyl,t16h,t16m,t16l);                               /* infty - 52/53 */

  /* Reconstruction:

     - Multiply by the inverted sign
     - Add Pi/2 in triple-double
     - Renormalize

  */

  zw1h = -sign * polyh;
  zw1m = -sign * polym;
  zw1l = -sign * polyl;

  Add33(&acoshover,&acosmover,&acoslover,PIHALFH,PIHALFM,PIHALFL,zw1h,zw1m,zw1l);

  Renormalize3(acosh,acosm,acosl,acoshover,acosmover,acoslover);

}


void acos_accurate_higher(double *acosh, double *acosm, double *acosl, double z, double sign) {
  double highPoly;
  double tt1h, tt1l;
  double t1h, t1l, t2h, t2l, t3h, t3l, t4h, t4l, t5h, t5l, t6h, t6l, t7h, t7l, t8h, t8l;
  double tt10h, tt10m, tt10l, tt11h, tt11m, tt11l, tt12h, tt12m, tt12l, tt13h, tt13m, tt13l;
  double tt14h, tt14m, tt14l, tt15h, tt15m, tt15l, tt16h, tt16m, tt16l, tt17h, tt17m, tt17l;
  double t9h, t9l, t10h, t10m, t10l, t11h, t11m, t11l, t12h, t12m, t12l, t13h, t13m, t13l;
  double t14h, t14m, t14l, t15h, t15m, t15l, t16h, t16m, t16l, t17h, t17m, t17l;
  double tt18h, tt18m, tt18l, polyh, polym, polyl;
  double sqrtzh, sqrtzm, sqrtzl, twoZ, pTimesSh, pTimesSm, pTimesSl;
  double allh, allm, alll;
  double tt13hover, tt13mover, tt13lover, tt16hover, tt16mover, tt16lover;
  double polyhover, polymover, polylover;

#if EVAL_PERF
  crlibm_second_step_taken++;
#endif

  /* We evaluate acos(x) with x > 0 as

     acos(x) = -1 * f(z) * sqrt(2*z)

     with z = 1 - x and

     f(z) = (asin(z) - Pi/2) / sqrt(2*z)

     f(z) is approximated by p(z)

     The polynomial p(z) is of degree 29
     Its coefficients start at tbl[TBLIDX10]
     Coefficients for degrees 29 to 18 are in double precision,
     for degrees 17 to 9 in double-double precision and
     finally for degrees 8 to 1 in triple-double.
     The constant coefficient (-1) is not stored in the table,
     the computations are nevertheless in triple-double
     We evaluate the monomials in the precision in which
     the correspondant coefficients are stored
     The coefficients' values decrease very quickly
     so even with |z| < 2^-2.18 we can compute degree 18
     already in double precision

     Compute than sqrt(2*z) as a triple-double
     multiply in triple-double.

  */

  /* Double computations */

#if defined(PROCESSOR_HAS_FMA) && !defined(AVOID_FMA)
  highPoly = FMA(FMA(FMA(FMA(FMA(FMA(FMA(FMA(FMA(FMA(FMA(
             tbl[TBLIDX10+53] ,z,tbl[TBLIDX10+52]),z,tbl[TBLIDX10+51]),z,
             tbl[TBLIDX10+50]),z,tbl[TBLIDX10+49]),z,tbl[TBLIDX10+48]),z,
             tbl[TBLIDX10+47]),z,tbl[TBLIDX10+46]),z,tbl[TBLIDX10+45]),z,
             tbl[TBLIDX10+44]),z,tbl[TBLIDX10+43]),z,tbl[TBLIDX10+42]);
#else
  highPoly = tbl[TBLIDX10+42] + z * (tbl[TBLIDX10+43] + z * (tbl[TBLIDX10+44] + z * (
             tbl[TBLIDX10+45] + z * (tbl[TBLIDX10+46] + z * (tbl[TBLIDX10+47] + z * (
             tbl[TBLIDX10+48] + z * (tbl[TBLIDX10+49] + z * (tbl[TBLIDX10+50] + z * (
             tbl[TBLIDX10+51] + z * (tbl[TBLIDX10+52] + z *  tbl[TBLIDX10+53]))))))))));
#endif

  /* Double-double computations */

  Mul12(&tt1h,&tt1l,z,highPoly);
  Add22(&t1h,&t1l,tbl[TBLIDX10+40],tbl[TBLIDX10+41],tt1h,tt1l);

  MulAdd212(&t2h,&t2l,tbl[TBLIDX10+38],tbl[TBLIDX10+39],z,t1h,t1l);
  MulAdd212(&t3h,&t3l,tbl[TBLIDX10+36],tbl[TBLIDX10+37],z,t2h,t2l);
  MulAdd212(&t4h,&t4l,tbl[TBLIDX10+34],tbl[TBLIDX10+35],z,t3h,t3l);
  MulAdd212(&t5h,&t5l,tbl[TBLIDX10+32],tbl[TBLIDX10+33],z,t4h,t4l);
  MulAdd212(&t6h,&t6l,tbl[TBLIDX10+30],tbl[TBLIDX10+31],z,t5h,t5l);
  MulAdd212(&t7h,&t7l,tbl[TBLIDX10+28],tbl[TBLIDX10+29],z,t6h,t6l);
  MulAdd212(&t8h,&t8l,tbl[TBLIDX10+26],tbl[TBLIDX10+27],z,t7h,t7l);
  MulAdd212(&t9h,&t9l,tbl[TBLIDX10+24],tbl[TBLIDX10+25],z,t8h,t8l);

  /* Triple-double computations */

  Mul123(&tt10h,&tt10m,&tt10l,z,t9h,t9l);                                                         /* 154 - 47/53 */
  Add33(&t10h,&t10m,&t10l,tbl[TBLIDX10+21],tbl[TBLIDX10+22],tbl[TBLIDX10+23],tt10h,tt10m,tt10l);  /* 144 - 42/53 */
  Mul133(&tt11h,&tt11m,&tt11l,z,t10h,t10m,t10l);                                                  /* 142 - 37/53 */
  Add33(&t11h,&t11m,&t11l,tbl[TBLIDX10+18],tbl[TBLIDX10+19],tbl[TBLIDX10+20],tt11h,tt11m,tt11l);  /* 134 - 32/53 */
  Mul133(&tt12h,&tt12m,&tt12l,z,t11h,t11m,t11l);                                                  /* 132 - 27/53 */
  Add33(&t12h,&t12m,&t12l,tbl[TBLIDX10+15],tbl[TBLIDX10+16],tbl[TBLIDX10+17],tt12h,tt12m,tt12l);  /* 124 - 22/53 */
  Mul133(&tt13hover,&tt13mover,&tt13lover,z,t12h,t12m,t12l);                                      /* 122 - 17/53 */

  Renormalize3(&tt13h,&tt13m,&tt13l,tt13hover,tt13mover,tt13lover);                               /* infty - 52/53 */

  Add33(&t13h,&t13m,&t13l,tbl[TBLIDX10+12],tbl[TBLIDX10+13],tbl[TBLIDX10+14],tt13h,tt13m,tt13l);  /* 149 - 47/53 */
  Mul133(&tt14h,&tt14m,&tt14l,z,t13h,t13m,t13l);                                                  /* 147 - 42/53 */
  Add33(&t14h,&t14m,&t14l,tbl[TBLIDX10+9],tbl[TBLIDX10+10],tbl[TBLIDX10+11],tt14h,tt14m,tt14l);   /* 139 - 37/53 */
  Mul133(&tt15h,&tt15m,&tt15l,z,t14h,t14m,t14l);                                                  /* 137 - 32/53 */
  Add33(&t15h,&t15m,&t15l,tbl[TBLIDX10+6],tbl[TBLIDX10+7],tbl[TBLIDX10+8],tt15h,tt15m,tt15l);     /* 129 - 27/53 */
  Mul133(&tt16hover,&tt16mover,&tt16lover,z,t15h,t15m,t15l);                                      /* 127 - 22/53 */

  Renormalize3(&tt16h,&tt16m,&tt16l,tt16hover,tt16mover,tt16lover);                               /* infty - 52/53 */

  Add33(&t16h,&t16m,&t16l,tbl[TBLIDX10+3],tbl[TBLIDX10+4],tbl[TBLIDX10+5],tt16h,tt16m,tt16l);     /* 149 - 47/53 */
  Mul133(&tt17h,&tt17m,&tt17l,z,t16h,t16m,t16l);                                                  /* 147 - 42/53 */
  Add33(&t17h,&t17m,&t17l,tbl[TBLIDX10+0],tbl[TBLIDX10+1],tbl[TBLIDX10+2],tt17h,tt17m,tt17l);     /* 139 - 37/53 */
  Mul133(&tt18h,&tt18m,&tt18l,z,t17h,t17m,t17l);                                                  /* 137 - 32/53 */
  Add133(&polyhover,&polymover,&polylover,-1,tt18h,tt18m,tt18l);                                  /* 136 - 30/53 */

  Renormalize3(&polyh,&polym,&polyl,polyhover,polymover,polylover);                               /* infty - 52/53 */

  /* Compute sqrt(2*z) as a triple-double */

  twoZ = 2 * z;
  Sqrt13(&sqrtzh,&sqrtzm,&sqrtzl,twoZ);                                                           /* 146 - 52/53 */

  /* Multiply p(z) by sqrt(2*z) */

  Mul33(&pTimesSh,&pTimesSm,&pTimesSl,polyh,polym,polyl,sqrtzh,sqrtzm,sqrtzl);                    /* 139 - 48/53 */

  /* Reconstruction:

     If the sign of x in acos(x) was positive:
      - Multiply pTimesSh + pTimesSm + pTimesSl approx f(x) * sqrt(2 * z) by -1
      - Renormalize
      - Return

     Otherwise:
      - Add Pi in triple-double to pTimesSh + pTimesSm + pTimesSl approx f(x) * sqrt(2 * z)
      - Renormalize
      - Return

  */

  if (sign > 0) {

    allh = -1.0 * pTimesSh;
    allm = -1.0 * pTimesSm;
    alll = -1.0 * pTimesSl;                                                                       /* 139 - 48/53 */

  } else {

    Add33(&allh,&allm,&alll,PIH,PIM,PIL,pTimesSh,pTimesSm,pTimesSl);                              /* 130 - 43/53 */

  }

  Renormalize3(acosh,acosm,acosl,allh,allm,alll);                                                 /* infty - 52/53 */

}


double acos_rn(double x) {
  db_number xdb;
  double sign, z, acosh, acosm, acosl;
  int i;
  double xSqh, xSql;
  double tt1h, tt1l;
  double tt6h, tt6l;
  double t1h, t1l, t2h, t2l, t3h, t3l, t4h, t4l, t5h, t5l, t6h, t6l;
  double t7h, t7l, t8h, t8l, polyh, polyl, twoZ, sqrtzh, sqrtzl;
  double pTimesSh, pTimesSl, highPoly, xCubeh, xCubel;
  double tmp1, tmp2, tmp3, tmp4, tmp5;
  double zw1h, zw1l;

  /*
#if CRLIBM_REQUIRES_ROUNDING_MODE_CHANGE
  SAVE_STATE_AND_SET_RNDOUBLE
#endif
  */

  /* Transform the argument into integer */
  xdb.d = x;

  /* Special case handling */

  /* Exact algebraic case x = 1, acos(1) = 0 */

  if (x == 1.0) return 0.0;

  /* Strip off the sign of argument x */
  if (xdb.i[HI] & 0x80000000) sign = -1; else sign = 1;
  xdb.i[HI] &= 0x7fffffff;

  /* acos is defined on -1 <= x <= 1, elsewhere it is NaN */
  if ((xdb.i[HI] > 0x3ff00000) || ((xdb.i[HI] == 0x3ff00000) && (xdb.i[LO] != 0x00000000))) {
    return (x-x)/0.0;    /* return NaN */
  }

  /* If |x| < 2^(-120) we have

     round(acos(x)) = round(pi/2)

     So we can decide the rounding without any computation
  */
  if (xdb.i[HI] < 0x38700000) {
    return PIHALFDOUBLERN;
  }

  /* Recast x */
  x = xdb.d;

  /* Find correspondant interval and compute index to the table
     We start by filtering the two special cases around 0 and 1
  */

  if (xdb.i[HI] < BOUND1) {
    /* Special interval 0..BOUND1
       The polynomial has no even monomials
       We must prove extra accuracy in the interval 0..sin(2^(-18))
    */

    /* Quick phase starts */

    /* Compute square of x for both quick and accurate phases */
    Mul12(&xSqh,&xSql,x,x);

    tmp4 = tbl[3];
    tmp5 = tbl[4];
    t4h = tmp4;
    t4l = tmp5;
    if (xdb.i[HI] > EXTRABOUND) {
      /* Double precision evaluation */
#if defined(PROCESSOR_HAS_FMA) && !defined(AVOID_FMA)
      highPoly = FMA(FMA(FMA(FMA(tbl[23],xSqh,tbl[21]),xSqh,tbl[19]),xSqh,tbl[17]),xSqh,tbl[15]);
#else
      highPoly = tbl[15] + xSqh * (tbl[17] + xSqh * (tbl[19] + xSqh * (tbl[21] + xSqh * tbl[23])));
#endif

      /* Double-double precision evaluation */
      Mul12(&tt1h,&tt1l,xSqh,highPoly);
      Add22(&t1h,&t1l,tbl[12],tbl[13],tt1h,tt1l);

      MulAdd212(&t2h,&t2l,tbl[9],tbl[10],xSqh,t1h,t1l);
      MulAdd212(&t3h,&t3l,tbl[6],tbl[7],xSqh,t2h,t2l);
      MulAdd22(&t4h,&t4l,tmp4,tmp5,xSqh,xSql,t3h,t3l);
    }

    MulAdd22(&t5h,&t5l,tbl[0],tbl[1],xSqh,xSql,t4h,t4l);

    Mul122(&xCubeh,&xCubel,x,xSqh,xSql);
    Mul22(&tt6h,&tt6l,xCubeh,xCubel,t5h,t5l);

    Add12(tmp1,tmp2,x,tt6h);
    tmp3 = tmp2 + tt6l;
    Add12(polyh,polyl,tmp1,tmp3);

    /* Reconstruction:

       - Multiply by the inverted sign
       - Add Pi/2 in double-double precision

    */

    zw1h = -sign * polyh;
    zw1l = -sign * polyl;

    Add22(&acosh,&acosm,PIHALFH,PIHALFM,zw1h,zw1l);

    /* Rounding test
       The RN rounding constant is at tbl[34]
    */
    if(acosh == (acosh + (acosm * tbl[34])))
      return acosh;

    /* Launch accurate phase */

    acos_accurate_lower(&acosh,&acosm,&acosl,x,xSqh,xSql,sign);

    ReturnRoundToNearest3(acosh,acosm,acosl);
  }

  if (xdb.i[HI] >= BOUND9) {
    /* Special interval BOUND9..1
       We use an asymptotic development of arcsin in sqrt(1 - x)
    */

    /* Argument reduction for quick and accurate phase
       z = 1 - x
       The operation is exact as per Sterbenz' lemma
    */

    z = 1 - x;

    /* Quick phase starts */

    /* Double precision evaluation */
#if defined(PROCESSOR_HAS_FMA) && !defined(AVOID_FMA)
    highPoly = FMA(FMA(FMA(FMA(FMA(FMA(FMA(FMA(FMA(
	       tbl[TBLIDX10+42] ,z,tbl[TBLIDX10+40]),z,tbl[TBLIDX10+38]),z,
               tbl[TBLIDX10+36]),z,tbl[TBLIDX10+34]),z,tbl[TBLIDX10+32]),z,
               tbl[TBLIDX10+30]),z,tbl[TBLIDX10+28]),z,tbl[TBLIDX10+26]),z,
               tbl[TBLIDX10+24]);
#else
    highPoly = tbl[TBLIDX10+24] + z * (tbl[TBLIDX10+26] + z * (tbl[TBLIDX10+28] + z * (
	       tbl[TBLIDX10+30] + z * (tbl[TBLIDX10+32] + z * (tbl[TBLIDX10+34] + z * (
	       tbl[TBLIDX10+36] + z * (tbl[TBLIDX10+38] + z * (tbl[TBLIDX10+40] + z *
               tbl[TBLIDX10+42]))))))));
#endif

    /* Double-double precision evaluation */
    Mul12(&tt1h,&tt1l,z,highPoly);
    Add22(&t1h,&t1l,tbl[TBLIDX10+21],tbl[TBLIDX10+22],tt1h,tt1l);

    MulAdd212(&t2h,&t2l,tbl[TBLIDX10+18],tbl[TBLIDX10+19],z,t1h,t1l);
    MulAdd212(&t3h,&t3l,tbl[TBLIDX10+15],tbl[TBLIDX10+16],z,t2h,t2l);
    MulAdd212(&t4h,&t4l,tbl[TBLIDX10+12],tbl[TBLIDX10+13],z,t3h,t3l);
    MulAdd212(&t5h,&t5l,tbl[TBLIDX10+9],tbl[TBLIDX10+10],z,t4h,t4l);
    MulAdd212(&t6h,&t6l,tbl[TBLIDX10+6],tbl[TBLIDX10+7],z,t5h,t5l);
    MulAdd212(&t7h,&t7l,tbl[TBLIDX10+3],tbl[TBLIDX10+4],z,t6h,t6l);
    MulAdd212(&t8h,&t8l,tbl[TBLIDX10+0],tbl[TBLIDX10+1],z,t7h,t7l);
    MulAdd212(&polyh,&polyl,-1,0,z,t8h,t8l);

    /* Compute sqrt(2*z) as a double-double */

    twoZ = 2 * z;
    sqrt12(&sqrtzh,&sqrtzl,twoZ);

    /* Multiply p(z) by sqrt(2*z) and add Pi/2 */

    Mul22(&pTimesSh,&pTimesSl,polyh,polyl,sqrtzh,sqrtzl);


    /* Reconstruction:

     If the sign of x in acos(x) was positive:
      - Multiply pTimesSh + pTimesSl approx f(x) * sqrt(2 * z) by -1
      - Return

     Otherwise:
      - Add Pi in triple-double to pTimesSh + pTimesSl approx f(x) * sqrt(2 * z)
      - Return

    */

    if (sign > 0) {

      acosh = -1.0 * pTimesSh;
      acosm = -1.0 * pTimesSl;

    } else {

      Add22(&acosh,&acosm,PIH,PIM,pTimesSh,pTimesSl);

    }

    /* Rounding test
       The RN rounding constant is at tbl[TBLIDX10+54]
    */

    if(acosh == (acosh + (acosm * tbl[TBLIDX10+54])))
      return acosh;

    /* Launch accurate phase */

    acos_accurate_higher(&acosh,&acosm,&acosl,z,sign);

    ReturnRoundToNearest3(acosh,acosm,acosl);
  }

  /* General 8 main intervals
     We can already suppose that BOUND1 <= x <= BOUND9
  */

  if (xdb.i[HI] < BOUND5) {
    if (xdb.i[HI] < BOUND3) {
      if (xdb.i[HI] < BOUND2) i = TBLIDX2; else i = TBLIDX3;
    } else {
      if (xdb.i[HI] < BOUND4) i = TBLIDX4; else i = TBLIDX5;
    }
  } else {
    if (xdb.i[HI] < BOUND7) {
      if (xdb.i[HI] < BOUND6) i = TBLIDX6; else i = TBLIDX7;
    } else {
      if (xdb.i[HI] < BOUND8) i = TBLIDX8; else i = TBLIDX9;
    }
  }

  /* Argument reduction
     i points to the interval midpoint value in the table
  */
  z = x - tbl[i];

  /* Quick phase starts */

  /* Double precision evaluation */

#if defined(PROCESSOR_HAS_FMA) && !defined(AVOID_FMA)
  highPoly = FMA(FMA(FMA(FMA(FMA(FMA(FMA(
	     tbl[i+35] ,z,tbl[i+33]),z,tbl[i+31]),z,tbl[i+29]),z,
             tbl[i+27]),z,tbl[i+25]),z,tbl[i+23]),z,tbl[i+21]);
#else
  highPoly = tbl[i+21] + z * (tbl[i+23] + z * (tbl[i+25] + z * (
             tbl[i+27] + z * (tbl[i+29] + z * (tbl[i+31] + z * (
             tbl[i+33] + z *  tbl[i+35]))))));
#endif

  /* Double-double precision evaluation */

  Mul12(&tt1h,&tt1l,z,highPoly);
  Add22(&t1h,&t1l,tbl[i+18],tbl[i+19],tt1h,tt1l);

  MulAdd212(&t2h,&t2l,tbl[i+15],tbl[i+16],z,t1h,t1l);
  MulAdd212(&t3h,&t3l,tbl[i+12],tbl[i+13],z,t2h,t2l);
  MulAdd212(&t4h,&t4l,tbl[i+9],tbl[i+10],z,t3h,t3l);
  MulAdd212(&t5h,&t5l,tbl[i+6],tbl[i+7],z,t4h,t4l);
  MulAdd212(&t6h,&t6l,tbl[i+3],tbl[i+4],z,t5h,t5l);
  MulAdd212(&polyh,&polyl,tbl[i+1],tbl[i+2],z,t6h,t6l);

  /* Reconstruction:

     - Multiply by the inverted sign
     - Add Pi/2 in double-double precision

  */

  zw1h = -sign * polyh;
  zw1l = -sign * polyl;

  Add22(&acosh,&acosm,PIHALFH,PIHALFM,zw1h,zw1l);

  /* Rounding test
     The RN rounding constant is at tbl[i+59]
  */
  if(acosh == (acosh + (acosm * tbl[i+59])))
    return acosh;

  /* Launch accurate phase */

  acos_accurate_middle(&acosh,&acosm,&acosl,z,i,sign);

  ReturnRoundToNearest3(acosh,acosm,acosl);
}

double acos_ru(double x) {
  db_number xdb;
  double sign, z, acosh, acosm, acosl;
  int i;
  double xSqh, xSql;
  double tt1h, tt1l;
  double tt6h, tt6l;
  double t1h, t1l, t2h, t2l, t3h, t3l, t4h, t4l, t5h, t5l, t6h, t6l;
  double t7h, t7l, t8h, t8l, polyh, polyl, twoZ, sqrtzh, sqrtzl;
  double pTimesSh, pTimesSl, highPoly, xCubeh, xCubel;
  double tmp1, tmp2, tmp3, tmp4, tmp5;
  double zw1h, zw1l;

  /* Transform the argument into integer */
  xdb.d = x;

  /* Special case handling */

  /* Exact algebraic case x = 1, acos(1) = 0 */

  if (x == 1.0) return 0.0;

  /* Strip off the sign of argument x */
  if (xdb.i[HI] & 0x80000000) sign = -1; else sign = 1;
  xdb.i[HI] &= 0x7fffffff;

  /* acos is defined on -1 <= x <= 1, elsewhere it is NaN */
  if ((xdb.i[HI] > 0x3ff00000) || ((xdb.i[HI] == 0x3ff00000) && (xdb.i[LO] != 0x00000000))) {
    return (x-x)/0.0;    /* return NaN */
  }

  /* If |x| < 2^(-120) we have

     round(acos(x)) = round(pi/2)

     So we can decide the rounding without any computation
  */
  if (xdb.i[HI] < 0x38700000) {
    return PIHALFDOUBLERU;
  }

  /* Recast x */
  x = xdb.d;

  /* Find correspondant interval and compute index to the table
     We start by filtering the two special cases around 0 and 1
  */

  if (xdb.i[HI] < BOUND1) {
    /* Special interval 0..BOUND1
       The polynomial has no even monomials
       We must prove extra accuracy in the interval 0..sin(2^(-18))
    */

    /* Quick phase starts */

    /* Compute square of x for both quick and accurate phases */
    Mul12(&xSqh,&xSql,x,x);

    tmp4 = tbl[3];
    tmp5 = tbl[4];
    t4h = tmp4;
    t4l = tmp5;
    if (xdb.i[HI] > EXTRABOUND) {
      /* Double precision evaluation */
#if defined(PROCESSOR_HAS_FMA) && !defined(AVOID_FMA)
      highPoly = FMA(FMA(FMA(FMA(tbl[23],xSqh,tbl[21]),xSqh,tbl[19]),xSqh,tbl[17]),xSqh,tbl[15]);
#else
      highPoly = tbl[15] + xSqh * (tbl[17] + xSqh * (tbl[19] + xSqh * (tbl[21] + xSqh * tbl[23])));
#endif

      /* Double-double precision evaluation */
      Mul12(&tt1h,&tt1l,xSqh,highPoly);
      Add22(&t1h,&t1l,tbl[12],tbl[13],tt1h,tt1l);

      MulAdd212(&t2h,&t2l,tbl[9],tbl[10],xSqh,t1h,t1l);
      MulAdd212(&t3h,&t3l,tbl[6],tbl[7],xSqh,t2h,t2l);
      MulAdd22(&t4h,&t4l,tmp4,tmp5,xSqh,xSql,t3h,t3l);
    }

    MulAdd22(&t5h,&t5l,tbl[0],tbl[1],xSqh,xSql,t4h,t4l);

    Mul122(&xCubeh,&xCubel,x,xSqh,xSql);
    Mul22(&tt6h,&tt6l,xCubeh,xCubel,t5h,t5l);

    Add12(tmp1,tmp2,x,tt6h);
    tmp3 = tmp2 + tt6l;
    Add12(polyh,polyl,tmp1,tmp3);

    /* Reconstruction:

       - Multiply by the inverted sign
       - Add Pi/2 in double-double precision

    */

    zw1h = -sign * polyh;
    zw1l = -sign * polyl;

    Add22(&acosh,&acosm,PIHALFH,PIHALFM,zw1h,zw1l);

    /* Rounding test
       The RU rounding constant is at tbl[35]
    */
    TEST_AND_RETURN_RU(acosh, acosm, tbl[35]);

    /* Launch accurate phase */

    acos_accurate_lower(&acosh,&acosm,&acosl,x,xSqh,xSql,sign);

    ReturnRoundUpwards3(acosh,acosm,acosl);
  }

  if (xdb.i[HI] >= BOUND9) {
    /* Special interval BOUND9..1
       We use an asymptotic development of arcsin in sqrt(1 - x)
    */

    /* Argument reduction for quick and accurate phase
       z = 1 - x
       The operation is exact as per Sterbenz' lemma
    */

    z = 1 - x;

    /* Quick phase starts */

    /* Double precision evaluation */
#if defined(PROCESSOR_HAS_FMA) && !defined(AVOID_FMA)
    highPoly = FMA(FMA(FMA(FMA(FMA(FMA(FMA(FMA(FMA(
	       tbl[TBLIDX10+42] ,z,tbl[TBLIDX10+40]),z,tbl[TBLIDX10+38]),z,
               tbl[TBLIDX10+36]),z,tbl[TBLIDX10+34]),z,tbl[TBLIDX10+32]),z,
               tbl[TBLIDX10+30]),z,tbl[TBLIDX10+28]),z,tbl[TBLIDX10+26]),z,
               tbl[TBLIDX10+24]);
#else
    highPoly = tbl[TBLIDX10+24] + z * (tbl[TBLIDX10+26] + z * (tbl[TBLIDX10+28] + z * (
	       tbl[TBLIDX10+30] + z * (tbl[TBLIDX10+32] + z * (tbl[TBLIDX10+34] + z * (
	       tbl[TBLIDX10+36] + z * (tbl[TBLIDX10+38] + z * (tbl[TBLIDX10+40] + z *
               tbl[TBLIDX10+42]))))))));
#endif

    /* Double-double precision evaluation */
    Mul12(&tt1h,&tt1l,z,highPoly);
    Add22(&t1h,&t1l,tbl[TBLIDX10+21],tbl[TBLIDX10+22],tt1h,tt1l);

    MulAdd212(&t2h,&t2l,tbl[TBLIDX10+18],tbl[TBLIDX10+19],z,t1h,t1l);
    MulAdd212(&t3h,&t3l,tbl[TBLIDX10+15],tbl[TBLIDX10+16],z,t2h,t2l);
    MulAdd212(&t4h,&t4l,tbl[TBLIDX10+12],tbl[TBLIDX10+13],z,t3h,t3l);
    MulAdd212(&t5h,&t5l,tbl[TBLIDX10+9],tbl[TBLIDX10+10],z,t4h,t4l);
    MulAdd212(&t6h,&t6l,tbl[TBLIDX10+6],tbl[TBLIDX10+7],z,t5h,t5l);
    MulAdd212(&t7h,&t7l,tbl[TBLIDX10+3],tbl[TBLIDX10+4],z,t6h,t6l);
    MulAdd212(&t8h,&t8l,tbl[TBLIDX10+0],tbl[TBLIDX10+1],z,t7h,t7l);
    MulAdd212(&polyh,&polyl,-1,0,z,t8h,t8l);

    /* Compute sqrt(2*z) as a double-double */

    twoZ = 2 * z;
    sqrt12(&sqrtzh,&sqrtzl,twoZ);

    /* Multiply p(z) by sqrt(2*z) and add Pi/2 */

    Mul22(&pTimesSh,&pTimesSl,polyh,polyl,sqrtzh,sqrtzl);


    /* Reconstruction:

     If the sign of x in acos(x) was positive:
      - Multiply pTimesSh + pTimesSl approx f(x) * sqrt(2 * z) by -1
      - Return

     Otherwise:
      - Add Pi in triple-double to pTimesSh + pTimesSl approx f(x) * sqrt(2 * z)
      - Return

    */

    if (sign > 0) {

      acosh = -1.0 * pTimesSh;
      acosm = -1.0 * pTimesSl;

    } else {

      Add22(&acosh,&acosm,PIH,PIM,pTimesSh,pTimesSl);

    }

    /* Rounding test
       The RU rounding constant is at tbl[TBLIDX10+55]
    */
    TEST_AND_RETURN_RU(acosh, acosm, tbl[TBLIDX10+55]);

    /* Launch accurate phase */

    acos_accurate_higher(&acosh,&acosm,&acosl,z,sign);

    ReturnRoundUpwards3(acosh,acosm,acosl);
  }

  /* General 8 main intervals
     We can already suppose that BOUND1 <= x <= BOUND9
  */

  if (xdb.i[HI] < BOUND5) {
    if (xdb.i[HI] < BOUND3) {
      if (xdb.i[HI] < BOUND2) i = TBLIDX2; else i = TBLIDX3;
    } else {
      if (xdb.i[HI] < BOUND4) i = TBLIDX4; else i = TBLIDX5;
    }
  } else {
    if (xdb.i[HI] < BOUND7) {
      if (xdb.i[HI] < BOUND6) i = TBLIDX6; else i = TBLIDX7;
    } else {
      if (xdb.i[HI] < BOUND8) i = TBLIDX8; else i = TBLIDX9;
    }
  }

  /* Argument reduction
     i points to the interval midpoint value in the table
  */
  z = x - tbl[i];

  /* Quick phase starts */

  /* Double precision evaluation */

#if defined(PROCESSOR_HAS_FMA) && !defined(AVOID_FMA)
  highPoly = FMA(FMA(FMA(FMA(FMA(FMA(FMA(
	     tbl[i+35] ,z,tbl[i+33]),z,tbl[i+31]),z,tbl[i+29]),z,
             tbl[i+27]),z,tbl[i+25]),z,tbl[i+23]),z,tbl[i+21]);
#else
  highPoly = tbl[i+21] + z * (tbl[i+23] + z * (tbl[i+25] + z * (
             tbl[i+27] + z * (tbl[i+29] + z * (tbl[i+31] + z * (
             tbl[i+33] + z *  tbl[i+35]))))));
#endif

  /* Double-double precision evaluation */

  Mul12(&tt1h,&tt1l,z,highPoly);
  Add22(&t1h,&t1l,tbl[i+18],tbl[i+19],tt1h,tt1l);

  MulAdd212(&t2h,&t2l,tbl[i+15],tbl[i+16],z,t1h,t1l);
  MulAdd212(&t3h,&t3l,tbl[i+12],tbl[i+13],z,t2h,t2l);
  MulAdd212(&t4h,&t4l,tbl[i+9],tbl[i+10],z,t3h,t3l);
  MulAdd212(&t5h,&t5l,tbl[i+6],tbl[i+7],z,t4h,t4l);
  MulAdd212(&t6h,&t6l,tbl[i+3],tbl[i+4],z,t5h,t5l);
  MulAdd212(&polyh,&polyl,tbl[i+1],tbl[i+2],z,t6h,t6l);

  /* Reconstruction:

     - Multiply by the inverted sign
     - Add Pi/2 in double-double precision

  */

  zw1h = -sign * polyh;
  zw1l = -sign * polyl;

  Add22(&acosh,&acosm,PIHALFH,PIHALFM,zw1h,zw1l);

  /* Rounding test
     The RU rounding constant is at tbl[i+60]
  */
  TEST_AND_RETURN_RU(acosh, acosm, tbl[i+60]);

  /* Launch accurate phase */

  acos_accurate_middle(&acosh,&acosm,&acosl,z,i,sign);

  ReturnRoundUpwards3(acosh,acosm,acosl);
}



double acos_rd(double x) {
  db_number xdb;
  double sign, z, acosh, acosm, acosl;
  int i;
  double xSqh, xSql;
  double tt1h, tt1l;
  double tt6h, tt6l;
  double t1h, t1l, t2h, t2l, t3h, t3l, t4h, t4l, t5h, t5l, t6h, t6l;
  double t7h, t7l, t8h, t8l, polyh, polyl, twoZ, sqrtzh, sqrtzl;
  double pTimesSh, pTimesSl, highPoly, xCubeh, xCubel;
  double tmp1, tmp2, tmp3, tmp4, tmp5;
  double zw1h, zw1l;

  /* Transform the argument into integer */
  xdb.d = x;

  /* Special case handling */

  /* Exact algebraic case x = 1, acos(1) = 0 */

  if (x == 1.0) return 0.0;

  /* Strip off the sign of argument x */
  if (xdb.i[HI] & 0x80000000) sign = -1; else sign = 1;
  xdb.i[HI] &= 0x7fffffff;

  /* acos is defined on -1 <= x <= 1, elsewhere it is NaN */
  if ((xdb.i[HI] > 0x3ff00000) || ((xdb.i[HI] == 0x3ff00000) && (xdb.i[LO] != 0x00000000))) {
    return (x-x)/0.0;    /* return NaN */
  }

  /* If |x| < 2^(-120) we have

     round(acos(x)) = round(pi/2)

     So we can decide the rounding without any computation
  */
  if (xdb.i[HI] < 0x38700000) {
    return PIHALFDOUBLERD;
  }

  /* Recast x */
  x = xdb.d;

  /* Find correspondant interval and compute index to the table
     We start by filtering the two special cases around 0 and 1
  */

  if (xdb.i[HI] < BOUND1) {
    /* Special interval 0..BOUND1
       The polynomial has no even monomials
       We must prove extra accuracy in the interval 0..sin(2^(-18))
    */

    /* Quick phase starts */

    /* Compute square of x for both quick and accurate phases */
    Mul12(&xSqh,&xSql,x,x);

    tmp4 = tbl[3];
    tmp5 = tbl[4];
    t4h = tmp4;
    t4l = tmp5;
    if (xdb.i[HI] > EXTRABOUND) {
      /* Double precision evaluation */
#if defined(PROCESSOR_HAS_FMA) && !defined(AVOID_FMA)
      highPoly = FMA(FMA(FMA(FMA(tbl[23],xSqh,tbl[21]),xSqh,tbl[19]),xSqh,tbl[17]),xSqh,tbl[15]);
#else
      highPoly = tbl[15] + xSqh * (tbl[17] + xSqh * (tbl[19] + xSqh * (tbl[21] + xSqh * tbl[23])));
#endif

      /* Double-double precision evaluation */
      Mul12(&tt1h,&tt1l,xSqh,highPoly);
      Add22(&t1h,&t1l,tbl[12],tbl[13],tt1h,tt1l);

      MulAdd212(&t2h,&t2l,tbl[9],tbl[10],xSqh,t1h,t1l);
      MulAdd212(&t3h,&t3l,tbl[6],tbl[7],xSqh,t2h,t2l);
      MulAdd22(&t4h,&t4l,tmp4,tmp5,xSqh,xSql,t3h,t3l);
    }

    MulAdd22(&t5h,&t5l,tbl[0],tbl[1],xSqh,xSql,t4h,t4l);

    Mul122(&xCubeh,&xCubel,x,xSqh,xSql);
    Mul22(&tt6h,&tt6l,xCubeh,xCubel,t5h,t5l);

    Add12(tmp1,tmp2,x,tt6h);
    tmp3 = tmp2 + tt6l;
    Add12(polyh,polyl,tmp1,tmp3);

    /* Reconstruction:

       - Multiply by the inverted sign
       - Add Pi/2 in double-double precision

    */

    zw1h = -sign * polyh;
    zw1l = -sign * polyl;

    Add22(&acosh,&acosm,PIHALFH,PIHALFM,zw1h,zw1l);

    /* Rounding test
       The RD rounding constant is at tbl[35]
    */
    TEST_AND_RETURN_RD(acosh, acosm, tbl[35]);

    /* Launch accurate phase */

    acos_accurate_lower(&acosh,&acosm,&acosl,x,xSqh,xSql,sign);

    ReturnRoundDownwards3(acosh,acosm,acosl);
  }

  if (xdb.i[HI] >= BOUND9) {
    /* Special interval BOUND9..1
       We use an asymptotic development of arcsin in sqrt(1 - x)
    */

    /* Argument reduction for quick and accurate phase
       z = 1 - x
       The operation is exact as per Sterbenz' lemma
    */

    z = 1 - x;

    /* Quick phase starts */

    /* Double precision evaluation */
#if defined(PROCESSOR_HAS_FMA) && !defined(AVOID_FMA)
    highPoly = FMA(FMA(FMA(FMA(FMA(FMA(FMA(FMA(FMA(
	       tbl[TBLIDX10+42] ,z,tbl[TBLIDX10+40]),z,tbl[TBLIDX10+38]),z,
               tbl[TBLIDX10+36]),z,tbl[TBLIDX10+34]),z,tbl[TBLIDX10+32]),z,
               tbl[TBLIDX10+30]),z,tbl[TBLIDX10+28]),z,tbl[TBLIDX10+26]),z,
               tbl[TBLIDX10+24]);
#else
    highPoly = tbl[TBLIDX10+24] + z * (tbl[TBLIDX10+26] + z * (tbl[TBLIDX10+28] + z * (
	       tbl[TBLIDX10+30] + z * (tbl[TBLIDX10+32] + z * (tbl[TBLIDX10+34] + z * (
	       tbl[TBLIDX10+36] + z * (tbl[TBLIDX10+38] + z * (tbl[TBLIDX10+40] + z *
               tbl[TBLIDX10+42]))))))));
#endif

    /* Double-double precision evaluation */
    Mul12(&tt1h,&tt1l,z,highPoly);
    Add22(&t1h,&t1l,tbl[TBLIDX10+21],tbl[TBLIDX10+22],tt1h,tt1l);

    MulAdd212(&t2h,&t2l,tbl[TBLIDX10+18],tbl[TBLIDX10+19],z,t1h,t1l);
    MulAdd212(&t3h,&t3l,tbl[TBLIDX10+15],tbl[TBLIDX10+16],z,t2h,t2l);
    MulAdd212(&t4h,&t4l,tbl[TBLIDX10+12],tbl[TBLIDX10+13],z,t3h,t3l);
    MulAdd212(&t5h,&t5l,tbl[TBLIDX10+9],tbl[TBLIDX10+10],z,t4h,t4l);
    MulAdd212(&t6h,&t6l,tbl[TBLIDX10+6],tbl[TBLIDX10+7],z,t5h,t5l);
    MulAdd212(&t7h,&t7l,tbl[TBLIDX10+3],tbl[TBLIDX10+4],z,t6h,t6l);
    MulAdd212(&t8h,&t8l,tbl[TBLIDX10+0],tbl[TBLIDX10+1],z,t7h,t7l);
    MulAdd212(&polyh,&polyl,-1,0,z,t8h,t8l);

    /* Compute sqrt(2*z) as a double-double */

    twoZ = 2 * z;
    sqrt12(&sqrtzh,&sqrtzl,twoZ);

    /* Multiply p(z) by sqrt(2*z) and add Pi/2 */

    Mul22(&pTimesSh,&pTimesSl,polyh,polyl,sqrtzh,sqrtzl);


    /* Reconstruction:

     If the sign of x in acos(x) was positive:
      - Multiply pTimesSh + pTimesSl approx f(x) * sqrt(2 * z) by -1
      - Return

     Otherwise:
      - Add Pi in triple-double to pTimesSh + pTimesSl approx f(x) * sqrt(2 * z)
      - Return

    */

    if (sign > 0) {

      acosh = -1.0 * pTimesSh;
      acosm = -1.0 * pTimesSl;

    } else {

      Add22(&acosh,&acosm,PIH,PIM,pTimesSh,pTimesSl);

    }

    /* Rounding test
       The RD rounding constant is at tbl[TBLIDX10+55]
    */
    TEST_AND_RETURN_RD(acosh, acosm, tbl[TBLIDX10+55]);

    /* Launch accurate phase */

    acos_accurate_higher(&acosh,&acosm,&acosl,z,sign);

    ReturnRoundDownwards3(acosh,acosm,acosl);
  }

  /* General 8 main intervals
     We can already suppose that BOUND1 <= x <= BOUND9
  */

  if (xdb.i[HI] < BOUND5) {
    if (xdb.i[HI] < BOUND3) {
      if (xdb.i[HI] < BOUND2) i = TBLIDX2; else i = TBLIDX3;
    } else {
      if (xdb.i[HI] < BOUND4) i = TBLIDX4; else i = TBLIDX5;
    }
  } else {
    if (xdb.i[HI] < BOUND7) {
      if (xdb.i[HI] < BOUND6) i = TBLIDX6; else i = TBLIDX7;
    } else {
      if (xdb.i[HI] < BOUND8) i = TBLIDX8; else i = TBLIDX9;
    }
  }

  /* Argument reduction
     i points to the interval midpoint value in the table
  */
  z = x - tbl[i];

  /* Quick phase starts */

  /* Double precision evaluation */

#if defined(PROCESSOR_HAS_FMA) && !defined(AVOID_FMA)
  highPoly = FMA(FMA(FMA(FMA(FMA(FMA(FMA(
	     tbl[i+35] ,z,tbl[i+33]),z,tbl[i+31]),z,tbl[i+29]),z,
             tbl[i+27]),z,tbl[i+25]),z,tbl[i+23]),z,tbl[i+21]);
#else
  highPoly = tbl[i+21] + z * (tbl[i+23] + z * (tbl[i+25] + z * (
             tbl[i+27] + z * (tbl[i+29] + z * (tbl[i+31] + z * (
             tbl[i+33] + z *  tbl[i+35]))))));
#endif

  /* Double-double precision evaluation */

  Mul12(&tt1h,&tt1l,z,highPoly);
  Add22(&t1h,&t1l,tbl[i+18],tbl[i+19],tt1h,tt1l);

  MulAdd212(&t2h,&t2l,tbl[i+15],tbl[i+16],z,t1h,t1l);
  MulAdd212(&t3h,&t3l,tbl[i+12],tbl[i+13],z,t2h,t2l);
  MulAdd212(&t4h,&t4l,tbl[i+9],tbl[i+10],z,t3h,t3l);
  MulAdd212(&t5h,&t5l,tbl[i+6],tbl[i+7],z,t4h,t4l);
  MulAdd212(&t6h,&t6l,tbl[i+3],tbl[i+4],z,t5h,t5l);
  MulAdd212(&polyh,&polyl,tbl[i+1],tbl[i+2],z,t6h,t6l);

  /* Reconstruction:

     - Multiply by the inverted sign
     - Add Pi/2 in double-double precision

  */

  zw1h = -sign * polyh;
  zw1l = -sign * polyl;

  Add22(&acosh,&acosm,PIHALFH,PIHALFM,zw1h,zw1l);

  /* Rounding test
     The RD rounding constant is at tbl[i+60]
  */
  TEST_AND_RETURN_RD(acosh, acosm, tbl[i+60]);

  /* Launch accurate phase */

  acos_accurate_middle(&acosh,&acosm,&acosl,z,i,sign);

  ReturnRoundDownwards3(acosh,acosm,acosl);
}

