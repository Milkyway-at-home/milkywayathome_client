/*
 *  Copyright (c) 2010-2011 Rensselaer Polytechnic Institute
 *  Copyright (c) 2010-2011 Ben Willett
 *  Copyright (c) 2010-2011 Matthew Arsenault
 *
 *  This file is part of Milkway@Home.
 *
 *  Milkway@Home is free software: you may copy, redistribute and/or modify it
 *  under the terms of the GNU General Public License as published by the
 *  Free Software Foundation, either version 3 of the License, or (at your
 *  option) any later version.
 *
 *  This file is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "nbody_config.h"

#include "nbody_priv.h"
#include "nbody_chisq.h"
#include "milkyway_util.h"
#include "nbody_mass.h"
#include "nbody_defaults.h"
#include "nbody_histogram.h"

static real nbSahaTerm(real* m, real* s)
{
    /*
      W = product_{i = 1 .. B }

                  (m_i + s_i)!
                 --------------
                   m_i! s_i!

       Using Stirling's approximation with an additional term
       ln(n!) ~= ln(2 *pi * n) / 2 + n * ln(n) - n
       ln(n!) ~= ln(2 pi) / 2 + (n + 1/2) * ln(n) - n

       This reduces to

       ln(W) = sum_{i = 1 .. B}


     */

    /* Stirling's approximation fails miserably at n = 0

       If s_i == 0 or m_i == 0
       the sum term reduces to 0.

       Suppose m_i == 0
       ln(0!) == 0

       ln(0!s_i!) - ln(0!) - ln(s_i) == ln(s_i!) - ln(s_i!) == 0

       If both == 0
       ln(0!0!) - ln(0!) - ln(0!) == 0
     */

    if (mw_fabs_0(showRealValue(m)) < 1.0e-10 || mw_fabs_0(showRealValue(s)) < 1.0e-10)
    {
        return ZERO_REAL;
    }
    else
    {
        real part1 = mw_real_const(-0.5 * mw_log_0(M_2PI));

        real add1 = mw_add(m, s);
        real add2 = mw_add_s(&add1, 0.5);
        real log1 = mw_log(&add1);
        real part2 = mw_mul(&add2, &log1);

        real add3 = mw_add_s(m, 0.5);
        real log2 = mw_log(m);
        real part3 = mw_mul(&add3, &log2);

        real add4 = mw_add_s(s, 0.5);
        real log3 = mw_log(s);
        real part4 = mw_mul(&add4, &log3);

        real part12 = mw_add(&part1, &part2);
        real part34 = mw_add(&part3, &part4);
        return mw_sub(&part12, &part34);
    }
}

static real nbPoissonTerm(real* f, real* y)
{
    /*
      Fitting a data set y(y1, y2, .. yn)  to a model function f(f1, f2, .. fn)
      sum = i .. N bins,
      2 sum(f_i - y_i) - sum(i != 1, y_i != 0) y_i * ln(f_i / y_i))
     */

    if (mw_fabs_0(showRealValue(f)) < 1.0e-10 || mw_fabs_0(showRealValue(y)) < 1.0e-10)
    {
        return mw_mul_s(f,2.0);
    }
    else
    {
        real fmy = mw_sub(f, y);
        real fdy = mw_div(f, y);
        real log1 = mw_log(&fdy);
        real mul1 = mw_mul(y, &log1);
        real sub1 = mw_sub(&fmy, &mul1);
        return mw_mul_s(&sub1 , 2.0);
    }
}

static real nbKullbackLeiblerTerm(real* h, real* k)
{
    /* Symmetrized version. (Jeffrey divergence?) */
    real m;

    if (mw_fabs_0(showRealValue(h)) < 1.0e-10 || mw_fabs_0(showRealValue(k)) < 1.0e-10)
    {
        return ZERO_REAL;
    }
    else
    {
        m = mw_add(h, k);
        m = mw_mul_s(&m, 0.5);
        real div1 = mw_div(h, &m);
        div1 = mw_log(&div1);
        div1 = mw_mul(h, &div1);
        real div2 = mw_div(k, &m);
        div2 = mw_log(&div2);
        div2 = mw_mul(k, &div2);
        return mw_add(&div1, &div2);
    }


#if 0
    real p = h;
    real q = k;
    /* Non-symmetrized version */
    if (fabs(p) < 1.0e-10 || fabs(q) < 1.0e-10)
    {
        /* Not sure this is really correct for q == 0 */
        return 0.0;
    }
    else
    {
        return p * (nbLog2(p) / nbLog2(q));
    }
#endif
}

static real nbChisqAlt(real* p, real* q)
{
    real ppq = mw_add(p, q);
    real pmq = mw_sub(p, q);
    real tmp = sqr(&pmq);
    tmp = mw_div(&tmp, &ppq);
    return mw_mul_s(&tmp,0.5);
}

/* Calculate chisq from read data histogram and the generated histogram */
real nbCalcChisq(const NBodyHistogram* data,        /* Data histogram */
                   const NBodyHistogram* histogram,   /* Generated histogram */
                   NBodyLikelihoodMethod method)
{
    unsigned int i;
    real tmp;
    real effTotalNum;
    real chiSq = ZERO_REAL;
    real n;
    real err;
    real simErr;
    real scale = mw_real_const(1.0);
    unsigned int nBin = data->lambdaBins * histogram->lambdaBins;

    assert(data->lambdaBins == histogram->lambdaBins);
    assert(data->betaBins == histogram->betaBins);


    if (!histogram->hasRawCounts)
    {
        mw_printf("FIXME: other likelihoods need raw count on generated histogram\n");
        return mw_real_const(NAN);
    }

    if (method == NBODY_SAHA)
    {
        /* We need to have the total number to scale to the correct
         * numbers for Saha likelihood */
        scale = (real) data->totalNum;
        if (showRealValue(&data->totalNum) == 0 || showRealValue(&histogram->totalNum) == 0)
        {
            mw_printf("Histogram scales required for Saha likelihood but missing\n");
            return mw_real_const(NAN);
        }
    }

    if (showRealValue(&histogram->totalNum) == 0)
    {
        return mw_real_const(INFINITY);
    }

    effTotalNum = (real) nbCorrectTotalNumberInHistogram(histogram, data);

    for (i = 0; i < nBin; ++i)
    {
        if (data->data[i].useBin)  /* Skip bins with missing data */
        {
            n = (real) histogram->data[i].variable;
            err = data->data[i].err;

            switch (method)
            {
                case NBODY_ORIG_CHISQ:
                    tmp = mw_div(&n, &effTotalNum);
                    tmp = mw_sub(&(data->data[i].variable), &tmp);
                    tmp = mw_div(&tmp, &err);
                    tmp = sqr(&tmp);
                    chiSq = mw_add(&chiSq, &tmp);
                    break;

                case NBODY_ORIG_ALT:
                    /* We already have errors from the simulation, but
                     * we need to correct the errors in case there
                     * were any bins we are skipping for matching to
                     * the data */
                    simErr = nbNormalizedHistogramError(&(histogram->data[i].variable), &effTotalNum);

                    /* effective error = sqrt( (data error)^2 + (sim count error)^2 ) */
                    err = mw_hypot(&err, &simErr);
                    tmp = mw_div(&n, &effTotalNum);
                    tmp = mw_sub(&(data->data[i].variable), &tmp);
                    tmp = mw_div(&tmp, &err);
                    tmp = sqr(&tmp);
                    chiSq = mw_add(&chiSq, &tmp);
                    break;

                case NBODY_CHISQ_ALT:
                    tmp = mw_div(&n, &effTotalNum);
                    tmp = nbChisqAlt(&(data->data[i].variable), &tmp);
                    chiSq = mw_add(&chiSq, &tmp);
                    break;

                case NBODY_POISSON:
                    /* Poisson one */
                    tmp = mw_div(&n, &effTotalNum);
                    tmp = nbPoissonTerm(&(data->data[i].variable), &tmp);
                    chiSq = mw_add(&chiSq, &tmp);
                    break;

                case NBODY_KOLMOGOROV:
                    tmp = mw_div(&n, &effTotalNum);
                    tmp = mw_sub(&(data->data[i].variable), &tmp);
                    tmp = mw_fabs(&tmp);
                    chiSq = mw_fmax(&chiSq, &tmp);
                    break;

                case NBODY_KULLBACK_LEIBLER:
                    /* "Relative entropy" */
                    tmp = mw_div(&n, &effTotalNum);
                    tmp = nbKullbackLeiblerTerm(&(data->data[i].variable), &tmp);
                    chiSq = mw_add(&chiSq, &tmp);
                    break;

                case NBODY_SAHA:
                    /* This will actually find ln(W). W is an unreasonably large number. */
                    tmp = mw_mul(&scale, &(data->data[i].variable));
                    tmp = nbSahaTerm(&n, &tmp);
                    chiSq = mw_add(&chiSq, &tmp);
                    break;

                case NBODY_INVALID_METHOD:
                case NBODY_EMD:
                default:
                    mw_fail("Invalid likelihood method\n");
            }
        }
    }

    return chiSq;
}


