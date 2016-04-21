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

static real nbSahaTerm(real m, real s)
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

    if (mw_fabs(m) < 1.0e-10 || mw_fabs(s) < 1.0e-10)
    {
        return 0.0;
    }
    else
    {
        return -0.5 * mw_log(M_2PI) + (m + s + 0.5) * mw_log(m + s) - (m + 0.5) * mw_log(m) - (s + 0.5) * mw_log(s);
    }
}

static real nbPoissonTerm(real f, real y)
{
    /*
      Fitting a data set y(y1, y2, .. yn)  to a model function f(f1, f2, .. fn)
      sum = i .. N bins,
      2 sum(f_i - y_i) - sum(i != 1, y_i != 0) y_i * ln(f_i / y_i))
     */

    if (mw_fabs(f) < 1.0e-10 || mw_fabs(y) < 1.0e-10)
    {
        return 2.0 * f;
    }
    else
    {
        return 2.0 * ((f - y) - y * mw_log(f / y));
    }
}

static real nbKullbackLeiblerTerm(real h, real k)
{
    /* Symmetrized version. (Jeffrey divergence?) */
    real m;

    if (mw_fabs(h) < 1.0e-10 || mw_fabs(k) < 1.0e-10)
    {
        return 0.0;
    }
    else
    {
        m = (h + k) / 2.0;
        return h * mw_log(h / m) + k * mw_log(k / m);
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

static real nbChisqAlt(real p, real q)
{
    return 0.5 * sqr(p - q) / (p + q);
}

/* Calculate chisq from read data histogram and the generated histogram */
real nbCalcChisq(const NBodyHistogram* data,        /* Data histogram */
                   const NBodyHistogram* histogram,   /* Generated histogram */
                   NBodyLikelihoodMethod method)
{
    unsigned int i;
    real tmp;
    real effTotalNum;
    real chiSq = 0.0;
    real n;
    real err;
    real simErr;
    real scale = 1.0;
    unsigned int nBin = data->lambdaBins * histogram->lambdaBins;

    assert(data->lambdaBins == histogram->lambdaBins);
    assert(data->betaBins == histogram->betaBins);


    if (!histogram->hasRawCounts)
    {
        mw_printf("FIXME: other likelihoods need raw count on generated histogram\n");
        return NAN;
    }

    if (method == NBODY_SAHA)
    {
        /* We need to have the total number to scale to the correct
         * numbers for Saha likelihood */
        scale = (real) data->totalNum;
        if (data->totalNum == 0 || histogram->totalNum == 0)
        {
            mw_printf("Histogram scales required for Saha likelihood but missing\n");
            return NAN;
        }
    }

    if (histogram->totalNum == 0)
    {
        return INFINITY;
    }

    effTotalNum = (real) nbCorrectTotalNumberInHistogram(histogram, data);

    for (i = 0; i < nBin; ++i)
    {
        if (data->data[i].useBin)  /* Skip bins with missing data */
        {
            n = (real) histogram->data[i].rawCount;
            err = data->data[i].err;

            switch (method)
            {
                case NBODY_ORIG_CHISQ:
                    tmp = (data->data[i].count - (n / effTotalNum)) / err;
                    chiSq += sqr(tmp);
                    break;

                case NBODY_ORIG_ALT:
                    /* We already have errors from the simulation, but
                     * we need to correct the errors in case there
                     * were any bins we are skipping for matching to
                     * the data */
                    simErr = nbNormalizedHistogramError(histogram->data[i].rawCount, effTotalNum);

                    /* effective error = sqrt( (data error)^2 + (sim count error)^2 ) */
                    err = sqrt(sqr(err) + sqr(simErr));
                    tmp = (data->data[i].count - (n / effTotalNum)) / err;
                    chiSq += sqr(tmp);
                    break;

                case NBODY_CHISQ_ALT:
                    chiSq += nbChisqAlt(data->data[i].count, n / effTotalNum);
                    break;

                case NBODY_POISSON:
                    /* Poisson one */
                    chiSq += nbPoissonTerm(data->data[i].count, n / effTotalNum);
                    break;

                case NBODY_KOLMOGOROV:
                    chiSq = mw_fmax(chiSq, fabs(data->data[i].count - (n / effTotalNum)));
                    break;

                case NBODY_KULLBACK_LEIBLER:
                    /* "Relative entropy" */
                    chiSq += nbKullbackLeiblerTerm(data->data[i].count, n / effTotalNum);
                    break;

                case NBODY_SAHA:
                    /* This will actually find ln(W). W is an unreasonably large number. */
                    chiSq += nbSahaTerm(n, scale * data->data[i].count);
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


