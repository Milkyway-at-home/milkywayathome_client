/*
 *  Copyright (c) 2008-2010 Travis Desell, Nathan Cole, Dave Przybylo
 *  Copyright (c) 2008-2010 Boleslaw Szymanski, Heidi Newberg
 *  Copyright (c) 2008-2010 Carlos Varela, Malik Magdon-Ismail
 *  Copyright (c) 2008-2011 Rensselaer Polytechnic Institute
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

#include "evaluation_state.h"
#include "integrals.h"
#include "coordinates.h"
#include "r_points.h"
#include "milkyway_util.h"
#include "calculated_constants.h"
#include "evaluation.h"
#include "probabilities_dispatch.h"

#include <time.h>

/* Marshaling into split r_points and qw_r3_N which helps with vectorization */
static RConsts* initRPoints(const AstronomyParameters* ap,
                            const IntegralArea* ia,
                            const StreamGauss sg,
                            real* RESTRICT rPoints,
                            real* RESTRICT qw_r3_N)
{
    unsigned int i, idx;
    int j;
    RPoints* rPts;
    RConsts* rc;

    rPts = precalculateRPts(ap, ia, sg, &rc, FALSE);

    for (i = 0; i < ia->r_steps; ++i)
    {
        for (j = 0; j < ap->convolve; ++j)
        {
            idx = i * ap->convolve + j;
            rPoints[idx] = rPts[idx].r_point;
            qw_r3_N[idx] = rPts[idx].qw_r3_N;
        }
    }

    mwFreeA(rPts);
    return rc;
}


#ifdef MILKYWAY_IPHONE_APP
double _milkywaySeparationGlobalProgress = 0.0;
#endif


static real progress(const AstronomyParameters* ap, const EvaluationState* es, const IntegralArea* ia, real totalCalcProbs)
{
    /* This integral's progress */
    /* When checkpointing is done, ia->mu_step would always be 0 */
    unsigned int i_prog =  ((uint64_t) es->nu_step * ia->mu_steps * ia->r_steps)
                        + ((uint64_t) es->mu_step * ia->r_steps); /* + es->r_step */

    return ((real)(i_prog + es->current_calc_probs) + (es->currentWU * totalCalcProbs)) / (ap->totalWUs * totalCalcProbs);
}


#if MILKYWAY_IPHONE_APP

static inline void doBoincCheckpoint(const AstronomyParameters* ap,
                                     const EvaluationState* es,
                                     const IntegralArea* ia,
                                     real total_calc_probs)
{
    static time_t lastCheckpoint = 0;
    static const time_t checkpointPeriod = 60;
    time_t now;

    if ((now = time(NULL)) - lastCheckpoint > checkpointPeriod)
    {
        lastCheckpoint = now;
        if (writeCheckpoint(es))
            mw_fail("Write checkpoint failed\n");
    }

    _milkywaySeparationGlobalProgress = progress(ap, es, ia, total_calc_probs);
}

#else /* Plain */

static inline void doBoincCheckpoint(const AstronomyParameters* ap,
                                     EvaluationState* es,
                                     const IntegralArea* ia,
                                     real total_calc_probs)
{
    if (mw_time_to_checkpoint())
    {
        if (writeCheckpoint(es))
        {
            mw_fail("Write checkpoint failed\n");
        }

        mw_checkpoint_completed();
    }

    mw_fraction_done(progress(ap, es, ia, total_calc_probs));
}

#endif /* BOINC_APPLICATION */

HOT
static inline void sumProbs(EvaluationState* es)
{
    int i;

    KAHAN_ADD(es->bgSum, es->bgTmp);
    for (i = 0; i < es->numberStreams; ++i)
        KAHAN_ADD(es->streamSums[i], es->streamTmps[i]);
}


HOT
static inline void r_sum(const AstronomyParameters* ap,
                         const StreamConstants* sc,
                         const real* RESTRICT sg_dx,
                         const real* RESTRICT rPoints,
                         const real* RESTRICT qw_r3_N,
                         LBTrig lbt,
                         real id,
                         EvaluationState* es,
                         const RConsts* rc,
                         unsigned int r_steps)
{
    unsigned int r_step;
    real reff_xr_rp3;

    for (r_step = 0; r_step < r_steps; ++r_step)
    {
        reff_xr_rp3 = id * rc[r_step].irv_reff_xr_rp3;
        es->bgTmp = probabilityFunc(ap,
                                    sc,
                                    sg_dx,
                                    &rPoints[r_step * ap->convolve],
                                    &qw_r3_N[r_step * ap->convolve],
                                    lbt,
                                    rc[r_step].gPrime,
                                    reff_xr_rp3,
                                    es->streamTmps);
        sumProbs(es);
    }
}

HOT
inline LBTrig lb_trig(LB lb)
{
    LBTrig lbt;
    real bCos;

    mw_sincos(d2r(LB_L(lb)), &lbt.lSinBCos, &lbt.lCosBCos);
    mw_sincos(d2r(LB_B(lb)), &lbt.bSin, &bCos);

    lbt.lCosBCos *= bCos;
    lbt.lSinBCos *= bCos;

    return lbt;
}

HOT
static inline void mu_sum(const AstronomyParameters* ap,
                          const IntegralArea* ia,
                          const StreamConstants* sc,
                          const RConsts* rc,
                          const real* RESTRICT sg_dx,
                          const real* RESTRICT rPoints,
                          const real* RESTRICT qw_r3_N,
                          const NuId nuid,
                          EvaluationState* es)
{
    real mu;
    LB lb;
    LBTrig lbt;

    const real mu_step_size = ia->mu_step_size;
    const real mu_min = ia->mu_min;

    for (; es->mu_step < ia->mu_steps; es->mu_step++)
    {
        doBoincCheckpoint(ap, es, ia, ap->total_calc_probs);

        mu = mu_min + (((real) es->mu_step + 0.5) * mu_step_size);

        lb = gc2lb(ap->wedge, mu, nuid.nu); /* integral point */
        lbt = lb_trig(lb);

        r_sum(ap, sc, sg_dx, rPoints, qw_r3_N, lbt, nuid.id, es, rc, ia->r_steps);
    }

    es->mu_step = 0;
}

static void nuSum(const AstronomyParameters* ap,
                  const IntegralArea* ia,
                  const StreamConstants* sc,
                  const RConsts* rc,
                  const real* RESTRICT sg_dx,
                  const real* RESTRICT rPoints,
                  const real* RESTRICT qw_r3_N,
                  EvaluationState* es)
{
    NuId nuid;

    for ( ; es->nu_step < ia->nu_steps; es->nu_step++)
    {
        nuid = calcNuStep(ia, es->nu_step);

        mu_sum(ap, ia, sc, rc, sg_dx, rPoints, qw_r3_N, nuid, es);
    }

    es->nu_step = 0;
}

void separationIntegralGetSums(EvaluationState* es)
{
    int i;

    /* Take the final sums */
    es->cut->bgIntegral = es->bgSum.sum;
    for (i  = 0; i < es->numberStreams; ++i)
        es->cut->streamIntegrals[i] = es->streamSums[i].sum;
}


/* returns background integral */
int integrate(const AstronomyParameters* ap,
              const IntegralArea* ia,
              const StreamConstants* sc,
              const StreamGauss sg,
              EvaluationState* es,
              const CLRequest* clr,
              const CLInfo* _ci)
{
    RConsts* rc;
    real* RESTRICT rPoints;
    real* RESTRICT qw_r3_N;

    (void) clr, (void) _ci;

    if (ap->q == 0.0)
    {
        /* if q is 0, there is no probability */
        /* Short circuit the entire integral rather than add up -1 many times. */
        mw_printf("q is 0.0\n");
        es->cut->bgIntegral = -1.0 * ia->nu_steps * ia->mu_steps * ia->r_steps;
        return 1;
    }

    rPoints = mwMallocA(sizeof(real) * ia->r_steps * ap->convolve);
    qw_r3_N = mwMallocA(sizeof(real) * ia->r_steps * ap->convolve);
    rc = initRPoints(ap, ia, sg, rPoints, qw_r3_N);

    nuSum(ap, ia, sc, rc, sg.dx, rPoints, qw_r3_N, es);
    separationIntegralGetSums(es);

    mwFreeA(rc);
    mwFreeA(rPoints);
    mwFreeA(qw_r3_N);

  #ifdef MILKYWAY_IPHONE_APP
    _milkywaySeparationGlobalProgress = 1.0;
  #endif

    return 0;
}

