/*
Copyright 2008-2010 Travis Desell, Dave Przybylo, Nathan Cole, Matthew
Arsenault, Boleslaw Szymanski, Heidi Newberg, Carlos Varela, Malik
Magdon-Ismail and Rensselaer Polytechnic Institute.

This file is part of Milkway@Home.

Milkyway@Home is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Milkyway@Home is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Milkyway@Home.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "evaluation_state.h"
#include "integrals.h"
#include "coordinates.h"
#include "r_points.h"
#include "milkyway_util.h"
#include "calculated_constants.h"

#include <time.h>

#ifdef MILKYWAY_IPHONE_APP
double _milkywaySeparationGlobalProgress = 0.0;
#endif


static inline real progress(const EvaluationState* es,
                            const IntegralArea* ia,
                            real total_calc_probs)
{
    /* This integral's progress */
    /* When checkpointing is done, ia->mu_step would always be 0 */
    unsigned int i_prog =  (es->nu_step * ia->mu_steps * ia->r_steps)
                         + (es->mu_step * ia->r_steps); /* + es->r_step */

    return (real)(i_prog + es->current_calc_probs) / total_calc_probs;
}


#if BOINC_APPLICATION

static inline void doBoincCheckpoint(const EvaluationState* es,
                                     const IntegralArea* ia,
                                     real total_calc_probs)
{
    if (boinc_time_to_checkpoint())
    {
        if (writeCheckpoint(es))
            fail("Write checkpoint failed\n");
        boinc_checkpoint_completed();
    }

    boinc_fraction_done(progress(es, ia, total_calc_probs));
}

#elif MILKYWAY_IPHONE_APP

static inline void doBoincCheckpoint(const EvaluationState* es,
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
            fail("Write checkpoint failed\n");
    }

    _milkywaySeparationGlobalProgress = progress(es, ia, total_calc_probs);
}

#else /* Plain */

#define doBoincCheckpoint(es, ia, total_calc_probs)

#endif /* BOINC_APPLICATION */

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

static inline mwvector lbr2xyz_2(const AstronomyParameters* ap,
                                 const real rPoint,
                                 const LBTrig lbt)
{
    mwvector xyz;

    // This mad for some reason increases GPR usage by 1 pushing into next level of unhappy
    xyz.x = mw_mad(rPoint, LCOS_BCOS(lbt), ap->m_sun_r0);
    xyz.y = rPoint * LSIN_BCOS(lbt);
    xyz.z = rPoint * BSIN(lbt);

    return xyz;
}

static inline real streamIncrement(const StreamConstants* sc, mwvector xyz)
{
    real xyz_norm, dotted;
    mwvector xyzs;

    xyzs = mw_subv(xyz, sc->c);
    dotted = mw_dotv(sc->a, xyzs);
    mw_incsubv_s(xyzs, sc->a, dotted);

    xyz_norm = mw_sqrv(xyzs);

    return mw_exp(-xyz_norm * sc->sigma_sq2_inv);
}

HOT
static inline void streamSums(real* st_probs,
                              const StreamConstants* sc,
                              const mwvector xyz,
                              const real qw_r3_N,
                              const unsigned int nstreams)
{
    unsigned int i;

    for (i = 0; i < nstreams; ++i)
        st_probs[i] += qw_r3_N * streamIncrement(&sc[i], xyz);
}

HOT
static inline real h_prob_fast(const AstronomyParameters* ap, real qw_r3_N, real rg)
{
    const real rs = rg + ap->r0;
    return qw_r3_N / (rg * cube(rs));
}

HOT
static inline real h_prob_slow(const AstronomyParameters* ap, real qw_r3_N, real rg)
{
    const real rs = rg + ap->r0;
    return qw_r3_N / (mw_powr(rg, ap->alpha) * mw_powr(rs, ap->alpha_delta3));
}

HOT
static inline real rg_calc(const AstronomyParameters* ap, const mwvector xyz)
{
    /* sqrt(x^2 + y^2 + q_inv_sqr * z^2) */

    real tmp;

    tmp = sqr(X(xyz));
    tmp = mw_mad(Y(xyz), Y(xyz), tmp);               /* x^2 + y^2 */
    tmp = mw_mad(ap->q_inv_sqr, sqr(Z(xyz)), tmp);   /* (q_invsqr * z^2) + (x^2 + y^2) */

    return mw_sqrt(tmp);
}

ALWAYS_INLINE OLD_GCC_EXTERNINLINE
inline void zero_st_probs(real* st_probs, const unsigned int nstream)
{
    unsigned int i;

    for (i = 0; i < nstream; ++i)
        st_probs[i] = 0.0;
}

ALWAYS_INLINE OLD_GCC_EXTERNINLINE
inline void sum_probs(Kahan* probs,
                      const real* st_probs,
                      const real V_reff_xr_rp3,
                      const unsigned int nstream)
{
    unsigned int i;

    for (i = 0; i < nstream; ++i)
        KAHAN_ADD(probs[i], V_reff_xr_rp3 * st_probs[i]);
}

static inline real aux_prob(const AstronomyParameters* ap,
                            const real qw_r3_N,
                            const real r_in_mag)
{
    return qw_r3_N * (ap->bg_a * sqr(r_in_mag) + ap->bg_b * r_in_mag + ap->bg_c);
}

HOT
static inline void multProbs(EvaluationState* es, real V_reff_xr_rp3)
{
    unsigned int i;

    es->bgTmp *= V_reff_xr_rp3;
    for (i = 0; i < es->numberStreams; ++i)
        es->streamTmps[i] *= V_reff_xr_rp3;
}

HOT
static inline real bg_probability_fast_hprob(const AstronomyParameters* ap,
                                             const StreamConstants* sc,
                                             const RPoints* r_pts,
                                             const real* sg_dx,
                                             LBTrig lbt,
                                             real gPrime,
                                             real* streamTmps)
{
    unsigned int i;
    real h_prob, g, rg;
    mwvector xyz;
    real bg_prob = 0.0;
    unsigned int convolve = ap->convolve;
    int aux_bg_profile = ap->aux_bg_profile;

    zero_st_probs(streamTmps, ap->number_streams);
    for (i = 0; i < convolve; ++i)
    {
        xyz = lbr2xyz_2(ap, r_pts[i].r_point, lbt);
        rg = rg_calc(ap, xyz);

        h_prob = h_prob_fast(ap, r_pts[i].qw_r3_N, rg);

        /* Add a quadratic term in g to the the Hernquist profile */
        if (aux_bg_profile)
        {
            g = gPrime + sg_dx[i];
            h_prob += aux_prob(ap, r_pts[i].qw_r3_N, g);
        }

        bg_prob += h_prob;
        streamSums(streamTmps, sc, xyz, r_pts[i].qw_r3_N, ap->number_streams);
    }

    return bg_prob;
}

HOT
static inline real bg_probability_slow_hprob(const AstronomyParameters* ap,
                                             const StreamConstants* sc,
                                             const RPoints* r_pts,
                                             const real* sg_dx,
                                             LBTrig lbt,
                                             real gPrime,
                                             real* streamTmps)
{
    unsigned int i;
    real rg, g;
    mwvector xyz;
    real bg_prob = 0.0;
    unsigned int convolve = ap->convolve;
    int aux_bg_profile = ap->aux_bg_profile;

    zero_st_probs(streamTmps, ap->number_streams);

    for (i = 0; i < convolve; ++i)
    {
        xyz = lbr2xyz_2(ap, r_pts[i].r_point, lbt);

        rg = rg_calc(ap, xyz);

        bg_prob += h_prob_slow(ap, r_pts[i].qw_r3_N, rg);
        if (aux_bg_profile)
        {
            g = gPrime + sg_dx[i];
            bg_prob += aux_prob(ap, r_pts[i].qw_r3_N, g);
        }

        streamSums(streamTmps, sc, xyz, r_pts[i].qw_r3_N, ap->number_streams);
    }

    return bg_prob;
}

HOT
static inline void sumProbs(EvaluationState* es)
{
    unsigned int i;

    KAHAN_ADD(es->bgSum, es->bgTmp);
    for (i = 0; i < es->numberStreams; ++i)
        KAHAN_ADD(es->streamSums[i], es->streamTmps[i]);
}


HOT ALWAYS_INLINE
inline void bg_probability(const AstronomyParameters* ap,
                           const StreamConstants* sc,
                           const RPoints* r_pts,
                           const real* sg_dx,
                           real gPrime,
                           real reff_xr_rp3,
                           LBTrig lbt, /* integral point */
                           EvaluationState* es)
{
    if (ap->fast_h_prob)
        es->bgTmp = bg_probability_fast_hprob(ap, sc, r_pts, sg_dx, lbt, gPrime, es->streamTmps);
    else
        es->bgTmp = bg_probability_slow_hprob(ap, sc, r_pts, sg_dx, lbt, gPrime, es->streamTmps);

    multProbs(es, reff_xr_rp3);
}


HOT
static inline void r_sum(const AstronomyParameters* ap,
                         const StreamConstants* sc,
                         const real* sg_dx,
                         LBTrig lbt,
                         real id,
                         EvaluationState* es,
                         const RPoints* r_pts,
                         const RConsts* rc,
                         unsigned int r_steps)
{
    unsigned int r_step;
    real reff_xr_rp3;

    for (r_step = 0; r_step < r_steps; ++r_step)
    {
        reff_xr_rp3 = id * rc[r_step].irv_reff_xr_rp3;
        bg_probability(ap, sc,
                       &r_pts[r_step * ap->convolve],
                       sg_dx,
                       rc[r_step].gPrime, reff_xr_rp3, lbt, es);
        sumProbs(es);
    }
}

HOT
static inline void mu_sum(const AstronomyParameters* ap,
                          const IntegralArea* ia,
                          const StreamConstants* sc,
                          const RConsts* rc,
                          const RPoints* r_pts,
                          const real* sg_dx,
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
        doBoincCheckpoint(es, ia, ap->total_calc_probs);

        mu = mu_min + (((real) es->mu_step + 0.5) * mu_step_size);

        lb = gc2lb(ap->wedge, mu, nuid.nu); /* integral point */
        lbt = lb_trig(lb);

        r_sum(ap, sc, sg_dx, lbt, nuid.id, es, r_pts, rc, ia->r_steps);
    }

    es->mu_step = 0;
}

static void nuSum(const AstronomyParameters* ap,
                  const IntegralArea* ia,
                  const StreamConstants* sc,
                  const RConsts* rc,
                  const RPoints* r_pts,
                  const real* sg_dx,
                  EvaluationState* es)
{
    NuId nuid;

    for ( ; es->nu_step < ia->nu_steps; es->nu_step++)
    {
        nuid = calcNuStep(ia, es->nu_step);

        mu_sum(ap, ia, sc, rc, r_pts, sg_dx, nuid, es);
    }

    es->nu_step = 0;
}

void separationIntegralApplyCorrection(EvaluationState* es)
{
    unsigned int i;

    es->cut->bgIntegral = es->bgSum.sum + es->bgSum.correction;
    for (i  = 0; i < es->numberStreams; ++i)
        es->cut->streamIntegrals[i] = es->streamSums[i].sum + es->streamSums[i].correction;
}


/* returns background integral */
int integrate(const AstronomyParameters* ap,
              const IntegralArea* ia,
              const StreamConstants* sc,
              const StreamGauss sg,
              EvaluationState* es)
{
    RPoints* r_pts;
    RConsts* rc;

    if (ap->q == 0.0)
    {
        /* if q is 0, there is no probability */
        /* Short circuit the entire integral rather than add up -1 many times. */
        warn("q is 0.0\n");
        es->cut->bgIntegral = -1.0 * ia->nu_steps * ia->mu_steps * ia->r_steps;
        return 1;
    }

    r_pts = precalculateRPts(ap, ia, sg, &rc, 0);

    nuSum(ap, ia, sc, rc, r_pts, sg.dx, es);
    separationIntegralApplyCorrection(es);

    mwFreeA(r_pts);
    mwFreeA(rc);

  #ifdef MILKYWAY_IPHONE_APP
    _milkywaySeparationGlobalProgress = 1.0;
  #endif

    return 0;
}


