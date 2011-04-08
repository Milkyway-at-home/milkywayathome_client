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

#include <stdio.h>
#include <time.h>

#include "separation_types.h"
#include "likelihood.h"
#include "integrals.h"
#include "integrals_common.h"
#include "r_points.h"
#include "calculated_constants.h"
#include "milkyway_util.h"
#include "separation_utils.h"
#include "evaluation_state.h"

/* CHECKME: What is this? */
static real probability_log(real bg, real sum_exp_weights)
{
    return mw_cmpzero_eps(bg, SEPARATION_EPS) ? -238.0 : mw_log10(bg / sum_exp_weights);
}

static real stream_sum(const SeparationResults* results,
                       const unsigned int number_streams,
                       const real* st_prob,
                       Kahan* st_only_sum,
                       const real* exp_stream_weights,
                       const real sum_exp_weights,
                       real bg_only)
{
    unsigned int i;
    real st_only;
    real star_prob = bg_only;

    for (i = 0; i < number_streams; ++i)
    {
        st_only = st_prob[i] / results->streamIntegrals[i] * exp_stream_weights[i];
        star_prob += st_only;
        st_only = probability_log(st_only, sum_exp_weights);
        KAHAN_ADD(st_only_sum[i], st_only);
    }
    star_prob /= sum_exp_weights;

    return star_prob;
}

static const int calculateSeparation = 1;
static const int twoPanel = 1;

/* get stream & background weight constants */
static real get_stream_bg_weight_consts(StreamStats* ss, const Streams* streams)
{
    unsigned int i;
    real epsilon_b;
    real denom = 1.0;

    for (i = 0; i < streams->number_streams; i++)
        denom += mw_exp(streams->stream_weight[i].weight);

    for (i = 0; i < streams->number_streams; i++)
    {
        ss[i].epsilon_s = mw_exp(streams->stream_weight[i].weight) / denom;
        printf("epsilon_s[%d]: %lf\n", i, ss[i].epsilon_s);
    }

    epsilon_b = 1.0 / denom;
    printf("epsilon_b:    %lf\n", epsilon_b);
    return epsilon_b;
}

static void twoPanelSeparation(const AstronomyParameters* ap,
                               const SeparationResults* results,
                               StreamStats* ss,
                               const real* st_probs,
                               real bg_prob,
                               real epsilon_b)
{
    unsigned int i;
    real pbx, psgSum;

    pbx = epsilon_b * bg_prob / results->backgroundIntegral;

    for (i = 0; i < ap->number_streams; i++)
        ss[i].psg = ss[i].epsilon_s * st_probs[i] / results->streamIntegrals[i];

    psgSum = 0;
    for (i = 0; i < ap->number_streams; i++)
        psgSum += ss[i].psg;

    for (i = 0; i < ap->number_streams; i++)
        ss[i].sprob = ss[i].psg / (psgSum + pbx);

    for (i = 0; i < ap->number_streams; i++)
        ss[i].nstars += ss[i].sprob;
}

static void nonTwoPanelSeparation(StreamStats* ss, unsigned int number_streams)
{
    unsigned int i;

    for (i = 0; i < number_streams; i++)
    {
        ss[i].sprob = 1.0;
        ss[i].nstars += 1.0;
    }
}

static void separation(FILE* f,
                       const AstronomyParameters* ap,
                       const SeparationResults* results,
                       const mwmatrix cmatrix,
                       StreamStats* ss,
                       const real* st_probs,
                       real bg_prob,
                       real epsilon_b,
                       mwvector current_star_point)
{
    int s_ok;
    mwvector starxyz;
    mwvector starxyzTransform;

    mwvector xsun = ZERO_VECTOR;
    X(xsun) = ap->m_sun_r0;

    if (twoPanel)
        twoPanelSeparation(ap, results, ss, st_probs, bg_prob, epsilon_b);
    else
        nonTwoPanelSeparation(ss, ap->number_streams);

    /* determine if star with sprob should be put into stream */
    s_ok = prob_ok(ss, ap->number_streams);
    if (s_ok >= 1)
        ss[s_ok-1].q++;

    starxyz = lbr2xyz(ap, current_star_point);
    starxyzTransform = transform_point(ap, starxyz, cmatrix, xsun);

    if (f)
    {
        fprintf(f,
                "%d %lf %lf %lf\n",
                s_ok,
                X(starxyzTransform), Y(starxyzTransform), Z(starxyzTransform));
    }
}

static real likelihood_probability(const AstronomyParameters* ap,
                                   const StreamConstants* sc,
                                   const Streams* streams,
                                   const RPoints* r_pts,
                                   const real* sg_dx,
                                   const LBTrig lbt,
                                   const RConsts rc,
                                   const real reff_xr_rp3,
                                   const SeparationResults* results,
                                   EvaluationState* es,

                                   const int do_separation,
                                   FILE* f,
                                   StreamStats* ss,
                                   const real epsilon_b,
                                   mwvector point,
                                   const mwmatrix cmatrix)
{
    real starProb;

    /* if q is 0, there is no probability */
    if (ap->q == 0.0)
    {
        es->bgTmp = -1.0;
    }
    else
    {
        bg_probability(ap, sc, r_pts, sg_dx, rc.gPrime, lbt, es);
        //multSumProbs(es, reff_xr_rp3);

#if 1
        es->bgTmp *= reff_xr_rp3;  /* bg only */
        for (unsigned int i = 0; i < ap->number_streams; ++i)
            es->streamTmps[i] *= reff_xr_rp3;
#endif

    }

    es->bgTmp = (es->bgTmp / results->backgroundIntegral) * ap->exp_background_weight;

    starProb = stream_sum(results,
                          streams->number_streams,
                          es->streamTmps,
                          es->streamSums,
                          streams->expStreamWeights,
                          streams->sumExpWeights,
                          es->bgTmp);

    es->bgTmp = probability_log(es->bgTmp, streams->sumExpWeights);
    KAHAN_ADD(es->bgSum, es->bgTmp);


    if (do_separation)
        separation(f, ap, results, cmatrix, ss, es->streamTmps, es->bgTmp, epsilon_b, point);


    return starProb;
}

static real calculateLikelihood(const Kahan* ksum, unsigned int nStars, unsigned int badJacobians)
{
    real sum;

    /*  log10(x * 0.001) = log10(x) - 3.0 */
    sum = ksum->sum + ksum->correction;
    sum /= (nStars - badJacobians);
    sum -= 3.0;

    return sum;
}

static void calculateLikelihoods(SeparationResults* results,
                                 Kahan* prob,
                                 Kahan* bgOnly,
                                 Kahan* streamOnly,
                                 unsigned int nStars,
                                 unsigned int nStreams,
                                 unsigned badJacobians)
{
    unsigned int i;

    /* CHECKME: badJacobians supposed to only be for final? */
    results->backgroundLikelihood = calculateLikelihood(bgOnly, nStars, 0);
    for (i = 0; i < nStreams; i++)
    {
        results->streamLikelihoods[i] = calculateLikelihood(&streamOnly[i], nStars, 0);
    }

    results->likelihood = calculateLikelihood(prob, nStars, badJacobians);
}

/* separation init stuffs */
static void setSeparationConstants(const AstronomyParameters* ap,
                                   const SeparationResults* results,
                                   mwmatrix cmatrix)
{
    unsigned int i;
    mwvector dnormal = ZERO_VECTOR;
    const mwvector dortho = mw_vec(0.0, 0.0, 1.0);

    if (!ap->sgr_coordinates)
        dnormal = stripe_normal(ap->wedge);
    else
        fail("sgr coordinates not implemented\n");

    get_transform(cmatrix, dnormal, dortho);

    printf("\nTransformation matrix:\n"
           "\t%lf %lf %lf\n"
           "\t%lf %lf %lf\n"
           "\t%lf %lf %lf\n",
           X(cmatrix[0]), Y(cmatrix[0]), Z(cmatrix[0]),
           X(cmatrix[1]), Y(cmatrix[1]), Z(cmatrix[1]),
           X(cmatrix[2]), Y(cmatrix[2]), Z(cmatrix[2]));

    printf("==============================================\n");
    printf("bint: %lf\n", results->backgroundIntegral);
    for (i = 0; i < ap->number_streams; i++)
        printf("sint[%d]: %lf\n", i, results->streamIntegrals[i]);
}

static void printSeparationStats(const StreamStats* ss,
                                 const unsigned int number_stars,
                                 const unsigned int number_streams)
{
    unsigned int i;
    real percent;

    printf("%d total stars\n", number_stars);
    for (i = 0; i < number_streams; ++i)
    {
        percent = 100.0 * (ss[i].nstars / (real) number_stars);
        printf("%lf in stream[%d] (%lf%%)\n", ss[i].nstars, i, percent);
    }

    for (i = 0; i < number_streams; ++i)
        printf("%d stars separated into stream\n", ss[i].q);
}

static int likelihood_sum(SeparationResults* results,
                          const AstronomyParameters* ap,
                          const StarPoints* sp,
                          const StreamConstants* sc,
                          const Streams* streams,
                          const StreamGauss sg,
                          RPoints* r_pts,
                          EvaluationState* es,

                          const int do_separation,
                          StreamStats* ss,
                          FILE* f)
{
    Kahan prob = ZERO_KAHAN;

    unsigned int current_star_point;
    mwvector point;
    real star_prob;
    LB lb;
    LBTrig lbt;
    real reff_xr_rp3;
    RConsts rc = { 0.0, 0.0 };

    real epsilon_b = 0.0;
    mwmatrix cmatrix;
    unsigned int num_zero = 0;
    unsigned int badJacobians = 0;  /* CHECKME: Seems like this never changes */

    if (do_separation)
    {
        setSeparationConstants(ap, results, cmatrix);
        epsilon_b = get_stream_bg_weight_consts(ss, streams);
    }

    for (current_star_point = 0; current_star_point < sp->number_stars; ++current_star_point)
    {
        point = sp->stars[current_star_point];
        rc.gPrime = calcG(Z(point));
        setRPoints(ap, sg, ap->convolve, rc.gPrime, r_pts);
        reff_xr_rp3 = calcReffXrRp3(Z(point), rc.gPrime);

        LB_L(lb) = L(point);
        LB_B(lb) = B(point);

        lbt = lb_trig(lb);

        star_prob = likelihood_probability(ap, sc, streams, r_pts, sg.dx, lbt, rc,
                                           reff_xr_rp3, results, es,
                                           do_separation, f, ss, epsilon_b, point, cmatrix);
        if (mw_cmpnzero_muleps(star_prob, SEPARATION_EPS))
        {
            star_prob = mw_log10(star_prob);
            KAHAN_ADD(prob, star_prob);
        }
        else
        {
            ++num_zero;
            prob.sum -= 238.0;
        }

    }

    calculateLikelihoods(results, &prob, &es->bgSum, es->streamSums,
                         sp->number_stars, streams->number_streams, badJacobians);


    if (do_separation)
        printSeparationStats(ss, sp->number_stars, ap->number_streams);

    return 0;
}

StreamStats* newStreamStats(const unsigned int number_streams)
{
    return (StreamStats*) mwCallocA(number_streams, sizeof(StreamStats));
}

int likelihood(SeparationResults* results,
               const AstronomyParameters* ap,
               const StarPoints* sp,
               const StreamConstants* sc,
               const Streams* streams,
               const StreamGauss sg,
               const int do_separation,
               const char* separation_outfile)
{
    RPoints* r_pts;
    EvaluationState* es;
    StreamStats* ss = NULL;
    FILE* f = NULL;

    int rc = 0;
    double t1, t2;


    if (do_separation)
    {
        f = mw_fopen(separation_outfile, "w");
        if (!f)
        {
            perror("Opening separation output file");
            return 1;
        }

        ss = newStreamStats(streams->number_streams);
    }

    /* New state for this sum */
    es = newEvaluationState(ap);

    r_pts = (RPoints*) mwMallocA(sizeof(RPoints) * ap->convolve);

    t1 = mwGetTime();
    rc = likelihood_sum(results,
                        ap, sp, sc, streams,
                        sg, r_pts,
                        es,

                        do_separation,
                        ss,
                        f);
    t2 = mwGetTime();
    warn("Likelihood time = %f s\n", t2 - t1);

    mwFreeA(r_pts);
    mwFreeA(ss);
    freeEvaluationState(es);

    if (f && fclose(f))
        perror("Closing separation output file");

    return rc;
}

