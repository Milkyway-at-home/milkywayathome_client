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
#include "integrals_common.h"
#include "r_points.h"
#include "calculated_constants.h"
#include "milkyway_util.h"
#include "separation_utils.h"

/* FIXME: Excessive duplication with stuff used in integrals which I
 * was too lazy to also fix here */

static inline real likelihood_bg_probability_main(const ASTRONOMY_PARAMETERS* ap,
                                                  const STREAM_CONSTANTS* sc,
                                                  const R_POINTS* r_pts,
                                                  const real* sg_dx,
                                                  const LB_TRIG lbt,
                                                  const R_CONSTS rc,
                                                  const unsigned int convolve,
                                                  real* st_probs)
{
    unsigned int i;
    real h_prob, g, rg;
    mwvector xyz;
    real bg_prob = 0.0;

    zero_st_probs(st_probs, ap->number_streams);

    for (i = 0; i < convolve; ++i)
    {
        xyz = lbr2xyz_2(ap, r_pts[i].r_point, lbt);
        rg = rg_calc(ap, xyz);

        /* CHECKME: Not having quadratic term on slow one looks like a bug but I'm not sure */
        if (ap->fast_h_prob)
        {
            h_prob = h_prob_fast(ap, r_pts[i].qw_r3_N, rg);
            /* the Hernquist profile includes a quadratic term in g */
            if (ap->aux_bg_profile)
            {
                g = rc.gPrime + sg_dx[i];
                h_prob += aux_prob(ap, r_pts[i].qw_r3_N, g);
            }
        }
        else
        {
            h_prob = h_prob_slow(ap, r_pts[i].qw_r3_N, rg);
        }

        stream_sums(st_probs, sc, xyz, r_pts[i].qw_r3_N, ap->number_streams);

        bg_prob += h_prob;
    }

    return bg_prob;
}

real likelihood_bg_probability(const ASTRONOMY_PARAMETERS* ap,
                               const STREAM_CONSTANTS* sc,
                               const R_POINTS* r_pts,
                               const real* sg_dx,
                               const LB_TRIG lbt,
                               const R_CONSTS rc,
                               const real reff_xr_rp3,
                               real* st_probs)
{
    real bg_prob;

    /* if q is 0, there is no probability */
    if (ap->zero_q)
        return -1.0;

    bg_prob = likelihood_bg_probability_main(ap, sc, r_pts, sg_dx, lbt, rc, ap->convolve, st_probs);
    bg_prob *= reff_xr_rp3;

    mult_probs(st_probs, reff_xr_rp3, ap->number_streams);



    return bg_prob;
}

/* CHECKME: What is this? */
static real probability_log(real bg, real sum_exp_weights)
{
    return (bg == 0.0) ? -238.0 : mw_log10(bg / sum_exp_weights);
}

static real stream_sum(const unsigned int number_streams,
                       const FINAL_STREAM_INTEGRALS* fsi,
                       const real* st_prob,
                       KAHAN* st_only_sum,
                       const real* exp_stream_weights,
                       const real sum_exp_weights,
                       real bg_only)
{
    unsigned int i;
    real st_only;
    real star_prob = bg_only;

    for (i = 0; i < number_streams; ++i)
    {
        st_only = st_prob[i] / fsi->stream_integrals[i] * exp_stream_weights[i];
        star_prob += st_only;

        st_only = probability_log(st_only, sum_exp_weights);
        KAHAN_ADD(st_only_sum[i], st_only);
    }
    star_prob /= sum_exp_weights;

    return star_prob;
}

/* Populates exp_stream_weights, and returns the sum */
real get_exp_stream_weights(real* exp_stream_weights,
                            const STREAMS* streams,
                            real exp_background_weight)
{
    unsigned int i;
    real sum_exp_weights = exp_background_weight;
    for (i = 0; i < streams->number_streams; i++)
    {
        exp_stream_weights[i] = mw_exp(streams->stream_weight[i].weight);
        sum_exp_weights += exp_stream_weights[i];
    }

    sum_exp_weights *= 0.001;

    return sum_exp_weights;
}

void get_stream_only_likelihood(KAHAN* st_only_sum,
                                const unsigned int number_stars,
                                const unsigned int number_streams)
{
    unsigned int i;
    fprintf(stderr, "<stream_only_likelihood>");
    for (i = 0; i < number_streams; i++)
    {
        st_only_sum[i].sum += st_only_sum[i].correction;
        st_only_sum[i].sum /= number_stars;
        st_only_sum[i].sum -= 3.0;

        fprintf(stderr, " %.20lf", st_only_sum[i].sum);
    }
    fprintf(stderr, " </stream_only_likelihood>\n");
}

const int calculateSeparation = 1;
const int twoPanel = 1;

/* get stream & background weight constants */
static real get_stream_bg_weight_consts(StreamStats* ss, const STREAMS* streams)
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

static void twoPanelSeparation(const ASTRONOMY_PARAMETERS* ap,
                               const FINAL_STREAM_INTEGRALS* fsi,
                               StreamStats* ss,
                               const real* st_probs,
                               real bg_prob,
                               real epsilon_b)
{
    unsigned int i;
    real pbx, psgSum;

    pbx = epsilon_b * bg_prob / fsi->background_integral;

    for (i = 0; i < ap->number_streams; i++)
        ss[i].psg = ss[i].epsilon_s * st_probs[i] / fsi->stream_integrals[i];

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
                       const ASTRONOMY_PARAMETERS* ap,
                       const FINAL_STREAM_INTEGRALS* fsi,
                       const mwmatrix cmatrix,
                       StreamStats* ss,
                       const real* st_probs,
                       real bg_prob,
                       real epsilon_b,
                       mwvector current_star_point)
{
    mwvector starxyz;
    mwvector starxyzTransform;

    mwvector xsun = ZERO_VECTOR;
    X(xsun) = ap->m_sun_r0;

    if (twoPanel)
        twoPanelSeparation(ap, fsi, ss, st_probs, bg_prob, epsilon_b);
    else
        nonTwoPanelSeparation(ss, ap->number_streams);

    /* determine if star with sprob should be put into stream */
    int s_ok = prob_ok(ss, ap->number_streams);
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

/* separation init stuffs */
static void setSeparationConstants(const ASTRONOMY_PARAMETERS* ap,
                                   const FINAL_STREAM_INTEGRALS* fsi,
                                   mwmatrix cmatrix)
{
    unsigned int i;
    real d;
    mwvector dnormal;
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
    printf("bint: %lf\n", fsi->background_integral);
    for (i = 0; i < ap->number_streams; i++)
        printf("sint[%d]: %lf\n", i, fsi->stream_integrals[i]);
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

static real likelihood_sum(const ASTRONOMY_PARAMETERS* ap,
                           const STAR_POINTS* sp,
                           const STREAM_CONSTANTS* sc,
                           const STREAMS* streams,
                           const FINAL_STREAM_INTEGRALS* fsi,
                           const STREAM_GAUSS sg,
                           R_POINTS* r_pts,
                           KAHAN* st_sum,
                           real* st_prob,
                           const real* exp_stream_weights,
                           const real sum_exp_weights,
                           const real exp_background_weight,
                           StreamStats* ss,
                           const int do_separation,
                           FILE* f)
{
    KAHAN prob = ZERO_KAHAN;
    KAHAN bg_only = ZERO_KAHAN;

    unsigned int current_star_point;
    mwvector point;
    real star_prob;
    real bg_prob, bg;
    LB lb;
    LB_TRIG lbt;
    real reff_xr_rp3;
    R_CONSTS rc = { 0.0, 0.0 };

    real epsilon_b;
    mwmatrix cmatrix;
    unsigned int num_zero = 0;
    unsigned int bad_jacobians = 0;

    if (do_separation)
    {
        setSeparationConstants(ap, fsi, cmatrix);
        epsilon_b = get_stream_bg_weight_consts(ss, streams);
    }

    for (current_star_point = 0; current_star_point < sp->number_stars; ++current_star_point)
    {
        point = sp->stars[current_star_point];
        rc.gPrime = calcG(Z(point));
        set_r_points(ap, sg, ap->convolve, rc.gPrime, r_pts);
        reff_xr_rp3 = calcReffXrRp3(Z(point), rc.gPrime);

        LB_L(lb) = L(point);
        LB_B(lb) = B(point);

        lbt = lb_trig(lb);

        bg_prob = likelihood_bg_probability(ap, sc, r_pts, sg.dx, lbt, rc, reff_xr_rp3, st_prob);

        bg = (bg_prob / fsi->background_integral) * exp_background_weight;

        star_prob = stream_sum(streams->number_streams,
                               fsi,
                               st_prob,
                               st_sum,
                               exp_stream_weights,
                               sum_exp_weights,
                               bg);

        if (star_prob != 0.0)
        {
            star_prob = mw_log10(star_prob);
            KAHAN_ADD(prob, star_prob);
        }
        else
        {
            ++num_zero;
            prob.sum -= 238.0;
        }

        bg = probability_log(bg, sum_exp_weights);
        KAHAN_ADD(bg_only, bg);

        if (do_separation)
            separation(f, ap, fsi, cmatrix, ss, st_prob, bg_prob, epsilon_b, point);
    }

    prob.sum += prob.correction;
    bg_only.sum += bg_only.correction;
    bg_only.sum /= sp->number_stars;
    bg_only.sum -= 3.0;

    fprintf(stderr, "<background_only_likelihood> %.20lf </background_only_likelihood>\n", bg_only.sum);

    if (do_separation)
        printSeparationStats(ss, sp->number_stars, ap->number_streams);

    /*  log10(x * 0.001) = log10(x) - 3.0 */
    return (prob.sum / (sp->number_stars - bad_jacobians)) - 3.0;
}

StreamStats* newStreamStats(const unsigned int number_streams)
{
    return (StreamStats*) callocSafe(number_streams, sizeof(StreamStats));
}

real likelihood(const ASTRONOMY_PARAMETERS* ap,
                const STAR_POINTS* sp,
                const STREAM_CONSTANTS* sc,
                const STREAMS* streams,
                const FINAL_STREAM_INTEGRALS* fsi,
                const STREAM_GAUSS sg,
                const int do_separation,
                const char* separation_outfile)
{
    real* st_prob;
    R_POINTS* r_pts;
    KAHAN* st_sum;
    StreamStats* ss = NULL;
    real* exp_stream_weights;
    real sum_exp_weights;
    real exp_background_weight;
    real likelihood_val;
    FILE* f = NULL;

    if (do_separation)
    {
        f = mw_fopen(separation_outfile, "w");
        if (!f)
        {
            perror("Opening separation output file");
            return NAN;
        }

        ss = newStreamStats(streams->number_streams);
    }

    st_prob = (real*) mwMallocAligned(sizeof(real) * streams->number_streams, 2 * sizeof(real));
    r_pts = (R_POINTS*) mwMallocAligned(sizeof(R_POINTS) * ap->convolve, sizeof(R_POINTS));
    st_sum = (KAHAN*) mwCallocAligned(sizeof(KAHAN), streams->number_streams, sizeof(KAHAN));
    exp_stream_weights = (real*) mwMallocAligned(sizeof(real) * streams->number_streams, 2 * sizeof(real));

    exp_background_weight = mw_exp(ap->background_weight);
    sum_exp_weights = get_exp_stream_weights(exp_stream_weights, streams, exp_background_weight);

    likelihood_val = likelihood_sum(ap, sp, sc, streams, fsi,
                                    sg, r_pts,
                                    st_sum, st_prob,
                                    exp_stream_weights,
                                    sum_exp_weights,
                                    exp_background_weight,
                                    ss,
                                    do_separation,
                                    f);

    get_stream_only_likelihood(st_sum, sp->number_stars, streams->number_streams);

    mwAlignedFree(st_prob);
    mwAlignedFree(r_pts);
    mwAlignedFree(st_sum);
    mwAlignedFree(exp_stream_weights);
    free(ss);

    if (f && fclose(f))
        perror("Closing separation output file");

    return likelihood_val;
}

