/*
 *  Copyright (c) 2008-2010 Travis Desell, Nathan Cole
 *  Copyright (c) 2008-2010 Boleslaw Szymanski, Heidi Newbergb
 *  Copyright (c) 2008-2010 Carlos Varela, Malik Magdon-Ismail
 *  Copyright (c) 2008-2011 Rensselaer Polytechnic Institute
 *  Copyright (c) 2010-2011 Matthew Arsenault
 *  Copyright (c) 2013      Jake Weiss
 *  Copyright (c) 1991-2000 University of Groningen, The Netherlands
 *  Copyright (c) 2001-2009 The GROMACS Development Team
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
 *
 *
 *  This file incorporates work covered by the following copyright and
 *  permission notice:
 *
 *    Fast exp(x) computation (with SSE2 optimizations).
 *
 *  Copyright (c) 2010, Naoaki Okazaki
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *      * Redistributions of source code must retain the above copyright
 *        notice, this list of conditions and the following disclaimer.
 *      * Redistributions in binary form must reproduce the above copyright
 *        notice, this list of conditions and the following disclaimer in the
 *        documentation and/or other materials provided with the distribution.
 *      * Neither the names of the authors nor the names of its contributors
 *        may be used to endorse or promote products derived from this
 *        software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 *  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 *  OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 *  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 *  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 *  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "milkyway_util.h"
#include "probabilities.h"
#include "coordinates.h"


static inline mwvector lbr2xyz_2(const AstronomyParameters* ap, real rPoint, LBTrig lbt)
{
    mwvector xyz;

    xyz.x = mw_mad(rPoint, lbt.lCosBCos, ap->m_sun_r0);
    xyz.y = rPoint * lbt.lSinBCos;
    xyz.z = rPoint * lbt.bSin;

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
static inline real hernquist_prob_fast(const AstronomyParameters* ap, real qw_r3_N, real rg)
{
    const real rs = rg + ap->r0;
    return qw_r3_N / (rg * cube(rs));
}

HOT
static inline real hernquist_prob_slow(const AstronomyParameters* ap, real qw_r3_N, real rg)
{
    const real rs = rg + ap->r0;
    /* This is safe to do because we know if we are using slow hernquist the innerPower will be equal to the old alpha as it is checked elsewhere in the program */
    return qw_r3_N / (mw_powr(rg, ap->innerPower) * mw_powr(rs, ap->alpha_delta3));
}


/* Takes cylindical r and z outputs thick disk probability from Xu et al. (2015) */
HOT
static inline real disk_prob_thick(const AstronomyParameters* ap, real qw_r3_N, mwvector rtz)
{
    const real ls = -0.285714286; /*scale length, -1/(3.5kpc) from Xu et al. (2015) */
    const real hs = -1.428571429; /*scale height, -1/(0.7kpc) from Xu et al. (2015) */
    return qw_r3_N * exp(CR(rtz)*ls + fabs(CZ(rtz))*hs);
}

HOT
static inline real broken_power_law_prob(const AstronomyParameters* ap, real qw_r3_N, real rg)  //p(R)=p0(R/R0)^-n  Power Law Equation (From a paper)
{
	  const real n = ap->innerPower + (real)(rg >= ap->r0) *  ap->alpha_delta3; /* Basically an if statement for checking if we are past R0 ensuring proper cpu pipelining for efficiency */
    return qw_r3_N * mw_powr(ap->sun_r0 / rg, n);
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

static inline void zero_st_probs(real* st_probs, const unsigned int nstream)
{
    unsigned int i;

    for (i = 0; i < nstream; ++i)
        st_probs[i] = 0.0;
}

static inline real aux_prob(const AstronomyParameters* ap,
                            const real qw_r3_N,
                            const real r_in_mag)
{
    return qw_r3_N * (ap->bg_a * sqr(r_in_mag) + ap->bg_b * r_in_mag + ap->bg_c);
}


void calculate_background_weights(AstronomyParameters* ap)
{
/*
    mwvector rtz;
    CR(rtz) = ap->sun_r0;
    CT(rtz) = 0.0;
    CZ(rtz) = .027; /* Sun is 27pc above the midplane Xu et al. 2015 */

    
    /*Precalculate this stuff in the future to make faster*/
 /*   real solar_thick_density = disk_prob_thick(ap, 1.0, rtz);
    real solar_thin_density = disk_prob_thin(ap, 1.0, rtz);
    real solar_halo_density = hernquist_prob_slow(ap, 1.0, ap->sun_r0);
    
    /* 11.484375 = .91875 / .08 */
/*    ap->thin_disk_weight =  11.484375 * (solar_thick_density) / solar_thin_density;

    /* 11.484375 = .00125 / .99875 */
/*    ap->bg_weight = 0.001251564 * (solar_thick_density + ap->thin_disk_weight * solar_thin_density) / solar_halo_density;*/
}

/*These functions are rarely every used.  Only if SSE2 and OpenCL are not working/not compiled*/
HOT
real probabilities_fast_hprob(const AstronomyParameters* ap,
                              const StreamConstants* sc,
                              const real* RESTRICT sg_dx,
                              const real* RESTRICT r_point,
                              const real* RESTRICT qw_r3_N,
                              LBTrig lbt,
                              real gPrime,
                              real reff_xr_rp3,
                              real* RESTRICT streamTmps)
{
    int i;
    real h_prob, g, rg;
    mwvector xyz, rtz;
    real bg_prob = 0.0;
    int convolve = ap->convolve;
    int aux_bg_profile = ap->aux_bg_profile;
    real bg_weight = ap->background_weight;
    real thick_disk_weight = ap->thick_disk_weight;

    zero_st_probs(streamTmps, ap->number_streams);
    for (i = 0; i < convolve; ++i)
    {
        xyz = lbr2xyz_2(ap, r_point[i], lbt);
        rtz = xyz2cyl(ap, xyz);

        rg = rg_calc(ap, xyz);

        h_prob = bg_weight * hernquist_prob_fast(ap, qw_r3_N[i], rg) + thick_disk_weight * disk_prob_thick(ap, qw_r3_N[i], rtz);

        /* Add a quadratic term in g to the the Hernquist profile */
        if (aux_bg_profile)
        {
            g = gPrime + sg_dx[i];
            h_prob += bg_weight * aux_prob(ap, qw_r3_N[i], g);
        }

        bg_prob += h_prob;
        streamSums(streamTmps, sc, xyz, qw_r3_N[i], ap->number_streams);
    }

    bg_prob *= reff_xr_rp3;
    for (i = 0; i < ap->number_streams; ++i)
        streamTmps[i] *= reff_xr_rp3;

    return bg_prob;
}

HOT
real probabilities_slow_hprob(const AstronomyParameters* ap,
                              const StreamConstants* sc,
                              const real* RESTRICT sg_dx,
                              const real* RESTRICT r_point,
                              const real* RESTRICT qw_r3_N,
                              LBTrig lbt,
                              real gPrime,
                              real reff_xr_rp3,
                              real* RESTRICT streamTmps)
{
    int i;
    real rg, g;
    mwvector xyz, rtz;
    real bg_prob = 0.0;
    int convolve = ap->convolve;
    int aux_bg_profile = ap->aux_bg_profile;
    real bg_weight = ap->background_weight;
    real thick_disk_weight = ap->thick_disk_weight;

    zero_st_probs(streamTmps, ap->number_streams);

    for (i = 0; i < convolve; ++i)
    {
        xyz = lbr2xyz_2(ap, r_point[i], lbt);
        rtz = xyz2cyl(ap, xyz);
        
        rg = rg_calc(ap, xyz);

        bg_prob += bg_weight * hernquist_prob_slow(ap, qw_r3_N[i], rg) + thick_disk_weight * disk_prob_thick(ap, qw_r3_N[i], rtz);
        if (aux_bg_profile)
        {
            g = gPrime + sg_dx[i];
            bg_prob += bg_weight * aux_prob(ap, qw_r3_N[i], g);
        }

        streamSums(streamTmps, sc, xyz, qw_r3_N[i], ap->number_streams);
    }

    bg_prob *= reff_xr_rp3;
    for (i = 0; i < ap->number_streams; ++i)
        streamTmps[i] *= reff_xr_rp3;

    return bg_prob;
}

HOT
real probabilities_broken_power_law(const AstronomyParameters* ap,
                              const StreamConstants* sc,
                              const real* RESTRICT sg_dx,
                              const real* RESTRICT r_point,
                              const real* RESTRICT qw_r3_N,
                              LBTrig lbt,
                              real gPrime,
                              real reff_xr_rp3,
                              real* RESTRICT streamTmps)
{
    int i;
    real rg, g;
    mwvector xyz;
    real bg_prob = 0.0;
    int convolve = ap->convolve;
    int aux_bg_profile = ap->aux_bg_profile;

    zero_st_probs(streamTmps, ap->number_streams);

    for (i = 0; i < convolve; ++i)
    {
        xyz = lbr2xyz_2(ap, r_point[i], lbt);

        rg = rg_calc(ap, xyz);

        bg_prob += broken_power_law_prob(ap, qw_r3_N[i], rg);

        streamSums(streamTmps, sc, xyz, qw_r3_N[i], ap->number_streams);
    }

    bg_prob *= reff_xr_rp3;
    for (i = 0; i < ap->number_streams; ++i)
        streamTmps[i] *= reff_xr_rp3;

    return bg_prob;
}
