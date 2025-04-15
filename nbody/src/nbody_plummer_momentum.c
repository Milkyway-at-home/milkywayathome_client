/*****************************************************\
* Calculates the angular momentum of a basic plummer  *
* dwarf accounting for differential effects to aid in *
* getting better reverse orbit parameters             *
*                                                     *
* NOTE: Currently assumes ∫ρ(rxv) dτ = (∫ρr dτ)xv     *
*       THIS MAY BE INVALID - DO DOUBLE CHECK MY WORK *
*                                                     *
* Aurora Chen - Apr. 2025 - Written for MilkyWay@Home *
\*****************************************************/

#include "nbody_priv.h"
//#include "milkyway_util.h"
//#include "milkyway_math.h"
//#include "milkyway_lua.h"
//#include "nbody_lua_types.h"
#include "nbody_plummer_momentum.h"

/*  r = radius from center of dwarf (integration parameter)     *
*   d = CoM dist to center of MW                                *
*   a = dwarf scale radius                                     */
static inline real dist_density_distr(real r, real d, real a)
{
    real plusR, minusR, denom;
    plusR = mw_pow((r+d), 2.0);
    minusR = mw_pow((r-d), 2.0);
    denom = mw_pow(r*r + a*a, 5.0);
    return r * (plusR * mw_sqrt(plusR/denom) - minusR * mw_sqrt(minusR/denom));
}

/* does a numerical integration to some distance k*a from the dwarf CoM        *
* accurate to ~4 decimal points with current delta r (upper bound ka / 50000) */
static inline mwvector mass_dist_vec(real scaleRadius, 
                                     mwvector pos,
                                     real dwarfMass,
                                     real k)
{
    real upperBound = k * scaleRadius;
    real delta = 0.00002 * upperBound;
    real dist = mw_length(pos);
    real mass_distr = 0;
    for (real r = 0; r <= upperBound; r += delta)
    {
        mass_distr += dist_density_distr(r, dist, scaleRadius) * delta;
    }
    real mag = mass_distr * dwarfMass * mw_pow(scaleRadius, 2.0) / (2.0 * dist);
    return mw_mulvs(pos, (mag / dist));
}

/* plummer dwarf enclosed mass as a proportion of k*a */
static inline real enclosedMass (real M0, real k)
{
    return M0 * mw_pow( (k / mw_sqrt(k*k + 1)), 3.0);
}

/* outputs velocity scaling factor based on calculated L versus. Mrxv */
real velocityAdj (real a, 
                  mwvector pos,
                  mwvector vel,
                  real M0,
                  real k)
{
    real massEnclosed = enclosedMass(M0, k);
    real L_base = massEnclosed * mw_length(mw_crossv(pos, vel));
    mwvector massDist = mass_dist_vec(a, pos, M0, k);
    real L_adj = mw_length(mw_crossv(massDist, vel));
    return L_adj/L_base;
}
