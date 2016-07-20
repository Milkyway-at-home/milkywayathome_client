/*
 * Copyright (c) 2012 Rensselaer Polytechnic Institute
 *
 * This file is part of Milkway@Home.
 *
 * Milkyway@Home is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Milkyway@Home is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Milkyway@Home.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "nbody_dwarf_potential.h"
#include "milkyway_math.h"
#include "nbody_types.h"



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*                             PLUMMER                                                                                   */
 static real plummer_pot(real r, real mass, real rscale)                                                                 //
{                                                                                                                        //
    return mass / mw_sqrt(sqr(r) + sqr(rscale));                                                                         //
}                                                                                                                        //
                                                                                                                         //
 static real plummer_den(real r, real mass, real rscale)                                                                 //
{                                                                                                                        //
    return  (3.0 / (4.0 * M_PI)) * (mass / cube(rscale)) * minusfivehalves( (1.0 + sqr(r)/sqr(rscale)) ) ;               //
}                                                                                                                        //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*                            NFW                                                                                        */
 static real nfw_den(real r, real mass, real rscale)                                                                     //
{                                                                                                                        //
    return (1.0 / (4.0 * M_PI)) * (mass * rscale / r) * (1.0 / sqr(1.0 + r / rscale));                                   //
}                                                                                                                        //
                                                                                                                         //
 static real nfw_pot(real r, real mass, real rscale)                                                                     //
{                                                                                                                        //
    return (mass / r) * mw_log(1.0 + r / rscale);                                                                        //
}                                                                                                                        //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                            /* GENERAL HERNQUIST */
static real gen_hern_den(real r, real mass, real rscale)
{
    return inv(2.0 * M_PI) * mass * rscale / ( r * cube(r + rscale));
}

static real gen_hern_pot(real r, real mass, real rscale)
{
    return mass / (r + rscale);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                            /* EINASTO */

// static real einasto_den(real r, real A, real alpha)
// {
//     return mw_exp(-A * mw_pow(r, alpha));
// }
// 
// static real einasto_pot(real r, real A, real alpha)
// {
// }



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

real get_potential(real r, real * args, const int type)
{
    const int plum      = 1;
    const int nfw       = 2;
    const int gen_hern  = 3;
    const int einasto   = 4;
    if(type == plum)
    {
        real mass   = args[0];
        real rscale = args[1];
        return plummer_pot(r, mass, rscale);
    }
    else if(type == nfw)
    {
        real mass   = args[0];
        real rscale = args[1];
        return nfw_pot(r, mass, rscale);
    }
    else if(type == gen_hern)
    {
        real mass   = args[0];
        real rscale = args[1];
        return gen_hern_pot(r, mass, rscale);
    }
//     else if(type == einasto)
//     {
//         real A     = args[0];
//         real alpha = args[1];
//         return einasto_pot(r, mass, rscale);
//     }
    else
    {
        mw_fail("Invalid dwarf potential");
    }
}



real get_density(real r, real * args, const int type)
{
    const int plum      = 1;
    const int nfw       = 2;
    const int gen_hern  = 3;
    const int einasto   = 4;
    if(type == plum)
    {
        real mass   = args[0];
        real rscale = args[1];
        return plummer_den(r, mass, rscale);
    }
    else if(type == nfw)
    {
        real mass   = args[0];
        real rscale = args[1];
        return nfw_den(r, mass, rscale);
    }
        else if(type == gen_hern)
    {
        real mass   = args[0];
        real rscale = args[1];
        return gen_hern_den(r, mass, rscale);
    }
//         else if(type == einasto)
//     {
//         real A     = args[0];
//         real alpha = args[1];
//         return einasto_den(r, A, alpha);
//     }
    else
    {
        mw_fail("Invalid dwarf potential");
    }
}