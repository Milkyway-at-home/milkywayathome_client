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


 static real plummer_pot(real r, real mass, real rscale)
{
    return mass / mw_sqrt(sqr(r) + sqr(rscale));
}

 static real plummer_den(real r, real mass, real rscale)
{
    return  (3.0 / (4.0 * M_PI)) * (mass / cube(rscale)) * minusfivehalves( (1.0 + sqr(r)/sqr(rscale)) ) ;
}

 static real nfw_den(real r, real mass, real rscale)
{
    return (1.0 / (4.0 * M_PI)) * (mass * rscale / r) * (1.0 / sqr(1.0 + r / rscale));
}

 static real nfw_pot(real r, real mass, real rscale)
{
    return (mass / r) * mw_log(1.0 + r / rscale);
}


real get_potential(real r, real * args, int type)
{
    static int plum = 1;
    static int nfw  = 2;
    switch (type)
    {
        case plum:
            real mass   = args[0];
            real rscale = args[1];
            return plummer_pot(r, mass, rscale);
        case nfw:
            real mass   = args[0];
            real rscale = args[1];
            return nfw_pot(r, mass, rscale);
        default:
            mw_fail("Invalid dwarf potential");
    }
}



real get_density(real r, real * args, int type)
{
    static int plum = 1;
    static int nfw  = 2;
    switch (type)
    {
        case plum:
            real mass   = args[0];
            real rscale = args[1];
            return plummer_den(r, mass, rscale);
        case nfw:
            real mass   = args[0];
            real rscale = args[1];
            return nfw_den(r, mass, rscale);
        default:
            mw_fail("Invalid dwarf density");
    }
}