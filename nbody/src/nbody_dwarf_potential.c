/*
 * Copyright (c) 2012 Rensselaer Polytechnic Institute
 * Copyright (c) 2016 Siddhartha Shelton
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
#include "nbody_potential_types.h"
#include "nbody_mass.h"

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*                             PLUMMER                                                                                   */
 static real plummer_pot(const Dwarf* model, real r)                                                                     //
{                                                                                                                        //
    const real mass = model->mass;                                                                                       //
    const real rscale = model->scaleLength;                                                                              //
    return mass / mw_sqrt(sqr(r) + sqr(rscale));                                                                         //
}                                                                                                                        //
                                                                                                                         //
 static real plummer_den(const Dwarf* model, real r)                                                                     //
{                                                                                                                        //
    const real mass = model->mass;                                                                                       //
    const real rscale = model->scaleLength;                                                                              //
    return  (3.0 / (4.0 * M_PI)) * (mass / cube(rscale)) * minusfivehalves( (1.0 + sqr(r)/sqr(rscale)) ) ;               //
}                                                                                                                        //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*                            NFW                                                                                        */
 static real nfw_den(const Dwarf* model, real r)                                                                         //
{                                                                                                                        //
    const real mass = model->mass;                                                                                       //
    const real rscale = model->scaleLength;                                                                              //
    return (1.0 / (4.0 * M_PI)) * (mass * rscale / r) * (1.0 / sqr(1.0 + r / rscale));                                   //
}                                                                                                                        //
                                                                                                                         //
 static real nfw_pot(const Dwarf* model, real r)                                                                         //
{                                                                                                                        //
    const real mass = model->mass;                                                                                       //
    const real rscale = model->scaleLength;                                                                              //
    return (mass / r) * mw_log(1.0 + r / rscale);                                                                        //
}                                                                                                                        //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*                             GENERAL HERNQUIST                                                                         */
static real gen_hern_den(const Dwarf* model, real r)                                                                     //
{                                                                                                                        //
    const real mass = model->mass;                                                                                       //
    const real rscale = model->scaleLength;                                                                              //
    return inv(2.0 * M_PI) * mass * rscale / ( r * cube(r + rscale));                                                    //
}                                                                                                                        //
                                                                                                                         //
static real gen_hern_pot(const Dwarf* model, real r)                                                                     //
{                                                                                                                        //
    const real mass = model->mass;                                                                                       //
    const real rscale = model->scaleLength;                                                                              //
    return mass / (r + rscale);                                                                                          //
}                                                                                                                        //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                            /* EINASTO */

static real einasto_den(const Dwarf* model, real r)
{
    const real mass = model->mass;
    const real h = model->scaleLength;
    const real n = model->n;

    real coeff = mass / ( 4.0 * M_PI * cube(h) * n * GammaFunc(3.0 * n));
    real thing = mw_pow(r, inv(n));
    return coeff * mw_exp(-thing);
}

static real einasto_pot(const Dwarf* model, real r)
{
    const real mass = model->mass;
    const real h = model->scaleLength;
    const real n = model->n;
    
    real coeff = mass / (h * r);
    real thing = mw_pow(r, 1.0 / n);
    
    real term = 1.0 - (   IncompleteGammaFunc(3.0 * n, thing) +  r * IncompleteGammaFunc(2.0 * n, thing)   ) / GammaFunc(3.0 * n);
    return coeff * term;
}



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

real get_potential(const Dwarf* model, real r)
{
    real pot_temp;
    
    switch(model->type)
    {
        case Plummer:
            pot_temp = plummer_pot(model, r);
            break;
        case NFW:
            pot_temp = nfw_pot(model, r );
            break;
        case General_Hernquist:
            pot_temp = gen_hern_pot(model, r );
            break;
//         case Einasto:
//             einasto_pot(model, r);
//             break;
        case InvalidDwarf:
        default:
            mw_fail("Invalid dwarf type\n");
    }

    return pot_temp;
}



real get_density(const Dwarf* model, real r)
{
    real den_temp;
    
    switch(model->type)
    {
        case Plummer:
            den_temp = plummer_den(model, r);
            break;
        case NFW:
            den_temp = nfw_den(model, r );
            break;
        case General_Hernquist:
            den_temp = gen_hern_den(model, r );
            break;
//         case Einasto:
//             einasto_den(model, r);
//             break;
        case InvalidDwarf:
        default:
            mw_fail("Invalid dwarf type");
            
    }
    
    return den_temp;
}

