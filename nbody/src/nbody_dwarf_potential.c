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

static real GammaFunc(const real z) 
{
    //Alogrithm for the calculation of the Lanczos Approx of the complete Gamma function 
    //as implemented in Numerical Recipes 3rd ed, 2007.
    real g = 4.7421875; //g parameter for the gamma function
    real x, tmp, y, A_g;
    
    //these are the cn's
    static const real coeff[14] = {57.1562356658629235,-59.5979603554754912,
                                14.1360979747417471,-0.491913816097620199,.339946499848118887e-4,
                                .465236289270485756e-4,-.983744753048795646e-4,.158088703224912494e-3,
                                -.210264441724104883e-3,.217439618115212643e-3,-.164318106536763890e-3,
                                .844182239838527433e-4,-.261908384015814087e-4,.368991826595316234e-5};
    y = x = z;
    tmp = x + g + 0.5;
    tmp = (x + 0.5) * mw_log(tmp) - tmp;
    A_g = 0.999999999999997092; //this is c0
    
    for (int j = 0; j < 14; j++) 
    {
        A_g += coeff[j] / ++y;
    } //calculates the series approx sum
        
    //sqrt(2 * pi) = 2.5066282746310005
    tmp += mw_log(2.5066282746310005 * A_g / x);//returns the log of the gamma function
    
    return mw_exp(tmp);
}
                            
                            
                            
static real einasto_den(const Dwarf* model, real r)
{
    const real mass = model->mass;
    const real h = model->h;
    const real s = 
    return mw_exp(-A * mw_pow(r, alpha));
}
// 
// static real einasto_pot(real r, real A, real alpha)
// {
// }



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
        
        case InvalidDwarf:
        default:
            mw_fail("Invalid dwarf type");
            
    }
    
    return den_temp;
}