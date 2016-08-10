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



static real IncompleteGammaFunc1(real a, real x)
{
    //the series approx returns gamma from 0 to X but we want from X to INF
    //Therefore, we subtract it from GammaFunc which is from 0 to INF
    //The continued frac approx is already from X to INF
    
    static const real max_a = 100;
    real gamma = GammaFunc(a);
    if (x == 0.0) return 0.0;
    else if (a >= max_a) 
    {
        return gamm_approx_integral(a, x, 1); //Quadrature.
    }
    else if (x < a + 1.0) 
    {
        // Use the series representation. 
        return gamma - series_approx(a,x);
    }
    else 
    {
        //Use the continued fraction representation.
        return continued_frac_approx(a,x);
    }
}



// real IncompleteGammaFunc2(real a, real x)
// {
//     //THIS IS THE INC GAMMA FUNC WITH INT FROM x TO inf
//     static const real max_a = 100;
//     
//     if (x == 0.0) return 1.0;
//     else if (a >= max_a) 
//     {
//         return gamm_approx_integral(a, x, 0); //Quadrature.
//     }
//     else if (x < a + 1.0) 
//     {
//         // Use the series representation.
//         return 1.0 - series_approx(a,x);
//     }
//     else 
//     {
//         //Use the continued fraction representation.
//         return continued_frac_approx(a,x);
//     }
// }


static real series_approx(real a, real x)
{

    real sum, del, ap;
    ap = a;
    del = sum = 1.0 / a;//starting: gammma(a) / gamma(a+1) = 1/a
    for (;;) 
    {
        ++ap;
        del *= x / ap;
        sum += del;
        if (mw_fabs(del) < mw_fabs(sum) * 1.0e-15) 
        {
            return sum * exp(-x + a * log(x));
        }
    }
    
}



struct Gamma : Gauleg18 {
    // Object for incomplete gamma function. Gauleg18 provides coefficients for Gauss-Legendre
    // quadrature.
    static const int ASWITCH=100;
    // When to switch to quadrature method.
    static const real EPS;
    // See end of struct for initializations.
    static const real FPMIN;
    real gln;
//     real gammp(const real a, const real x) 
//     {
//         // Returns the incomplete gamma function P .a; x/.
//         if (x < 0.0 || a <= 0.0) throw("bad args in gammp");
//         if (x == 0.0) return 0.0;
//         else if ((int)a >= ASWITCH) return gammpapprox(a,x,1); //Quadrature.
//         else if (x < a+1.0) return gser(a,x);
//         // Use the series representation.
//         else return 1.0-gcf(a,x);
//         //Use the continued fraction representation.
//     }
// 
//     real gammq(const real a, const real x) {
//     //     Returns the incomplete gamma function Q.a; x/ Á 1 P .a; x/.
//         if (x < 0.0 || a <= 0.0) throw("bad args in gammq");
//         if (x == 0.0) return 1.0;
//         else if ((int)a >= ASWITCH) return gammpapprox(a,x,0); //Quadrature.
//         else if (x < a+1.0) return 1.0-gser(a,x);
//         //Use the series representation.
//         else return gcf(a,x);
//         //Use the continued fraction representation.
//     }
    real gser(const real a, const real x) 
    {
//         Returns the incomplete gamma function P .a; x/ evaluated by its series representation.
//         Also sets ln .a/ as gln. User should not call directly.
        real sum,del,ap;
        gln=gammln(a);
        ap=a;
        del=sum=1.0/a;
        for (;;) 
        {
            ++ap;
            del *= x/ap;
            sum += del;
            if (fabs(del) < fabs(sum)*EPS) 
            {
                return sum*exp(-x+a*log(x)-gln);
            }
        }
    }
    real gcf(const real a, const real x) 
    {
//         Returns the incomplete gamma function Q.a; x/ evaluated by its continued fraction rep-
//         resentation. Also sets ln .a/ as gln. User should not call directly.
        int i;
        real an,b,c,d,del,h;
        gln=gammln(a);
        b=x+1.0-a;
//         Set up for evaluating continued fraction
        c=1.0/FPMIN;
//         by modified Lentz’s method (5.2)
        d=1.0/b;
//         with b 0 D 0.
        h=d;
        for (i=1;;i++) 
        {
    //         Iterate to convergence.
            an = -i*(i-a);
            b += 2.0;
            d=an*d+b;
            Chapter 6. Special Functions
            if (fabs(d) < FPMIN) d=FPMIN;
            c=b+an/c;
            if (fabs(c) < FPMIN) c=FPMIN;
            d=1.0/d;
            del=d*c;
            h *= del;
            if (fabs(del-1.0) <= EPS) break;
        }
            return exp(-x+a*log(x)-gln)*h;
    //         Put factors in front.
    }
    real gammpapprox(real a, real x, int psig)
    {
//         Incomplete gamma by quadrature. Returns P .a; x/ or Q.a; x/, when psig is 1 or 0,
//         respectively. User should not call directly.
        int j;
        real xu,t,sum,ans;
        real a1 = a-1.0, lna1 = log(a1), sqrta1 = sqrt(a1);
        gln = gammln(a);
//         Set how far to integrate into the tail:
        if (x > a1) xu = MAX(a1 + 11.5*sqrta1, x + 6.0*sqrta1);
        else xu = MAX(0.,MIN(a1 - 7.5*sqrta1, x - 5.0*sqrta1));
        sum = 0;
        for (j=0;j<ngau;j++) 
        {
    //         Gauss-Legendre.
            t = x + (xu-x)*y[j];
            sum += w[j]*exp(-(t-a1)+a1*(log(t)-lna1));
        }
        ans = sum*(xu-x)*exp(a1*(lna1-1.)-gln);
        return (psig?(ans>0.0? 1.0-ans:-ans):(ans>=0.0? ans:1.0+ans));
    }
    real invgammp(real p, real a);
//     Inverse function on x of P .a; x/. See 6.2.1.
};
    
    const real Gamma::EPS = numeric_limits<real>::epsilon();
    const real Gamma::FPMIN = numeric_limits<real>::min()/EPS;
//     Remember that since Gamma is an object, you have to declare an instance of it
//     before you can use its member functions. We habitually write
    Gamma gam;
//     as a global declaration, and then call gam.gammp or gam.gammq as needed. The
//     structure Gauleg18 just contains the abscissas and weights for the Gauss-Legendre
//     quadrature.
    struct Gauleg18 
    {
    //     Abscissas and weights for Gauss-Legendre quadrature.
        static const int ngau = 18;
        static const real y[18];
        static const real w[18];
    };
    const real Gauleg18::y[18] = {0.0021695375159141994,
    0.011413521097787704,0.027972308950302116,0.051727015600492421,
    0.082502225484340941, 0.12007019910960293,0.16415283300752470,
    0.21442376986779355, 0.27051082840644336, 0.33199876341447887,
    0.39843234186401943, 0.46931971407375483, 0.54413605556657973,
    0.62232745288031077, 0.70331500465597174, 0.78649910768313447,
    0.87126389619061517, 0.95698180152629142};
    const real Gauleg18::w[18] = {0.0055657196642445571,
    0.012915947284065419,0.020181515297735382,0.027298621498568734,
    0.034213810770299537,0.040875750923643261,0.047235083490265582,
    0.053244713977759692,0.058860144245324798,0.064039797355015485,
    0.068745323835736408,0.072941885005653087,0.076598410645870640,
    0.079687828912071670,0.082187266704339706,0.084078218979661945,
    0.085346685739338721,0.085983275670394821};








































