/*
 * Copyright (c) 2012 Rensselaer Polytechnic Institute
 * Copyright (c) 2016-2018 Siddhartha Shelton
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
#include "nbody_potential_types.h"
#include "nbody_mass.h"

/* NOTE
 * we want the term nu which is the density per mass unit. However, these return just normal density.
 * In galactic dynamics 2nd edition, equation 4.48 defines nu which the distribution function is written 
 * in terms of. However, this mass is not the mass of each component but the total mass of both. Therefore,
 * this term can be pulled out of the integral. since we are rejection sampling, it cancels with the denom. 
 * It does not change the distribution so it would be ok if it had not canceled.
 * the potential functions return the negative version of the potential, psi, which is what is needed. 
 */


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*                             PLUMMER                                                                                   */
/* this potential and density are both taken from binney 2nd ed                                                          */
static real_0 plummer_den(const Dwarf* model, real_0 r)                                                                  //
{                                                                                                                        //
    const real_0 mass = model->mass;                                                                                     //
    const real_0 rscale = model->scaleLength;                                                                            //
    return  (3.0 / (4.0 * M_PI)) * (mass / cube_0(rscale)) * minusfivehalves_0( (1.0 + sqr_0(r / rscale)) ) ;            //
}                                                                                                                        //
                                                                                                                         //
static real plummer_den_real(real* mass, real* rscale, real* r)                                                          //
{                                                                                                                        //
    real tmp1, tmp2;                                                                                                     //
    tmp1 =  cube(rscale);                                                                                                //
    tmp1 =  mw_div(mass, &tmp1);                                                                                         //
    tmp2 =  mw_div(r, rscale);                                                                                           //
    tmp2 = sqr(&tmp2);                                                                                                   //
    tmp2 = mw_add_s(&tmp2, 1.0);                                                                                         //
    tmp2 = minusfivehalves(&tmp2);                                                                                       //
    tmp1 = mw_mul(&tmp1, &tmp2);                                                                                         //
    real den = mw_mul_s(&tmp1, (3.0 / (4.0 * M_PI)));                                                                     //
    return  den;                                                                                                         //
}                                                                                                                        //
                                                                                                                         //
static real_0 d_dr_plummer_den(const Dwarf* model, real_0 r)
{
    const real_0 mass = model->mass;
    const real_0 rscale = model->scaleLength;
    return (-15.0 / (4.0*M_PI)) * r * mass / fifth_0(rscale)
           * minusfivehalves_0( (1.0 + sqr_0(r / rscale)) ) / (1.0 + sqr_0(r / rscale));
}

static real d_dr_plummer_den_real(real* mass, real* rscale, real* r)
{
    real tmp1, tmp2;
    tmp1 = mw_mul_s(r, (-15.0 / (4.0*M_PI)));
    tmp1 = mw_mul(&tmp1, mass);
    tmp2 = fifth(rscale);
    tmp1 = mw_div(&tmp1, &tmp2);
    tmp2 = mw_div(r, rscale);
    tmp2 = sqr(&tmp2);
    tmp2 = mw_add_s(&tmp2, 1.0);
    tmp2 = minusfivehalves(&tmp2);
    tmp1 = mw_mul(&tmp1, &tmp2);
    tmp2 = mw_div(r, rscale);
    tmp2 = sqr(&tmp2);
    tmp2 = mw_add_s(&tmp2, 1.0);
    real d_dr = mw_div(&tmp1, &tmp2);
    return d_dr;
}

static real_0 d2_dr2_plummer_den(const Dwarf* model, real_0 r)
{
    const real_0 mass = model->mass;
    const real_0 rscale = model->scaleLength;
    return (-15.0 / (4.0*M_PI)) * mass * (sqr_0(rscale) - 6.0*sqr_0(r)) / fifth_0(rscale) / sqr_0(rscale)
           * minusfivehalves_0( (1.0 + sqr_0(r / rscale)) ) / sqr_0(1.0 + sqr_0(r / rscale));
}

static real d2_dr2_plummer_den_real(real* mass, real* rscale, real* r)
{
    real tmp1, tmp2, tmp3;
    tmp1 = mw_mul_s(mass, (-15.0 / (4.0*M_PI)));
    tmp2 = sqr(rscale);
    tmp3 = sqr(r);
    tmp3 = mw_mul_s(&tmp3, 6.0);
    tmp2 = mw_sub(&tmp2, &tmp3);
    tmp1 = mw_mul(&tmp1, &tmp2);
    tmp2 = fifth(rscale);
    tmp1 = mw_div(&tmp1, &tmp2);
    tmp2 = sqr(rscale);
    tmp1 = mw_div(&tmp1, &tmp2);
    tmp2 = mw_div(r, rscale);
    tmp2 = sqr(&tmp2);
    tmp2 = mw_add_s(&tmp2, 1.0);
    tmp2 = minusfivehalves(&tmp2);
    tmp1 = mw_mul(&tmp1, &tmp2);
    tmp2 = mw_div(r, rscale);
    tmp2 = sqr(&tmp2);
    tmp2 = mw_add_s(&tmp2, 1.0);
    tmp2 = sqr(&tmp2);
    real d2_dr2 = mw_div(&tmp1, &tmp2);
    return  d2_dr2;
}

 static real_0 plummer_pot(const Dwarf* model, real_0 r)                                                                 //
{                                                                                                                        //
    const real_0 mass = model->mass;                                                                                     //
    const real_0 rscale = model->scaleLength;                                                                            //
    return mass / mw_hypot_0(r, rscale);                                                                                 //
}                                                                                                                        //
                                                                                                                         //
 static real plummer_pot_real(real* mass, real* rscale, real* r)                                                       //
{                                                                                                                        //
    real tmp = mw_hypot(r, rscale);                                                                                      //
    real pot = mw_div(mass, &tmp);                                                                                       //
    return pot;                                                                                                          //
}                                                                                                                        //

static real_0 d_dr_plummer_pot(const Dwarf* model, real_0 r)
{
    const real_0 mass = model->mass;
    const real_0 rscale = model->scaleLength;
    return -mass * r / cube_0(mw_hypot_0(r, rscale));
}

static real d_dr_plummer_pot_real(real* mass, real* rscale, real* r)
{
    real tmp1, tmp2;
    tmp1 = mw_mul(mass, r);
    tmp1 = mw_neg(&tmp1);
    tmp2 = mw_hypot(r, rscale);
    tmp2 = cube(&tmp2);
    real d_dr = mw_div(&tmp1, &tmp2);
    return d_dr;
}

static real_0 d2_dr2_plummer_pot(const Dwarf* model, real_0 r)
{
    const real_0 mass = model->mass;
    const real_0 rscale = model->scaleLength;
    return mass * (2*sqr_0(r) - sqr_0(rscale)) / fifth_0(mw_hypot_0(r, rscale));
}

static real d2_dr2_plummer_pot_real(real* mass, real* rscale, real* r)
{
    real tmp1, tmp2;
    tmp1 = sqr(r);
    tmp1 = mw_mul_s(&tmp1, 2.0);
    tmp2 = sqr(rscale);
    tmp1 = mw_sub(&tmp1, &tmp2);
    tmp1 = mw_mul(mass, &tmp1);
    tmp2 = mw_hypot(r, rscale);
    tmp2 = fifth(&tmp2);
    real d2_dr2 = mw_div(&tmp1, &tmp2);
    return d2_dr2;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*                            NFW                                                                                        */
/* this density is taken from the 1997 paper by nfw. the potential is taken from binney 2nd ed                           */                         
 static real_0 nfw_den(const Dwarf* model, real_0 r)                                                                     //
{                                                                                                                        //                                                                                       //
    const real_0 rscale = model->scaleLength;                                                                            //
    const real_0 p0 = model->p0;                                                                                         //
    real_0 R = r / rscale;                                                                                               //
    /* at r = 0 the density goes to inf. however, the sampling is guarded against r = 0 anyway.*/                        //
    return p0 * inv_0(R) * inv_0(sqr_0(1.0 + R));                                                                        //
}                                                                                                                        //
                                                                                                                         //
static real nfw_den_real(real* p0, real* rscale, real* r)                                                                //
{                                                                                                                        //
    real tmp1, tmp2;                                                                                                     //
    tmp1 = mw_div(r, rscale);                                                                                            //
    tmp2 = mw_add_s(&tmp1, 1.0);                                                                                         //
    tmp2 = sqr(&tmp2);                                                                                                   //
    real den = mw_div(p0, &tmp1);                                                                                        //
    den = mw_div(&den, &tmp2);                                                                                           //
    return  den;                                                                                                         //                                                                                                        //
}                                                                                                                        //

static real_0 d_dr_nfw_den(const Dwarf* model, real_0 r)
{
    const real_0 rscale = model->scaleLength;
    const real_0 p0 = model->p0; 
    return -p0 * cube_0(rscale) * (rscale + 3*r) / sqr_0(r) / cube_0(rscale + r);
}

static real d_dr_nfw_den_real(real* p0, real* rscale, real* r)
{
    real tmp1, tmp2;
    tmp1 = cube(rscale);
    tmp1 = mw_mul(p0, &tmp1);
    tmp1 = mw_neg(&tmp1);
    tmp2 = mw_mul_s(r, 3.0);
    tmp2 = mw_add(rscale, &tmp2);
    tmp1 = mw_mul(&tmp1, &tmp2);
    tmp2 = sqr(r);
    tmp1 = mw_div(&tmp1, &tmp2);
    tmp2 = mw_add(rscale, r);
    tmp2 = cube(&tmp2);
    real d_dr = mw_div(&tmp1, &tmp2);
    return d_dr;
}

static real_0 d2_dr2_nfw_den(const Dwarf* model, real_0 r)
{
    const real_0 rscale = model->scaleLength;
    const real_0 p0 = model->p0; 
    return 2.0 * p0 * cube_0(rscale) * (sqr_0(rscale) + 4.0*rscale*r + 6.0*sqr_0(r)) / cube_0(r) / fourth_0(rscale + r);
}

static real d2_dr2_nfw_den_real(real* p0, real* rscale, real* r)
{
    real tmp1, tmp2, tmp3;
    tmp1 = mw_mul_s(p0, 2.0);
    tmp2 = cube(rscale);
    tmp1 = mw_mul(&tmp1, &tmp2);
    tmp2 = sqr(rscale);
    tmp3 = mw_mul(rscale, r);
    tmp3 = mw_mul_s(&tmp3, 4.0);
    tmp2 = mw_add(&tmp2, &tmp3);
    tmp3 = sqr(r);
    tmp3 = mw_mul_s(&tmp3, 6.0);
    tmp2 = mw_add(&tmp2, &tmp3);
    tmp1 = mw_mul(&tmp1, &tmp2);
    tmp2 = cube(r);
    tmp1 = mw_div(&tmp1, &tmp2);
    tmp2 = mw_add(rscale, r);
    tmp2 = fourth(&tmp2);
    real d2_dr2 = mw_div(&tmp1, &tmp2);
    return d2_dr2;
}
                                                                                                                         //
 static real_0 nfw_pot(const Dwarf* model, real_0 r)                                                                     //
{                                                                                                                        //                                                                                      //
    const real_0 rscale = model->scaleLength;                                                                            //
    const real_0 p0 = model->p0;                                                                                         //
    real_0 R = r / rscale;                                                                                               //
    /* at r = 0 the pot goes to inf. however, the sampling is guarded against r = 0 anyway. */                           //
    return  4.0 * M_PI * sqr_0(rscale) * p0 * inv_0(R) * mw_log_0(1.0 + R);                                              //
}                                                                                                                        //
                                                                                                                         //
static real nfw_pot_real(real* p0, real* rscale, real* r)                                                                //
{                                                                                                                        //
    real tmp1, tmp2;                                                                                                     //
    tmp1 = sqr(rscale);                                                                                                  //
    tmp2 = mw_div(r, rscale);                                                                                            //
    tmp2 = inv(&tmp2);                                                                                                   //
    tmp1 = mw_mul(&tmp1, p0);                                                                                            //
    tmp1 = mw_mul(&tmp1, &tmp2);                                                                                         //
    tmp2 = mw_div(r, rscale);                                                                                            //
    tmp2 = mw_add_s(&tmp2, 1.0);                                                                                         //
    tmp2 = mw_log(&tmp2);                                                                                                //
    tmp1 = mw_mul(&tmp1, &tmp2);                                                                                         //
    real pot = mw_mul_s(&tmp1, 4*M_PI);                                                                                  //
    return  pot;                                                                                                         //
}                                                                                                                        //

static real_0 d_dr_nfw_pot(const Dwarf* model, real_0 r)
{
    const real_0 rscale = model->scaleLength;
    const real_0 p0 = model->p0; 
    return -4.0*M_PI * p0 * cube_0(rscale) * ( (rscale + r)*mw_log_0(1 + (r/rscale)) - r ) / sqr_0(r) / (rscale + r);
}

static real d_dr_nfw_pot_real(real* p0, real* rscale, real* r)
{
    real tmp1, tmp2, tmp3;
    tmp1 = mw_mul_s(p0, -4.0*M_PI);
    tmp2 = cube(rscale);
    tmp1 = mw_mul(&tmp1, &tmp2);
    tmp2 = mw_add(rscale, r);
    tmp3 = mw_div(r, rscale);
    tmp3 = mw_add_s(&tmp3, 1.0);
    tmp3 = mw_log(&tmp3);
    tmp2 = mw_mul(&tmp2, &tmp3);
    tmp2 = mw_sub(&tmp2, r);
    tmp1 = mw_mul(&tmp1, &tmp2);
    tmp2 = sqr(r);
    tmp1 = mw_div(&tmp1, &tmp2);
    tmp2 = mw_add(rscale, r);
    real d_dr = mw_div(&tmp1, &tmp2);
    return d_dr;
}

static real_0 d2_dr2_nfw_pot(const Dwarf* model, real_0 r)
{
    const real_0 rscale = model->scaleLength;
    const real_0 p0 = model->p0; 
    return 4.0*M_PI * p0 * cube_0(rscale) * ( 2.0*sqr_0(rscale + r)*mw_log_0(1 + (r/rscale))
            - r * (2.0 * rscale + 3.0 * r)) / cube_0(r) / sqr_0(rscale + r);
}

static real d2_dr2_nfw_pot_real(real* p0, real* rscale, real* r)
{
    real tmp1, tmp2, tmp3, tmp4;
    tmp1 = mw_mul_s(p0, 4.0*M_PI);
    tmp2 = cube(rscale);
    tmp1 = mw_mul(&tmp1, &tmp2);
    tmp2 = mw_add(rscale, r);
    tmp2 = sqr(&tmp2);
    tmp2 = mw_mul_s(&tmp2, 2.0);
    tmp3 = mw_div(r, rscale);
    tmp3 = mw_log1p(&tmp3);
    tmp2 = mw_mul(&tmp2, &tmp3);
    tmp3 = mw_mul_s(rscale, 2.0);
    tmp4 = mw_mul_s(r, 3.0);
    tmp3 = mw_add(&tmp3, &tmp4);
    tmp3 = mw_mul(r, &tmp3);
    tmp2 = mw_sub(&tmp2, &tmp3);
    tmp1 = mw_mul(&tmp1, &tmp2);
    tmp2 = cube(r);
    tmp1 = mw_div(&tmp1, &tmp2);
    tmp2 = mw_add(rscale, r);
    tmp2 = sqr(&tmp2);
    real d2_dr2 = mw_div(&tmp1, &tmp2);
    return d2_dr2;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*                             GENERAL HERNQUIST                                                                         */
/* this potential and density are both taken from the 1990 paper by hernquist                                            */
static real_0 gen_hern_den(const Dwarf* model, real_0 r)                                                                 //
{                                                                                                                        //
    const real_0 mass = model->mass;                                                                                     //
    const real_0 rscale = model->scaleLength;                                                                            //
    return inv_0(2.0 * M_PI) * mass * rscale / ( r * cube_0(r + rscale));                                                //
}                                                                                                                        //

static real gen_hern_den_real(real* mass, real* rscale, real* r)                                                         //
{                                                                                                                        //
    real tmp = mw_add(r, rscale);
    tmp = cube(&tmp);
    tmp = mw_mul(r, &tmp);
    tmp = mw_div(rscale, &tmp);
    tmp = mw_mul(mass, &tmp);
    real den = mw_mul_s(&tmp, inv_0(2.0 * M_PI));
    return den;
}

static real_0 d_dr_gen_hern_den(const Dwarf* model, real_0 r)
{
    const real_0 mass = model->mass;
    const real_0 rscale = model->scaleLength;
    return -mass * rscale * (rscale + 4.0*r) / sqr_0(r) / fourth_0(rscale + r) / (2.0*M_PI);
}

static real d_dr_gen_hern_den_real(real* mass, real* rscale, real* r)
{
    real tmp1, tmp2;
    tmp1 = mw_neg(mass);
    tmp1 = mw_mul(&tmp1, rscale);
    tmp2 = mw_mul_s(r, 4.0);
    tmp2 = mw_add(rscale, &tmp2);
    tmp1 = mw_mul(&tmp1, &tmp2);
    tmp2 = sqr(r);
    tmp1 = mw_div(&tmp1, &tmp2);
    tmp2 = mw_add(rscale, r);
    tmp2 = fourth(&tmp2);
    tmp1 = mw_div(&tmp1, &tmp2);
    real d_dr = mw_mul_s(&tmp1, inv_0(2.0*M_PI));
    return d_dr;
}

static real_0 d2_dr2_gen_hern_den(const Dwarf* model, real_0 r)
{
    const real_0 mass = model->mass;
    const real_0 rscale = model->scaleLength;
    return mass * rscale * (sqr_0(rscale) + 5.0*rscale*r + 10.0*sqr_0(r)) / cube_0(r) / fifth_0(rscale + r) / M_PI;
}

static real d2_dr2_gen_hern_den_real(real* mass, real* rscale, real* r)
{
    real tmp1, tmp2;
    tmp1 = sqr(rscale);
    tmp2 = mw_mul(rscale, r);
    tmp2 = mw_mul_s(&tmp2, 5.0);
    tmp1 = mw_add(&tmp1, &tmp2);
    tmp2 = sqr(r);
    tmp2 = mw_mul_s(&tmp2, 10.0);
    tmp1 = mw_add(&tmp1, &tmp2);
    tmp1 = mw_mul(rscale, &tmp1);
    tmp1 = mw_mul(mass, &tmp1);
    tmp2 = cube(r);
    tmp1 = mw_div(&tmp1, &tmp2);
    tmp2 = mw_add(rscale, r);
    tmp2 = fifth(&tmp2);
    tmp1 = mw_div(&tmp1, &tmp2);
    real d2_dr2 = mw_mul_s(&tmp1, inv_0(M_PI));
    return d2_dr2;
}
                                                                                                                         //
static real_0 gen_hern_pot(const Dwarf* model, real_0 r)                                                                 //
{                                                                                                                        //
    const real_0 mass = model->mass;                                                                                     //
    const real_0 rscale = model->scaleLength;                                                                            //
    return mass / (r + rscale);                                                                                          //
}                                                                                                                        //

static real gen_hern_pot_real(real* mass, real* rscale, real* r)                                                                 //
{                                                                                                                        //
    real tmp = mw_add(r, rscale);
    real pot = mw_div(mass, &tmp);
    return pot;                                                                                          //
}                                                                                                                        //

static real_0 d_dr_gen_hern_pot(const Dwarf* model, real_0 r)
{
    const real_0 mass = model->mass;
    const real_0 rscale = model->scaleLength;
    return -mass / sqr_0(rscale + r);
}

static real d_dr_gen_hern_pot_real(real* mass, real* rscale, real* r)
{
    real tmp;
    tmp = mw_add(rscale, r);
    tmp = sqr(&tmp);
    tmp = mw_div(mass, &tmp);
    real d_dr = mw_neg(&tmp);
    return d_dr;
}

static real_0 d2_dr2_gen_hern_pot(const Dwarf* model, real_0 r)
{
    const real_0 mass = model->mass;
    const real_0 rscale = model->scaleLength;
    return 2.0 * mass / cube_0(rscale + r);
}

static real d2_dr2_gen_hern_pot_real(real* mass, real* rscale, real* r)
{
    real tmp;
    tmp = mw_add(rscale, r);
    tmp = cube(&tmp);
    tmp = mw_div(mass, &tmp);
    real d2_dr2 = mw_mul_s(&tmp, 2.0);
    return d2_dr2;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*                             EINASTO                                                                                   */
/* these are taken from the einasto paper. There are many problems with this, so it is currently unused.                 */
static real_0 einasto_den(const Dwarf* model, real_0 r)                                                                  //
{                                                                                                                        //
    const real_0 mass = model->mass;                                                                                     //
    const real_0 h = model->scaleLength;                                                                                 //
    const real_0 n = model->n;                                                                                           //
                                                                                                                         //
    real_0 coeff = 1.0 / ( 4.0 * M_PI * cube_0(h) * n * GammaFunc(3.0 * n));                                             //
    real_0 thing = mw_pow_0(r, inv_0(n));                                                                                //
    return coeff * mw_exp_0(-thing);                                                                                     //
}                                                                                                                        //
                                                                                                                         //
static real_0 einasto_pot(const Dwarf* model, real_0 r)                                                                  //
{                                                                                                                        //
    const real_0 mass = model->mass;                                                                                     //
    const real_0 h = model->scaleLength;                                                                                 //
    const real_0 n = model->n;                                                                                           //
                                                                                                                         //
    real_0 coeff = mass / (h * r);                                                                                       //
    real_0 thing = mw_pow_0(r, 1.0 / n);                                                                                 //
                                                                                                                         //
    real_0 term1 = IncompleteGammaFunc(3.0 * n, thing);                                                                  //
    real_0 term2 = r * IncompleteGammaFunc(2.0 * n, thing);                                                              //
    real_0 term = 1.0 - ( term1 + term2 ) / GammaFunc(3.0 * n);                                                          //
    return coeff * term;                                                                                                 //
}                                                                                                                        //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

static real get_p0(real* mass, real* scale)
{
    real tmp1, tmp2;

    tmp1 = mw_mul_s(mass, 1.0/(vol_pcrit));
    real r200 = mw_cbrt(&tmp1);

    real c = mw_div(&r200, scale);

    tmp1 = mw_add_s(&c, 1.0);
    tmp1 = mw_div(&c, &tmp1);
    tmp2 = mw_log1p(&c);
    real term = mw_sub(&tmp2, &tmp1);

    tmp1 = mw_mul_s(&term, 3.0);
    tmp1 = inv(&tmp1);
    tmp2 = cube(&c);
    tmp1 = mw_mul(&tmp1, &tmp2);
    real p0 = mw_mul_s(&tmp1, 200.0 * pcrit);
    return p0;
}

real_0 get_potential(const Dwarf* model, real_0 r)
{
    real_0 pot_temp;
    
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

real get_potential_real(const Dwarf* model_light, const Dwarf* model_dark, real* r, mwbool isLight)
{
    real pot_temp, tmp, rho_light, rho_dark;

    real_0 scale_b = model_light->scaleLength;
    real_0 mass_b = model_light->originmass;
    real_0 xi_scale = model_light->scaleLength / (model_light->scaleLength + model_dark->scaleLength);
    real_0 xi_mass = model_light->originmass / (model_light->originmass + model_dark->originmass);

    real scale_light = mw_real_var(scale_b, BARYON_RADIUS_POS);
    real mass_light = mw_real_var(mass_b, BARYON_MASS_POS);
    real xi_R = mw_real_var(xi_scale, RADIUS_RATIO_POS);
    real xi_M = mw_real_var(xi_mass, MASS_RATIO_POS);

    tmp = inv(&xi_R);
    tmp = mw_add_s(&tmp, -1.0);
    real scale_dark = mw_mul(&tmp, &scale_light);

    tmp = inv(&xi_M);
    tmp = mw_add_s(&tmp, -1.0);
    real mass_dark = mw_mul(&tmp, &mass_light);

    if(isLight)
    {
        switch(model_light->type)
        {
            case Plummer:
                pot_temp = plummer_pot_real(&mass_light, &scale_light, r);
                break;
            case NFW:
                rho_light = get_p0(&mass_light, &scale_light);
                pot_temp = nfw_pot_real(&rho_light, &scale_light, r);
                break;
            case General_Hernquist:
                pot_temp = gen_hern_pot_real(&mass_light, &scale_light, r);
                break;
    //         case Einasto:
    //             einasto_pot(model, r);
    //             break;
            case InvalidDwarf:
            default:
                mw_fail("Invalid dwarf type");     
        }
    }
    else
    {
        switch(model_dark->type)
        {
            case Plummer:
                pot_temp = plummer_pot_real(&mass_dark, &scale_dark, r);
                break;
            case NFW:
                rho_dark = get_p0(&mass_dark, &scale_dark);
                pot_temp = nfw_pot_real(&rho_dark, &scale_dark, r);
                break;
            case General_Hernquist:
                pot_temp = gen_hern_pot_real(&mass_dark, &scale_dark, r);
                break;
    //         case Einasto:
    //             einasto_pot(model, r);
    //             break;
            case InvalidDwarf:
            default:
                mw_fail("Invalid dwarf type");     
        }
    }
    
    return pot_temp;
}

real_0 get_first_derv_potential(const Dwarf* model, real_0 r)
{
    real_0 pot_temp;
    
    switch(model->type)
    {
        case Plummer:
            pot_temp = d_dr_plummer_pot(model, r);
            break;
        case NFW:
            pot_temp = d_dr_nfw_pot(model, r );
            break;
        case General_Hernquist:
            pot_temp = d_dr_gen_hern_pot(model, r );
            break;
        case InvalidDwarf:
        default:
            mw_fail("Invalid dwarf type\n");
    }

    return pot_temp;
}

real get_first_derv_potential_real(const Dwarf* model_light, const Dwarf* model_dark, real* r, mwbool isLight)
{
    real pot_temp, tmp, rho_light, rho_dark;

    real_0 scale_b = model_light->scaleLength;
    real_0 mass_b = model_light->originmass;
    real_0 xi_scale = model_light->scaleLength / (model_light->scaleLength + model_dark->scaleLength);
    real_0 xi_mass = model_light->originmass / (model_light->originmass + model_dark->originmass);

    real scale_light = mw_real_var(scale_b, BARYON_RADIUS_POS);
    real mass_light = mw_real_var(mass_b, BARYON_MASS_POS);
    real xi_R = mw_real_var(xi_scale, RADIUS_RATIO_POS);
    real xi_M = mw_real_var(xi_mass, MASS_RATIO_POS);

    tmp = inv(&xi_R);
    tmp = mw_add_s(&tmp, -1.0);
    real scale_dark = mw_mul(&tmp, &scale_light);

    tmp = inv(&xi_M);
    tmp = mw_add_s(&tmp, -1.0);
    real mass_dark = mw_mul(&tmp, &mass_light);

    if(isLight)
    {
        switch(model_light->type)
        {
            case Plummer:
                pot_temp = d_dr_plummer_pot_real(&mass_light, &scale_light, r);
                break;
            case NFW:
                rho_light = get_p0(&mass_light, &scale_light);
                pot_temp = d_dr_nfw_pot_real(&rho_light, &scale_light, r);
                break;
            case General_Hernquist:
                pot_temp = d_dr_gen_hern_pot_real(&mass_light, &scale_light, r);
                break;
            case InvalidDwarf:
            default:
                mw_fail("Invalid dwarf type");     
        }
    }
    else
    {
        switch(model_dark->type)
        {
            case Plummer:
                pot_temp = d_dr_plummer_pot_real(&mass_dark, &scale_dark, r);
                break;
            case NFW:
                rho_dark = get_p0(&mass_dark, &scale_dark);
                pot_temp = d_dr_nfw_pot_real(&rho_dark, &scale_dark, r);
                break;
            case General_Hernquist:
                pot_temp = d_dr_gen_hern_pot_real(&mass_dark, &scale_dark, r);
                break;
            case InvalidDwarf:
            default:
                mw_fail("Invalid dwarf type");     
        }
    }
    
    return pot_temp;
}

real_0 get_second_derv_potential(const Dwarf* model, real_0 r)
{
    real_0 pot_temp;
    
    switch(model->type)
    {
        case Plummer:
            pot_temp = d2_dr2_plummer_pot(model, r);
            break;
        case NFW:
            pot_temp = d2_dr2_nfw_pot(model, r );
            break;
        case General_Hernquist:
            pot_temp = d2_dr2_gen_hern_pot(model, r );
            break;
        case InvalidDwarf:
        default:
            mw_fail("Invalid dwarf type\n");
    }

    return pot_temp;
}

real get_second_derv_potential_real(const Dwarf* model_light, const Dwarf* model_dark, real* r, mwbool isLight)
{
    real pot_temp, tmp, rho_light, rho_dark;

    real_0 scale_b = model_light->scaleLength;
    real_0 mass_b = model_light->originmass;
    real_0 xi_scale = model_light->scaleLength / (model_light->scaleLength + model_dark->scaleLength);
    real_0 xi_mass = model_light->originmass / (model_light->originmass + model_dark->originmass);

    real scale_light = mw_real_var(scale_b, BARYON_RADIUS_POS);
    real mass_light = mw_real_var(mass_b, BARYON_MASS_POS);
    real xi_R = mw_real_var(xi_scale, RADIUS_RATIO_POS);
    real xi_M = mw_real_var(xi_mass, MASS_RATIO_POS);

    tmp = inv(&xi_R);
    tmp = mw_add_s(&tmp, -1.0);
    real scale_dark = mw_mul(&tmp, &scale_light);

    tmp = inv(&xi_M);
    tmp = mw_add_s(&tmp, -1.0);
    real mass_dark = mw_mul(&tmp, &mass_light);

    if(isLight)
    {
        switch(model_light->type)
        {
            case Plummer:
                pot_temp = d2_dr2_plummer_pot_real(&mass_light, &scale_light, r);
                break;
            case NFW:
                rho_light = get_p0(&mass_light, &scale_light);
                pot_temp = d2_dr2_nfw_pot_real(&rho_light, &scale_light, r);
                break;
            case General_Hernquist:
                pot_temp = d2_dr2_gen_hern_pot_real(&mass_light, &scale_light, r);
                break;
            case InvalidDwarf:
            default:
                mw_fail("Invalid dwarf type");     
        }
    }
    else
    {
        switch(model_dark->type)
        {
            case Plummer:
                pot_temp = d2_dr2_plummer_pot_real(&mass_dark, &scale_dark, r);
                break;
            case NFW:
                rho_dark = get_p0(&mass_dark, &scale_dark);
                pot_temp = d2_dr2_nfw_pot_real(&rho_dark, &scale_dark, r);
                break;
            case General_Hernquist:
                pot_temp = d2_dr2_gen_hern_pot_real(&mass_dark, &scale_dark, r);
                break;
            case InvalidDwarf:
            default:
                mw_fail("Invalid dwarf type");     
        }
    }
    
    return pot_temp;
}

real_0 get_density(const Dwarf* model, real_0 r)
{
    real_0 den_temp;
    
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


real get_density_real(const Dwarf* model_light, const Dwarf* model_dark, real* r, mwbool isLight)
{
    real den_temp, tmp;
    real rho_light, rho_dark;

    real_0 scale_b = model_light->scaleLength;
    real_0 mass_b = model_light->originmass;
    real_0 xi_scale = model_light->scaleLength / (model_light->scaleLength + model_dark->scaleLength);
    real_0 xi_mass = model_light->originmass / (model_light->originmass + model_dark->originmass);

    real scale_light = mw_real_var(scale_b, BARYON_RADIUS_POS);
    real mass_light = mw_real_var(mass_b, BARYON_MASS_POS);
    real xi_R = mw_real_var(xi_scale, RADIUS_RATIO_POS);
    real xi_M = mw_real_var(xi_mass, MASS_RATIO_POS);

    tmp = inv(&xi_R);
    tmp = mw_add_s(&tmp, -1.0);
    real scale_dark = mw_mul(&tmp, &scale_light);

    tmp = inv(&xi_M);
    tmp = mw_add_s(&tmp, -1.0);
    real mass_dark = mw_mul(&tmp, &mass_light);

    if(isLight)
    {
        switch(model_light->type)
        {
            case Plummer:
                den_temp = plummer_den_real(&mass_light, &scale_light, r);
                break;
            case NFW:
                rho_light = get_p0(&mass_light, &scale_light);
                den_temp = nfw_den_real(&rho_light, &scale_light, r);
                break;
            case General_Hernquist:
                den_temp = gen_hern_den_real(&mass_light, &scale_light, r);
                break;
    //         case Einasto:
    //             einasto_den(model, r);
    //             break;
            case InvalidDwarf:
            default:
                mw_fail("Invalid dwarf type");     
        }
    }
    else
    {
        switch(model_dark->type)
        {
            case Plummer:
                den_temp = plummer_den_real(&mass_dark, &scale_dark, r);
                break;
            case NFW:
                rho_dark = get_p0(&mass_dark, &scale_dark);
                den_temp = nfw_den_real(&rho_dark, &scale_dark, r);
                break;
            case General_Hernquist:
                den_temp = gen_hern_den_real(&mass_dark, &scale_dark, r);
                break;
    //         case Einasto:
    //             einasto_den(model, r);
    //             break;
            case InvalidDwarf:
            default:
                mw_fail("Invalid dwarf type");     
        }
    }
    
    return den_temp;
}

real_0 get_first_derv_density(const Dwarf* model, real_0 r)
{
    real_0 den_temp;
    
    switch(model->type)
    {
        case Plummer:
            den_temp = d_dr_plummer_den(model, r);
            break;
        case NFW:
            den_temp = d_dr_nfw_den(model, r );
            break;
        case General_Hernquist:
            den_temp = d_dr_gen_hern_den(model, r );
            break;
        case InvalidDwarf:
        default:
            mw_fail("Invalid dwarf type");
            
    }
    
    return den_temp;
}

real get_first_derv_density_real(const Dwarf* model_light, const Dwarf* model_dark, real* r, mwbool isLight)
{
    real den_temp, tmp;
    real rho_light, rho_dark;

    real_0 scale_b = model_light->scaleLength;
    real_0 mass_b = model_light->originmass;
    real_0 xi_scale = model_light->scaleLength / (model_light->scaleLength + model_dark->scaleLength);
    real_0 xi_mass = model_light->originmass / (model_light->originmass + model_dark->originmass);

    real scale_light = mw_real_var(scale_b, BARYON_RADIUS_POS);
    real mass_light = mw_real_var(mass_b, BARYON_MASS_POS);
    real xi_R = mw_real_var(xi_scale, RADIUS_RATIO_POS);
    real xi_M = mw_real_var(xi_mass, MASS_RATIO_POS);

    tmp = inv(&xi_R);
    tmp = mw_add_s(&tmp, -1.0);
    real scale_dark = mw_mul(&tmp, &scale_light);

    tmp = inv(&xi_M);
    tmp = mw_add_s(&tmp, -1.0);
    real mass_dark = mw_mul(&tmp, &mass_light);

    if(isLight)
    {
        switch(model_light->type)
        {
            case Plummer:
                den_temp = d_dr_plummer_den_real(&mass_light, &scale_light, r);
                break;
            case NFW:
                rho_light = get_p0(&mass_light, &scale_light);
                den_temp = d_dr_nfw_den_real(&rho_light, &scale_light, r);
                break;
            case General_Hernquist:
                den_temp = d_dr_gen_hern_den_real(&mass_light, &scale_light, r);
                break;
            case InvalidDwarf:
            default:
                mw_fail("Invalid dwarf type");     
        }
    }
    else
    {
        switch(model_dark->type)
        {
            case Plummer:
                den_temp = d_dr_plummer_den_real(&mass_dark, &scale_dark, r);
                break;
            case NFW:
                rho_dark = get_p0(&mass_dark, &scale_dark);
                den_temp = d_dr_nfw_den_real(&rho_dark, &scale_dark, r);
                break;
            case General_Hernquist:
                den_temp = d_dr_gen_hern_den_real(&mass_dark, &scale_dark, r);
                break;
            case InvalidDwarf:
            default:
                mw_fail("Invalid dwarf type");     
        }
    }
    
    return den_temp;
}

real_0 get_second_derv_density(const Dwarf* model, real_0 r)
{
    real_0 den_temp;
    
    switch(model->type)
    {
        case Plummer:
            den_temp = d2_dr2_plummer_den(model, r);
            break;
        case NFW:
            den_temp = d2_dr2_nfw_den(model, r );
            break;
        case General_Hernquist:
            den_temp = d2_dr2_gen_hern_den(model, r );
            break;
        case InvalidDwarf:
        default:
            mw_fail("Invalid dwarf type");
            
    }
    
    return den_temp;
}

real get_second_derv_density_real(const Dwarf* model_light, const Dwarf* model_dark, real* r, mwbool isLight)
{
    real den_temp, tmp;
    real rho_light, rho_dark;

    real_0 scale_b = model_light->scaleLength;
    real_0 mass_b = model_light->originmass;
    real_0 xi_scale = model_light->scaleLength / (model_light->scaleLength + model_dark->scaleLength);
    real_0 xi_mass = model_light->originmass / (model_light->originmass + model_dark->originmass);

    real scale_light = mw_real_var(scale_b, BARYON_RADIUS_POS);
    real mass_light = mw_real_var(mass_b, BARYON_MASS_POS);
    real xi_R = mw_real_var(xi_scale, RADIUS_RATIO_POS);
    real xi_M = mw_real_var(xi_mass, MASS_RATIO_POS);

    tmp = inv(&xi_R);
    tmp = mw_add_s(&tmp, -1.0);
    real scale_dark = mw_mul(&tmp, &scale_light);

    tmp = inv(&xi_M);
    tmp = mw_add_s(&tmp, -1.0);
    real mass_dark = mw_mul(&tmp, &mass_light);

    if(isLight)
    {
        switch(model_light->type)
        {
            case Plummer:
                den_temp = d2_dr2_plummer_den_real(&mass_light, &scale_light, r);
                break;
            case NFW:
                rho_light = get_p0(&mass_light, &scale_light);
                den_temp = d2_dr2_nfw_den_real(&rho_light, &scale_light, r);
                break;
            case General_Hernquist:
                den_temp = d2_dr2_gen_hern_den_real(&mass_light, &scale_light, r);
                break;
            case InvalidDwarf:
            default:
                mw_fail("Invalid dwarf type");     
        }
    }
    else
    {
        switch(model_dark->type)
        {
            case Plummer:
                den_temp = d2_dr2_plummer_den_real(&mass_dark, &scale_dark, r);
                break;
            case NFW:
                rho_dark = get_p0(&mass_dark, &scale_dark);
                den_temp = d2_dr2_nfw_den_real(&rho_dark, &scale_dark, r);
                break;
            case General_Hernquist:
                den_temp = d2_dr2_gen_hern_den_real(&mass_dark, &scale_dark, r);
                break;
            case InvalidDwarf:
            default:
                mw_fail("Invalid dwarf type");     
        }
    }
    
    return den_temp;
}
