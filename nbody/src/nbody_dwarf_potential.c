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
 static real plummer_den(const Dwarf* model, real r)                                                                     //
{                                                                                                                        //
    const real mass = model->mass;                                                                                       //
    const real rscale = model->scaleLength;                                                                              //
    return  (3.0 / (4.0 * M_PI)) * (mass / cube(rscale)) * minusfivehalves( (1.0 + sqr(r / rscale)) ) ;                  //
}                                                                                                                        //
                                                                                                                         //
 static real plummer_pot(const Dwarf* model, real r)                                                                     //
{                                                                                                                        //
    const real mass = model->mass;                                                                                       //
    const real rscale = model->scaleLength;                                                                              //
    return mass / mw_sqrt(sqr(r) + sqr(rscale));                                                                         //
}                                                                                                                        //
                                                                                                                         //
 static real plummer_vel_disp(const Dwarf* model, real r)                                                                //
{                                                                                                                        //
    const real mass = model->mass;                                                                                       //
    const real rscale = model->scaleLength;                                                                              //
    return mass / (6* mw_sqrt(sqr(r)+sqr(rscale)));                                                                      //
}                                                                                                                        //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*                            NFW                                                                                        */
/* this density is taken from the 1997 paper by nfw. the potential is taken from binney 2nd ed                           */
/* Cutoff for density is addapted from Zemp et al. 2008                                                                  */
 static real nfw_den(const Dwarf* model, real r)                                                                         //
{                                                                                                                        //
    const real rscale = model->scaleLength;                                                                              //
    const real p0 = model->p0;                                                                                           //
    const real rcut = model->rcut;                                                                                       //
    real R = r / rscale;                                                                                                 //
    if (rcut != 0.0) {                                                                                                   //
        const real rdecay = model->rdecay;                                                                               //
        const real pcut = model->pcut;                                                                                   //
        const real delta = model->delta;                                                                                 //
        if (r > rcut) {                                                                                                  //
            return pcut * mw_pow(r / rcut, delta) * mw_exp(-(r - rcut) / rdecay);                                        //
        }                                                                                                                //
        else {                                                                                                           //
            return p0 * inv(R) * inv(sqr(1.0 + R));                                                                      //
        }                                                                                                                //
    }                                                                                                                    //
    /* at r = 0 the density goes to inf. however, the sampling is guarded against r = 0 anyway.*/                        //
    return p0 * inv(R) * inv(sqr(1.0 + R));                                                                              //
}                                                                                                                        //
                                                                                                                         //
 static real nfw_pot(const Dwarf* model, real r)                                                                         //
{                                                                                                                        //
    const real rscale = model->scaleLength;                                                                              //
    const real p0 = model->p0;                                                                                           //
    const real rcut = model->rcut;                                                                                       //
    real R = r / rscale;                                                                                                 //
    if (rcut != 0.0) {                                                                                                   //
        const real rdecay = model->rdecay;                                                                               //
        const real pcut = model->pcut;                                                                                   //
        const real delta = model->delta;                                                                                 //
        const real m_nfw_cut = model->m_nfw_cut;                                                                         //
        const real gamma1 = model->gamma1;                                                                               //
        if (r > rcut) {                                                                                                  //
            return (                                                                                                     //
                4.0 * M_PI * pcut * mw_pow(rcut, -delta) * mw_exp(rcut / rdecay) * mw_pow(rdecay, delta + 3)             //
                * (((gamma1 - UpperIncompleteGammaFunc(delta + 3, r / rdecay)) / r)                                      //
                + (UpperIncompleteGammaFunc(delta + 2, r / rdecay) / rdecay)) + m_nfw_cut / r                            //
            );                                                                                                           //
        } else {                                                                                                         //
            const real psi_nfw_cut = model->psi_nfw_cut;                                                                 //
            const real psi_cut_cut = model->psi_cut_cut;                                                                 //
            const real m_nfw_cut = model->m_nfw_cut;                                                                     //
            return (4.0 * M_PI * p0 * cube(rscale) * mw_log(1.0 + R) * inv(r)                                            //
                - psi_nfw_cut + psi_cut_cut + m_nfw_cut / rcut);                                                         //
        }                                                                                                                //
    }                                                                                                                    //
    /* at r = 0 the pot goes to inf. however, the sampling is guarded against r = 0 anyway. */                           //
    return  4.0 * M_PI * sqr(rscale) * p0 * inv(R) * mw_log(1.0 + R);                                                    //
}                                                                                                                        //
                                                                                                                         //
 static real nfw_vel_disp(const Dwarf* model, real r)                                                                    //
{                                                                                                                        //
    printf("WARNING: currently using plummer velocity dispersion for NFW");                                              //
    const real mass = model->mass;                                                                                       //
    const real rscale = model->scaleLength;                                                                              //
    return mass / (6* mw_sqrt(sqr(r)+sqr(rscale)));                                                                      //
}                                                                                                                        //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*                             GENERAL HERNQUIST                                                                         */
/* this potential and density are both taken from the 1990 paper by hernquist                                            */
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
                                                                                                                         //
static real gen_hern_vel_disp(const Dwarf* model, real r)                                                                //
{                                                                                                                        //
    printf("WARNING: currently using plummer velocity dispersion for Hernquist");                                        //
    const real mass = model->mass;                                                                                       //
    const real rscale = model->scaleLength;                                                                              //
    return mass / (6* mw_sqrt(sqr(r)+sqr(rscale)));                                                                      //
}                                                                                                                        //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*                             EINASTO                                                                                   */
/* these are taken from the einasto paper. There are many problems with this, so it is currently unused.                 */
/* Might be fixed now that gamma functions are now working correctly.                                                    */
static real einasto_den(const Dwarf* model, real r)                                                                      //                                                                     
{                                                                                                                        //
    const real mass __attribute__((unused)) = model->mass;                                                               //
    const real h = model->scaleLength;                                                                                   //
    const real n = model->n;                                                                                             //
                                                                                                                         //
    real coeff = 1.0 / ( 4.0 * M_PI * cube(h) * n * GammaFunc(3.0 * n));                                                 //
    real thing = mw_pow(r, inv(n));                                                                                      //
    return coeff * mw_exp(-thing);                                                                                       //
}                                                                                                                        //
                                                                                                                         //
static real einasto_pot(const Dwarf* model, real r)                                                                      //
{                                                                                                                        //
    const real mass = model->mass;                                                                                       //
    const real h = model->scaleLength;                                                                                   //
    const real n = model->n;                                                                                             //
                                                                                                                         //
    real coeff = mass / (h * r);                                                                                         //
    real thing = mw_pow(r, 1.0 / n);                                                                                     //
                                                                                                                         //
    real term1 = UpperIncompleteGammaFunc(3.0 * n, thing);                                                               //
    real term2 = r * UpperIncompleteGammaFunc(2.0 * n, thing);                                                           //
    real term = 1.0 - ( term1 + term2 ) / GammaFunc(3.0 * n);                                                            //
    return coeff * term;                                                                                                 //
}                                                                                                                        //
                                                                                                                         //
static real einasto_vel_disp(const Dwarf* model, real r)                                                                 //
{                                                                                                                        //
    printf("WARNING: currently using plummer velocity dispersion for Einasto");                                          //
    const real mass = model->mass;                                                                                       //
    const real rscale = model->scaleLength;                                                                              //
    return mass / (6* mw_sqrt(sqr(r)+sqr(rscale)));                                                                      //
}                                                                                                                        //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*                             CORED                                                                                     */
/* this potential and density are cored NFW profiles to be used with SIDM.                                               */
static real cored_den(const Dwarf* model, real r)                                                                        //
{                                                                                                                        //
    const real r1 = model->r1;                                                                                           //
    const real rcut = model->rcut;                                                                                       //
                                                                                                                         //
    if (rcut != 0.0 && r > rcut)                                                                                         //
    {                                                                                                                    //
        const real pcut = model->pcut;                                                                                   //
        const real delta = model->delta;                                                                                 //
        const real rdecay = model->rdecay;                                                                               //
        return pcut * mw_pow(r / rcut, delta) * mw_exp(-(r - rcut) / rdecay);                                            //
    }                                                                                                                    //
    else if (r <= r1)                                                                                                    //
    {                                                                                                                    //
        const real p0 = model->p0;                                                                                       //
        const real rc = model->rc;                                                                                       //
        return p0 / (1.0 + sqr(r / rc));                                                                                 //
    }                                                                                                                    //
    else                                                                                                                 //
    {                                                                                                                    //
        const real ps = model->ps;                                                                                       //
        const real rs = model->scaleLength;                                                                              //
        return ps / ((r / rs) * sqr(1.0 + r / rs));                                                                      //
    }                                                                                                                    //
}                                                                                                                        //
                                                                                                                         //
static real cored_pot(const Dwarf* model, real r)                                                                        //
{                                                                                                                        //
    const real r1 = model->r1;                                                                                           //
    const real p0 = model->p0;                                                                                           //
    const real rc = model->rc;                                                                                           //
    const real ps = model->ps;                                                                                           //
    const real rs = model->scaleLength;                                                                                  //
    const real rcut = model->rcut;                                                                                       //
    const real m_iso_r1 = model->m_iso_r1;                                                                               //
    const real m_nfw_r1 = model->m_nfw_r1;                                                                               //
    const real m_nfw_cut = model->m_nfw_cut;                                                                             //
                                                                                                                         //
    if (rcut != 0.0 && r > rcut)                                                                                         //
    {                                                                                                                    //
        const real pcut = model->pcut;                                                                                   //
        const real delta = model->delta;                                                                                 //
        const real rdecay = model->rdecay;                                                                               //
        const real gamma1 = model->gamma1;                                                                               //
        return (                                                                                                         //
            4.0 * M_PI * pcut * mw_pow(rcut, -delta) * mw_exp(rcut / rdecay) * mw_pow(rdecay, delta + 3)                 //
            * (((gamma1 - UpperIncompleteGammaFunc(delta + 3, r / rdecay)) * inv(r))                                     //
            + (UpperIncompleteGammaFunc(delta + 2, r / rdecay) * inv(rdecay)))                                           //
            + ((m_nfw_cut + m_iso_r1 - m_nfw_r1) * inv(r))                                                               //
        );                                                                                                               //
    }                                                                                                                    //
    else if (r <= r1)                                                                                                    //
    {                                                                                                                    //
        const real p0 = model->p0;                                                                                       //
        const real rc = model->rc;                                                                                       //
        const real psi_iso_r1 = model->psi_iso_r1;                                                                       //
        const real psi_nfw_r1 = model->psi_nfw_r1;                                                                       //
        real psi = (                                                                                                     //
            -4.0 * M_PI * p0 * sqr(rc) * ((mw_log(sqr(rc) + sqr(r)) * inv(2.0)) + ((rc * mw_atan(r / rc) * inv(r))))     //
            - psi_iso_r1 + psi_nfw_r1 + ((m_iso_r1 - m_nfw_r1) * inv(r1))                                                //
        );                                                                                                               //
        if (rcut != 0.0) {                                                                                               //
            const real psi_nfw_cut = model->psi_nfw_cut;                                                                 //
            const real psi_cut_cut = model->psi_cut_cut;                                                                 //
            psi += -psi_nfw_cut - ((m_iso_r1 - m_nfw_r1) * inv(rcut))                                                    //
                + psi_cut_cut + ((m_nfw_cut + m_iso_r1 - m_nfw_r1) * inv(rcut));                                         //
        }                                                                                                                //
        return psi;                                                                                                      //
    }                                                                                                                    //
    else                                                                                                                 //
    {                                                                                                                    //
        const real ps = model->ps;                                                                                       //
        const real rs = model->scaleLength;                                                                              //
        real psi = (                                                                                                     //
            4.0 * M_PI * ps * cube(rs) * inv(r) * mw_log(1.0 + r / rs) + ((m_iso_r1 - m_nfw_r1) * inv(r))                //
        );                                                                                                               //
        if (rcut != 0.0) {                                                                                               //
            const real psi_nfw_cut = model->psi_nfw_cut;                                                                 //
            const real psi_cut_cut = model->psi_cut_cut;                                                                 //
            const real m_nfw_cut = model->m_nfw_cut;                                                                     //
            psi += -psi_nfw_cut - ((m_iso_r1 - m_nfw_r1) * inv(rcut))                                                    //
                + psi_cut_cut + ((m_nfw_cut + m_iso_r1 - m_nfw_r1) * inv(rcut));                                         //
        }                                                                                                                //
        return psi;                                                                                                      //
    }                                                                                                                    //
}                                                                                                                        //
                                                                                                                         //
static real cored_vel_disp(const Dwarf* model, real r)                                                                   //
{                                                                                                                        //
    printf("WARNING: currently using plummer velocity dispersion for Cored");                                            //
    const real mass = model->mass;                                                                                       //
    const real rscale = model->scaleLength;                                                                              //
    return mass / (6* mw_sqrt(sqr(r)+sqr(rscale)));                                                                      //
}                                                                                                                        //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

real get_potential(const Dwarf* model, real r)
{
    real pot_temp = 0.0;

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
        case Einasto:
            printf("WARNING: Einsato dwarf currently has problems and should not be used \n");
            pot_temp = einasto_pot(model, r);
            break;
        case Cored:
            pot_temp = cored_pot(model, r);
            break;
        case InvalidDwarf:
        default:
            mw_fail("Invalid dwarf type, %d\n", model->type);
    }

    return pot_temp;
}



real get_density(const Dwarf* model, real r)
{
    real den_temp = 0.0;

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
        case Einasto:
            printf("WARNING: Einsato dwarf currently has problems and should not be used \n");
            den_temp = einasto_den(model, r);
            break;
        case Cored:
            den_temp = cored_den(model, r);
            break;
        case InvalidDwarf:
        default:
            mw_fail("Invalid dwarf type, %d\n", model->type);

    }

    return den_temp;
}

real get_vel_disp(const Dwarf* model) //radii calculated here are for softening length calculation
{
    real vel_disp_temp = 0;
    real r = 0;

    switch(model->type)
    {
        case Plummer:
            r = 1.3*model->scaleLength;
            vel_disp_temp = plummer_vel_disp(model, r);
            break;
        case NFW:
            vel_disp_temp = nfw_vel_disp(model, r );
            break;
        case General_Hernquist:
            r = (1 + mw_sqrt(2))*model->scaleLength;
            vel_disp_temp = gen_hern_vel_disp(model, r );
            break;
        case Einasto:
            printf("WARNING: Einsato dwarf currently has problems and should not be used \n");
            vel_disp_temp=einasto_vel_disp(model, r);
            break;
        case Cored:
            vel_disp_temp = cored_vel_disp(model, r);
            break;
        case InvalidDwarf:
        default:
            mw_fail("Invalid dwarf type, %d\n", model->type);
    }

    return vel_disp_temp;
}
