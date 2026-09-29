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
#include "nbody_king_model.h"
#include "nbody_mixeddwarf.h"

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
static real plummer_den(const Dwarf* model, real r)                                                                      //
{                                                                                                                        //
    const real mass = model->mass;                                                                                       //
    const real rscale = model->scaleLength;                                                                              //
    return  (3.0 / (4.0 * M_PI)) * (mass / cube(rscale)) * minusfivehalves( (1.0 + sqr(r / rscale)) ) ;                  //
}                                                                                                                        //
                                                                                                                         //
static real plummer_pot(const Dwarf* model, real r)                                                                      //
{                                                                                                                        //
    const real mass = model->mass;                                                                                       //
    const real rscale = model->scaleLength;                                                                              //
    return mass / mw_sqrt(sqr(r) + sqr(rscale));                                                                         //
}                                                                                                                        //
                                                                                                                         //
__attribute__((unused)) static real plummer_vel_disp(const Dwarf* model, real r)                                         //
{                                                                                                                        //
    const real mass = model->mass;                                                                                       //
    const real rscale = model->scaleLength;                                                                              //
    return mass / (6* mw_sqrt(sqr(r)+sqr(rscale)));                                                                      //
}                                                                                                                        //
                                                                                                                         //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*                            NFW                                                                                        */
/* this density is taken from the 1997 paper by nfw. the potential is taken from binney 2nd ed                           */
/* Cutoff for density is addapted from Zemp et al. 2008                                                                  */
static real nfw_den(const Dwarf* model, real r)                                                                          //
{                                                                                                                        //
    const real rscale = model->scaleLength;                                                                              //
    const real p0 = model->p0;                                                                                           //
    const real rcut = model->rcut;                                                                                       //
    real R = r / rscale;                                                                                                 //
#pragma GCC diagnostic push                                                                                              //
#pragma GCC diagnostic ignored "-Wfloat-equal"                                                                           //
    if (rcut != 0.0) {                                                                                                   //
#pragma GCC diagnostic pop                                                                                               //
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
static real nfw_pot(const Dwarf* model, real r)                                                                          //
{                                                                                                                        //
    const real rscale = model->scaleLength;                                                                              //
    const real p0 = model->p0;                                                                                           //
    const real rcut = model->rcut;                                                                                       //
    real R = r / rscale;                                                                                                 //
#pragma GCC diagnostic push                                                                                              //
#pragma GCC diagnostic ignored "-Wfloat-equal"                                                                           //
    if (rcut != 0.0) {                                                                                                   //
#pragma GCC diagnostic pop                                                                                               //
        const real rdecay = model->rdecay;                                                                               //
        const real pcut = model->pcut;                                                                                   //
        const real delta = model->delta;                                                                                 //
        const real m_nfw_cut = model->m_nfw_cut;                                                                         //
        const real const_gamma_func = model->const_gamma_func; // UpperIncompleteGammaFunc(delta + 3, rcut / rdecay)     //
        if (r > rcut) {                                                                                                  //
            return (                                                                                                     //
                4.0 * M_PI * pcut * mw_pow(rcut, -delta) * mw_exp(rcut / rdecay) * mw_pow(rdecay, delta + 3)             //
                * (((const_gamma_func - UpperIncompleteGammaFunc(delta + 3, r / rdecay)) / r)                            //
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
__attribute__((unused)) static real gen_hern_vel_disp(const Dwarf* model, real r)                                        //
{                                                                                                                        //
    const real mass = model->mass;                                                                                       //
    const real rscale = model->scaleLength;                                                                              //
    const real term1 = 12 * r * mw_pow(r + rscale, 3.0) * mw_log((r + rscale) / r) / mw_pow(rscale, 4.0);                //
    const real term2 = 25 + 52 * r / rscale + 42 * sqr(r / rscale) + 12 * mw_pow(r / rscale, 3.0);                       //
    return mass / (12*rscale) * (term1 - r / (r + rscale) * term2);                                                      //
}                                                                                                                        //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*                             EINASTO                                                                                   */
/* Referenced Retana-Montenegro et al. 2012 for the formulas                                                             */
static real einasto_den(const Dwarf* model, real r)                                                                      //
{                                                                                                                        //
    const real M = model->mass;                                                                                          //
    const real rs = model->scaleLength; // half mass radius                                                              //
    const real n = model->n; // Einasto index                                                                            //
    const real h = model->h; // scale length                                                                             //
    const real d = model->d; // dimenstionless constant that depends on n and gurantees that rs is the half mass radius  //
    const real const_gamma_func = model->const_gamma_func;  // GammaFunc(3.0 * n)                                        //
                                                                                                                         //
    real coeff = M / ( 4.0 * M_PI * cube(h) * n * const_gamma_func);                                                     //
    real s = mw_pow(d, n) * r / rs;                                                                                      //
    real s_term = mw_pow(s, inv(n));                                                                                     // 
    return coeff * mw_exp(-s_term);                                                                                      //
}                                                                                                                        //
                                                                                                                         //
static real einasto_pot(const Dwarf* model, real r)                                                                      //
{                                                                                                                        //
    const real M = model->mass;                                                                                          //
    const real rs = model->scaleLength; // half mass radius                                                              //
    const real n = model->n; // Einasto index                                                                            //
    const real h = model->h; // scale length                                                                             //
    const real d = model->d; // dimenstionless constant that depends on n and gurantees that rs is the half mass radius  //
    const real const_gamma_func = model->const_gamma_func;  // GammaFunc(3.0 * n)                                        //
                                                                                                                         //
    real s = mw_pow(d, n) * r / rs;                                                                                      //
    real s_term = mw_pow(s, inv(n));                                                                                     //
    real coeff = M / (h * s);                                                                                            //
                                                                                                                         //
    real term1 = UpperIncompleteGammaFunc(3.0 * n, s_term);                                                              //
    real term2 = s * UpperIncompleteGammaFunc(2.0 * n, s_term);                                                          //
    real term = 1.0 - ( term1 - term2 ) / const_gamma_func;                                                              //
    return coeff * term;                                                                                                 //
}                                                                                                                        //
                                                                                                                         //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*                             CORED                                                                                     */
/* this potential and density are cored NFW profiles to be used with SIDM.                                               */
static real cored_den(const Dwarf* model, real r)                                                                        //
{                                                                                                                        //
    const real r1 = model->r1;                                                                                           //
    const real rcut = model->rcut;                                                                                       //
                                                                                                                         //
#pragma GCC diagnostic push                                                                                              //
#pragma GCC diagnostic ignored "-Wfloat-equal"                                                                           //                                           
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
#pragma GCC diagnostic pop                                                                                               //
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
#pragma GCC diagnostic push                                                                                              //
#pragma GCC diagnostic ignored "-Wfloat-equal"                                                                           //
    if (rcut != 0.0 && r > rcut)                                                                                         //
    {                                                                                                                    //
        const real pcut = model->pcut;                                                                                   //
        const real delta = model->delta;                                                                                 //
        const real rdecay = model->rdecay;                                                                               //
        const real const_gamma_func = model->const_gamma_func; // UpperIncompleteGammaFunc(delta + 3, rcut / rdecay)     //
        return (                                                                                                         //
            4.0 * M_PI * pcut * mw_pow(rcut, -delta) * mw_exp(rcut / rdecay) * mw_pow(rdecay, delta + 3)                 //
            * (((const_gamma_func - UpperIncompleteGammaFunc(delta + 3, r / rdecay)) * inv(r))                           //
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
#pragma GCC diagnostic pop                                                                                               //
}                                                                                                                        //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*                            KING                                                                                       */
/* Model is computed numerically, theory from Galactic Dynamics Binney & Tremaine 2nd ed.                                */
/* (lowered isothermal models sec. 4.3). See nbody_king_model.c for the full function content                            */
                                                                                                                         //
static real king_pot(Dwarf* model, real r)                                                                               //
{                                                                                                                        //
    real r_t = model->r_t;
    real r_0 = model->r_0;
    real Rt = r_t/r_0;
    real potential;

    if (r <= model->r_t) {
        // find the 'index' for r to get neighboring indices for R and W array
        real i_ref = (r/r_0)*((real)model->king_steps)/Rt;
        int i1 = mw_floor(i_ref);
        int i2 = mw_ceil(i_ref);

        // floor and ceil have weird behaviors if the value is an integer
        if (i1 == i2) {
            if (i1 != 0 && i2 != 0) {
                i1 -= 1;
            } else {
                i1 = 0;
                i2 = 1;
            }
        }
        //printf("\n\tfor potential: i1=%d, i2=%d, r=%lf", i1, i2, r);
        //fflush(stdout);
        if (i1 < 0 || i2 < 0) {
            printf("weird pot, i1=%d, i2=%d, r=%lf", i1, i2, r);
        }
        //interpolate between W(R[i1]) and W(R[i2]) to get W(R), solve for potential using this and other params
        real W_R = interpolateLinear(r/r_0, model->king_R[i1], model->king_R[i2], model->king_W[i1], model->king_W[i2]);

        potential = (model->sigma)*(model->sigma)*W_R - model->phi0;
    } else {
        potential = (model->mass)/r;
    }

    //printf("(i, i1, i2)=(%lf, %d, %d): R=%lf, W=%lf, potential=%lf\n", i_ref, i1, i2, r/r_0, W_R, potential);

    return potential;                                                                                                      //
}                                                                                                                        //
                                                                                                                         //
static real king_den(Dwarf* model, real r)                                                                               //
{                                                                                                                        //
    real r_t = model->r_t;
    real r_0 = model->r_0;
    real Rt = r_t/r_0;
    real density;
    if (r < model->r_t) {
        // find the 'index' for r to get neighboring indices for R and W array
        real i_ref = (r/r_0)*((real)model->king_steps)/Rt;
        int i1 = mw_floor(i_ref);
        int i2 = mw_ceil(i_ref);
        
        // floor and ceil have weird behaviors if the value is an integer
        if (i1 == i2) {
            if (i1 != 0 && i2 != 0) {
                i1 -= 1;
            } else {
                i1 = 0;
                i2 = 1;
            }
        }
        if (i1 < 0 || i2 < 0) {
            printf("weird dens, i1=%d, i2=%d", i1, i2);
        }
        //interpolate between W(R[i1]) and W(R[i2]) to get W(R), solve for density using this and rho(W) function
        real W_R = interpolateLinear(r/r_0, model->king_R[i1], model->king_R[i2], model->king_W[i1], model->king_W[i2]);
        real erfTerm = mw_exp(W_R) * mw_erf(mw_sqrt(W_R));

        density = (model->rho1) * (erfTerm - mw_sqrt(4.0*W_R/(M_PI))*(1.0 + (2.0*W_R)/3.0));
    } else {
        density = 0.0;
    }
    //printf("\n\tdensity results: (i, i1, i2)=(%lf, %d, %d): r=%lf, R=%lf, W=%lf, density=%lf", i_ref, i1, i2, r, r/r_0, W_R, density);
    
    return density;                                                                                                     //
}                                                                                                                        //
                                                                                                                         // 
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
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
            if (model->p0 == 0.0) {
#pragma GCC diagnostic pop
                set_model_params(model);
            }
            pot_temp = nfw_pot(model, r );
            break;
        case General_Hernquist:
            pot_temp = gen_hern_pot(model, r );
            break;
        case Einasto:
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
            if (model->h == 0.0) {
#pragma GCC diagnostic pop
                set_model_params(model);
            }
            pot_temp = einasto_pot(model, r);
            break;
        case Cored:
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
            if (model->p0 == 0.0) {
#pragma GCC diagnostic pop
                set_model_params(model);
            }
            pot_temp = cored_pot(model, r);
            break;
        case King:
            #pragma GCC diagnostic push
            #pragma GCC diagnostic ignored "-Wfloat-equal"
            if (model->r_0 == 0.0) {
            #pragma GCC diagnostic pop
                set_model_params(model);
            }

            pot_temp = king_pot(model, r);
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
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
            if (model->p0 == 0.0) {
#pragma GCC diagnostic pop
                set_model_params(model);
            }
            den_temp = nfw_den(model, r );
            break;
        case General_Hernquist:
            den_temp = gen_hern_den(model, r );
            break;
        case Einasto:
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
            if (model->h == 0.0) {
#pragma GCC diagnostic pop
                set_model_params(model);
            }
            den_temp = einasto_den(model, r);
            break;
        case Cored:
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
            if (model->p0 == 0.0) {
#pragma GCC diagnostic pop
                set_model_params(model);
            }
            den_temp = cored_den(model, r);
            break;
        case King:
            #pragma GCC diagnostic push
            #pragma GCC diagnostic ignored "-Wfloat-equal"
            if (model->r_0 == 0.0) {
            #pragma GCC diagnostic pop
                set_model_params(model);
            }
            
            den_temp = king_den(model, r);
            break;
        case InvalidDwarf:
        default:
            mw_fail("Invalid dwarf type, %d\n", model->type);

    }

    return den_temp;
}

real get_vel_disp_radius(const Dwarf* model)
//radii calculated here are for the velocity dispersion approximation within the softening length calculation
{
    real hmr_temp = 1e-4;

    switch(model->type)
    {
        case Plummer:
        // this is the half mass radius for a Plummer profile. The exact value is 1 / mw_sqrt(1 / mw_pow(5, 2 / 3) - 1) * rscale
            hmr_temp = 1.3*model->scaleLength;
            break;
        case NFW:
        // This is (allegedly) the point of maximum density, but it should be sufficient for softening length calculations
            hmr_temp = model->scaleLength;
            break;
        case General_Hernquist:
        // This is the half mass radius for a Hernquist profile.
            hmr_temp = (1 + mw_sqrt(2))*model->scaleLength;
            break;
        case Einasto:
        // This is the half mass radius for an Einasto profile (which is also the scale radius).
            hmr_temp = model->scaleLength;
            break;
        case Cored:
        // This is also the point of maximum density, but it should be sufficient for softening length calculations
            hmr_temp = (model->scaleLength > model->r1) ? model->scaleLength : model->r1;
            break;
        case King:
        // This is the radius at half of the central surface brightness, aka the King/Core radius r0
            #pragma GCC diagnostic push
            #pragma GCC diagnostic ignored "-Wfloat-equal"
            if (model->r_0 == 0.0) {
            #pragma GCC diagnostic pop
                set_model_params(model);
            }
            hmr_temp = model->r_0;
            break;
        case InvalidDwarf:
        default:
            mw_fail("Invalid dwarf type, %d\n", model->type);
    }

    return hmr_temp;
}