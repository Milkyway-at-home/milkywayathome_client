/*
 * Copyright (c) 2020-2021 Rensselaer Polytechnic Institute
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

/*
 * This code tests the Poisson Relation between the acceleration and density
 * functions for each component type. NOTE: This test still has problems
 * when the density is really close to zero.
 */

#include "milkyway_util.h"
#include "nbody_potential.h"
#include "nbody_potential_types.h"
#include "nbody_density.h"
#include "nbody_check_params.h"
#include "nbody_show.h"
#include <stdio.h> 
#include <stdlib.h>
#include <time.h>

static const real h = mw_pow(2.0,-12.0);
static const real thresh = 0.1;
static const real mass_thresh = mw_pow(10.0,-6.0);
static const real pi = 3.1415926535;

static const int nPotentials = 10;
static const int nPositions = 40;
static const real min_pos = -10.0;
static const real max_pos = 10.0;

static inline real randomValueReal(real min_val, real max_val)
{
    real x = ((real)rand()/(real)RAND_MAX)*(max_val - min_val) + min_val;
    return x;
}

static inline void generateRandomSphericalParams(real* mass,
                                                 real* scale)
{
    *mass  = randomValueReal(1000.0,100000.0);
    *scale = randomValueReal(0.1,5.0);
}

static inline void generateRandomDiskParams(real* mass,
                                            real* scaleLength,
                                            real* scaleHeight)
{
    *mass        = randomValueReal(1000.0,100000.0);
    *scaleLength = randomValueReal(0.1,5.0);
    *scaleHeight = randomValueReal(0.1,5.0);
}

static inline void generateRandomHaloParams(real* vhalo,
                                            real* scaleLength,
                                            real* flattenX,
                                            real* flattenY,
                                            real* flattenZ,
                                            real* triaxAngle,
                                            real* gamma,
                                            real* lambda,
                                            real* mass,
                                            real* rho0)
{
    *vhalo       = randomValueReal(1.0,200.0);
    *scaleLength = randomValueReal(0.1,30.0);
    *flattenX    = randomValueReal(1.0,1.8);
    *flattenY    = 1.0;
    *flattenZ    = randomValueReal(mw_sqrt(*flattenX * *flattenX / (*flattenX * *flattenX + 1.0)),(1.0/mw_pow(*flattenX,2.71828))+1);
    *triaxAngle  = randomValueReal(0.0,180.0);
    *gamma       = randomValueReal(2.0,5.0);
    *lambda      = randomValueReal(20.0,300.0);
    *mass        = randomValueReal(1000.0,100000.0);
    *rho0        = randomValueReal(5.0,50.0);
}

static inline mwbool checkPoisson(const Potential* pot, mwvector pos)//const pos
{

    mwvector pos_xp, pos_xm, pos_yp, pos_ym, pos_zp, pos_zm;
    mwvector acc_xp, acc_xm, acc_yp, acc_ym, acc_zp, acc_zm;
    real diff;
    real h_r = h * (mw_pow(X(pos),2.0)+mw_pow(Y(pos),2.0)+mw_pow(Z(pos),2.0));
    //real h_r = h * mw_pow((mw_pow(X(pos),2.0)+mw_pow(Y(pos),2.0)+mw_pow(Z(pos),2.0)),0.5);
    //real h_r = h;
    real time = 0.0;
    real flux = 0.0;
    real boxes = 4.0;
    //mw_printf("%s\n",showPotential(pot));
    //mw_printf("%.15f,%.15f,%.15f,%.15f\n",X(pos),Y(pos),Z(pos),mw_pow(mw_pow(X(nbExtAcceleration(pot, pos, time)),2.0)+mw_pow(Y(nbExtAcceleration(pot, pos, time)),2.0)+mw_pow(Z(nbExtAcceleration(pot, pos, time)),2.0),0.5));
    for (real i = -boxes; i < boxes + 0.5; i = i + 1.0) {
        for (real j = -boxes; j < boxes + 0.5; j = j + 1.0) {
            SET_VECTOR(pos_xp, X(pos)+(boxes + 0.5)*h_r, Y(pos)+i*h_r, Z(pos)+j*h_r);
            SET_VECTOR(pos_xm, X(pos)-(boxes + 0.5)*h_r, Y(pos)+i*h_r, Z(pos)+j*h_r);
            SET_VECTOR(pos_yp, X(pos)+i*h_r, Y(pos)+(boxes + 0.5)*h_r, Z(pos)+j*h_r);
            SET_VECTOR(pos_ym, X(pos)+i*h_r, Y(pos)-(boxes + 0.5)*h_r, Z(pos)+j*h_r);
            SET_VECTOR(pos_zp, X(pos)+i*h_r, Y(pos)+j*h_r, Z(pos)+(boxes + 0.5)*h_r);
            SET_VECTOR(pos_zm, X(pos)+i*h_r, Y(pos)+j*h_r, Z(pos)-(boxes + 0.5)*h_r);

            acc_xp = nbExtAcceleration(pot, pos_xp, time);
            acc_xm = nbExtAcceleration(pot, pos_xm, time);
            acc_yp = nbExtAcceleration(pot, pos_yp, time);
            acc_ym = nbExtAcceleration(pot, pos_ym, time);
            acc_zp = nbExtAcceleration(pot, pos_zp, time);
            acc_zm = nbExtAcceleration(pot, pos_zm, time);
        
            flux = flux + (X(acc_xp) - X(acc_xm) + Y(acc_yp) - Y(acc_ym) + Z(acc_zp) - Z(acc_zm))*h_r*h_r;
        }
    }
    mwvector pos_den_xyz;
    real M_enc = 0.0;
    for (real i = -boxes; i < boxes + 0.5; i = i + 1.0) {
        for (real j = -boxes; j < boxes + 0.5; j = j + 1.0) {
            for (real k = -boxes; k < boxes + 0.5; k = k + 1.0) {
                SET_VECTOR(pos_den_xyz, X(pos)+i*h_r, Y(pos)+j*h_r, Z(pos)+k*h_r);
                real density = nbExtDensity(pot, pos_den_xyz, time);
                M_enc = M_enc + density*h_r*h_r*h_r;
            }
        }
    }

    if (M_enc <= mass_thresh)
    {
        return 2;       /*For this case, we want poisson test to look at different random point*/
    }
    else
    {
        diff = mw_abs(1.0 + 4.0*pi*M_enc/flux);
    }
    //mw_printf("%.15f,", diff);
    if (diff > thresh)
    {
        mw_printf("|1 + 4*pi*M_enc/flux| = |1 + %.15f/%.15f| = %.15f\n", 4.0*pi*M_enc, flux, diff);
        mw_printf("  Potential = %s\n", showPotential(pot));
        return 1;
    }

    return 0;

}

static inline void createNullPotential(Potential* p)
{
    p->sphere[0].type = NoSpherical;
    p->disk.type = NoDisk;
    p->disk2.type = NoDisk;
    p->halo.type = NoHalo;
}

static inline mwbool testSphericalPotential(spherical_t t)
{
    Potential p = EMPTY_POTENTIAL;
    real mass, scale;
    int j;
    unsigned int poisson_check;
    mwbool pot_failed = FALSE;
    mwvector pos;
    createNullPotential(&p);

    generateRandomSphericalParams(&mass, &scale);
    p.sphere[0].type = t;
    p.sphere[0].mass = mass;
    p.sphere[0].scale = scale;
    if(checkSphericalConstants(&(p.sphere[0])))
    {
        mw_printf("Bad Spherical Constants in %s\n", showSphericalT(t));
        return TRUE;
    }
    for (j = 0; j < nPositions; j++)
    {
        poisson_check = 2;
        while (poisson_check == 2)
        {
            SET_VECTOR(pos, randomValueReal(min_pos,max_pos), randomValueReal(min_pos,max_pos), randomValueReal(min_pos,max_pos));
            poisson_check = checkPoisson(&p,pos);
        }
        if(poisson_check)
        {
            pot_failed = TRUE;
            mw_printf("    %s failed Poisson Test at pos %d = [%.15f,%.15f,%.15f]!\n", showSphericalT(t), j, X(pos), Y(pos), Z(pos));
            break;
        }
    }
    return pot_failed;
}

static inline mwbool testDiskPotential(disk_t t)
{
    Potential p = EMPTY_POTENTIAL;
    real mass, scaleLength, scaleHeight;
    int j;
    unsigned int poisson_check;
    mwbool pot_failed = FALSE;
    mwvector pos;
    createNullPotential(&p);

    generateRandomDiskParams(&mass, &scaleLength, &scaleHeight);
    p.disk.type = t;
    p.disk.mass = mass;
    p.disk.scaleLength = scaleLength;
    p.disk.scaleHeight = scaleHeight;
    if(p.disk.type == OrbitingBar){
        p.disk.startAngle = 0;
        p.disk.patternSpeed = 250;
    }
    if(checkDiskConstants(&(p.disk)))
    {
        mw_printf("Bad Disk Constants in %s\n", showDiskT(t));
        return TRUE;
    }
    for (j = 0; j < nPositions; j++)
    {
        poisson_check = 2;
        while (poisson_check == 2)
        {
            SET_VECTOR(pos, randomValueReal(min_pos,max_pos), randomValueReal(min_pos,max_pos), randomValueReal(min_pos,max_pos));
            poisson_check = checkPoisson(&p,pos);
        }
        if(poisson_check)
        {
            pot_failed = TRUE;
            mw_printf("    %s failed Poisson Test at pos %d = [%.15f,%.15f,%.15f]!\n", showDiskT(t), j, X(pos), Y(pos), Z(pos));
            break;
        }
    }
    return pot_failed;
}

static inline mwbool testBarPotential(disk_t t)
{
    Potential p = EMPTY_POTENTIAL;
    real mass, scaleLength, scaleHeight;
    int j;
    unsigned int poisson_check;
    mwbool pot_failed = FALSE;
    mwvector pos;
    createNullPotential(&p);

    generateRandomDiskParams(&mass, &scaleLength, &scaleHeight);
    p.disk2.type = t;
    p.disk2.mass = mass;
    p.disk2.scaleLength = scaleLength;
    p.disk2.scaleHeight = scaleHeight;
    if(p.disk2.type == OrbitingBar){
        p.disk2.startAngle = 0;
        p.disk2.patternSpeed = 250;
    }
    if(checkDiskConstants(&(p.disk2)))
    {
        mw_printf("Bad Disk Constants in %s\n", showDiskT(t));
        return TRUE;
    }
    for (j = 0; j < nPositions; j++)
    {
        poisson_check = 2;
        while (poisson_check == 2)
        {
            SET_VECTOR(pos, randomValueReal(-4,4), randomValueReal(-0.2,0.2), randomValueReal(-0.2,0.02));
            while(mw_abs(X(pos)) < 0.01){
                pos.x = randomValueReal(-4,4);
            }
            poisson_check = checkPoisson(&p,pos);
        }
        if(poisson_check)
        {
            pot_failed = TRUE;
            mw_printf("    %s failed Poisson Test at pos %d = [%.15f,%.15f,%.15f]!\n", showDiskT(t), j, X(pos), Y(pos), Z(pos));
            break;
        }
    }
    return pot_failed;
}

static inline mwbool testHaloPotential(halo_t t)
{
    Potential p = EMPTY_POTENTIAL;
    real vhalo, scaleLength, flattenX, flattenY, flattenZ, triaxAngle, gamma, lambda, mass, rho0, mag;
    int j;
    unsigned int poisson_check;
    mwbool pot_failed = FALSE;
    mwvector pos;
    createNullPotential(&p);

    generateRandomHaloParams(&vhalo, &scaleLength, &flattenX, &flattenY, &flattenZ, &triaxAngle, &gamma, &lambda, &mass, &rho0);
    p.halo.type = t;
    p.halo.vhalo = vhalo;
    p.halo.scaleLength = scaleLength;
    p.halo.flattenX = flattenX;
    p.halo.flattenY = flattenY;
    p.halo.flattenZ = flattenZ;
    p.halo.triaxAngle = triaxAngle;
    p.halo.gamma = gamma;
    p.halo.lambda = lambda;
    p.halo.mass = mass;
    p.halo.rho0 = rho0;
    if(checkHaloConstants(&(p.halo)))
    {
        mw_printf("Bad Halo Constants in %s\n", showHaloT(t));
        return TRUE;
    }
    for (j = 0; j < nPositions; j++)
    {
        poisson_check = 2;
        while (poisson_check == 2)
        {
            SET_VECTOR(pos, randomValueReal(min_pos,max_pos), randomValueReal(min_pos,max_pos), randomValueReal(min_pos,max_pos));
            mag = mw_pow(X(pos)*X(pos) + Y(pos)*Y(pos) + Z(pos)*Z(pos), 0.5);
            while (mag > lambda)  /** This is to prevent calculations at radii where the density is set to zero. **/
            {
                mw_incmulvs(pos,0.5);
                mag = mw_pow(X(pos)*X(pos) + Y(pos)*Y(pos) + Z(pos)*Z(pos), 0.5);
            }
            poisson_check = checkPoisson(&p,pos);
        }
        if(poisson_check)
        {
            pot_failed = TRUE;
            mw_printf("    %s failed Poisson Test at pos %d = [%.15f,%.15f,%.15f]!\n", showHaloT(t), j, X(pos), Y(pos), Z(pos));
            break;
        }
    }
    return pot_failed;
}

int main()
{
    srand(time(0));

    int i;
    int failed = 0;

    //TEST SPHERICAL COMPONENTS
    //Hernquist Spherical Potential
    
    for (i = 0; i < nPotentials; i++)
    {
        if (testSphericalPotential(HernquistSpherical))
        {
            failed = 1;
            break;
        }
    }

    //Plummer Spherical Potential  
    for (i = 0; i < nPotentials; i++)
    {
        if (testSphericalPotential(PlummerSpherical))
        {
            failed = 1;
            break;
        }
    }

    //TEST DISK COMPONENTS
//    //Freeman Disk Potential                  /** This potential is an infinitely thin disk, so it's difficult to test this one. May need more terms at lower z values. **/
//    for (i = 0; i < nPotentials; i++)
//    {
//        if (testDiskPotential(FreemanDisk))
//        {
//            failed = 1;
//            break;
//        }
//    }

    //Miyamoto-Nagai Disk Potential  
    for (i = 0; i < nPotentials; i++)
    {
        if (testDiskPotential(MiyamotoNagaiDisk))
        {
            failed = 1;
            break;
        }
    }
    
    //Orbiting Bar Potential
    for (i = 0; i < nPotentials; i++)
    {
        if (testBarPotential(OrbitingBar))
        {
            failed = 1;
            break;
        }
    }

    //Double Exponential Disk Potential
    for (i = 0; i < nPotentials; i++)
    {
        if (testDiskPotential(DoubleExponentialDisk))
        {
            failed = 1;
            break;
        }
    }
    
    //Hyperbolic Exponential Disk Potential  /** FIXME: This potential fails the poisson test. Do not use this potential until we can more accurately calculate it **/
//    for (i = 0; i < nPotentials; i++)
//    {
//        if (testDiskPotential(Sech2ExponentialDisk))
//        {
//            failed = 1;
//            break;
//        }
//    }

    //TEST HALO COMPONENTS
    //Logarithmic Halo Potential  
    for (i = 0; i < nPotentials; i++)
    {
        if (testHaloPotential(LogarithmicHalo))
        {
            failed = 1;
            break;
        }
    }

    //NFW Halo Potential
    for (i = 0; i < nPotentials; i++)
    {
        if (testHaloPotential(NFWHalo))
        {
            failed = 1;
            break;
        }
    }

    //Triaxial Halo Potential
    for (i = 0; i < nPotentials; i++)
    {
        if (testHaloPotential(TriaxialHalo))
        {
            failed = 1;
            break;
        }
    }

    //Allen-Santillan Halo Potential
    for (i = 0; i < nPotentials; i++)
    {
        if (testHaloPotential(AllenSantillanHalo))
        {
            failed = 1;
            break;
        }
    }

    //Wilkinson-Evans Halo Potential
    for (i = 0; i < nPotentials; i++)
    {
        if (testHaloPotential(WilkinsonEvansHalo))
        {
            failed = 1;
            break;
        }
    }

    //NFW (MASS) Halo Potential
    for (i = 0; i < nPotentials; i++)
    {
        if (testHaloPotential(NFWMassHalo))
        {
            failed = 1;
            break;
        }
    }

    //Plummer Halo Potential
    for (i = 0; i < nPotentials; i++)
    {
        if (testHaloPotential(PlummerHalo))
        {
            failed = 1;
            break;
        }
    }

    //Hernquist Halo Potential
    for (i = 0; i < nPotentials; i++)
    {
        if (testHaloPotential(HernquistHalo))
        {
            failed = 1;
            break;
        }
    }

    //Ninkovic Halo Potential
    for (i = 0; i < nPotentials; i++)
    {
        if (testHaloPotential(NinkovicHalo))
        {
            failed = 1;
            break;
        }
    }

    return failed;

}
