/* Copyright 2010 Matthew Arsenault, Travis Desell, Boleslaw
Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail and
Rensselaer Polytechnic Institute.

Copyright (c) 2016-2018 Siddhartha Shelton

This file is part of Milkway@Home.

Milkyway@Home is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Milkyway@Home is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Milkyway@Home.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "nbody_priv.h"
#include "nbody_orbit_integrator.h"
#include "nbody_potential.h"
#include "nbody_io.h"
#include "nbody_coordinates.h"
#include "nbody_defaults.h"
#include "milkyway_util.h"
#include "nbody.h"
#include "nbody_friction.h"
/* Simple orbit integrator in user-defined potential
    Written for BOINC Nbody
    willeb 10 May 2010 */
/* Altered to be consistent with nbody integrator.
 * shelton June 25 2018 */

mwvector* shiftByLMC = NULL; //Ptr to LMC Shift Array (default is NULL)
size_t nShiftLMC = 0;

mwvector LMCpos = ZERO_VECTOR; //Ptr to LMC position (default is NULL)
mwvector LMCvel = ZERO_VECTOR; //Ptr to LMC velocity (default is NULL)

void nbReverseOrbit(mwvector* finalPos,
                    mwvector* finalVel,
                    const Potential* pot,
                    mwvector* pos,
                    mwvector* vel,
                    real_0 tstop,
                    real_0 dt,
                    real_0 sun_dist)
{
    mw_printf("Performing Reverse Orbit Calculation...\n");
    mwvector v_var, x_var, x_lbr, tmp;
    mwvector acc, v, x;
    real_0 t;
    real_0 dt_half = dt / 2.0;
    int initialLArrayIndex = tstop/dt;

    // Set derivative information for AUTODIFF
    //mw_printf("OLD X = [%.15f, %.15f, %.15f]\n", showRealValue(&pos->x), showRealValue(&pos->y), showRealValue(&pos->z));
    x_lbr = cartesianToLbr(pos, sun_dist);
    B(&x_lbr) = mw_real_var(showRealValue(&B(&x_lbr)), B_COORD_POS);
    R(&x_lbr) = mw_real_var(showRealValue(&R(&x_lbr)), R_COORD_POS);
    x_var = lbrToCartesian(&x_lbr, sun_dist);
    //mw_printf("NEW X = [%.15f, %.15f, %.15f]\n", showRealValue(&x_var.x), showRealValue(&x_var.y), showRealValue(&x_var.z));

    //mw_printf("OLD V = [%.15f, %.15f, %.15f]\n", showRealValue(&vel->x), showRealValue(&vel->y), showRealValue(&vel->z));
    v_var.x = mw_real_var(showRealValue(&X(vel)), VX_COORD_POS);
    v_var.y = mw_real_var(showRealValue(&Y(vel)), VY_COORD_POS);
    v_var.z = mw_real_var(showRealValue(&Z(vel)), VZ_COORD_POS);
    //mw_printf("NEW V = [%.15f, %.15f, %.15f]\n", showRealValue(&v_var.x), showRealValue(&v_var.y), showRealValue(&v_var.z));

    // Set the initial conditions
    x = x_var;
    v = mw_negv(&v_var);

    // Get the initial acceleration
    acc = nbExtAcceleration(pot, &x, 0);
    //do this loop backward in order to get an accurate time for time-dependent potentials
    for (t = 0; t >= tstop*(-1); t -= dt)
    {
        //mw_printf("POS  = [%.15f, %.15f, %.15f]\n", showRealValue(&x.x), showRealValue(&x.y), showRealValue(&x.z));
        //mw_printf("VEL  = [%.15f, %.15f, %.15f]\n", showRealValue(&v.x), showRealValue(&v.y), showRealValue(&v.z));
        // Update the velocities and positions
        tmp.x = mw_mul_s(&X(&acc), dt_half);
        tmp.y = mw_mul_s(&Y(&acc), dt_half);
        tmp.z = mw_mul_s(&Z(&acc), dt_half);
        v = mw_addv(&v, &tmp);

        tmp.x = mw_mul_s(&X(&v), dt);
        tmp.y = mw_mul_s(&Y(&v), dt);
        tmp.z = mw_mul_s(&Z(&v), dt);
        x = mw_addv(&x, &tmp); 

        
        // Compute the new acceleration
        acc = nbExtAcceleration(pot, &x, t);
        
        tmp.x = mw_mul_s(&X(&acc), dt_half);
        tmp.y = mw_mul_s(&Y(&acc), dt_half);
        tmp.z = mw_mul_s(&Z(&acc), dt_half);
        v = mw_addv(&v, &tmp);
    }
    
    /* Report the final values (don't forget to reverse the velocities) */
    v = mw_negv(&v);
    //mw_printf("FINAL POS  = [%.15f, %.15f, %.15f]\n", showRealValue(&x.x), showRealValue(&x.y), showRealValue(&x.z));
    //mw_printf("FINAL VEL  = [%.15f, %.15f, %.15f]\n", showRealValue(&v.x), showRealValue(&v.y), showRealValue(&v.z));

    //FIXME: NEEDS EVOLUTION TIME DERIVATIVE INFORMATION!
    *finalPos = x;
    *finalVel = v;

    mw_printf("Initial Dwarf Position = [%.15f, %.15f, %.15f]\n", showRealValue(&x.x), showRealValue(&x.y), showRealValue(&x.z));
    mw_printf("Initial Dwarf Velocity = [%.15f, %.15f, %.15f]\n", showRealValue(&v.x), showRealValue(&v.y), showRealValue(&v.z));
}

void nbReverseOrbit_LMC(mwvector* finalPos,
                    mwvector* finalVel,
                    mwvector* LMCfinalPos,
                    mwvector* LMCfinalVel,
                    const Potential* pot,
                    mwvector* pos,
                    mwvector* vel,
                    mwvector* LMCposition,
                    mwvector* LMCvelocity,
                    mwbool LMCDynaFric,
                    real_0 ftime,
                    real_0 tstop,
                    real_0 dt,
                    real* LMCmass,
                    real* LMCscale,
                    real_0 sun_dist)
{
    mw_printf("Performing Reverse Orbit Calculation with LMC...\n");	
    unsigned int steps = mw_ceil_0((tstop)/(dt)) + 1;
    unsigned int exSteps = mw_abs_0(mw_ceil_0((ftime-tstop)/(dt)) + 1);
    unsigned int i = 0, j = 0, k = 0;
    real_0 dt_half = dt / 2.0;
    mwvector v_var, x_var, x_lbr;
    mwvector acc, v, x, mw_acc, LMC_acc, DF_acc, LMCv, LMCx, tmp, dv2, dx, dLMCv2, dLMCx;
    mwvector mw_x = ZERO_VECTOR;
    mwvector* bacArray = NULL;
    mwvector* forArray = NULL;

    unsigned int vectorSize = sizeof(mwvector);

    //Placeholder arrays for LMC acceleration corrections
    bacArray = (mwvector*)mwCallocA(steps + 1, vectorSize);
    forArray = (mwvector*)mwCallocA(exSteps + 1, vectorSize);

    real_0 t;

    // Set derivative information for AUTODIFF
    x_lbr = cartesianToLbr(pos, sun_dist);
    B(&x_lbr) = mw_real_var(showRealValue(&B(&x_lbr)), B_COORD_POS);
    R(&x_lbr) = mw_real_var(showRealValue(&R(&x_lbr)), R_COORD_POS);
    x_var = lbrToCartesian(&x_lbr, sun_dist);

    v_var.x = mw_real_var(showRealValue(&X(vel)), VX_COORD_POS);
    v_var.y = mw_real_var(showRealValue(&Y(vel)), VY_COORD_POS);
    v_var.z = mw_real_var(showRealValue(&Z(vel)), VZ_COORD_POS);

    // Check if forward time is larger than backward time. We will need to manually compute additional LMC accelerations in that case.
    if (ftime > tstop) {

        // Set the initial conditions for forward orbit
        mw_printf("    Calculating forward orbit...\n");
        x = x_var;
        v = v_var;
        LMCv = *LMCvelocity;
        LMCx = *LMCposition;

        // Get the initial acceleration
        mw_acc = plummerAccel(&mw_x, &LMCx, LMCmass, LMCscale);
        DF_acc = dynamicalFriction_LMC(pot, &LMCx, &LMCv, LMCmass, LMCscale, LMCDynaFric, 0.0);
        LMC_acc = nbExtAcceleration(pot, &LMCx, 0.0);
        LMC_acc = mw_addv(&LMC_acc, &DF_acc);
        acc = nbExtAcceleration(pot, &x, 0.0);
        tmp = plummerAccel(&x, &LMCx, LMCmass, LMCscale);
        acc = mw_addv(&acc, &tmp);

        // Shift the body
        LMC_acc = mw_subv(&LMC_acc, &mw_acc);
        acc = mw_subv(&acc, &mw_acc);

        for (t = 0; t <= (ftime-tstop); t += dt)
        {   
    	    exSteps = (int) mw_round_0(t/dt);
    	    if ((exSteps % 10 == 0)&&(t!=0)) { 
    	        forArray[k] = mw_negv(&mw_acc);
                k++;
    	    }

            // Update the velocities and positions
            dv2.x = mw_mul_s(&acc.x, dt_half);
            dv2.y = mw_mul_s(&acc.y, dt_half);
            dv2.z = mw_mul_s(&acc.z, dt_half);
            v = mw_addv(&v, &dv2);

            dx.x = mw_mul_s(&v.x, dt);
            dx.y = mw_mul_s(&v.y, dt);
            dx.z = mw_mul_s(&v.z, dt);
            x = mw_addv(&x, &dx);

            dLMCv2.x = mw_mul_s(&LMC_acc.x, dt_half);
            dLMCv2.y = mw_mul_s(&LMC_acc.y, dt_half);
            dLMCv2.z = mw_mul_s(&LMC_acc.z, dt_half);
            LMCv = mw_addv(&LMCv, &dLMCv2);

            dLMCx.x = mw_mul_s(&LMCv.x, dt);
            dLMCx.y = mw_mul_s(&LMCv.y, dt);
            dLMCx.z = mw_mul_s(&LMCv.z, dt);
            LMCx = mw_addv(&LMCx, &dLMCx);
        
            // Compute the new acceleration
            mw_acc = plummerAccel(&mw_x, &LMCx, LMCmass, LMCscale);
            DF_acc = dynamicalFriction_LMC(pot, &LMCx, &LMCv, LMCmass, LMCscale, LMCDynaFric, t);
            LMC_acc = nbExtAcceleration(pot, &LMCx, t);
            LMC_acc = mw_addv(&LMC_acc, &DF_acc);
            acc = nbExtAcceleration(pot, &x, t);
            tmp = plummerAccel(&x, &LMCx, LMCmass, LMCscale);
    	    acc = mw_addv(&acc, &tmp);

    	    // Shift the body
            LMC_acc = mw_subv(&LMC_acc, &mw_acc);
            acc = mw_subv(&acc, &mw_acc);
        
            dv2.x = mw_mul_s(&acc.x, dt_half);
            dv2.y = mw_mul_s(&acc.y, dt_half);
            dv2.z = mw_mul_s(&acc.z, dt_half);
            v = mw_addv(&v, &dv2);

            dLMCv2.x = mw_mul_s(&LMC_acc.x, dt_half);
            dLMCv2.y = mw_mul_s(&LMC_acc.y, dt_half);
            dLMCv2.z = mw_mul_s(&LMC_acc.z, dt_half);
            LMCv = mw_addv(&LMCv, &dLMCv2);

        }
        forArray[k] = mw_negv(&mw_acc); //set the last index after the loop ends
    }

    // Set the initial conditions for reverse orbit
    mw_printf("    Calculating backward orbit...\n");
    x = x_var;
    v = mw_negv(&v_var);
    LMCv = mw_negv(LMCvelocity);
    LMCx = *LMCposition;


    // Get the initial acceleration
    LMC_acc = nbExtAcceleration(pot, &LMCx, 0);
    if (LMCDynaFric) {
        DF_acc = dynamicalFriction_LMC(pot, &LMCx, &LMCv, LMCmass, LMCscale, TRUE, 0);
        LMC_acc = mw_subv(&LMC_acc, &DF_acc); /* Inverting drag force for reverse orbit */
     }
    acc = nbExtAcceleration(pot, &x, 0);
    tmp = plummerAccel(&x, &LMCx, LMCmass, LMCscale);
    acc = mw_addv(&acc, &tmp);

    // Shift the body
    mw_acc = plummerAccel(&mw_x, &LMCx, LMCmass, LMCscale);
    LMC_acc = mw_subv(&LMC_acc, &mw_acc);
    acc = mw_subv(&acc, &mw_acc);

    real_0 negT = 0;
    for (t = 0; t <= tstop; t += dt)
    {   
        //negate this time for use in time-dependent potentials
        negT = t*-1;
    	steps = t/dt;
    	if( steps % 10 == 0){ 
    		bacArray[i] = mw_negv(&mw_acc);
        	i++;
    	}

        // Update the velocities and positions
        dv2.x = mw_mul_s(&acc.x, dt_half);
        dv2.y = mw_mul_s(&acc.y, dt_half);
        dv2.z = mw_mul_s(&acc.z, dt_half);
        v = mw_addv(&v, &dv2);

        dx.x = mw_mul_s(&v.x, dt);
        dx.y = mw_mul_s(&v.y, dt);
        dx.z = mw_mul_s(&v.z, dt);
        x = mw_addv(&x, &dx);

        dLMCv2.x = mw_mul_s(&LMC_acc.x, dt_half);
        dLMCv2.y = mw_mul_s(&LMC_acc.y, dt_half);
        dLMCv2.z = mw_mul_s(&LMC_acc.z, dt_half);
        LMCv = mw_addv(&LMCv, &dLMCv2);

        dLMCx.x = mw_mul_s(&LMCv.x, dt);
        dLMCx.y = mw_mul_s(&LMCv.y, dt);
        dLMCx.z = mw_mul_s(&LMCv.z, dt);
        LMCx = mw_addv(&LMCx, &dLMCx);
        
        // Compute the new acceleration
        LMC_acc = nbExtAcceleration(pot, &LMCx, negT);
        if (LMCDynaFric) {
            DF_acc = dynamicalFriction_LMC(pot, &LMCx, &LMCv, LMCmass, LMCscale, TRUE, negT);
            LMC_acc = mw_subv(&LMC_acc, &DF_acc); /* Inverting drag force for reverse orbit */
        }
        acc = nbExtAcceleration(pot, &x, negT);
        tmp = plummerAccel(&x, &LMCx, LMCmass, LMCscale);
        acc = mw_addv(&acc, &tmp);

    	// Shift the body
        mw_acc = plummerAccel(&mw_x, &LMCx, LMCmass, LMCscale);
        LMC_acc = mw_subv(&LMC_acc, &mw_acc);
        acc = mw_subv(&acc, &mw_acc);
        
        dv2.x = mw_mul_s(&acc.x, dt_half);
        dv2.y = mw_mul_s(&acc.y, dt_half);
        dv2.z = mw_mul_s(&acc.z, dt_half);
        v = mw_addv(&v, &dv2);

        dLMCv2.x = mw_mul_s(&LMC_acc.x, dt_half);
        dLMCv2.y = mw_mul_s(&LMC_acc.y, dt_half);
        dLMCv2.z = mw_mul_s(&LMC_acc.z, dt_half);
        LMCv = mw_addv(&LMCv, &dLMCv2);

    }
    bacArray[i] = mw_negv(&mw_acc); //set the last index after the loop ends
    
    //Allocate memory for the shift array equal to (x,y,z) i times with extra wiggle room dependent on evolve time
    unsigned int size = i + k + 2;
    shiftByLMC = (mwvector*)mwCallocA(size, vectorSize); 

    //Fill reverse orbit of shift array
    for(j = 0; j < i+1; j++) {
        tmp = bacArray[i-j];
        shiftByLMC[j] = tmp;
    }

    //Fill forward orbit of shift array
    if (ftime > tstop) {
        for(j = 0; j < k+1; j++) {
            tmp = forArray[j];
            shiftByLMC[i+1+j] = tmp;
        }
    }

    //Free placeholder arrays
    mwFreeA(bacArray);
    mwFreeA(forArray);

    nShiftLMC = size;

    /* Report the final values (don't forget to reverse the velocities) */
    v = mw_negv(&v);
    LMCv = mw_negv(&LMCv);

    //FIXME: NEEDS EVOLUTION TIME DERIVATIVE INFORMATION!
    *finalPos = x;
    *finalVel = v;
    *LMCfinalPos = LMCx;
    *LMCfinalVel = LMCv;

    mw_printf("Initial Dwarf Position = [%.15f, %.15f, %.15f]\n", showRealValue(&x.x), showRealValue(&x.y), showRealValue(&x.z));
    mw_printf("Initial Dwarf Velocity = [%.15f, %.15f, %.15f]\n", showRealValue(&v.x), showRealValue(&v.y), showRealValue(&v.z));
    mw_printf("Initial LMC Position = [%.15f, %.15f, %.15f]\n", showRealValue(&LMCx.x), showRealValue(&LMCx.y), showRealValue(&LMCx.z));
    mw_printf("Initial LMC Velocity = [%.15f, %.15f, %.15f]\n", showRealValue(&LMCv.x), showRealValue(&LMCv.y), showRealValue(&LMCv.z));

    //Store LMC position and velocity
    LMCpos = LMCx;
    LMCvel = LMCv;
    mw_printf("Reverse Orbit Calculated!\n");
}

void getLMCArray(mwvector ** shiftArrayPtr, size_t * shiftSizePtr) {
    //Allows access to shift array
    *shiftArrayPtr = shiftByLMC;
    *shiftSizePtr = nShiftLMC;
}

void getLMCPosVel(mwvector * LMCposPtr, mwvector * LMCvelPtr) {
    //Allows access to LMC position and velocity
    *LMCposPtr = LMCpos;
    *LMCvelPtr = LMCvel;
}


void nbPrintReverseOrbit(mwvector* finalPos,
                         mwvector* finalVel,
                         const Potential* pot,
                         mwvector* pos,
                         mwvector* vel,
                         real_0 tstop,
                         real_0 tstopforward,
                         real_0 dt)
{
    mwvector acc, v, x;
    mwvector v_for, x_for, tmp;
    mwvector lbr;
    real_0 t;
    real_0 dt_half = dt / 2.0;

    // Set the initial conditions
    x = *pos;
    v = mw_negv(vel);
    x_for = *pos;
    v_for = *vel;

    // Get the initial acceleration
    acc = nbExtAcceleration(pot, &x, 0);

    FILE * fp;
    fp = fopen("reverse_orbit.out", "w");
    // Loop through time
    for (t = 0; t >= tstop*(-1); t -= dt)
    {
        // Update the velocities and positions
        tmp.x = mw_mul_s(&X(&acc), dt_half);
        tmp.y = mw_mul_s(&Y(&acc), dt_half);
        tmp.z = mw_mul_s(&Z(&acc), dt_half);
        v = mw_addv(&v, &tmp);

        tmp.x = mw_mul_s(&X(&v), dt);
        tmp.y = mw_mul_s(&Y(&v), dt);
        tmp.z = mw_mul_s(&Z(&v), dt);
        x = mw_addv(&x, &tmp);
        
        
        // Compute the new acceleration
        acc = nbExtAcceleration(pot, &x, t);
        tmp.x = mw_mul_s(&X(&acc), dt_half);
        tmp.y = mw_mul_s(&Y(&acc), dt_half);
        tmp.z = mw_mul_s(&Z(&acc), dt_half);
        v = mw_addv(&v, &tmp);
        
        lbr = cartesianToLbr(&x, DEFAULT_SUN_GC_DISTANCE);
        fprintf(fp, "%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\n", showRealValue(&X(&x)), showRealValue(&Y(&x)), showRealValue(&Z(&x)), showRealValue(&X(&lbr)), showRealValue(&Y(&lbr)), showRealValue(&Z(&lbr)), showRealValue(&X(&v)), showRealValue(&Y(&v)), showRealValue(&Z(&v)));
    }

    fclose(fp);
    fp = fopen("forward_orbit.out", "w");
    acc = nbExtAcceleration(pot, &x_for, 0);
    for (t = 0; t <= tstopforward; t += dt)
    {
        // Update the velocities and positions
        tmp.x = mw_mul_s(&X(&acc), dt_half);
        tmp.y = mw_mul_s(&Y(&acc), dt_half);
        tmp.z = mw_mul_s(&Z(&acc), dt_half);
        v_for = mw_addv(&v_for, &tmp);

        tmp.x = mw_mul_s(&X(&v_for), dt);
        tmp.y = mw_mul_s(&Y(&v_for), dt);
        tmp.z = mw_mul_s(&Z(&v_for), dt);
        x_for = mw_addv(&x_for, &tmp);
        
        
        // Compute the new acceleration
        acc = nbExtAcceleration(pot, &x_for, t);
        tmp.x = mw_mul_s(&X(&acc), dt_half);
        tmp.y = mw_mul_s(&Y(&acc), dt_half);
        tmp.z = mw_mul_s(&Z(&acc), dt_half);
        v_for = mw_addv(&v_for, &tmp);
        
        lbr = cartesianToLbr(&x_for, DEFAULT_SUN_GC_DISTANCE);
        fprintf(fp, "%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\n", showRealValue(&X(&x_for)), showRealValue(&Y(&x_for)), showRealValue(&Z(&x_for)), showRealValue(&X(&lbr)), showRealValue(&Y(&lbr)), showRealValue(&Z(&lbr)), showRealValue(&X(&v_for)), showRealValue(&Y(&v_for)), showRealValue(&Z(&v_for)));
    }
    fclose(fp);
    
    /* Report the final values (don't forget to reverse the velocities) */
    v = mw_negv(&v);

    *finalPos = x;
    *finalVel = v;
}

void nbPrintReverseOrbit_LMC(mwvector* finalPos,
                         mwvector* finalVel,
                         mwvector* LMCfinalPos,
                         mwvector* LMCfinalVel,
                         const Potential* pot,
                         mwvector* pos,
                         mwvector* vel,
                         mwvector* LMCposition,
                         mwvector* LMCvelocity,
                         mwbool LMCDynaFric,
                         real_0 tstop,
                         real_0 tstopforward,
                         real_0 dt,
                         real* LMCmass,
                         real* LMCscale)
{
    mwvector acc, mw_acc, LMC_acc, DF_acc, tmp, dv2, dx, dLMCv2, dLMCx;
    mwvector v, x, LMCv, LMCx;
    mwvector v_for, x_for, LMCv_for, LMCx_for;
    mwvector mw_x = mw_vec(ZERO_REAL, ZERO_REAL, ZERO_REAL);
    mwvector lbr;
    real_0 t;
    real_0 dt_half = dt/2.0;

    // Set the initial conditions
    x = *pos;
    v = *vel;
    LMCx = *LMCposition;
    LMCv = *LMCvelocity;
    x_for = *pos;
    v_for = *vel;
    LMCx_for = *LMCposition;
    LMCv_for = *LMCvelocity;
    
    v = mw_negv(&v);
    LMCv = mw_negv(&LMCv);

    // Get the initial acceleration (reverse)
    mw_acc = plummerAccel(&mw_x, &LMCx, LMCmass, LMCscale);
    LMC_acc = nbExtAcceleration(pot, &LMCx, 0);
    if (LMCDynaFric) {
        DF_acc = dynamicalFriction_LMC(pot, &LMCx, &LMCv, LMCmass, LMCscale, TRUE, 0);
        DF_acc = mw_negv(&DF_acc); /* Inverting drag force for reverse orbit */
        LMC_acc = mw_addv(&LMC_acc, &DF_acc);
     }
    acc = nbExtAcceleration(pot, &x, 0);
    tmp = plummerAccel(&x, &LMCx, LMCmass, LMCscale);
    acc = mw_addv(&acc, &tmp);
    LMC_acc = mw_subv(&LMC_acc,& mw_acc);
    acc = mw_subv(&acc, &mw_acc);

    FILE * fp;
    fp = fopen("reverse_orbit.out", "w");

    // Loop through time (reverse)
    for (t = 0; t >= tstop*(-1); t -= dt)
    {
        // Update the velocities and positions
        dv2.x = mw_mul_s(&acc.x, dt_half);
        dv2.y = mw_mul_s(&acc.y, dt_half);
        dv2.z = mw_mul_s(&acc.z, dt_half);
        v = mw_addv(&v, &dv2);

        dx.x = mw_mul_s(&v.x, dt);
        dx.y = mw_mul_s(&v.y, dt);
        dx.z = mw_mul_s(&v.z, dt);
        x = mw_addv(&x, &dx);

        dLMCv2.x = mw_mul_s(&LMC_acc.x, dt_half);
        dLMCv2.y = mw_mul_s(&LMC_acc.y, dt_half);
        dLMCv2.z = mw_mul_s(&LMC_acc.z, dt_half);
        LMCv = mw_addv(&LMCv, &dLMCv2);

        dLMCx.x = mw_mul_s(&LMCv.x, dt);
        dLMCx.y = mw_mul_s(&LMCv.y, dt);
        dLMCx.z = mw_mul_s(&LMCv.z, dt);
        LMCx = mw_addv(&LMCx, &dLMCx);
        
        // Compute the new acceleration
        mw_acc = plummerAccel(&mw_x, &LMCx, LMCmass, LMCscale);
        LMC_acc = nbExtAcceleration(pot, &LMCx, t);
        if (LMCDynaFric) {
            DF_acc = dynamicalFriction_LMC(pot, &LMCx, &LMCv, LMCmass, LMCscale, TRUE, t);
            LMC_acc = mw_subv(&LMC_acc, &DF_acc); /* Inverting drag force for reverse orbit */
        }
        acc = nbExtAcceleration(pot, &x, t);
        tmp = plummerAccel(&x, &LMCx, LMCmass, LMCscale);
    	acc = mw_addv(&acc, &tmp);

    	// Shift the body
    	mw_acc = mw_negv(&mw_acc);
        LMC_acc = mw_addv(&LMC_acc, &mw_acc);
        acc = mw_addv(&acc, &mw_acc);
        
        dv2.x = mw_mul_s(&acc.x, dt_half);
        dv2.y = mw_mul_s(&acc.y, dt_half);
        dv2.z = mw_mul_s(&acc.z, dt_half);
        v = mw_addv(&v, &dv2);

        dLMCv2.x = mw_mul_s(&LMC_acc.x, dt_half);
        dLMCv2.y = mw_mul_s(&LMC_acc.y, dt_half);
        dLMCv2.z = mw_mul_s(&LMC_acc.z, dt_half);
        LMCv = mw_addv(&LMCv, &dLMCv2);

        
        lbr = cartesianToLbr(&x, DEFAULT_SUN_GC_DISTANCE);
        fprintf(fp, "%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\n", showRealValue(&X(&x)), showRealValue(&Y(&x)), showRealValue(&Z(&x)), showRealValue(&X(&lbr)), showRealValue(&Y(&lbr)), showRealValue(&Z(&lbr)), showRealValue(&X(&v)), showRealValue(&Y(&v)), showRealValue(&Z(&v)));
    }
    fclose(fp);

    fp = fopen("forward_orbit.out", "w");

    //Get the initial acceleration (forward)
    mw_acc = plummerAccel(&mw_x, &LMCx_for, LMCmass, LMCscale);
    DF_acc = dynamicalFriction_LMC(pot, &LMCx_for, &LMCv_for, LMCmass, LMCscale, LMCDynaFric, 0);
    LMC_acc = nbExtAcceleration(pot, &LMCx_for, 0);
    LMC_acc = mw_addv(&LMC_acc, &DF_acc);
    acc = nbExtAcceleration(pot, &x_for, 0);
    tmp = plummerAccel(&x_for, &LMCx_for, LMCmass, LMCscale);
    acc = mw_addv(&acc, &tmp);
    LMC_acc = mw_subv(&LMC_acc, &mw_acc);
    acc = mw_subv(&acc, &mw_acc);

    //Loop through time (forward)
    for (t = 0; t <= tstopforward; t += dt)
    {
        // Update the velocities and positions
        dv2.x = mw_mul_s(&acc.x, dt_half);
        dv2.y = mw_mul_s(&acc.y, dt_half);
        dv2.z = mw_mul_s(&acc.z, dt_half);
        v_for = mw_addv(&v_for, &dv2);

        dx.x = mw_mul_s(&v_for.x, dt);
        dx.y = mw_mul_s(&v_for.y, dt);
        dx.z = mw_mul_s(&v_for.z, dt);
        x_for = mw_addv(&x_for, &dx);

        dLMCv2.x = mw_mul_s(&LMC_acc.x, dt_half);
        dLMCv2.y = mw_mul_s(&LMC_acc.y, dt_half);
        dLMCv2.z = mw_mul_s(&LMC_acc.z, dt_half);
        LMCv_for = mw_addv(&LMCv_for, &dLMCv2);

        dLMCx.x = mw_mul_s(&LMCv_for.x, dt);
        dLMCx.y = mw_mul_s(&LMCv_for.y, dt);
        dLMCx.z = mw_mul_s(&LMCv_for.z, dt);
        LMCx_for = mw_addv(&LMCx_for, &dLMCx);
    
        // Compute the new acceleration
        mw_acc = plummerAccel(&mw_x, &LMCx_for, LMCmass, LMCscale);
        DF_acc = dynamicalFriction_LMC(pot, &LMCx_for, &LMCv_for, LMCmass, LMCscale, LMCDynaFric, t);
        LMC_acc = nbExtAcceleration(pot, &LMCx_for, t);
        LMC_acc = mw_addv(&LMC_acc, &DF_acc);
        acc = nbExtAcceleration(pot, &x_for, t);
        tmp = plummerAccel(&x_for, &LMCx_for, LMCmass, LMCscale);
        acc = mw_addv(&acc, &tmp);
        // Shift the body
        mw_acc = mw_negv(&mw_acc);
        LMC_acc = mw_addv(&LMC_acc, &mw_acc);
        acc = mw_addv(&acc, &mw_acc);
        
        dv2.x = mw_mul_s(&acc.x, dt_half);
        dv2.y = mw_mul_s(&acc.y, dt_half);
        dv2.z = mw_mul_s(&acc.z, dt_half);
        v_for = mw_addv(&v_for, &dv2);

        dLMCv2.x = mw_mul_s(&LMC_acc.x, dt_half);
        dLMCv2.y = mw_mul_s(&LMC_acc.y, dt_half);
        dLMCv2.z = mw_mul_s(&LMC_acc.z, dt_half);
        LMCv_for = mw_addv(&LMCv_for, &dLMCv2);
        
        lbr = cartesianToLbr(&x_for, DEFAULT_SUN_GC_DISTANCE);
        fprintf(fp, "%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\n", showRealValue(&X(&x_for)), showRealValue(&Y(&x_for)), showRealValue(&Z(&x_for)), showRealValue(&X(&lbr)), showRealValue(&Y(&lbr)), showRealValue(&Z(&lbr)), showRealValue(&X(&v_for)), showRealValue(&Y(&v_for)), showRealValue(&Z(&v_for)));
    }
    fclose(fp);
    
    /* Report the final values (don't forget to reverse the velocities) */
    *finalPos = x;
    *finalVel = mw_negv(&v);
    *LMCfinalPos = LMCx;
    *LMCfinalVel = mw_negv(&LMCv);
}
