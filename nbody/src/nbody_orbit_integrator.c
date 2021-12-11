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
                    mwvector pos,
                    mwvector vel,
                    real_0 tstop,
                    real_0 dt,
                    real_0 sun_dist)
{
    mwvector v_var, x_var, x_lbr;
    mwvector acc, v, x;
    real_0 t;
    real_0 dt_half = dt / 2.0;
    int initialLArrayIndex = tstop/dt;

    // Set derivative information for AUTODIFF
    x_lbr = cartesianToLbr(pos, sun_dist);
    B(x_lbr) = mw_real_var(showRealValue(B(x_lbr)), 6);
    R(x_lbr) = mw_real_var(showRealValue(R(x_lbr)), 7);
    x_var = lbrToCartesian(x_lbr, sun_dist);

    X(v_var) = mw_real_var(showRealValue(X(vel)), 8);
    Y(v_var) = mw_real_var(showRealValue(Y(vel)), 9);
    Z(v_var) = mw_real_var(showRealValue(Z(vel)), 10);

    // Set the initial conditions
    x = x_var;
    v = v_var;
    mw_incnegv(v);

    // Get the initial acceleration
    acc = nbExtAcceleration(pot, x, 0);
    //do this loop backward in order to get an accurate time for time-dependent potentials
    for (t = 0; t >= tstop*(-1); t -= dt)
    {
        // Update the velocities and positions
        mw_incaddv_s(v, acc, mw_real_const(dt_half));
        mw_incaddv_s(x, v, mw_real_const(dt)); 

        
        // Compute the new acceleration
        acc = nbExtAcceleration(pot, x, t);
        
        mw_incaddv_s(v, acc, mw_real_const(dt_half));
    }
    
    /* Report the final values (don't forget to reverse the velocities) */
    mw_incnegv(v);

    //FIXME: NEEDS EVOLUTION TIME DERIVATIVE INFORMATION!
    *finalPos = x;
    *finalVel = v;
}

void nbReverseOrbit_LMC(mwvector* finalPos,
                    mwvector* finalVel,
                    mwvector* LMCfinalPos,
                    mwvector* LMCfinalVel,
                    const Potential* pot,
                    mwvector pos,
                    mwvector vel,
                    mwvector LMCposition,
                    mwvector LMCvelocity,
                    mwbool LMCDynaFric,
                    real_0 ftime,
                    real_0 tstop,
                    real_0 dt,
                    real LMCmass,
                    real LMCscale,
                    real_0 sun_dist)
{	
    unsigned int steps = mw_ceil((tstop)/(dt)) + 1;
    unsigned int exSteps = mw_abs(mw_ceil((ftime-tstop)/(dt)) + 1);
    unsigned int i = 0, j = 0, k = 0;
    mwvector v_var, x_var, x_lbr;
    mwvector acc, v, x, mw_acc, LMC_acc, DF_acc, LMCv, LMCx, tmp;
    mwvector mw_x = mw_vec(ZERO_REAL, ZERO_REAL, ZERO_REAL);
    mwvector* bacArray = NULL;
    mwvector* forArray = NULL;

    //Placeholder arrays for LMC acceleration corrections
    bacArray = (mwvector*)mwCallocA(steps + 1, sizeof(mwvector));
    forArray = (mwvector*)mwCallocA(exSteps + 1, sizeof(mwvector));

    real_0 t;
    real_0 dt_half = dt / 2.0;

    // Set derivative information for AUTODIFF
    x_lbr = cartesianToLbr(pos, sun_dist);
    B(x_lbr) = mw_real_var(showRealValue(B(x_lbr)), 6);
    R(x_lbr) = mw_real_var(showRealValue(R(x_lbr)), 7);
    x_var = lbrToCartesian(x_lbr, sun_dist);

    X(v_var) = mw_real_var(showRealValue(X(vel)), 8);
    Y(v_var) = mw_real_var(showRealValue(Y(vel)), 9);
    Z(v_var) = mw_real_var(showRealValue(Z(vel)), 10);

    // Check if forward time is larger than backward time. We will need to manually compute additional LMC accelerations in that case.
    if (ftime > tstop) {

        // Set the initial conditions for forward orbit
        x = x_var;
        v = v_var;
        LMCv = LMCvelocity;
        LMCx = LMCposition;

        // Get the initial acceleration
        mw_acc = plummerAccel(mw_x, LMCx, LMCmass, LMCscale);
        LMC_acc = mw_addv(nbExtAcceleration(pot, LMCx, 0), dynamicalFriction_LMC(pot, LMCx, LMCv, LMCmass, LMCscale, LMCDynaFric, 0));
        acc = nbExtAcceleration(pot, x, 0);
        tmp = plummerAccel(x, LMCx, LMCmass, LMCscale);
        mw_incaddv(acc, tmp);

        // Shift the body
        mw_incnegv(mw_acc);
        mw_incaddv(LMC_acc, mw_acc);
        mw_incaddv(acc, mw_acc);

        for (t = 0; t <= (ftime-tstop); t += dt)
        {   
    	    exSteps = (int) mw_round_0(t/dt);
    	    if ((exSteps % 10 == 0)&&(t!=0)) { 
    	        forArray[k] = mw_acc;
                k++;
    	    }

            // Update the velocities and positions
            mw_incaddv_s(v, acc, mw_real_const(dt_half));
            mw_incaddv_s(x, v, mw_real_const(dt));
            mw_incaddv_s(LMCv, LMC_acc, mw_real_const(dt_half));
            mw_incaddv_s(LMCx, LMCv, mw_real_const(dt));
        
            // Compute the new acceleration
            mw_acc = plummerAccel(mw_x, LMCx, LMCmass, LMCscale);
            LMC_acc = mw_addv(nbExtAcceleration(pot, LMCx, t), dynamicalFriction_LMC(pot, LMCx, LMCv, LMCmass, LMCscale, LMCDynaFric, t));
            acc = nbExtAcceleration(pot, x, t);
            tmp = plummerAccel(x, LMCx, LMCmass, LMCscale);
    	    mw_incaddv(acc, tmp);

    	    // Shift the body
    	    mw_incnegv(mw_acc);
            mw_incaddv(LMC_acc, mw_acc);
            mw_incaddv(acc, mw_acc);
        
            mw_incaddv_s(v, acc, mw_real_const(dt_half));
            mw_incaddv_s(LMCv, LMC_acc, mw_real_const(dt_half));

        }
        forArray[k] = mw_acc; //set the last index after the loop ends
    }

    // Set the initial conditions for reverse orbit
    x = x_var;
    v = v_var;
    LMCv = LMCvelocity;
    LMCx = LMCposition;
    mw_incnegv(v);
    mw_incnegv(LMCv);


    // Get the initial acceleration
    mw_acc = plummerAccel(mw_x, LMCx, LMCmass, LMCscale);
    LMC_acc = nbExtAcceleration(pot, LMCx, 0);
    if (LMCDynaFric) {
        DF_acc = dynamicalFriction_LMC(pot, LMCx, LMCv, LMCmass, LMCscale, TRUE, 0);
        mw_incnegv(DF_acc); /* Inverting drag force for reverse orbit */
        mw_incaddv(LMC_acc, DF_acc)
     }
    acc = nbExtAcceleration(pot, x, 0);
    tmp = plummerAccel(x, LMCx, LMCmass, LMCscale);
    mw_incaddv(acc, tmp);

    // Shift the body
    mw_incnegv(mw_acc);
    mw_incaddv(LMC_acc, mw_acc);
    mw_incaddv(acc, mw_acc);

    real_0 negT = 0;
    for (t = 0; t <= tstop; t += dt)
    {   
        //negate this time for use in time-dependent potentials
        negT = t*-1;
    	steps = t/dt;
    	if( steps % 10 == 0){ 
    		bacArray[i] = mw_acc;
        	i++;
    	}

        // Update the velocities and positions
        mw_incaddv_s(v, acc, mw_real_const(dt_half));
        mw_incaddv_s(x, v, mw_real_const(dt));
        mw_incaddv_s(LMCv, LMC_acc, mw_real_const(dt_half));
        mw_incaddv_s(LMCx, LMCv, mw_real_const(dt));
        
        // Compute the new acceleration
        mw_acc = plummerAccel(mw_x, LMCx, LMCmass, LMCscale);
        LMC_acc = nbExtAcceleration(pot, LMCx, negT);
        if (LMCDynaFric) {
            DF_acc = dynamicalFriction_LMC(pot, LMCx, LMCv, LMCmass, LMCscale, TRUE, negT);
            //mw_printf("DF: [%.15f,%.15f,%.15f]\n",X(DF_acc),Y(DF_acc),Z(DF_acc));
            mw_incnegv(DF_acc); /* Inverting drag force for reverse orbit */
            mw_incaddv(LMC_acc, DF_acc)
        }
        acc = nbExtAcceleration(pot, x, negT);
        tmp = plummerAccel(x, LMCx, LMCmass, LMCscale);
    	mw_incaddv(acc, tmp);

    	// Shift the body
    	mw_incnegv(mw_acc);
        mw_incaddv(LMC_acc, mw_acc);
        mw_incaddv(acc, mw_acc);
        
        mw_incaddv_s(v, acc, mw_real_const(dt_half));
        mw_incaddv_s(LMCv, LMC_acc, mw_real_const(dt_half));

        //mw_printf("LMCx: [%.15f,%.15f,%.15f] | ",X(LMCx),Y(LMCx),Z(LMCx));
        //mw_printf("LMCv: [%.15f,%.15f,%.15f]\n",X(LMCv),Y(LMCv),Z(LMCv));

    }
    bacArray[i] = mw_acc; //set the last index after the loop ends
    
    //Allocate memory for the shift array equal to (x,y,z) i times with extra wiggle room dependent on evolve time
    unsigned int size = i + k + 2;
    shiftByLMC = (mwvector*)mwCallocA(size, sizeof(mwvector)); 

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
    mw_incnegv(v);
    mw_incnegv(LMCv);

    //FIXME: NEEDS EVOLUTION TIME DERIVATIVE INFORMATION!
    *finalPos = x;
    *finalVel = v;
    *LMCfinalPos = LMCx;
    *LMCfinalVel = LMCv;

    //mw_printf("Initial LMC position: [%.15f,%.15f,%.15f]\n",X(LMCx),Y(LMCx),Z(LMCx));
    //mw_printf("Initial LMC velocity: [%.15f,%.15f,%.15f]\n",X(LMCv),Y(LMCv),Z(LMCv));

    //Store LMC position and velocity
    LMCpos = LMCx;
    LMCvel = LMCv;
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
                         mwvector pos,
                         mwvector vel,
                         real_0 tstop,
                         real_0 tstopforward,
                         real_0 dt)
{
    mwvector acc, v, x;
    mwvector v_for, x_for;
    mwvector lbr;
    real_0 t;
    real_0 dt_half = dt / 2.0;

    // Set the initial conditions
    x = pos;
    v = vel;
    x_for = pos;
    v_for = vel;
    
    mw_incnegv(v);

    // Get the initial acceleration
    acc = nbExtAcceleration(pot, x, 0);

    FILE * fp;
    fp = fopen("reverse_orbit.out", "w");
    // Loop through time
    for (t = 0; t >= tstop*(-1); t -= dt)
    {
        // Update the velocities and positions
        mw_incaddv_s(v, acc, dt_half);
        mw_incaddv_s(x, v, dt);
        
        
        // Compute the new acceleration
        acc = nbExtAcceleration(pot, x, t);
        mw_incaddv_s(v, acc, dt_half);
        
        lbr = cartesianToLbr(x, DEFAULT_SUN_GC_DISTANCE);
        fprintf(fp, "%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\n", X(x), Y(x), Z(x), X(lbr), Y(lbr), Z(lbr), X(v), Y(v), Z(v));
    }

    fclose(fp);
    fp = fopen("forward_orbit.out", "w");
    acc = nbExtAcceleration(pot, x_for, 0);
    for (t = 0; t <= tstopforward; t += dt)
    {
        // Update the velocities and positions
        mw_incaddv_s(v_for, acc, dt_half);
        mw_incaddv_s(x_for, v_for, dt);
        
        
        // Compute the new acceleration
        acc = nbExtAcceleration(pot, x_for, t);
        mw_incaddv_s(v_for, acc, dt_half);
        
        lbr = cartesianToLbr(x_for, DEFAULT_SUN_GC_DISTANCE);
        fprintf(fp, "%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\n", X(x_for), Y(x_for), Z(x_for), X(lbr), Y(lbr), Z(lbr), X(v_for), Y(v_for), Z(v_for));
    }
    fclose(fp);
    
    /* Report the final values (don't forget to reverse the velocities) */
    mw_incnegv(v);

    *finalPos = x;
    *finalVel = v;
}

void nbPrintReverseOrbit_LMC(mwvector* finalPos,
                         mwvector* finalVel,
                         mwvector* LMCfinalPos,
                         mwvector* LMCfinalVel,
                         const Potential* pot,
                         mwvector pos,
                         mwvector vel,
                         mwvector LMCposition,
                         mwvector LMCvelocity,
                         mwbool LMCDynaFric,
                         real_0 tstop,
                         real_0 tstopforward,
                         real_0 dt,
                         real LMCmass,
                         real LMCscale)
{
    mwvector acc, mw_acc, LMC_acc, DF_acc, tmp;
    mwvector v, x, LMCv, LMCx;
    mwvector v_for, x_for, LMCv_for, LMCx_for;
    mwvector mw_x = mw_vec(ZERO_REAL, ZERO_REAL, ZERO_REAL);
    mwvector lbr;
    real_0 t;
    real_0 dt_half = dt / 2.0;

    // Set the initial conditions
    x = pos;
    v = vel;
    LMCx = LMCposition;
    LMCv = LMCvelocity;
    x_for = pos;
    v_for = vel;
    LMCx_for = LMCposition;
    LMCv_for = LMCvelocity;
    
    mw_incnegv(v);
    mw_incnegv(LMCv);

    // Get the initial acceleration (reverse)
    mw_acc = plummerAccel(mw_x, LMCx, LMCmass, LMCscale);
    LMC_acc = nbExtAcceleration(pot, LMCx, 0);
    if (LMCDynaFric) {
        DF_acc = dynamicalFriction_LMC(pot, LMCx, LMCv, LMCmass, LMCscale, TRUE, 0);
        mw_incnegv(DF_acc); /* Inverting drag force for reverse orbit */
        mw_incaddv(LMC_acc, DF_acc)
     }
    acc = nbExtAcceleration(pot, x, 0);
    tmp = plummerAccel(x, LMCx, LMCmass, LMCscale);
    mw_incaddv(acc, tmp);
    mw_incnegv(mw_acc);
    mw_incaddv(LMC_acc, mw_acc);
    mw_incaddv(acc, mw_acc);

    FILE * fp;
    fp = fopen("reverse_orbit.out", "w");

    // Loop through time (reverse)
    for (t = 0; t >= tstop*(-1); t -= dt)
    {
        // Update the velocities and positions
        mw_incaddv_s(v, acc, mw_real_const(dt_half));
        mw_incaddv_s(x, v, mw_real_const(dt));
        mw_incaddv_s(LMCv, LMC_acc, mw_real_const(dt_half));
        mw_incaddv_s(LMCx, LMCv, mw_real_const(dt));
        
        // Compute the new acceleration
        mw_acc = plummerAccel(mw_x, LMCx, LMCmass, LMCscale);
        LMC_acc = nbExtAcceleration(pot, LMCx, t);
        if (LMCDynaFric) {
            DF_acc = dynamicalFriction_LMC(pot, LMCx, LMCv, LMCmass, LMCscale, TRUE, t);
            mw_incnegv(DF_acc); /* Inverting drag force for reverse orbit */
            mw_incaddv(LMC_acc, DF_acc)
        }
        acc = nbExtAcceleration(pot, x, t);
        tmp = plummerAccel(x, LMCx, LMCmass, LMCscale);
    	mw_incaddv(acc, tmp);

    	// Shift the body
    	mw_incnegv(mw_acc);
        mw_incaddv(LMC_acc, mw_acc);
        mw_incaddv(acc, mw_acc);
        
        mw_incaddv_s(v, acc, mw_real_const(dt_half));
        mw_incaddv_s(LMCv, LMC_acc, mw_real_const(dt_half));

        
        lbr = cartesianToLbr(x, DEFAULT_SUN_GC_DISTANCE);
        fprintf(fp, "%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\n", showRealValue(X(x)), showRealValue(Y(x)), showRealValue(Z(x)), showRealValue(X(lbr)), showRealValue(Y(lbr)), showRealValue(Z(lbr)), showRealValue(X(v)), showRealValue(Y(v)), showRealValue(Z(v)));
    }
    fclose(fp);

    fp = fopen("forward_orbit.out", "w");

    //Get the initial acceleration (forward)
    mw_acc = plummerAccel(mw_x, LMCx_for, LMCmass, LMCscale);
    LMC_acc = mw_addv(nbExtAcceleration(pot, LMCx_for, 0), dynamicalFriction_LMC(pot, LMCx_for, LMCv_for, LMCmass, LMCscale, LMCDynaFric, 0));
    acc = nbExtAcceleration(pot, x_for, 0);
    tmp = plummerAccel(x_for, LMCx_for, LMCmass, LMCscale);
    mw_incaddv(acc, tmp);
    mw_incnegv(mw_acc);
    mw_incaddv(LMC_acc, mw_acc);
    mw_incaddv(acc, mw_acc);

    //Loop through time (forward)
    for (t = 0; t <= tstopforward; t += dt)
    {
        // Update the velocities and positions
        mw_incaddv_s(v_for, acc, mw_real_const(dt_half));
        mw_incaddv_s(x_for, v_for, mw_real_const(dt));
        mw_incaddv_s(LMCv_for, LMC_acc, mw_real_const(dt_half));
        mw_incaddv_s(LMCx_for, LMCv_for, mw_real_const(dt));
    
        // Compute the new acceleration
        mw_acc = plummerAccel(mw_x, LMCx_for, LMCmass, LMCscale);
        LMC_acc = mw_addv(nbExtAcceleration(pot, LMCx_for, t), dynamicalFriction_LMC(pot, LMCx_for, LMCv_for, LMCmass, LMCscale, LMCDynaFric, t));
        acc = nbExtAcceleration(pot, x_for, t);
        tmp = plummerAccel(x_for, LMCx_for, LMCmass, LMCscale);
        mw_incaddv(acc, tmp);
        // Shift the body
        mw_incnegv(mw_acc);
        mw_incaddv(LMC_acc, mw_acc);
        mw_incaddv(acc, mw_acc);
        
        mw_incaddv_s(v_for, acc, mw_real_const(dt_half));
        mw_incaddv_s(LMCv_for, LMC_acc, mw_real_const(dt_half));
        
        lbr = cartesianToLbr(x_for, DEFAULT_SUN_GC_DISTANCE);
        fprintf(fp, "%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\n", showRealValue(X(x_for)), showRealValue(Y(x_for)), showRealValue(Z(x_for)), showRealValue(X(lbr)), showRealValue(Y(lbr)), showRealValue(Z(lbr)), showRealValue(X(v_for)), showRealValue(Y(v_for)), showRealValue(Z(v_for)));
    }
    fclose(fp);
    
    /* Report the final values (don't forget to reverse the velocities) */
    mw_incnegv(v);
    mw_incnegv(LMCv);

    *finalPos = x;
    *finalVel = v;
    *LMCfinalPos = LMCx;
    *LMCfinalVel = LMCv;
}
