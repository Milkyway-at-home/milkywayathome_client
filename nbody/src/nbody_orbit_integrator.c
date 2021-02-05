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
#include "nbody_friction.h"
/* Simple orbit integrator in user-defined potential
    Written for BOINC Nbody
    willeb 10 May 2010 */
/* Altered to be consistent with nbody integrator.
 * shelton June 25 2018 */



mwvector* shiftByLMC = NULL; //Ptr to LMC Shift Array (default is NULL)
size_t nShiftLMC = 0;

mwvector* LMCpos = NULL; //Ptr to LMC position (default is NULL)
mwvector* LMCvel = NULL; //Ptr to LMC velocity (default is NULL)

void nbReverseOrbit(mwvector* finalPos,
                    mwvector* finalVel,
                    const Potential* pot,
                    mwvector pos,
                    mwvector vel,
                    real tstop,
                    real dt)
{
    mwvector acc, v, x;
    real t;
    real dt_half = dt / 2.0;
    
    // Set the initial conditions
    x = pos;
    v = vel;
    mw_incnegv(v);

    // Get the initial acceleration
    acc = nbExtAcceleration(pot, x);

    for (t = 0; t <= tstop; t += dt)
    {
        // Update the velocities and positions
        mw_incaddv_s(v, acc, dt_half);
        mw_incaddv_s(x, v, dt);
        
        // Compute the new acceleration
        acc = nbExtAcceleration(pot, x);
        
        mw_incaddv_s(v, acc, dt_half);
    }
    
    /* Report the final values (don't forget to reverse the velocities) */
    mw_incnegv(v);
    
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
                    real ftime,
                    real tstop,
                    real dt,
                    real LMCmass,
                    real LMCscale
                    )
{	
    unsigned int steps = mw_ceil((tstop)/(dt)) + 1;
    unsigned int exSteps = mw_abs(mw_ceil((ftime-tstop)/(dt)) + 1);
    unsigned int i = 0, j = 0, k = 0;
    mwvector acc, v, x, mw_acc, LMC_acc, DF_acc, LMCv, LMCx, tmp;
    mwvector mw_x = mw_vec(0, 0, 0);
    mwvector* bacArray = NULL;
    mwvector* forArray = NULL;

    //Placeholder arrays for LMC acceleration corrections
    bacArray = (mwvector*)mwCallocA(steps + 1, sizeof(mwvector));
    forArray = (mwvector*)mwCallocA(exSteps + 1, sizeof(mwvector));

    real t;
    real dt_half = dt / 2.0;

    // Check if forward time is larger than backward time. We will need to manually compute additional LMC accelerations in that case.
    if (ftime > tstop) {

        // Set the initial conditions for forward orbit
        x = pos;
        v = vel;
        LMCv = LMCvelocity;
        LMCx = LMCposition;

        // Get the initial acceleration
        mw_acc = plummerAccel(mw_x, LMCx, LMCmass, LMCscale);
        LMC_acc = mw_addv(nbExtAcceleration(pot, LMCx), dynamicalFriction_LMC(pot, LMCx, LMCv, LMCmass, LMCscale, LMCDynaFric));
        acc = nbExtAcceleration(pot, x);
        tmp = plummerAccel(x, LMCx, LMCmass, LMCscale);
        mw_incaddv(acc, tmp);

        // Shift the body
        mw_incnegv(mw_acc);
        mw_incaddv(LMC_acc, mw_acc);
        mw_incaddv(acc, mw_acc);

        for (t = 0; t <= (ftime-tstop); t += dt)
        {   
    	    exSteps = t/dt;
    	    if ((exSteps % 10 == 0)&&(t!=0)) { 
    	        forArray[k] = mw_acc;
                k++;
    	    }

            // Update the velocities and positions
            mw_incaddv_s(v, acc, dt_half);
            mw_incaddv_s(x, v, dt);
            mw_incaddv_s(LMCv, LMC_acc, dt_half);
            mw_incaddv_s(LMCx, LMCv, dt);
        
            // Compute the new acceleration
            mw_acc = plummerAccel(mw_x, LMCx, LMCmass, LMCscale);
            LMC_acc = mw_addv(nbExtAcceleration(pot, LMCx), dynamicalFriction_LMC(pot, LMCx, LMCv, LMCmass, LMCscale, LMCDynaFric));
            acc = nbExtAcceleration(pot, x);
            tmp = plummerAccel(x, LMCx, LMCmass, LMCscale);
    	    mw_incaddv(acc, tmp);

    	    // Shift the body
    	    mw_incnegv(mw_acc);
            mw_incaddv(LMC_acc, mw_acc);
            mw_incaddv(acc, mw_acc);
        
            mw_incaddv_s(v, acc, dt_half);
            mw_incaddv_s(LMCv, LMC_acc, dt_half);

        }
        forArray[k] = mw_acc; //set the last index after the loop ends
    }

    // Set the initial conditions for reverse orbit
    x = pos;
    v = vel;
    LMCv = LMCvelocity;
    LMCx = LMCposition;
    mw_incnegv(v);
    mw_incnegv(LMCv);


    // Get the initial acceleration
    mw_acc = plummerAccel(mw_x, LMCx, LMCmass, LMCscale);
    LMC_acc = nbExtAcceleration(pot, LMCx);
    if (LMCDynaFric) {
        DF_acc = dynamicalFriction_LMC(pot, LMCx, LMCv, LMCmass, LMCscale, TRUE);
        mw_incnegv(DF_acc); /* Inverting drag force for reverse orbit */
        mw_incaddv(LMC_acc, DF_acc)
     }
    acc = nbExtAcceleration(pot, x);
    tmp = plummerAccel(x, LMCx, LMCmass, LMCscale);
    mw_incaddv(acc, tmp);

    // Shift the body
    mw_incnegv(mw_acc);
    mw_incaddv(LMC_acc, mw_acc);
    mw_incaddv(acc, mw_acc);

    for (t = 0; t <= tstop; t += dt)
    {   
    	steps = t/dt;
    	if( steps % 10 == 0){ 
    		bacArray[i] = mw_acc;
        	i++;
    	}

        // Update the velocities and positions
        mw_incaddv_s(v, acc, dt_half);
        mw_incaddv_s(x, v, dt);
        mw_incaddv_s(LMCv, LMC_acc, dt_half);
        mw_incaddv_s(LMCx, LMCv, dt);
        
        // Compute the new acceleration
        mw_acc = plummerAccel(mw_x, LMCx, LMCmass, LMCscale);
        LMC_acc = nbExtAcceleration(pot, LMCx);
        if (LMCDynaFric) {
            DF_acc = dynamicalFriction_LMC(pot, LMCx, LMCv, LMCmass, LMCscale, TRUE);
            //mw_printf("DF: [%.15f,%.15f,%.15f]\n",X(DF_acc),Y(DF_acc),Z(DF_acc));
            mw_incnegv(DF_acc); /* Inverting drag force for reverse orbit */
            mw_incaddv(LMC_acc, DF_acc)
        }
        acc = nbExtAcceleration(pot, x);
        tmp = plummerAccel(x, LMCx, LMCmass, LMCscale);
    	mw_incaddv(acc, tmp);

    	// Shift the body
    	mw_incnegv(mw_acc);
        mw_incaddv(LMC_acc, mw_acc);
        mw_incaddv(acc, mw_acc);
        
        mw_incaddv_s(v, acc, dt_half);
        mw_incaddv_s(LMCv, LMC_acc, dt_half);

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
    *finalPos = x;
    *finalVel = v;
    *LMCfinalPos = LMCx;
    *LMCfinalVel = LMCv;

    //mw_printf("Initial LMC position: [%.15f,%.15f,%.15f]\n",X(LMCx),Y(LMCx),Z(LMCx));
    //mw_printf("Initial LMC velocity: [%.15f,%.15f,%.15f]\n",X(LMCv),Y(LMCv),Z(LMCv));

    //Allocate memory for LMC position and velocity
    LMCpos = (mwvector*)mwMalloc(sizeof(mwvector));
    LMCvel = (mwvector*)mwMalloc(sizeof(mwvector));

    //Store LMC position and velocity
    LMCpos[0] = LMCx;
    LMCvel[0] = LMCv;
}

void getLMCArray(mwvector ** shiftArrayPtr, size_t * shiftSizePtr) {
    //Allows access to shift array
    *shiftArrayPtr = shiftByLMC;
    *shiftSizePtr = nShiftLMC;
}

void getLMCPosVel(mwvector ** LMCposPtr, mwvector ** LMCvelPtr) {
    //Allows access to LMC position and velocity
    *LMCposPtr = LMCpos;
    *LMCvelPtr = LMCvel;
}



void nbPrintReverseOrbit(mwvector* finalPos,
                         mwvector* finalVel,
                         const Potential* pot,
                         mwvector pos,
                         mwvector vel,
                         real tstop,
                         real tstopforward,
                         real dt)
{
    mwvector acc, v, x;
    mwvector v_for, x_for;
    mwvector lbr;
    real t;
    real dt_half = dt / 2.0;

    // Set the initial conditions
    x = pos;
    v = vel;
    x_for = pos;
    v_for = vel;
    
    mw_incnegv(v);

    // Get the initial acceleration
    acc = nbExtAcceleration(pot, x);

    FILE * fp;
    fp = fopen("reverse_orbit.out", "w");
    // Loop through time
    for (t = 0; t <= tstop; t += dt)
    {
        // Update the velocities and positions
        mw_incaddv_s(v, acc, dt_half);
        mw_incaddv_s(x, v, dt);
        
        
        // Compute the new acceleration
        acc = nbExtAcceleration(pot, x);
        mw_incaddv_s(v, acc, dt_half);
        
        lbr = cartesianToLbr(x, DEFAULT_SUN_GC_DISTANCE);
        fprintf(fp, "%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\n", X(x), Y(x), Z(x), X(lbr), Y(lbr), Z(lbr), X(v), Y(v), Z(v));
    }

    fclose(fp);
    fp = fopen("forward_orbit.out", "w");
    acc = nbExtAcceleration(pot, x_for);
    for (t = 0; t <= tstopforward; t += dt)
    {
        // Update the velocities and positions
        mw_incaddv_s(v_for, acc, dt_half);
        mw_incaddv_s(x_for, v_for, dt);
        
        
        // Compute the new acceleration
        acc = nbExtAcceleration(pot, x_for);
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
