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
/* Simple orbit integrator in user-defined potential
    Written for BOINC Nbody
    willeb 10 May 2010 */
/* Altered to be consistent with nbody integrator.
 * shelton June 25 2018 */

mwvector* backwardOrbitPositions;//this is the l coordinate of the backwards
//orbit point mass at each timestep (timestep being the index)
//index 0 is the end of the backward orbit
//note: make sure to free it
int backwardOrbitArraySize;


mwvector** shiftByLMC = NULL; //Ptr to LMC Shift Array (default is NULL)

void nbReverseOrbit(mwvector* finalPos,
                    mwvector* finalVel,
                    const Potential* pot,
                    mwvector pos,
                    mwvector vel,
                    real tstop,
                    real dt)
{
    mw_printf("start pos: x: %f y: %f z: %f\n", pos.x, pos.y, pos.z);
    mw_printf("start vel: x: %f y: %f z: %f\n", vel.x, vel.y, vel.z);
    mwvector acc, v, x;
    real t;
    real dt_half = dt / 2.0;
    int LArrayIndex = tstop / dt;
    backwardOrbitArraySize = LArrayIndex;
    mw_printf("Array Length = %d\n", backwardOrbitArraySize);
    // Set the initial conditions
    x = pos;
    v = vel;
    mw_incnegv(v);

    // Get the initial acceleration
    acc = nbExtAcceleration(pot, x, 0);
    //mw_printf("tstop: %f\n", tstop);
    
    //allocate backwardsOrbitL
    backwardOrbitPositions = (mwvector*)mwMalloc(LArrayIndex * sizeof(mwvector));

    for (t = 0; t >= tstop*(-1); t -= dt)
    {
        // Update the velocities and positions
        mw_incaddv_s(v, acc, dt_half);
        mw_incaddv_s(x, v, dt); 

        
        // Compute the new acceleration
        acc = nbExtAcceleration(pot, x, t);
        
        mw_incaddv_s(v, acc, dt_half);

        backwardOrbitPositions[LArrayIndex] = x;
        //record the current l coordinate
        LArrayIndex--;
    }
    
    /* Report the final values (don't forget to reverse the velocities) */
    mw_incnegv(v);
    
    *finalPos = x;
    *finalVel = v;
}

void getBackwardOrbitArray(mwvector** ptr) {
    //call : real* tmpOrbitPtr;
    //getBackwardOrbitArray(&tmpOrbitPtr);
    *ptr = backwardOrbitPositions;
}

int getOrbitArraySize(){
    return backwardOrbitArraySize;
}

void nbReverseOrbit_LMC(mwvector* finalPos,
                    mwvector* finalVel,
                    mwvector* LMCfinalPos,
                    mwvector* LMCfinalVel,
                    const Potential* pot,
                    mwvector pos,
                    mwvector vel,
                    mwvector LMCpos,
                    mwvector LMCvel,
                    real tstop,
                    real dt,
                    real LMCmass
                    )
{	
   
	unsigned int steps = (tstop)/ (dt) + 1;
    unsigned int i = 0, j = 0;
    mwvector acc, v, x, mw_acc, LMC_acc, LMCv, LMCx, tmp;
    mwvector mw_x = mw_vec(0, 0, 0);
    mwvector array[steps + 1];
    real t;
    real dt_half = dt / 2.0;
    // Set the initial conditions
    x = pos;
    v = vel;
    LMCv = LMCvel;
    LMCx = LMCpos;
    mw_incnegv(v);
    mw_incnegv(LMCv);


    // Get the initial acceleration
    mw_acc = pointAccel(mw_x, LMCx, LMCmass);
    LMC_acc = nbExtAcceleration(pot, LMCx);
    acc = nbExtAcceleration(pot, x);
    tmp = pointAccel(x, LMCx, LMCmass);
    mw_incaddv(acc, tmp);

    // Shift the body
    mw_incnegv(mw_acc);
    mw_incaddv(LMC_acc, mw_acc);
    mw_incaddv(acc, mw_acc);

    for (t = 0; t <= tstop; t += dt)
    {   
    	steps = t/dt;
    	if( steps % 10 == 0){ 
    		array[i] = mw_acc;
        	i++;
    	}

        // Update the velocities and positions
        mw_incaddv_s(v, acc, dt_half);
        mw_incaddv_s(x, v, dt);
        mw_incaddv_s(LMCv, LMC_acc, dt_half);
        mw_incaddv_s(LMCx, LMCv, dt);
        
        // Compute the new acceleration
        mw_acc = pointAccel(mw_x, LMCx, LMCmass);
        LMC_acc = nbExtAcceleration(pot, LMCx);
        acc = nbExtAcceleration(pot, x);
    	tmp = pointAccel(x, LMCx, LMCmass);
    	mw_incaddv(acc, tmp);

    	// Shift the body
    	mw_incnegv(mw_acc);
        mw_incaddv(LMC_acc, mw_acc);
        mw_incaddv(acc, mw_acc);
        
        mw_incaddv_s(v, acc, dt_half);
        mw_incaddv_s(LMCv, LMC_acc, dt_half);

    }
    array[i] = mw_acc; //set the last index after the loop ends
    
    //Allocate memory for the shift array equal to (x,y,z) i times with extra wiggle room dependent on evolve time
    unsigned int size = i+ 40*tstop;
    shiftByLMC = (mwvector**)mwMalloc((size)* sizeof(mwvector*)); 
    for(j = 0; j < (size); j++) {
        shiftByLMC[j] = (mwvector*)mwMalloc(sizeof(mwvector));
    }
    
    //Fill the shift array forward
    for(j = 0; j < i+1; j++) {
        tmp = array[i-j];
        shiftByLMC[j][0] = tmp;
    }

    /* Report the final values (don't forget to reverse the velocities) */
    mw_incnegv(v);
    mw_incnegv(LMCv);
    *finalPos = x;
    *finalVel = v;
    *LMCfinalPos = LMCx;
    *LMCfinalVel = LMCv;
}

void getLMCArray(mwvector *** shiftArrayPtr) {
    //Allows access to shift array
    *shiftArrayPtr = shiftByLMC;
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
