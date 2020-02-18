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
/* Simple orbit integrator in user-defined potential
    Written for BOINC Nbody
    willeb 10 May 2010 */
/* Altered to be consistent with nbody integrator.
 * shelton June 25 2018 */

//Ptr to LMC Shift Array (default is NULL)

real** shiftByLMC = NULL;

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
                    mwvector LMCpos,
                    mwvector LMCvel,
                    real tstop,
                    real dt,
                    real LMCmass
                    )
{	
   
	int steps = (tstop/dt/10+1);
    mwvector acc, v, x, mw_acc, LMC_acc, LMCv, LMCx, tmp;
    mwvector mw_x = mw_vec(0, 0, 0);
    mwvector array[steps];
    real t;
    real dt_half = dt / 2.0;
    int i = 0;
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

    array[i] = mw_acc;
    /*FILE *fA = fopen("testArray.txt", "w");
    fprintf(fA,"i: %d\n\n", i);*/ //print code in place to confirm shift array is working--will remove for final version

    //Allocate memory for the shift array equal to (x,y,z) i times
    shiftByLMC = (real**)mwMalloc( (i+1)* sizeof(real*)); 
    int idx;
    for(idx = 0; idx < (i+1); idx++) {
        shiftByLMC[idx] = (real*)mwMalloc(3*sizeof(real));
    }
    
    //Fill the shift array with the calculated values from "array" from index=i to index=2
    int j;
    for(j = i; j > 1; j--) {
        tmp = array[j];
        shiftByLMC[j][0] = tmp.x;
        shiftByLMC[j][1] = tmp.y;
        shiftByLMC[j][2] = tmp.z;
        //fprintf(fA, "j: %d   %f %f %f\n", j, tmp.x, tmp.y, tmp.z);
    }
    //fclose(fA);

    /* Report the final values (don't forget to reverse the velocities) */
    mw_incnegv(v);
    mw_incnegv(LMCv);
    *finalPos = x;
    *finalVel = v;
    *LMCfinalPos = LMCx;
    *LMCfinalVel = LMCv;
}

void getLMCArray(real *** shiftArrayPtr) {
    //call : real** tmpShiftPtr;
    //getLMCArray(&tmpShiftPtr);
    *shiftArrayPtr = shiftByLMC;
}
void freeLMCArray(unsigned int size) {
    int i;
    for(i = 0; i < size; i++) {
        mwFreeA(shiftByLMC[i]);
    }
    mwFreeA(shiftByLMC);
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
