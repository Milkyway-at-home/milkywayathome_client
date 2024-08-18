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
#include <math.h>
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
                    real tstop,
                    real dt)
{
    mwvector acc, v, x;
    real t;
    real dt_half = dt / 2.0;
    int initialLArrayIndex = tstop/dt;

    // Set the initial conditions
    x = pos;
    v = vel;
    mw_incnegv(v);

    // Get the initial acceleration
    acc = nbExtAcceleration(pot, x, 0);
    //do this loop backward in order to get an accurate time for time-dependent potentials
    for (t = 0; t >= tstop*(-1); t -= dt)
    {
        // Update the velocities and positions
        mw_incaddv_s(v, acc, dt_half);
        mw_incaddv_s(x, v, dt); 

        
        // Compute the new acceleration
        acc = nbExtAcceleration(pot, x, t);
        
        mw_incaddv_s(v, acc, dt_half);
    }
    
    /* Report the final values (don't forget to reverse the velocities) */
    mw_incnegv(v);
    
    *finalPos = x;
    *finalVel = v;

    mw_printf("Dwarf Initial Position: [%.15f,%.15f,%.15f]\n", X(x), Y(x), Z(x));
    mw_printf("Dwarf Initial Velocity: [%.15f,%.15f,%.15f]\n", X(v), Y(v), Z(v));
}

void nbReverseOrbitS(mwvector* finalPos,
                    mwvector* finalVel,
                    const Potential* pot,
                    mwvector* pos,
                    mwvector* vel,
                    size_t len,
                    real tstop,
                    real dt,
                    real* masses)
{    
    mwvector acc[len], v[len], x[len];  // Set Initial Empty List for acc, v and x.
    mwvector track[200][len];

    for (size_t i = 0; i < len; ++i)  // Input the dwarfs data
    {
        x[i] = pos[i];
        v[i] = vel[i]; 
        mw_incnegv(v[i]);
    }
    int step = 0;
    for (real t = 0; t >= tstop*(-1); t -= dt)
    {   
        real dt_half = dt / 2.0;
        for (size_t i = 0; i < len; ++i)
        {
            acc[i] = nbExtAcceleration(pot, x[i], t);
            //if(step % 1000 == 0){mw_printf("Dwarf %d Acc : [%.15f,%.15f,%f] at t = %.15f\n", i+1, X(acc[i]), Y(acc[i]), Z(acc[i]), t);}
            // Total Accleration From other Dwarfs
            for (size_t j = 0; j < len; ++j)
            {
                if (i != j)
                {   
                    mwvector accFromJ = gravity(x[i], x[j], masses[j]);
                    acc[i] = mw_addv(acc[i], accFromJ);                   
                }
            }
            
            // Update velocity for first time at dt_half and update position using whole t
            mw_incaddv_s(v[i], acc[i], dt_half);
            mw_incaddv_s(x[i], v[i], dt); 

            acc[i] = nbExtAcceleration(pot, x[i], t-dt_half);
            //mwvector test = nbExtAcceleration(pot, x[i], t+dt_half);
            //if(step % 1000 == 0){mw_printf("Dwarf %d half Acc : [%.15f,%.15f,%f] at t = %.15f\n", i+1, X(acc[i]), Y(acc[i]), Z(acc[i]), t+dt_half);
            //                    mw_printf("Dwarf %d test half Acc : [%.15f,%.15f,%f] at t = %.15f\n", i+1, X(test), Y(test), Z(test), t-dt_half);}

            for (size_t j = 0; j < len; ++j)
            {
                if (i != j)
                {
                    mwvector accFromJ = gravity(x[i], x[j], masses[j]);
                    acc[i] = mw_addv(acc[i], accFromJ);
                }
            }
            mw_incaddv_s(v[i], acc[i], dt_half);// Update velocity for second time by dt_half step
        }
                // Track the coordinate during the calculation for 30000 steps
        if (step % 150 == 0) {
            //mwvector temp_vectors[len];
            for (size_t i = 0; i < len; ++i) {
                mwvector temp_vector;
                temp_vector.x = X(x[i]);
                temp_vector.y = Y(x[i]);
                temp_vector.z = Z(x[i]);
                track[step/150-1][i] = temp_vector;
            }
        }
        ++step;
    }


    for (size_t i = 0; i < len; ++i)// Store the Final velocity and position
    {   
        mw_incnegv(v[i]);
        finalPos[i] = x[i];        
        track[199][i] = x[i];
        finalVel[i] = v[i];
        mw_printf("Dwarf %d Initial Position: [%.15f,%.15f,%.15f]\n", i+1, X(x[i]), Y(x[i]), Z(x[i]));
        mw_printf("Dwarf %d Initial Velocity: [%.15f,%.15f,%.15f]\n", i+1, X(v[i]), Y(v[i]), Z(v[i]));
    }

    for (size_t i = 0; i < 200; ++i) {
        mw_printf("[");
        for (size_t j = 0; j < len; ++j){
            mw_printf("(%f, %f, %f)",track[i][j].x, track[i][j].y, track[i][j].z);
        }
        mw_printf("]\n");
    }
}

void nbReverseOrbitS_LMC(mwvector* finalPos,
                    mwvector* finalVel,
                    mwvector* LMCfinalPos,
                    mwvector* LMCfinalVel,
                    const Potential* pot,
                    mwvector* pos,
                    mwvector* vel,
                    size_t len,
                    mwvector LMCposition,
                    mwvector LMCvelocity,
                    mwbool LMCDynaFric,
                    real ftime,
                    real tstop,
                    real dt,
                    real LMCmass,
                    real LMCscale,
                    real coulomb_log,
                    real* masses,
                    real* scales
                    )
{	        
    unsigned int steps = mw_ceil((tstop)/(dt)) + 1;
    unsigned int exSteps = mw_abs(mw_ceil((ftime-tstop)/(dt)) + 1);
    unsigned int a = 0, b = 0, c = 0;
    mwvector acc[len], v[len], x[len];  // Set Initial Empty List for acc, v and x.    
    mwvector mw_acc, LMC_acc, DF_acc, LMCv, LMCx, tmp, LMC_acci;
    mwvector mw_x = mw_vec(0, 0, 0);
    mwvector* bacArray = NULL;
    mwvector* forArray = NULL;
    mwvector track[200][len];

    //Placeholder arrays for LMC acceleration corrections
    bacArray = (mwvector*)mwCallocA(steps + 1, sizeof(mwvector));
    forArray = (mwvector*)mwCallocA(exSteps + 1, sizeof(mwvector));
    
    // Set the initial conditions for reverse orbit
    for (size_t i = 0; i < len; ++i) // Input the dwarfs data
    {
        x[i] = pos[i];
        v[i] = vel[i]; 
        mw_incnegv(v[i]);
    }

    real t;
    real dt_half = dt / 2.0;

    LMCv = LMCvelocity;
    LMCx = LMCposition;
    mw_incnegv(LMCv);

    int step = 0;
    for (real t = 0; t >= tstop*(-1); t -= dt)
    {   
        mw_acc = plummerAccel(mw_x, LMCx, LMCmass, LMCscale); // Relative Accleration from Milkyway accleration from LMC
        steps = t/dt;

        // Shift the body
        mw_incnegv(mw_acc);

    	if( steps % 10 == 0){ 
    		bacArray[a] = mw_acc;
        	a++;
    	}   

        LMC_acc = nbExtAcceleration(pot, LMCx, t);
        if (LMCDynaFric) 
        {
            DF_acc = dynamicalFriction_LMC(pot, LMCx, LMCv, LMCmass, LMCscale, TRUE, 0, coulomb_log);
            mw_incnegv(DF_acc); /* Inverting drag force for reverse orbit */
            mw_incaddv(LMC_acc, DF_acc);
        }
        mw_incaddv(LMC_acc, mw_acc);
        
        //if (step % 10 == 0 & step <= 500) {mw_printf("LMC_acc1 from MW(both) and DF = [%.15f,%.15f,%.15f] at t = %f\n", X(LMC_acc), Y(LMC_acc), Z(LMC_acc), t);}

        for (size_t i = 0; i < len; ++i)
        {
            // Get the Starting accelerations
            acc[i] = nbExtAcceleration(pot, x[i], t);
            mw_incaddv(acc[i], mw_acc);
            mw_incaddv(acc[i], plummerAccel(x[i], LMCx, LMCmass, LMCscale));//Aceleration from LMC to dwarfs
            LMC_acci = plummerAccel(LMCx, x[i], masses[i], scales[i]);//Acceleration from dwarfs to LMC
            //LMC_acci = pointAccel(LMCx, x[i], masses[i]);
            mw_incaddv(LMC_acc, LMC_acci);
            // if (step % 10 == 0 & step <= 500) {mw_printf("LMC_acc1 from Dwarfs %d = [%.15f,%.15f,%.15f] at t = %f\n", i+1 , X(LMC_acci), Y(LMC_acci), Z(LMC_acci), t);}
            

            //if (step % 1000 == 0) // For Testing
            // Total Accleration From other Dwarfs
            for (size_t j = 0; j < len; ++j)
            {
                if (i != j)
                {   
                    mwvector accFromJ = plummerAccel(x[i], x[j], masses[j], scales[j]);
                    //mwvector accFromJ = pointAccel(x[i], x[j], masses[j]);
                    acc[i] = mw_addv(acc[i], accFromJ);                   
                }
            } 

            // Update Velocities for first time by dt_half step
            mw_incaddv_s(v[i], acc[i], dt_half);
            // Update the Dwarfs Positions for whole dt step
            mw_incaddv_s(x[i], v[i], dt); 
        }
        
        //if (step % 10 == 0 & step <= 500) {mw_printf("LMC_acc1 from MW(both) and DF = [%.15f,%.15f,%.15f] at t = %f\n", X(LMC_acc), Y(LMC_acc), Z(LMC_acc), t);}


        // //if (steps % 10 == 0 & steps <= 500) {
        //     mw_printf("LMC_acc = [%.15f,%.15f,%.15f] at t = %f\n", X(LMC_acc), Y(LMC_acc), Z(LMC_acc), t);
        //     mw_printf("LMCv = [%.15f,%.15f,%.15f] at t = %f\n",X(LMCv),Y(LMCv),Z(LMCv), t);}// For Testing
        mw_incaddv_s(LMCv, LMC_acc, dt_half);

        // Update the LMC Position for whole dt step
        mw_incaddv_s(LMCx, LMCv, dt);

        // Track the coordinate during the calculation for 30000 steps
        if (step % 150 == 0) {
            //mwvector temp_vectors[len];
            for (size_t i = 0; i < len; ++i) {
                mwvector temp_vector;
                temp_vector.x = X(x[i]);
                temp_vector.y = Y(x[i]);
                temp_vector.z = Z(x[i]);
                track[step/150-1][i] = temp_vector;
            }
        }

        //** Get the Second Half acceleration    
        mw_acc = plummerAccel(mw_x, LMCx, LMCmass, LMCscale); // Relative Accleration from Milkyway
        LMC_acc = nbExtAcceleration(pot, LMCx, t-dt_half);
        // Shift the body
        mw_incnegv(mw_acc);
        
        if (LMCDynaFric) {
            DF_acc = dynamicalFriction_LMC(pot, LMCx, LMCv, LMCmass, LMCscale, TRUE, 0, coulomb_log);
            mw_incnegv(DF_acc); /* Inverting drag force for reverse orbit */
            mw_incaddv(LMC_acc, DF_acc);
        }
        mw_incaddv(LMC_acc, mw_acc);

        //if (step % 10 == 0 & step <= 500) {mw_printf("LMC_acc2 from MW(both) and DF = [%.15f,%.15f,%.15f] at t = %f\n", X(LMC_acc), Y(LMC_acc), Z(LMC_acc), t-dt_half);}

        for (size_t i = 0; i < len; ++i)
        {
            //Get the second Starting accelerations
            acc[i] = nbExtAcceleration(pot, x[i], t-dt_half);
            mw_incaddv(acc[i], mw_acc);
            mw_incaddv(acc[i], plummerAccel(x[i], LMCx, LMCmass, LMCscale));//Aceleration from LMC to dwarfs

            LMC_acci = plummerAccel(LMCx, x[i], masses[i], scales[i]);
            //LMC_acci = pointAccel(LMCx, x[i], masses[i]);
            mw_incaddv(LMC_acc, LMC_acci);

            for (size_t j = 0; j < len; ++j)
            {
                if (i != j)
                {
                    mwvector accFromJ = plummerAccel(x[i], x[j], masses[j], scales[j]);
                    //mwvector accFromJ = pointAccel(x[i], x[j], masses[j]);
                    acc[i]= mw_addv(acc[i], accFromJ);
                }
            }
            
            // Update Velocities for second time by dt_half step
            mw_incaddv_s(v[i], acc[i], dt_half);
        }

        //if (step % 10 == 0 & step <= 500) {mw_printf("LMC_acc2 from MW and DF and Dwarfs= [%.15f,%.15f,%.15f] at t = %f\n\n", X(LMC_acc), Y(LMC_acc), Z(LMC_acc), t-dt_half);}

        // //if (steps % 10 == 0 & steps <= 500) {
        //     mw_printf("LMC_acc = [%.15f,%.15f,%.15f] at t = %f\n", X(LMC_acc), Y(LMC_acc), Z(LMC_acc), t);
        //     mw_printf("LMCv = [%.15f,%.15f,%.15f] at t = %f\n",X(LMCv),Y(LMCv),Z(LMCv), t);
        // }// For Testing
        
        mw_incaddv_s(LMCv, LMC_acc, dt_half);
        ++step;
    }
    
            //mw_printf("DF: [%.15f,%.15f,%.15f]\n",X(DF_acc),Y(DF_acc),Z(DF_acc));
            //mw_printf("LMCx: [%.15f,%.15f,%.15f] | ",X(LMCx),Y(LMCx),Z(LMCx));
            //mw_printf("LMCv: [%.15f,%.15f,%.15f]\n",X(LMCv),Y(LMCv),Z(LMCv));    
    bacArray[a] = mw_acc; //set the last index after the loop ends
    
    //Allocate memory for the shift array equal to (x,y,z) i times with extra wiggle room dependent on evolve time
    unsigned int size = a + c + 2;
    shiftByLMC = (mwvector*)mwCallocA(size, sizeof(mwvector)); 

    //Fill reverse orbit of shift array
    for(b = 0; b < a+1; b++) {
        tmp = bacArray[a-b];
        shiftByLMC[b] = tmp;
    }

    //Fill forward orbit of shift array
    if (ftime > tstop) {
        for(b = 0; b < c+1; b++) {
            tmp = forArray[b];
            shiftByLMC[a+1+b] = tmp;
        }
    }

    //Free placeholder arrays
    mwFreeA(bacArray);
    mwFreeA(forArray);

    /* Report the final values (don't forget to reverse the velocities) */
    for (size_t i = 0; i < len; ++i)// Store the Final velocity and position
    {   
        mw_incnegv(v[i]);
        finalPos[i] = x[i];
        track[199][i] = x[i];
        finalVel[i] = v[i];
        mw_printf("Dwarf %d Initial Position: %.15f,%.15f,%.15f\n", i+1, X(x[i]), Y(x[i]), Z(x[i]));
        mw_printf("Dwarf %d Initial Velocity: %.15f,%.15f,%.15f\n", i+1, X(v[i]), Y(v[i]), Z(v[i]));
    }

    mw_incnegv(LMCv);
    *LMCfinalPos = LMCx;
    *LMCfinalVel = LMCv;

    mw_printf("Initial LMC position: %.15f,%.15f,%.15f\n",X(LMCx),Y(LMCx),Z(LMCx));
    mw_printf("Initial LMC velocity: %.15f,%.15f,%.15f\n",X(LMCv),Y(LMCv),Z(LMCv));

    //Store LMC position and velocity
    LMCpos = LMCx;
    LMCvel = LMCv;

    for (size_t i = 0; i < 200; ++i) {
        mw_printf("[");
        for (size_t j = 0; j < len; ++j){
            mw_printf("(%f, %f, %f)",track[i][j].x, track[i][j].y, track[i][j].z);
        }
        mw_printf("]\n");
    }
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
                    real LMCscale,
                    real coulomb_log
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
        LMC_acc = mw_addv(nbExtAcceleration(pot, LMCx, 0), dynamicalFriction_LMC(pot, LMCx, LMCv, LMCmass, LMCscale, LMCDynaFric, 0, coulomb_log));
        acc = nbExtAcceleration(pot, x, 0);
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
            LMC_acc = mw_addv(nbExtAcceleration(pot, LMCx, t), dynamicalFriction_LMC(pot, LMCx, LMCv, LMCmass, LMCscale, LMCDynaFric, t, coulomb_log));
            acc = nbExtAcceleration(pot, x, t);
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
    LMC_acc = nbExtAcceleration(pot, LMCx, 0);
    if (LMCDynaFric) {
        DF_acc = dynamicalFriction_LMC(pot, LMCx, LMCv, LMCmass, LMCscale, TRUE, 0, coulomb_log);
        mw_incnegv(DF_acc); /* Inverting drag force for reverse orbit */
        mw_incaddv(LMC_acc, DF_acc);
     }
    acc = nbExtAcceleration(pot, x, 0);
    tmp = plummerAccel(x, LMCx, LMCmass, LMCscale);
    mw_incaddv(acc, tmp);

    // Shift the body
    mw_incnegv(mw_acc);
    mw_incaddv(LMC_acc, mw_acc);
    mw_incaddv(acc, mw_acc);

    real negT = 0;
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
        mw_incaddv_s(v, acc, dt_half);
        mw_incaddv_s(x, v, dt);
        mw_incaddv_s(LMCv, LMC_acc, dt_half);
        mw_incaddv_s(LMCx, LMCv, dt);
        
        // Compute the new acceleration
        mw_acc = plummerAccel(mw_x, LMCx, LMCmass, LMCscale);
        LMC_acc = nbExtAcceleration(pot, LMCx, negT);
        if (LMCDynaFric) {
            DF_acc = dynamicalFriction_LMC(pot, LMCx, LMCv, LMCmass, LMCscale, TRUE, negT, coulomb_log);
            //mw_printf("DF: [%.15f,%.15f,%.15f]\n",X(DF_acc),Y(DF_acc),Z(DF_acc));
            mw_incnegv(DF_acc); /* Inverting drag force for reverse orbit */
            mw_incaddv(LMC_acc, DF_acc);
        }
        acc = nbExtAcceleration(pot, x, negT);
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

    mw_printf("Dwarf Initial Position: [%.15f,%.15f,%.15f]\n", X(x), Y(x), Z(x));
    mw_printf("Dwarf Initial Velocity: [%.15f,%.15f,%.15f]\n", X(v), Y(v), Z(v));
    mw_printf("Initial LMC position: [%.15f,%.15f,%.15f]\n",X(LMCx),Y(LMCx),Z(LMCx));
    mw_printf("Initial LMC velocity: [%.15f,%.15f,%.15f]\n",X(LMCv),Y(LMCv),Z(LMCv));

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

mwvector gravity(mwvector v1, mwvector v2, real m2) {
    mwvector acc = ZERO_VECTOR; 
    real rx, ry, rz;
    rx = X(v2) - X(v1);
    ry = Y(v2) - Y(v1);
    rz = Z(v2) - Z(v1);
    real r = sqrt(rx * rx + ry * ry + rz * rz);

    if (r == 0.0) {
        mw_printf("Detected dwarfs too close");
        return acc;
    }

    real factor = m2 / pow(r, 3);
    X(acc) = factor * rx;
    Y(acc) = factor * ry;
    Z(acc) = factor * rz;
    
    return acc;
}