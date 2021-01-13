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
#include "nbody_bar_time.h"
#include "nbody_friction.h"
/* Simple orbit integrator in user-defined potential
    Written for BOINC Nbody
    willeb 10 May 2010 */
/* Altered to be consistent with nbody integrator.
 * shelton June 25 2018 */

mwvector* backwardOrbitPositions;//this is the l coordinate of the backwards
//orbit point mass at each timestep (timestep being the index)
//index 0 is the end of the backward orbit
//note: make sure to free it
int backwardOrbitArraySize, oldOrbitArraySize, prevForwardTime;

const int fitReverseOrbit = 0;

SingleParticleOrbitParams* reverseOrbitParams;


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
                    real dt,
                    real forwardTime)
{
    mw_printf("start pos: x: %f y: %f z: %f\n", pos.x, pos.y, pos.z);
    mw_printf("start vel: x: %f y: %f z: %f\n", vel.x, vel.y, vel.z);
    mw_printf("tstop: %lf dt: %lf\n", tstop, dt);
    mwvector acc, v, x;
    real t;
    real dt_half = dt / 2.0;
    int initialLArrayIndex = tstop/dt;
    backwardOrbitArraySize = ((forwardTime > tstop) ? forwardTime : tstop) / dt;
    oldOrbitArraySize = tstop/dt;
    int LArrayIndex = initialLArrayIndex;
    int barTime;
    mw_printf("Array Length = %d\n", backwardOrbitArraySize);
    
    //allocate backwardsOrbitL
    backwardOrbitPositions = (mwvector*)mwMalloc(backwardOrbitArraySize * sizeof(mwvector));
    reverseOrbitParams = (SingleParticleOrbitParams*)mwMalloc(sizeof(SingleParticleOrbitParams));
    
    //set initial conditions (forward orbit)
    //this is only to fill
    x = pos;
    v = vel;

    //this loop only occurs when forwardTime > tstop
    //for(t = tstop - forwardTime; t <= 0)//for bar time
    //starting at end of previous iteration
    for(t = 0; t <= forwardTime - tstop; t += dt){
        // Update the velocities and positions
        mw_incaddv_s(v, acc, dt_half);
        mw_incaddv_s(x, v, dt); 
        
        // Compute the new acceleration
        acc = nbExtAcceleration(pot, x, t);
        
        mw_incaddv_s(v, acc, dt_half);

        backwardOrbitPositions[LArrayIndex] = x;
        //record the current l coordinate
        LArrayIndex++;
    }

    // Set the initial conditions (reverse orbit)
    x = pos;
    v = vel;
    mw_incnegv(v);
    LArrayIndex = initialLArrayIndex;

    // Get the initial acceleration
    acc = nbExtAcceleration(pot, x, 0);

    for (t = 0; t >= tstop*(-1); t -= dt)
    {
        //barTime = t + (tstop - forwardTime);
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

    mw_printf("end of backward orbit: (%f, %f, %f)\n", x.x, x.y, x.z);
    mw_printf("vel: (%f, %f, %f)\n", v.x, v.y, v.z);
    //set the state variables
    reverseOrbitParams->revOrbitPos = pos;
    reverseOrbitParams->revOrbitVel = vel;
    reverseOrbitParams->revOrbitLMCPos = (mwvector)ZERO_VECTOR;
    reverseOrbitParams->revOrbitLMCVel = (mwvector)ZERO_VECTOR;
    reverseOrbitParams->revOrbitdt = dt;
    reverseOrbitParams->LMCmass = -1;
    reverseOrbitParams->revOrbitTstop = tstop;
    //_st.pot = memcpy(&_st.pot, pot, sizeof(Potential));
    reverseOrbitParams->pot.sphere[0]  = pot->sphere[0];
    reverseOrbitParams->pot.disk  = pot->disk;
    reverseOrbitParams->pot.disk2  = pot->disk2;
    reverseOrbitParams->pot.halo  = pot->halo;
    reverseOrbitParams->previousForwardTime = forwardTime;
    prevForwardTime = forwardTime;
    
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

int getPrevForwardTime(){
    return prevForwardTime;
}

//this one is only based off tstop
int getOldOrbitArraySize(){
    return oldOrbitArraySize;
}

SingleParticleOrbitParams* getOrbitParams(){
    return reverseOrbitParams;
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
    const mwvector backwardOrbitTarget = mw_vec(0,0,0);

    //mwvector mw_x = mw_vec(-35.849678, -57.758362, 51.630098);
    real t, bestTime, bestDist = -1;
    int LArrayIndex, bestTimeStep, timeStep = 0;
    backwardOrbitArraySize = 0;

    // Set the initial conditions
    unsigned int steps = mw_ceil((tstop)/(dt)) + 1;
    unsigned int exSteps = mw_abs(mw_ceil((ftime-tstop)/(dt)) + 1);
    unsigned int i = 0, j = 0, k = 0;
    mwvector acc, v, x, mw_acc, LMC_acc, DF_acc, LMCv, LMCx, tmp;
    mwvector mw_x = mw_vec(0, 0, 0);
    mwvector bacArray[steps + 1];
    mwvector forArray[exSteps + 1];
    real dt_half = dt / 2.0;

    if(!fitReverseOrbit){
        LArrayIndex = tstop / dt;
        backwardOrbitArraySize = LArrayIndex;
        mw_printf("Array Length = %d\n", backwardOrbitArraySize);
        //allocate backwardsOrbitL
        backwardOrbitPositions = (mwvector*)mwMalloc(LArrayIndex * sizeof(mwvector));
    }

    // Check if forward time is larger than backward time. We will need to manually compute additional LMC accelerations in that case.
    if (ftime > tstop) {

        // Set the initial conditions for forward orbit
        x = pos;
        v = vel;
        LMCv = LMCvelocity;
        LMCx = LMCposition;

        // Get the initial acceleration
        mw_acc = plummerAccel(mw_x, LMCx, LMCmass, LMCscale);
        LMC_acc = mw_addv(nbExtAcceleration(pot, LMCx, 0), dynamicalFriction_LMC(pot, LMCx, LMCv, LMCmass, LMCscale, LMCDynaFric));
        acc = nbExtAcceleration(pot, x, 0); //check time argument
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
            LMC_acc = mw_addv(nbExtAcceleration(pot, LMCx, t), dynamicalFriction_LMC(pot, LMCx, LMCv, LMCmass, LMCscale, LMCDynaFric));
            acc = nbExtAcceleration(pot, x, t); //check time arg here
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


    
    //set the state variables
    reverseOrbitParams->revOrbitPos = pos;
    reverseOrbitParams->revOrbitVel = vel;
    reverseOrbitParams->revOrbitLMCPos = *LMCpos;
    reverseOrbitParams->revOrbitLMCVel = *LMCvel;
    reverseOrbitParams->revOrbitdt = dt;
    reverseOrbitParams->LMCmass = LMCmass;
    reverseOrbitParams->revOrbitTstop = tstop;
    //_st.pot = memcpy(&_st.pot, pot, sizeof(Potential));
    reverseOrbitParams->pot.sphere[0]  = pot->sphere[0];
    reverseOrbitParams->pot.disk  = pot->disk;
    reverseOrbitParams->pot.disk2  = pot->disk2;
    reverseOrbitParams->pot.halo  = pot->halo;


    // Get the initial acceleration
    mw_acc = plummerAccel(mw_x, LMCx, LMCmass, LMCscale);
    LMC_acc = nbExtAcceleration(pot, LMCx, 0);
    if (LMCDynaFric) {
        DF_acc = dynamicalFriction_LMC(pot, LMCx, LMCv, LMCmass, LMCscale, TRUE);
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


    if(fitReverseOrbit){
        //I am trying to calculate tstop so that it makes
        //the particle as close to the target value as possible
        //once fully implemented this should be optional!
        //get as close to possible to target backward orbit
        for (t = 0; t >= 1.1*tstop*(-1); t -= dt)
        {
            // Update the velocities and positions
            mw_incaddv_s(v, acc, dt_half);
            mw_incaddv_s(x, v, dt);
            mw_incaddv_s(LMCv, LMC_acc, dt_half);
            mw_incaddv_s(LMCx, LMCv, dt);
            
            // Compute the new acceleration
            mw_acc = pointAccel(mw_x, LMCx, LMCmass);
            LMC_acc = nbExtAcceleration(pot, LMCx, 0);
            acc = nbExtAcceleration(pot, x, t);
            tmp = pointAccel(x, LMCx, LMCmass);
            mw_incaddv(acc, tmp);

            // Shift the body
            mw_incnegv(mw_acc);
            mw_incaddv(LMC_acc, mw_acc);
            mw_incaddv(acc, mw_acc);
            
            mw_incaddv_s(v, acc, dt_half);
            mw_incaddv_s(LMCv, LMC_acc, dt_half);

            if(bestDist == -1 || bestDist > mw_distv(x, backwardOrbitTarget)){
                bestDist = mw_distv(x, backwardOrbitTarget);
                bestTime = t;
                bestTimeStep = timeStep;
            }
            ++timeStep;
        }

        //reset initial conditions
        x = pos;
        v = vel;
        acc = nbExtAcceleration(pot, x, 0);
        
        LMCv = *LMCvel;
        LMCx = *LMCpos;
        mw_incnegv(v);
        mw_incnegv(LMCv);

        // Get the initial acceleration
        mw_acc = pointAccel(mw_x, LMCx, LMCmass);
        LMC_acc = nbExtAcceleration(pot, LMCx, 0);
        acc = nbExtAcceleration(pot, x, 0);
        tmp = pointAccel(x, LMCx, LMCmass);
        mw_incaddv(acc, tmp);

        // Shift the body
        mw_incnegv(mw_acc);
        mw_incaddv(LMC_acc, mw_acc);
        mw_incaddv(acc, mw_acc);

        //allocate backwardsOrbitL
        backwardOrbitPositions = (mwvector*)mwMalloc((bestTimeStep+1) * sizeof(mwvector));
        tstop = (-1) * bestTime;
        LArrayIndex = bestTimeStep;
        backwardOrbitArraySize = LArrayIndex;
        mw_printf("bestTimeStep: %d tstop: %f\nbestDist: %f\n", bestTimeStep, tstop, bestDist);
    }


    for (t = 0; t >= tstop*(-1); t -= dt)
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
        LMC_acc = nbExtAcceleration(pot, LMCx, t);
        if (LMCDynaFric) {
            DF_acc = dynamicalFriction_LMC(pot, LMCx, LMCv, LMCmass, LMCscale, TRUE);
            //mw_printf("DF: [%.15f,%.15f,%.15f]\n",X(DF_acc),Y(DF_acc),Z(DF_acc));
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
        
        mw_incaddv_s(v, acc, dt_half);
        mw_incaddv_s(LMCv, LMC_acc, dt_half);

        //record the current xyz coordinate
        backwardOrbitPositions[LArrayIndex] = x;
        LArrayIndex--;
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

//using a target lambda value and the initial orbit params,
//go backward and then forward to find the point on the single body
//orbit that is closest to the result of the stream
//useful for the first ("calibration") run to get the previous forward time correct
void fitOrbitStart(mwvector* finalPos,
                    mwvector* finalVel,
                    real* finalTime,
                    real* dt,
                    NBodyState* st,
                    NBodyCtx* ctx,
                    real targetLambda,
                    mwvector targetPos
                    )
{	
    const real MAX_DIST_FROM_TARGET = 0.4;
    //only arccept single orbit positions this close
    //to the target pos

	unsigned int steps = (reverseOrbitParams->revOrbitTstop)/ (reverseOrbitParams->revOrbitdt) + 1;
    unsigned int i = 0, j = 0;
    mwvector acc, v, x, mw_acc, LMC_acc, LMCv, LMCx, tmp, bestX, bestV;
    mwvector mw_x = mw_vec(-35.849678, -57.758362, 51.630098);
    mwvector array[steps + 1];
    real t, bestTime, lambda, bestLambda = 10;
    real bestDist = -1;
    real dt_half = reverseOrbitParams->revOrbitdt / 2.0;
    int LArrayIndex, bestTimeStep, timeStep = 0;

    // Set the initial conditions
    x = reverseOrbitParams->revOrbitPos;
    v = reverseOrbitParams->revOrbitVel;
    LMCv = reverseOrbitParams->revOrbitLMCVel;
    LMCx = reverseOrbitParams->revOrbitLMCPos;
    mw_incnegv(v);
    mw_incnegv(LMCv);

    // Get the initial acceleration
    if(reverseOrbitParams->LMCmass != -1){
        mw_acc = pointAccel(mw_x, LMCx, reverseOrbitParams->LMCmass);
        LMC_acc = nbExtAcceleration(&reverseOrbitParams->pot, LMCx, 0);
    }
    acc = nbExtAcceleration(&reverseOrbitParams->pot, x, 0);
    tmp = pointAccel(x, LMCx, reverseOrbitParams->LMCmass);
    mw_incaddv(acc, tmp);


    if(reverseOrbitParams->LMCmass != -1){
        // Shift the body
        mw_incnegv(mw_acc);
        mw_incaddv(LMC_acc, mw_acc);
        mw_incaddv(acc, mw_acc);
    }else
        mw_printf("LMC not in fit orbid\n");

    for (t = 0; t >= 0.03*reverseOrbitParams->revOrbitTstop*(-1); t -= reverseOrbitParams->revOrbitdt)
    {
        // Update the velocities and positions
        mw_incaddv_s(v, acc, dt_half);
        mw_incaddv_s(x, v, reverseOrbitParams->revOrbitdt);

        if(reverseOrbitParams->LMCmass != -1){
            mw_incaddv_s(LMCv, LMC_acc, dt_half);
            mw_incaddv_s(LMCx, LMCv, reverseOrbitParams->revOrbitdt);
        
            // Compute the new acceleration
            mw_acc = pointAccel(mw_x, LMCx, reverseOrbitParams->LMCmass);
            LMC_acc = nbExtAcceleration(&reverseOrbitParams->pot, LMCx, 0);

            tmp = pointAccel(x, LMCx, reverseOrbitParams->LMCmass);
            mw_incaddv(acc, tmp);

            // Shift the body
            mw_incnegv(mw_acc);
            mw_incaddv(LMC_acc, mw_acc);
            mw_incaddv(acc, mw_acc);

            mw_incaddv_s(LMCv, LMC_acc, dt_half);
        }

        acc = nbExtAcceleration(&reverseOrbitParams->pot, x, t);
        
        mw_incaddv_s(v, acc, dt_half);

        lambda = getLambda(x, st);

        real tmpDiff = getAngleDiff(lambda, targetLambda);
        //check to see if this is better than our current bestDist
        if(bestDist == -1 || (bestDist > tmpDiff && mw_distv(targetPos, x) < MAX_DIST_FROM_TARGET)){
            bestDist = tmpDiff;
            bestTime = t;
            bestTimeStep = timeStep;
            bestLambda = lambda;
            bestX = x;
            bestV = v;
            mw_incnegv(bestV);
        }
        ++timeStep;
    }

    //mw_printf("best dist 1: %f\n", bestDist);

    // Reset the initial conditions but dont negate the velocity
    x = reverseOrbitParams->revOrbitPos;
    v = reverseOrbitParams->revOrbitVel;
    LMCv = reverseOrbitParams->revOrbitLMCVel;
    LMCx = reverseOrbitParams->revOrbitLMCPos;

    // Get the initial acceleration
    mw_acc = pointAccel(mw_x, LMCx, reverseOrbitParams->LMCmass);
    LMC_acc = nbExtAcceleration(&reverseOrbitParams->pot, LMCx, 0);
    acc = nbExtAcceleration(&reverseOrbitParams->pot, x, 0);
    tmp = pointAccel(x, LMCx, reverseOrbitParams->LMCmass);
    mw_incaddv(acc, tmp);

    // Shift the body
    mw_incnegv(mw_acc);
    mw_incaddv(LMC_acc, mw_acc);
    mw_incaddv(acc, mw_acc);

    for (t = 0; t <= 0.03*reverseOrbitParams->revOrbitTstop; t += reverseOrbitParams->revOrbitdt)
    {
        // Update the velocities and positions
        mw_incaddv_s(v, acc, dt_half);
        mw_incaddv_s(x, v, reverseOrbitParams->revOrbitdt);

        if(reverseOrbitParams->LMCmass != -1){
            mw_incaddv_s(LMCv, LMC_acc, dt_half);
            mw_incaddv_s(LMCx, LMCv, reverseOrbitParams->revOrbitdt);
        
            // Compute the new acceleration
            mw_acc = pointAccel(mw_x, LMCx, reverseOrbitParams->LMCmass);
            LMC_acc = nbExtAcceleration(&reverseOrbitParams->pot, LMCx, 0);

            tmp = pointAccel(x, LMCx, reverseOrbitParams->LMCmass);
            mw_incaddv(acc, tmp);

            // Shift the body
            mw_incnegv(mw_acc);
            mw_incaddv(LMC_acc, mw_acc);
            mw_incaddv(acc, mw_acc);

            mw_incaddv_s(LMCv, LMC_acc, dt_half);
        }
        
        acc = nbExtAcceleration(&reverseOrbitParams->pot, x, t);
        
        mw_incaddv_s(v, acc, dt_half);

        lambda = getLambda(x, st);

        if(bestDist == -1 || bestDist > getAngleDiff(lambda, targetLambda)&& mw_distv(targetPos, x) < MAX_DIST_FROM_TARGET){
            bestDist = getAngleDiff(lambda, targetLambda);
            bestTime = t;
            bestTimeStep = timeStep;
            bestLambda = lambda;
            bestX = x;
            bestV = v;

        }
        ++timeStep;
    }
    
    //mw_printf("best dist 2: %f\n", bestDist);

    /* Report the final values (don't forget to reverse the velocities) */
    *finalPos = bestX;
    *finalVel = bestV;
    *finalTime = bestTime;
    *dt = reverseOrbitParams->revOrbitdt;
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
