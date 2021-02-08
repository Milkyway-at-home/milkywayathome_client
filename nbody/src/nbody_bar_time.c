#include "milkyway_util.h"
#include "nbody_bar_time.h"
#include "nbody.h"


void fillBackwardOrbitAngles(NBodyState* st){
    mwvector* xyPositions;
    mwvector tmpPos;
    getBackwardOrbitArray(&xyPositions);
    int backwardArrayLength = getOrbitArraySize();
    int i;
    int totalZ = 0;
    st->backwardOrbitAngles = (real*)mwMalloc(backwardArrayLength * sizeof(real));
    //mw_printf("len: %d\n", backwardArrayLength);
    for(i = 0; i < backwardArrayLength; i++){
        st->backwardOrbitAngles[i] = getLambda(xyPositions[i], st);
    }
}

//calculates the orbit rotation needed to align our coordinates
//with the plane of the star stream
void setBackwardOrbitRotation(NBodyState* st){
    mwvector* xyPositions;
    mwvector rotationAngles, position1, position2;
    int backwardArrayLength = getOrbitArraySize();
    int oldBackwardOrbitLength = getOldOrbitArraySize();
    int time1 = backwardArrayLength/5;
    int tmpTime = 0;
    int time2 = 0;
    real bestAngle = 0;
    real tmpAngle = 0;
    rotationAngles.x = 0;
    rotationAngles.y = 0;
    rotationAngles.z = 0;
    getBackwardOrbitArray(&xyPositions);
    position1 = xyPositions[time1];
    position1.w = 1;
    //iterate through xyPositions to get a position 90 degrees
    //from position 1
    while(tmpTime < oldBackwardOrbitLength){
        position2 = xyPositions[tmpTime];
        position2.w = 1;
        tmpAngle = mw_vecangle(position1, position2);
        if(M_PI_2-tmpAngle < M_PI_2-bestAngle){
            time2 = tmpTime;
            bestAngle = tmpAngle;
        }
        //might as well break once we get really close
        if(M_PI_2 - bestAngle < 0.5)
            break;
        tmpTime += 10;
    }
    position2 = xyPositions[time2];

    //first, set particle 1 to lie on the x axis
    rotationAngles.z = -mw_atan(position1.y/position1.x);
    position1 = mw_rotateZ(position1, rotationAngles.z);
    rotationAngles.y = mw_atan(position1.z/position1.x);
    position1 = mw_rotateY(position1, rotationAngles.y);

    //set particle 2 to lie on the x-y plane
    position2 = mw_rotateZ(position2, rotationAngles.z);
    position2 = mw_rotateY(position2, rotationAngles.y);
    rotationAngles.x = -mw_atan(position2.z/position2.y);
    position2 = mw_rotateX(position2, rotationAngles.x);
    position1 = mw_rotateX(position1, rotationAngles.x);
    
    st->backwardOrbitRotationAngles = rotationAngles;
    st->backwardOrbitArrayLength = backwardArrayLength;
}

//returns the mean of the lambda angles of these bodies
real meanBodyAngle(Body *bodies, int nbody, NBodyState* st)
{
    mwvector tmpPos;
    real y_part = 0, x_part = 0, angle;
    int i;
    
    for (i = 0; i < nbody; i++)
        {
            angle = getLambda(bodies[i].bodynode.pos, st);
            x_part += cos (angle);
            y_part += sin (angle);
        }
    
    angle = mw_atan2 (y_part / nbody, x_part / nbody);
    return angle;
}

//returns the standard deviation of the lambda angles of the bodies
real stdevBodyAngle(Body *bodies, int nbody, real meanAngle, NBodyState* st){
    mwvector tmpPos;
    real sum = 0, angle;
    int i;
    
    for (i = 0; i < nbody; i++)
        {
            angle = getLambda(bodies[i].bodynode.pos, st);
            sum += (angle - meanAngle)*(angle - meanAngle);
        }
    
    return mw_sqrt(sum/(nbody - 1));
}

//performs the rotation specified in st->backwardOrbitRotationAngles
//and returns the lambda valuen of pos in radians
real getLambda(mwvector pos, NBodyState* st){
    mwvector tmpPos = mw_3DRotation(pos, st->backwardOrbitRotationAngles);
    real angle = mw_atan(tmpPos.y/tmpPos.x);
    if(tmpPos.x < 0 && tmpPos.y > 0)
        angle += M_PI;
    if(tmpPos.x < 0 && tmpPos.y < 0)
        angle += M_PI;
    if(tmpPos.x > 0 && tmpPos.y < 0)
        angle += 2*M_PI;
    return angle;
}

//creates a histogram of bodies (using backward orbit rotations in st) and 
//returns angle of highest histogram bin. If returnBinCenter returns center coordinate
//in XYZ
real highestHistPeak(Body *bodies, int nbody, NBodyState* st, mwbool returnBinCenter, 
mwvector* histCenter, mwvector* meanBinCenter, mwvector* histCenterVelocity, mwvector* meanBinVelocity){
    //this value is in radians and determines how far from the meanAngle we
    //search for the highest histogram peak. Useful for when stream splits into
    //multiple peaks
    const real SEARCH_AREA = 0.075;
    const int MAX_BINS = 5000;

    real meanAngle;
    meanAngle = meanBodyAngle(bodies, nbody, st);
    if(meanAngle < 0){
        meanAngle += 2*M_PI;
    }
    real angleStdev = stdevBodyAngle(bodies, nbody, meanAngle, st)*2;
    if(angleStdev > 2*M_PI)
        angleStdev = 2*M_PI;
    real histSpread = 2*M_PI;

    //cap numBins at MAX_BINS
    int numBins = (st->numBarBins * (histSpread/angleStdev) <= MAX_BINS) ? 
        st->numBarBins * (histSpread/angleStdev) : MAX_BINS;
    unsigned int* histogram = (unsigned int*)mwMalloc(numBins * sizeof(unsigned int));
    mwvector* coordinateTotals;
    real* lambdaTotals;
    mwvector* velocityTotals;
    mwvector tmpPos;
    real angle, angleDiff;
    real binSize = histSpread/numBins;
    int binNum;
    int highestBinNum=0;
    int highestBinVal=0;
    int i;
    int meanBinNum = (numBins*(meanAngle))/histSpread;
    int leftSearch = (int)((numBins*(meanAngle-SEARCH_AREA))/histSpread);
    int rightSearch = (int)((numBins*(meanAngle+SEARCH_AREA))/histSpread);
    if(returnBinCenter){
        coordinateTotals = (mwvector*)mwMalloc(numBins * sizeof(mwvector));
        velocityTotals = (mwvector*)mwMalloc(numBins * sizeof(mwvector));
        lambdaTotals = (mwvector*)mwMalloc(numBins * sizeof(real));
        for(i = 0; i < numBins; i++){
            coordinateTotals[i].x = 0;
            coordinateTotals[i].y = 0;
            coordinateTotals[i].z = 0;
            velocityTotals[i].x = 0;
            velocityTotals[i].y = 0;
            velocityTotals[i].z = 0;
        }
    }

    //initialize histogram to 0
    for(i = 0; i < numBins; i++){
        histogram[i] = 0;
    }
    //fill histogram
    for(i = 0; i < nbody; i++){
        angle = getLambda(bodies[i].bodynode.pos, st);
        
        binNum = (int)((numBins*(angle))/histSpread);
        histogram[binNum] += 1;
        if(returnBinCenter){
            coordinateTotals[binNum].x += bodies[i].bodynode.pos.x;
            coordinateTotals[binNum].y += bodies[i].bodynode.pos.y;
            coordinateTotals[binNum].z += bodies[i].bodynode.pos.z;
            velocityTotals[binNum].x += bodies[i].vel.x;
            velocityTotals[binNum].y += bodies[i].vel.y;
            velocityTotals[binNum].z += bodies[i].vel.z;
        }
    }

    //find highest bin *could be a way to go out of bounds
    for(i = leftSearch; i < rightSearch; i++){
        if((int)highestBinVal < (int)histogram[i]){
            highestBinNum = i;
            highestBinVal = histogram[i];
        }
    }
    /*if(st->step % 100 == 0)
        mw_printf("numBins: %d\n", numBins);
        for(i = 0; i < numBins; i++){
            mw_printf("bin %d: %d\n", i, histogram[i]);
        }*/
    //return average of coordinates in given bin
    if(returnBinCenter){
        histCenter->x = coordinateTotals[highestBinNum].x/highestBinVal;
        histCenter->y = coordinateTotals[highestBinNum].y/highestBinVal;
        histCenter->z = coordinateTotals[highestBinNum].z/highestBinVal;

        meanBinCenter->x = coordinateTotals[meanBinNum].x/histogram[meanBinNum];
        meanBinCenter->y = coordinateTotals[meanBinNum].y/histogram[meanBinNum];
        meanBinCenter->z = coordinateTotals[meanBinNum].z/histogram[meanBinNum];
        
        histCenterVelocity->x = velocityTotals[highestBinNum].x/highestBinVal;
        histCenterVelocity->y = velocityTotals[highestBinNum].y/highestBinVal;
        histCenterVelocity->z = velocityTotals[highestBinNum].z/highestBinVal;

        meanBinVelocity->x = velocityTotals[meanBinNum].x/histogram[meanBinNum];
        meanBinVelocity->y = velocityTotals[meanBinNum].y/histogram[meanBinNum];
        meanBinVelocity->z = velocityTotals[meanBinNum].z/histogram[meanBinNum];

    }

    //calculate angle to return 
    angle = ((highestBinNum)/(float)numBins) * histSpread;
    if(angle < 0){
        angle += 2*M_PI;
    }else if(angle > 2*M_PI){
        angle -= 2*M_PI;
    }
    /*if(st->step % 10 == 0){
        mw_printf("mean angle: %f angle: %f highest bin: %d\n", meanAngle, angle, highestBinNum);
        mw_printf("highest binval: %d LS: %d RS: %d Bodies: %d\n", highestBinVal, leftSearch, rightSearch, nbody);
    }*/
    
    //freeing here causes the program to crash. Not sure why.
    /*mwFreeA(histogram);
    if(returnBinCenter){
        mwFreeA(coordinateTotals);
        mwFreeA(velocityTotals);
        mwFreeA(lambdaTotals);
    }*/
    return angle;
}

real getAngleDiff(real a1, real a2){
    return fabs(M_PI - fabs(fabs(a1 - a2) - M_PI));
}

real getAngleDiffDegrees(real a1, real a2){
    return fabs(180 - fabs(fabs(a1 - a2) - 180));
}

mwbool angleIsBetween(real start, real end, real mid) {     
    end = (end - start) < 0.0 ? end - start + 360.0 : end - start;    
    mid = (mid - start) < 0.0 ? mid - start + 360.0 : mid - start; 
    return (mid < end); 
}

//also updates the state with the new bar timestep
int getBarTime(Body* bodies, int nbody, NBodyState* st, NBodyCtx* ctx){
    int MAX_STEP = 50;
    mwvector tmp;
    real oldBackwardOrbitTheta, newBackwardOrbitTheta, oldDiff, newDiff;
    int oldTime, newTime, arraySize, defaultReturnTime;
    arraySize = getOrbitArraySize();
    real streamCenter;

    //uncomment this if/else to use meanBodyAngle for the majority
    //of the run. Could be useful if the stream has 2+ histogram peaks
    //at times.
    //if(st->barTimeStep/(float)arraySize > 0.9)
    streamCenter = highestHistPeak(bodies, nbody, st, FALSE, &tmp, &tmp, &tmp, &tmp);
    //else{
        //streamCenter = meanBodyAngle(bodies, nbody, st);
    //    highestHistPeak(bodies, nbody, st, FALSE, &tmp, &tmp, &tmp);
    //}
    

    /*if(st->step % 10 == 0)
        mw_printf("stream center %f:\n", streamCenter);*/
    oldTime = st->lastFittedBarTimeStep;
    oldBackwardOrbitTheta = st->backwardOrbitAngles[oldTime];
    defaultReturnTime = st->barTimeStep + 10;

    /*if(st->step % 10 == 0)
        mw_printf("theta %f:\n", oldBackwardOrbitTheta);*/

    if(st->barTimeStep + 1 < arraySize){
        newTime = st->barTimeStep + 1;
        newBackwardOrbitTheta = st->backwardOrbitAngles[newTime];
    }else{
        st->barTimeStep = defaultReturnTime;
        mw_printf("array bounds exceeded\n");
        return defaultReturnTime;
    }


    oldDiff = getAngleDiff(oldBackwardOrbitTheta, streamCenter);
    newDiff = getAngleDiff(newBackwardOrbitTheta, streamCenter);
    if(oldDiff <= newDiff){
        st->barTimeStep = defaultReturnTime;
        return defaultReturnTime;
    }
    while(newDiff < oldDiff){
        if(newTime + 1 < arraySize){
            oldTime = newTime;
            newTime++;
            oldBackwardOrbitTheta = newBackwardOrbitTheta;
            newBackwardOrbitTheta = st->backwardOrbitAngles[newTime];
        }else{//if we have gone out of array bounds
            st->barTimeStep = oldTime;
            return oldTime;
        }
        oldDiff = newDiff;
        newDiff = getAngleDiff(newBackwardOrbitTheta, streamCenter);
    }
    if(oldTime - defaultReturnTime > MAX_STEP){
        st->barTimeStep = defaultReturnTime;
        st->lastFittedBarTimeStep = defaultReturnTime;
        return defaultReturnTime;
    }
    st->barTimeStep = oldTime;
    st->lastFittedBarTimeStep = oldTime;
    /*if(st->step % 100 == 0 && st->step > 10){
        for(int i = oldTime - 10; i < oldTime + 10; i++){
            mw_printf("theta: %f\n", st->backwardOrbitAngles[i]);
        }
    }*/
    return oldTime;
}

//returns hist center, calculates and printfs a new starting position and back time
//that is on the original single particle orbit. I'm planning on deprecating this
mwvector getStreamCenter(NBodyState* st, NBodyCtx* ctx, mwvector* meanBinCenter,
mwvector* histCenterVelocity, mwvector* meanBinVelocity){
    mwvector histCenter, finalPos, finalVel;
    real finalTime;
    real streamAngle = highestHistPeak(st->bodytab, st->nbody, st, TRUE, &histCenter, meanBinCenter, histCenterVelocity, meanBinVelocity);
    real dt;
    //This function did a single body orbit to find out how far off the forward time was
    //from the backtime. I'm replacing this with simply checking how far off the bar time
    //is from 0.
    fitOrbitStart(&finalPos, &finalVel, &finalTime, &dt, st, ctx, streamAngle, histCenter);
    mw_printf("final pos (on single orbit): (%f, %f, %f)\n", finalPos.x, finalPos.y, finalPos.z);
    mw_printf("final vel: (%f, %f, %f)\n", finalVel.x, finalVel.y, finalVel.z);
    mw_printf("final time: %f\n", finalTime);
    mw_printf("lastFitted: %d\n", st->lastFittedBarTimeStep);
    mw_printf("dt: %f\n", dt);
    return histCenter;
}