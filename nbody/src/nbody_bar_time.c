#include "milkyway_util.h"
#include "nbody_bar_time.h"

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
        tmpPos = mw_3DRotation(xyPositions[i], st->backwardOrbitRotationAngles);
        st->backwardOrbitAngles[i] = mw_atan(tmpPos.y/tmpPos.x);
        if(tmpPos.x < 0 && tmpPos.y > 0)
            st->backwardOrbitAngles[i] += M_PI;
        if(tmpPos.x < 0 && tmpPos.y < 0)
            st->backwardOrbitAngles[i] += M_PI;
        if(tmpPos.x > 0 && tmpPos.y < 0)
            st->backwardOrbitAngles[i] += 2*M_PI;
    }
}

//calculates the orbit rotation needed to align our coordinates
//with the plane of the star stream
void setBackwardOrbitRotation(NBodyState* st){
    mwvector* xyPositions;
    mwvector rotationAngles, position1, position2;
    int backwardArrayLength = getOrbitArraySize();
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
    while(tmpTime < backwardArrayLength){
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

real meanBodyAngle(Body *bodies, int nbody, NBodyState* st)
{
    mwvector tmpPos;
    real y_part = 0, x_part = 0, angle;
    int i;
    
    for (i = 0; i < nbody; i++)
        {
            tmpPos = mw_3DRotation(bodies[i].bodynode.pos, st->backwardOrbitRotationAngles);
            angle = mw_atan(tmpPos.y/tmpPos.x);
            if(tmpPos.x < 0 && tmpPos.y > 0)
                angle += M_PI;
            if(tmpPos.x < 0 && tmpPos.y < 0)
                angle += M_PI;
            if(tmpPos.x > 0 && tmpPos.y < 0)
                angle += 2*M_PI;
            x_part += cos (angle);
            y_part += sin (angle);
        }
    
    angle = mw_atan2 (y_part / nbody, x_part / nbody);
    return angle;
}

real stdevBodyAngle(Body *bodies, int nbody, real meanAngle, NBodyState* st){
    mwvector tmpPos;
    real sum = 0, angle;
    int i;
    
    for (i = 0; i < nbody; i++)
        {
            tmpPos = mw_3DRotation(bodies[i].bodynode.pos, st->backwardOrbitRotationAngles);
            angle = mw_atan(tmpPos.y/tmpPos.x);
            if(tmpPos.x < 0 && tmpPos.y > 0)
                angle += M_PI;
            if(tmpPos.x < 0 && tmpPos.y < 0)
                angle += M_PI;
            if(tmpPos.x > 0 && tmpPos.y < 0)
                angle += 2*M_PI;
            sum += (angle - meanAngle)*(angle - meanAngle);
        }
    
    return mw_sqrt(sum/(nbody - 1));
}

//creates a histogram of bodies (using backward orbit rotations in st) and 
//returns angle of highest histogram bin. If returnBinCenter returns center coordinate
//in XYZ
real highestHistPeak(Body *bodies, int nbody, NBodyState* st, mwbool returnBinCenter, 
mwvector* histCenter, mwvector* meanBinCenter, mwvector* histCenterVelocity, mwvector* meanBinVelocity){
    real meanAngle;
    meanAngle = meanBodyAngle(bodies, nbody, st);
    if(meanAngle < 0){
        meanAngle += 2*M_PI;
    }
    /*real histSpread = stdevBodyAngle(bodies, nbody, meanAngle, st)*2;
    if(histSpread > 2*M_PI)
        histSpread = 2*M_PI;*/
    real histSpread = 2*M_PI;
    unsigned int* histogram = mwMalloc(st->numBarBins * sizeof(unsigned int));
    mwvector* coordinateTotals;
    mwvector* velocityTotals;
    mwvector tmpPos;
    real angle, angleDiff, binSize = histSpread/st->numBarBins;
    int binNum, highestBinNum=0, highestBinVal=0, i, meanBinNum=(st->numBarBins*(meanAngle))/histSpread;
    int leftSearch = (st->numBarBins*(meanAngle-0.2))/histSpread;
    int rightSearch = (st->numBarBins*(meanAngle+0.2))/histSpread;
    if(returnBinCenter){
        coordinateTotals = (mwvector*)mwMalloc(st->numBarBins * sizeof(mwvector));
        velocityTotals = (mwvector*)mwMalloc(st->numBarBins * sizeof(mwvector));
        for(i = 0; i < st->numBarBins; i++){
            coordinateTotals[i].x = 0;
            coordinateTotals[i].y = 0;
            coordinateTotals[i].z = 0;
            velocityTotals[i].x = 0;
            velocityTotals[i].y = 0;
            velocityTotals[i].z = 0;

        }
    }

    //initialize histogram to 0
    for(i = 0; i < st->numBarBins; i++){
        histogram[i] = 0;
    }
    //fill histogram
    for(i = 0; i < nbody; i++){
        tmpPos = mw_3DRotation(bodies[i].bodynode.pos, st->backwardOrbitRotationAngles);
        angle = mw_atan(tmpPos.y/tmpPos.x);
        if(tmpPos.x < 0 && tmpPos.y > 0)
            angle += M_PI;
        if(tmpPos.x < 0 && tmpPos.y < 0)
            angle += M_PI;
        if(tmpPos.x > 0 && tmpPos.y < 0)
            angle += 2*M_PI;
        
        binNum = (st->numBarBins*(angle))/histSpread;
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

    //find highest bin
    for(i = leftSearch; i < rightSearch; i++){
        if(histogram[i] > highestBinVal){
            highestBinNum = i;
            highestBinVal = histogram[i];
        }
    }
    if(st->step % 100 == 0)
        for(i = 0; i < st->numBarBins; i++){
            mw_printf("bin %d: %d\n", i, histogram[i]);
        }
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
    angle = ((highestBinNum)/(float)st->numBarBins) * histSpread;
    if(angle < 0){
        angle += 2*M_PI;
    }else if(angle > 2*M_PI){
        angle -= 2*M_PI;
    }
    if(st->step % 10 == 0)
        mw_printf("mean angle: %f angle: %f highest bin: %d\n", meanAngle, angle, highestBinNum);
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
    mwvector tmp;
    real oldBackwardOrbitTheta, newBackwardOrbitTheta, oldDiff, newDiff;
    int oldTime, newTime, arraySize;
    arraySize = getOrbitArraySize();
    real streamCenter;

    //if(st->barTimeStep/(float)arraySize > 0.9)
    streamCenter = highestHistPeak(bodies, nbody, st, FALSE, &tmp, &tmp, &tmp, &tmp);
    //else{
    //    streamCenter = meanBodyAngle(bodies, nbody, st);
    //    highestHistPeak(bodies, nbody, st, FALSE, &tmp, &tmp, &tmp);
    //}

    if(st->step % 10 == 0)
        mw_printf("stream center %f:\n", streamCenter);
    oldTime = st->barTimeStep;
    oldBackwardOrbitTheta = st->backwardOrbitAngles[oldTime];
    if(st->barTimeStep + 1 < arraySize){
        newTime = st->barTimeStep + 1;
        newBackwardOrbitTheta = st->backwardOrbitAngles[newTime];
    }else
        return oldTime;


    oldDiff = getAngleDiff(oldBackwardOrbitTheta, streamCenter);
    newDiff = getAngleDiff(newBackwardOrbitTheta, streamCenter);
    while(newDiff < oldDiff){
        if(newTime + 1 < arraySize){
            oldTime = newTime;
            newTime++;
            oldBackwardOrbitTheta = newBackwardOrbitTheta;
            newBackwardOrbitTheta = st->backwardOrbitAngles[newTime];
        }else{
            st->barTimeStep = oldTime;
            return oldTime;
        }
        oldDiff = newDiff;
        newDiff = getAngleDiff(newBackwardOrbitTheta, streamCenter);
    }
    st->barTimeStep = oldTime;
    if(st->step % 100 == 0 && st->step > 10){
        for(int i = oldTime - 10; i < oldTime + 10; i++){
            mw_printf("theta: %f\n", st->backwardOrbitAngles[i]);
        }
    }
    return oldTime;
}

//returns hist center
mwvector getStreamCenter(NBodyState* st, mwvector* meanBinCenter, 
mwvector* histCenterVelocity, mwvector* meanBinVelocity){
    mwvector ans;
    highestHistPeak(st->bodytab, st->nbody, st, TRUE, &ans, meanBinCenter, histCenterVelocity, meanBinVelocity);
    return ans;
}