/*****************************************************************************
 *                                                                           *
 *  Copyright (C) 2010 Shane Reilly, Heidi Newberg, Malik Magdon-Ismail,     *
 *  Carlos Varela, Boleslaw Szymanski, and Rensselaer Polytechnic Institute  *
 *                                                                           *
 *  This file is part of the Light Modeling Library (LModL).                 *
 *                                                                           *
 *  This library is free software: you can redistribute it and/or modify     *
 *  it under the terms of the GNU General Public License as published by     *
 *  the Free Software Foundation, either version 3 of the License, or        *
 *  (at your option) any later version.                                      *
 *                                                                           *
 *  This library is distributed in the hope that it will be useful,          *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of           *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the             *
 *  GNU General Public License for more details.                             *
 *                                                                           *
 *  You should have received a copy of the GNU General Public License        *
 *  along with this library. If not, see <http://www.gnu.org/licenses/>.     *
 *                                                                           *
 *  Shane Reilly                                                             *
 *  reills2@cs.rpi.edu                                                       *
 *                                                                           *
 *****************************************************************************/

#include "3dproj.hpp"


AccelTimer::AccelTimer()
{
    tStart = 0.;
    cHalf = 0.;
    cDuration = 0.;
    totChange = 0.;
    cAccel = 0.;
    hPos = 0.;
    hVel = 0.;
    accelType = MAX_ACCEL;
}

void AccelTimer::set( float duration, float totalChange, AccelType type )
{
    cDuration = duration;
    tStart = SDL_GetTicks();
    cHalf = cDuration/2.;
    totChange = totalChange;
    accelType = type;
    if( accelType==MAX_ACCEL ) {
        cAccel = totChange/(cHalf*cHalf);
        hPos = cAccel/2.*cHalf*cHalf;
        hVel = cHalf*cAccel;
    }

}

bool AccelTimer::getPosition( float &position )
{

    position = 0.;  // Avoids compiler warnings
    float t = SDL_GetTicks()-tStart;

    // Check if time is up
    if( t>=cDuration ) {
        position = totChange;
        return false;
    }

    switch( accelType ){

    case MAX_ACCEL:

        if( t<cHalf ) {
            // Accelerating
            position = cAccel/2.*t*t;
        }
        else {
            // Decelerating
            position = hPos + hVel*(t-cHalf) - cAccel/2.*(t-cHalf)*(t-cHalf);
        }
        break;

    case CONST_VEL:
        position = t/cDuration*totChange;
        break;

    }

    return true;

}

void Camera::setVUp()
{
    vUp = vRight;
    vUp *= vEye;
    vUp.normalize();
}

void Camera::setVRight()
{
    vRight = vEye;
    vRight *= vUp;
    vRight.normalize();
}

void Camera::setVEye()
{
    vEye = vUp;
    vEye *= vRight;
    vEye.normalize();
}

void Camera::setMatrix()
{
    // X transform component
    rMat[0][0] = vRight.x;
    rMat[0][1] = vRight.y;
    rMat[0][2] = vRight.z;
    // Y component
    rMat[1][0] = vUp.x;
    rMat[1][1] = vUp.y;
    rMat[1][2] = vUp.z;
    // Z component
    rMat[2][0] = vEye.x;
    rMat[2][1] = vEye.y;
    rMat[2][2] = vEye.z;
}

Camera::Camera( int xSize, int ySize, float fps )
{
    setSize(xSize, ySize);
    setFps(fps);
    setPosition(Vector3d(0., 0., 0.));
    setFocusPosition(0.);
    storeCamera(camStore);
}

void Camera::setFps( float fps )
{
    this->fps = fps;
    setPanSpeed();
    setRotateSpeed();
    setRevolveSpeed();
    setAcceleration();
}

void Camera::setPanSpeed( float radiansPerSecond )
{
    panSpeed = radiansPerSecond/fps;
}

void Camera::setRotateSpeed( float radiansPerSecond )
{
    rotSpeed = radiansPerSecond/fps;
}

void Camera::setRevolveSpeed( float radiansPerSecondUnit )
{
    revSpeed = radiansPerSecondUnit/fps;
}

void Camera::setAcceleration( float unitsPerSecond, float multiplier )
{
    float unitsPerFrame = unitsPerSecond/fps;
    for( int i = 0; i<(int)(sizeof(accel)/sizeof(*accel)); i++ ) {
        accel[i] = unitsPerFrame;
        unitsPerFrame *= multiplier;
    }
    accInc = accel[0];
}

void Camera::setSize( int xSize, int ySize )
{
    this->xSize = xSize;
    this->ySize = ySize;
    xCenter = float(xSize)/2.;
    yCenter = float(ySize)/2.;
    setZoom();
}

void Camera::setZoom( float viewAngle )
{
    /// TODO /// Determine if viewAngle uses an ideal default
    /// TODO /// float-check this equation;
    zMult = xSize/tan(viewAngle);
}

void Camera::setFocusPosition( float distance, float longitude, float latitude, float rotation )
{
    // Set Camera to l=0, b=0, rotation=0
    vEye.set(1., 0., 0.);
    vUp.set(0., 0., 1.);
    vRight.set(0., -1., 0.);

    // Change values according to parameters
    panHorizontal(longitude*TRIG_DEG_TO_RAD);
    panVertical(latitude*TRIG_DEG_TO_RAD);
    rotate(rotation*TRIG_DEG_TO_RAD);

    // Set position relative to distance
    position = vEye*-distance;
    position.x -= 8;

    stop();
    setMatrix();
}

void Camera::setFocusPosition( float distance, const Vector3d& focusVector, float rotation )
{
    // Set Camera eye to specified direction
    vEye = focusVector;
    vEye.normalize();

    // Determine the appropriate vectors representing 'up' and 'right'
    if( vEye.x==0. && vEye.y==0. )
        vUp.set(0, -vEye.z, 0.);
    else
        vUp.set(0., 0., 1.);

    setVRight();
    setVUp();

    // Rotate
    rotate(rotation*TRIG_DEG_TO_RAD);

    // Set position relative to distance
    position = vEye*-distance;
    position.x -= 8;

    stop();
    setMatrix();
}

void Camera::setFocusDirection( Vector3d focusVector, float seconds )
{

    focusVector.normalize();

    // Check for indeterminate focus
    if( focusVector!=focusVector )
        return;

    // Set orientation by rotating view directly toward destination
    float vAngle = acos(vEye.dotProduct(focusVector));
    if( vAngle<.001 ) {
        vEye = focusVector;
        setVRight();
        setVUp();
        setMatrix();
        return;
    }

    if( vAngle>3.141 && vAngle<3.142 )
        rotV = vUp;
    else {
        rotV = vEye * focusVector;
        rotV.normalize();
    }

    finalV = focusVector;
    cType = ROTATE;
    storeCamera(tCam);
    cAccel.set(seconds*1000., vAngle);

    step();

}

void Camera::setFocusPoint( const Vector3d &focusCoordinate, float seconds )
{

    Vector3d direction = focusCoordinate-position;
    setFocusDirection(direction, seconds);

}

void Camera::moveToPoint( const Vector3d &newPosition, float finalDist, float seconds )
{

    if( position==newPosition )
        return;
    moveV = newPosition-position;
    float distance = sqrt(moveV.dotProduct(moveV)) - finalDist;
    moveV.normalize();
    finalPos = newPosition-moveV*finalDist;

    cType = MOVE;
    storeCamera(tCam);
    cAccel.set(seconds*1000., distance);

    step();

}

bool Camera::getAutoStatus()
{
    return cType!=NONE;
}

Vector3d Camera::getFocusPoint( float distance, float longitude, float latitude )
{

    // Get point coordinates
    storeCamera(tCam);
    setFocusPosition(-distance, longitude, latitude, 0.);
    Vector3d focus = position - tCam.position;
    if( focus==Vector3d(0., 0., 0.) ) {
        restoreCamera(tCam);
        return tCam.eye;
    }
    focus.normalize();
    restoreCamera(tCam);
    return focus;

}

Vector3d Camera::getFocusPoint( float distance, const Vector3d &focusVector )
{

    // Get point coordinates
    storeCamera(tCam);
    setFocusPosition(-distance, focusVector, 0.);
    Vector3d focus = position - tCam.position;
    if( focus==Vector3d(0., 0., 0.) ) {
        restoreCamera(tCam);
        return tCam.eye;
    }
    focus.normalize();
    restoreCamera(tCam);

    return focus;

}

void Camera::storeCamera( CameraStore &camStore )
{
    camStore.position = position;
    camStore.eye = vEye;
    camStore.up = vUp;
}

void Camera::restoreCamera( CameraStore &camStore )
{
    position = camStore.position;
    vEye = camStore.eye;
    vUp = camStore.up;
    setVRight();
    setMatrix();
}

void Camera::getCameraDebug( Vector3d& eye, Vector3d& up, Vector3d& right, float mult )
{
    eye.set(mult*vEye.x, mult*vEye.y, mult*vEye.z);
    up.set(mult*vUp.x, mult*vUp.y, mult*vUp.z);
    right.set(mult*vRight.x, mult*vRight.y, mult*vRight.z);
}

void Camera::rotate( float angle )
{
    vUp.rotate(vEye, 3.*angle);
    vUp.normalize();
    setVRight();
    setMatrix();
}

void Camera::panHorizontal( float angle )
{
    vEye.rotate(vUp, angle);
    vEye.normalize();
    setVRight();
    setMatrix();
}

void Camera::panVertical( float angle )
{
    vEye.rotate(vRight, angle);
    vEye.normalize();
    setVUp();
    setMatrix();
}

void Camera::revHorizontal( float angle )
{
    position.rotate(vUp, angle);
    momentum.rotate(vUp, angle);
    vEye.rotate(vUp, angle);
    vEye.normalize();
    setVRight();
    setMatrix();
}

void Camera::revVertical( float angle )
{
    position.rotate(vRight, angle);
    momentum.rotate(vRight, angle);
    vEye.rotate(vRight, angle);
    vEye.normalize();
    setVUp();
    setMatrix();
}

void Camera::rotateCW()
{
    rotate(-rotSpeed);
}

void Camera::rotateCCW()
{
    rotate(rotSpeed);
}

void Camera::panUp()
{
    panVertical(-panSpeed);
}

void Camera::panDown()
{
    panVertical(panSpeed);
}

void Camera::panRight()
{
    panHorizontal(-panSpeed);
}

void Camera::panLeft()
{
    panHorizontal(panSpeed);
}

void Camera::revolveUp()
{
    revVertical(revSpeed);
}

void Camera::revolveDown()
{
    revVertical(-revSpeed);
}

void Camera::revolveRight()
{
    revHorizontal(revSpeed);
}

void Camera::revolveLeft()
{
    revHorizontal(-revSpeed);
}

void Camera::accelerateUp( float acc )
{
    accelerateVertical(acc);
}

void Camera::accelerateDown( float acc )
{
    accelerateVertical(-acc);
}

void Camera::accelerateRight( float acc )
{
    accelerateHorizontal(acc);
}

void Camera::accelerateLeft( float acc )
{
    accelerateHorizontal(-acc);
}

void Camera::setPosition( Vector3d position )
{
    this->position = position;
    stop();
}

void Camera::stop()
{
    momentum.set(0., 0., 0.);
    cType = NONE;
}

void Camera::accelerate( float acc )
{
    momentum += vEye*acc;
}

void Camera::accelerateVertical( float acc )
{
    momentum += vUp*acc;
}

void Camera::accelerateHorizontal( float acc )
{
    momentum += vRight*acc;
}

void Camera::step()
{

    switch( cType ) {

    case ROTATE:

        // Apply angle change with constant acceleration/deceleration
        float cAngle;
        if( !cAccel.getPosition(cAngle) ) {
            stop();
            // Ensure final direction is accurate
            vEye = finalV;
            setVUp();
            setVRight();
            setMatrix();
        }
        else {
            restoreCamera(tCam);
            vEye.rotate(rotV, cAngle);
            vEye.normalize();
            vUp.rotate(rotV, cAngle);
            vUp.normalize();
            setVRight();
            setMatrix();
        }

        break;

    case REVOLVE:
        /// TODO /// STUB
        break;

    case MOVE:
        float distance;
        if( !cAccel.getPosition(distance) ) {
            stop();
            position = finalPos;
        }
        else
            position = tCam.position + moveV*distance;
        break;

    default:
        position += momentum;
        break;

    }

}

bool Camera::pollEvent( SDL_Event *event )
    // This should be called at regular intervals to avoid key-response lag
{

    if( !SDL_PollEvent(event) ) {

        if( event->type==SDL_QUIT )
            exit(0);
        if( event->type!=SDL_KEYDOWN )
        {
            step();
            return false;
        }

    }

    switch( event->key.keysym.sym ) {

    case SDLK_LEFT:
        if( event->key.keysym.mod & KMOD_SHIFT )
            revolveLeft();
        else if( event->key.keysym.mod & KMOD_CTRL )
            rotateCW();
        else if( event->key.keysym.mod & KMOD_ALT )
            accelerateLeft(accInc);
        else
            panLeft();
        break;
    case SDLK_RIGHT:
        if( event->key.keysym.mod & KMOD_SHIFT )
            revolveRight();
        else if( event->key.keysym.mod & KMOD_CTRL )
            rotateCCW();
        else if( event->key.keysym.mod & KMOD_ALT )
            accelerateRight(accInc);
        else
            panRight();
        break;
    case SDLK_UP:
        if( event->key.keysym.mod & KMOD_SHIFT )
            revolveUp();
        else if( event->key.keysym.mod & KMOD_ALT )
            accelerateUp(accInc);
        else
            panUp();
        break;
    case SDLK_DOWN:
        if( event->key.keysym.mod & KMOD_SHIFT )
            revolveDown();
        else if( event->key.keysym.mod & KMOD_ALT )
            accelerateDown(accInc);
        else
            panDown();
        break;
	case SDLK_f:
	case SDLK_PERIOD:
        accelerate(accInc);
        break;
	case SDLK_b:
    case SDLK_COMMA:
        accelerate(-accInc);
        break;
    case SDLK_SPACE:
        stop();
        break;
    case SDLK_F1:
        accInc = accel[0];
        break;
    case SDLK_F2:
        accInc = accel[1];
        break;
    case SDLK_F3:
        accInc = accel[2];
        break;
    case SDLK_F4:
        accInc = accel[3];
        break;
    case SDLK_F5:
        accInc = accel[4];
        break;
    case SDLK_F6:
        accInc = accel[5];
        break;
    case SDLK_F7:
        accInc = accel[6];
        break;
    case SDLK_F8:
        accInc = accel[7];
        break;
    case SDLK_F9:
        accInc = accel[8];
        break;
    case SDLK_ESCAPE:
        exit(0);
        break;
    case SDLK_TAB:
        if( event->key.keysym.mod & KMOD_SHIFT )
            storeCamera(camStore);
        else
            restoreCamera(camStore);
            stop();
        break;

    default:
        break;

    }

    step();
    return true;

}
