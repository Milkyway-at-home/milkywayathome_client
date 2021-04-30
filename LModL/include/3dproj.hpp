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

#ifndef _3DPROJ_HPP_
#define _3DPROJ_HPP_

#include <cmath>
#include "drawcore.hpp"
#include "vectmath.hpp"
#include "trigtbl.hpp"

using namespace std;


class AccelTimer
{

    float tStart, cHalf, cDuration, totChange, cAccel, hPos, hVel;
    enum AccelType {MAX_ACCEL, CONST_VEL} accelType;

public:

    AccelTimer();

    void set( float duration, float totalChange, AccelType type = MAX_ACCEL );

    bool getPosition( float &position );

};

struct CameraStore
{
    Vector3d position;
    Vector3d eye;
    Vector3d up;
};

class Camera
{

    int xSize, ySize, xCenter, yCenter;
    float zMult, accInc, panSpeed, rotSpeed, revSpeed;
    float accel[9];
    float fps;
    Vector3d position;
    Vector3d momentum;
    Vector3d vEye, vUp, vRight;
    CameraStore camStore, tCam;
    float rMat[3][3];

    enum { ROTATE, MOVE, REVOLVE, NONE } cType;
    Vector3d rotV;
    Vector3d finalV;

    Vector3d moveV;
    Vector3d finalPos;

    AccelTimer cAccel;

    void setVUp();

    void setVRight();

    void setVEye();

    void setMatrix();

public:

    Camera( int xSize, int ySize, float fps = 30. );

    void setFps( float fps );

    void setPanSpeed( float radiansPerSecond = TRIG_2PI/25. );

    void setRotateSpeed( float radiansPerSecond = TRIG_2PI/15. );

    void setRevolveSpeed( float radiansPerSecondUnit = TRIG_2PI/10. );

    void setAcceleration( float unitsPerSecond = .25, float multiplier = 10. );

    void setSize( int xSize, int ySize );

    void setZoom( float viewAngle = TRIG_2PI/16. );

    void setFocusPosition( float distance, float longitude = TRIG_2PI/4., float latitude = 0., float rotation = 0. );

    void setFocusPosition( float distance, const Vector3d& focusVector, float rotation = 0. );

    void setFocusDirection( const Vector3d focusVector, float seconds = 1. );

    void setFocusPoint( const Vector3d &focusVector, float seconds = 1. );

    void moveToPoint( const Vector3d &newPosition, float finalDist = 0., float seconds = 1. );

    bool getAutoStatus();

    Vector3d getFocusPoint( float distance, float longitude = TRIG_2PI/4, float latitude = 0. );

    Vector3d getFocusPoint( float distance, const Vector3d &focusVector );

    void storeCamera( CameraStore &camStore );

    void restoreCamera( CameraStore &camStore );

    void getCameraDebug( Vector3d& eye, Vector3d& up, Vector3d& right, float mult );

    void rotate( float angle );

    void panHorizontal( float angle );

    void panVertical( float angle );

    void revHorizontal( float angle );

    void revVertical( float angle );

    void rotateCW();

    void rotateCCW();

    void panUp();

    void panDown();

    void panRight();

    void panLeft();

    void revolveUp();

    void revolveDown();

    void revolveRight();

    void revolveLeft();

    void accelerateUp( float acc );

    void accelerateDown( float acc );

    void accelerateRight( float acc );

    void accelerateLeft( float acc );

    void setPosition( Vector3d position );

    void stop();

    void accelerate( float acc );

    void accelerateVertical( float acc );

    void accelerateHorizontal( float acc );

    void step();

    bool pollEvent( SDL_Event *event );
        // This should be called at regular intervals to avoid key-response lag

    inline void get3dTransform( const Vector3d &p, Vector3d &map ) const;

    inline bool get3dProjection( Vector3d &map ) const;

    inline bool getCameraProjection( const Vector3d &p, Vector3d &map ) const;

    inline void getDisplayOffset( Vector3d &map ) const;

};


inline void Camera::get3dTransform( const Vector3d &p, Vector3d &map ) const
    // Returns true if xm, ym, and lum represent the screen coordinates of the point for the given parameters
    // Returns false if the point is not visible on the screen
{

    Vector3d vt = p;
    vt -= position;

    map.x = rMat[0][0]*vt.x + rMat[0][1]*vt.y + rMat[0][2]*vt.z;
    map.y = rMat[1][0]*vt.x + rMat[1][1]*vt.y + rMat[1][2]*vt.z;
    map.z = rMat[2][0]*vt.x + rMat[2][1]*vt.y + rMat[2][2]*vt.z;

}

inline bool Camera::get3dProjection( Vector3d &map ) const
    // Returns true if xm, ym, and lum represent the screen coordinates of the point for the given parameters
    // Returns false if the point is not visible on the screen
{

    /// TODO ///  Clip x and y field as well if it enhances speed - preferably before cos/sin calculation
    // Clipping
    bool r;
    if( map.z==0 ) {
        map.x *= INFINITY;
        map.y *= INFINITY;
        r = false;
    }
    else {
        if( map.z>0 )
            r = true;
        else
            r = false;
        float t = zMult/map.z;
        map.x *= t;
        map.y *= t;
    }

    return r;

}

inline bool Camera::getCameraProjection( const Vector3d &p, Vector3d &map ) const
    // Returns true if xm, ym, and lum represent the screen coordinates of the point for the given parameters
    // Returns false if the point is not visible on the screen
{
    get3dTransform(p, map);
    return get3dProjection(map);
}

inline void Camera::getDisplayOffset( Vector3d &map ) const
    // Returns true if xm, ym, and lum represent the screen coordinates of the point for the given parameters
    // Returns false if the point is not visible on the screen
{
    map.x += xCenter;
    map.y = yCenter - map.y;
}


#endif /* _3DPROJ_HPP_ */
