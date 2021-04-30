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

#ifndef _QUATMATH_HPP_
#define _QUATMATH_HPP_

#include <cmath>
#include "vectmath.hpp"

using namespace std;


class Quaternion
{

public:

    float x;
    float y;
    float z;
    float w;

    inline Quaternion( float nx, float ny, float nz, float nw )
    {
        x = nx;
        y = ny;
        z = nz;
        w = nw;
    }

    inline Quaternion()
    {
        x = 0.;
        y = 0.;
        z = 1.;
        w = 0.;
    }

    inline Quaternion( const Vector3d v, float angle )
    {
        setVectorAngle(v, angle);
    }

    inline void set( float nx, float ny, float nz, float nw )
    {
        x = nx;
        y = ny;
        z = nz;
        w = nw;
    }

    inline Quaternion& operator=( const Quaternion& q )
    {
        x = q.x;
        y = q.y;
        z = q.z;
        w = q.w;

        return *this;
    }

    inline Quaternion& operator*=( const Quaternion& q )
    {
        Quaternion l = *this;
        x = l.w*q.x + l.x*q.w + l.y*q.z - l.z*q.y;
        y = l.w*q.y + l.y*q.w + l.z*q.x - l.x*q.z;
        z = l.w*q.z + l.z*q.w + l.x*q.y - l.y*q.x;
        w = l.w*q.w - l.x*q.x - l.y*q.y - l.z*q.z;
        normalize();
        return *this;
    }

    inline Quaternion& operator+=( const Quaternion& q )
    {
        x += q.x;
        y += q.y;
        z += q.z;
        w += q.w;
        return *this;
    }

    inline Quaternion& operator-=( const Quaternion& q )
    {
        x -= q.x;
        y -= q.y;
        z -= q.z;
        w -= q.w;
        return *this;
    }

    inline Quaternion conjugate()
    {
        return Quaternion(-x, -y, -z, w);
    }

    inline void normalize()
    {
        float d = sqrt(x*x+y*y+z*z+w*w);
        x /= d;
        y /= d;
        z /= d;
        w /= d;
    }

    inline void setVectorAngle( const Vector3d v, float angle )
    {
        Vector3d nv = v;
        nv.normalize();
        float sina = sin(angle/2.);
        float cosa = cos(angle/2.);
        x = nv.x*sina;
        y = nv.y*sina;
        z = nv.z*sina;
        w = cosa;
        normalize();
    }

    inline void getVectorAngle( Vector3d &v, float &angle )
    {
        angle = acos(w)*2.;
        float sina = sqrt(1.-w*w);
        v.x = x/sina;
        v.y = y/sina;
        v.z = z/sina;
    }

    inline void getRotMat( float m[3][3] )
    {
        float yy = y*y;
        float zz = z*z;
        m[0][0] = 1.-2.*(yy+zz);
        float xy = x*y;
        float zw = z*w;
        m[0][1] = 2.*(xy-zw);
        float xz = x*z;
        float yw = y*w;
        m[0][2] = 2.*(xz+yw);
        m[1][0] = 2.*(xy+zw);
        float xx = x*x;
        m[1][1] = 1.-2.*(xx+zz);
        float yz = y*z;
        float xw = x*w;
        m[1][2] = 2.*(yz-xw);
        m[2][0] = 2.*(xz-yw);
        m[2][1] = 2.*(yz+xw);
        m[2][2] = 1.-2.*(xx+yy);
/*cout << "***SMAT***\n";
cout << m[0][0] << "\t" << m[0][1] << "\t" << m[0][2] << endl;
cout << m[1][0] << "\t" << m[1][1] << "\t" << m[1][2] << endl;
cout << m[2][0] << "\t" << m[2][1] << "\t" << m[2][2] << endl;
cout << "***EMAT***\n\n";
*/
    }

};


#endif /* _QUATMATH_HPP_ */
