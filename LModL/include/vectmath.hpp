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

#ifndef _VECTMATH_HPP_
#define _VECTMATH_HPP_

#include <iostream>
#include <cmath>

using namespace std;


class Vector3d
{

public:

    float x;
    float y;
    float z;

    inline Vector3d()
    {
        x = 0.;
        y = 0.;
        z = 1.;
    }

    inline Vector3d( float x, float y, float z )
    {
        this->x = x;
        this->y = y;
        this->z = z;
    }

    inline void set( float x, float y, float z )
    {
        this->x = x;
        this->y = y;
        this->z = z;
    }

    inline void normalize()
    {
        float d = sqrt(x*x+y*y+z*z);
        x /= d;
        y /= d;
        z /= d;
    }

    inline bool operator==( const Vector3d& v ) const
    {
        return ( x==v.x && y==v.y && z==v.z );
    }

    inline bool operator!=( const Vector3d& v ) const
    {
        return !( (*this) == v );
    }

    inline Vector3d& operator=( const Vector3d& v )
    {
        x = v.x;
        y = v.y;
        z = v.z;

        return *this;
    }

    inline Vector3d& operator+=( const Vector3d& v )
    {
        x += v.x;
        y += v.y;
        z += v.z;

        return *this;
    }

    inline Vector3d& operator-=( const Vector3d& v )
    {
        x -= v.x;
        y -= v.y;
        z -= v.z;

        return *this;
    }

    inline Vector3d& operator*=( const Vector3d& v )
    {
        Vector3d l = *this;
        x = l.y*v.z-l.z*v.y;
        y = l.z*v.x-l.x*v.z;
        z = l.x*v.y-l.y*v.x;
        return *this;
    }

    inline Vector3d& operator*=( float s )
    {
        x *= s;
        y *= s;
        z *= s;
        return *this;
    }

    inline Vector3d& operator/=( float s )
    {
        x /= s;
        y /= s;
        z /= s;
        return *this;
    }

    inline Vector3d operator*( const Vector3d v ) const
    {
        Vector3d r;
        r.x = y*v.z-z*v.y;
        r.y = z*v.x-x*v.z;
        r.z = x*v.y-y*v.x;
        return r;
    }

    inline Vector3d operator+( const Vector3d v ) const
    {
        Vector3d r = *this;
        r += v;
        return r;
    }

    inline Vector3d operator-( const Vector3d v ) const
    {
        Vector3d r = *this;
        r -= v;
        return r;
    }

    inline Vector3d operator*( float s ) const
    {
        Vector3d r(x*s, y*s, z*s);
        return r;
    }

    inline Vector3d operator/( float s ) const
    {
        Vector3d r(x/s, y/s, z/s);
        return r;
    }

    inline Vector3d operator+( float s ) const
    {
        Vector3d r(x+s, y+s, z+s);
        return r;
    }

    inline Vector3d operator-( float s ) const
    {
        Vector3d r(x-s, y-s, z-s);
        return r;
    }

    inline float dotProduct( const Vector3d& r ) const
    {
        return x*r.x + y*r.y + z*r.z;
    }

    inline Vector3d crossProduct( const Vector3d& v ) const
    {
        Vector3d r;
        r.x = y*v.z-z*v.y;
        r.y = z*v.x-x*v.z;
        r.z = x*v.y-y*v.x;
        return r;
    }

    inline void rotate( const Vector3d& r, float a )
        // Rotate vector about another vector
    {
        float u = r.x;
        float v = r.y;
        float w = r.z;

        float sa = sin(a);
        float ca = cos(a);

        float ux = u*x;
        float vv = v*v;
        float ww = w*w;
        float vy = v*y;
        float wz = w*z;
        float t = (ux+vy+wz);

        float xn = u*t + ( x*(vv+ww) - u*(vy+wz) )*ca + (-w*y+v*z)*sa;
        float uu = u*u;
        float yn = v*(ux+vy+wz) + (y*(u*u+w*w) - v*(ux+wz))*ca + (w*x-u*z)*sa;
        float zn = w*t + ( z*(uu+vv) - w*(ux+vy) )*ca + (-v*x+u*y)*sa;

        x = xn;
        y = yn;
        z = zn;

    }

    inline void display() const
    {
        cout << x << ", " << y << ", " << z << endl;
    }

};


#endif /* _VECTMATH_HPP_ */
