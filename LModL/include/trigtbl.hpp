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

#ifndef _TRIGTBL_HPP_
#define _TRIGTBL_HPP_

#include <iostream>
#include <cstdlib>

using namespace std;

const float TRIG_PI           = 3.1415926535897932e+000;
const float TRIG_2PI          = 6.2831853071795864e+000;
const float TRIG_DEG_TO_RAD   = 1.7453292519943296e-002;
const float TRIG_RAD_TO_DEG   = 5.7295779513082321e+001;
const float TRIG_SQRT2_OVER_2 = 7.0710678118654752e-001;

// Lookup table accurate to 4 places

#define TRIGTBL_SIZE        126
#define TRIGTBL_COS_OFFSET  25
#define TRIGTBL_GRAN        1.5915494309189533e+01

/*
#define TRIGTBL_SIZE        12501
#define TRIGTBL_COS_OFFSET  2500
#define TRIGTBL_GRAN        1.5915494309189535e+003
*/

extern const float TRIG_TABLE[TRIGTBL_SIZE];


inline float sin_lu( float x );
    // 'sin' with lookup tables - many times faster than cos
    // should only be called with x between 0. and 2*pi

inline float cos_lu( float x );
    // 'cos' with lookup tables - many times faster than sin
    // should only be called with x between 0. and 2*pi


inline float sin_lu( float x )
{

    int ix = int(x*TRIGTBL_GRAN);

#ifndef NDEBUG
    if( (unsigned int) ix >= TRIGTBL_SIZE )
    {
        cerr << "Sine table access out of range (0-2*pi)\n";
        exit(1);
    }
#endif

   return TRIG_TABLE[ix];

}

inline float cos_lu( float x )
{

    int ix = int(x*TRIGTBL_GRAN) + TRIGTBL_COS_OFFSET;

#ifndef NDEBUG
    if( (unsigned int) ix >= TRIGTBL_SIZE )
    {
        cerr << "Cosine table access out of range (0-2*pi)\n";
        exit(1);
    }
#endif

    return TRIG_TABLE[ix];

}


#endif /* _TRIGTBL_HPP_ */
