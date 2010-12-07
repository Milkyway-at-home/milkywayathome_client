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

#ifndef _RGBCONV_HPP_
#define _RGBCONV_HPP_

#include <cmath>
#include "SDL.h"
#include "trigtbl.hpp"


using namespace std;


inline Uint32 rgbToColor( const SDL_Surface* surface, Uint8 r, Uint8 g, Uint8 b );

inline void colorToRgb( const SDL_Surface* surface, Uint32 color, Uint8& r, Uint8& g, Uint8& b );

inline Uint32 lToColor32( Uint8 l );

inline Uint8 colorToL( const SDL_Surface* surface, Uint32 color );

inline Uint8 rgbToL( Uint8 r, Uint8 g, Uint8 b );

inline Uint32 rgbToColor( const SDL_Surface* surface, Uint8 r, Uint8 g, Uint8 b )
{
    return SDL_MapRGB(surface->format, r, g, b);
}

inline void colorToRgb( const SDL_Surface* surface, Uint32 color, Uint8& r, Uint8& g, Uint8& b )
{
	SDL_GetRGB(color, surface->format, &r, &g, &b);
}

inline Uint32 lToColor32( const SDL_Surface* surface, Uint8 l )
{
   return rgbToColor(surface, l, l, l);
}

inline Uint8 colorToL( const SDL_Surface* surface, Uint32 color )
{
   Uint8 r, g, b;
   colorToRgb(surface, color, r, g, b);
   return rgbToL(r, g, b);
}

inline Uint8 rgbToL( Uint8 r, Uint8 g, Uint8 b )
{
   return (r+g+b)/3;
}

void rgbToXyl( Uint8 r, Uint8 g, Uint8 b, double &x, double &y, double &l );

void xylToRgb( double x, double y, double l, Uint8 &r, Uint8 &g, Uint8 &b );

void xylToHsl( double x, double y, double li, double &h, double &s, double &l );

void hslToXyl( double h, double s, double li, double &x, double &y, double &l );

void rgbToHsl( Uint8 r, Uint8 g, Uint8 b, double &h, double &s, double &l );
    // Returns hue, saturation and lighting
    //   h is a degree in radians
    //   s values are from 0. - 1.
    //   l values are from 0. - 1.
    //   r, g, b, are integers from 0-255

void hslToRgb( double h, double s, double l, Uint8 &r, Uint8 &g, Uint8 &b );
    // Returns machine R, G, B values
    //   h is a degree in radians
    //   s values are from 0. - 1.
    //   l values are from 0. - 1.
    //   r, g, b, are integers from 0-255

void hslToRgb( double h, double s, double l, double &r, double &g, double &b );
    // Returns machine R, G, B values
    //   h is a degree in radians
    //   s values are from 0. - 1.
    //   l values are from 0. - 1.
    //   r, g, b, values are from 0.-1.

#endif /* _RGBCONV_HPP_ */
