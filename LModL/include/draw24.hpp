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

#ifndef _DRAW24_HPP_
#define _DRAW24_HPP_

#include <iostream>
#include <cstdlib>
#include "drawcore.hpp"

using namespace std;


inline void putPixel24( SDL_Surface *surface, int x, int y, Uint32 color )
{

#ifndef NDEBUG
    if( (unsigned int) x >= (unsigned int) surface->w || (unsigned int) y >= (unsigned int) surface->h ) {
        cerr << "Putpixel access out of range (0-" << surface->w-1 << ", 0-" << surface->h-1 << ")\n";
        exit(1);
    }
#endif

    Uint8 *p = (Uint8*) surface->pixels + y * surface->pitch + x*3;

#if SDL_BYTEORDER == SDL_BIG_ENDIAN
        p[0] = (color >> 16) & 0xff;
        p[1] = (color >> 8) & 0xff;
        p[2] = color & 0xff;
#else
        p[0] = color & 0xff;
        p[1] = (color >> 8) & 0xff;
        p[2] = (color >> 16) & 0xff;
#endif

}

inline Uint32 getPixel24( const SDL_Surface *surface, int x, int y )
{

#ifndef NDEBUG
    if( (unsigned int) x >= (unsigned int) surface->w || (unsigned int) y >= (unsigned int) surface->h ) {
        cerr << "Putpixel access out of range (0-" << surface->w-1 << ", 0-" << surface->h-1 << ")\n";
        exit(1);
    }
#endif

    Uint8 *p = (Uint8*) surface->pixels + y * surface->pitch + x*3;
    Uint32 color;

#if SDL_BYTEORDER == SDL_BIG_ENDIAN
        color = p[0]<<16;
        color |= p[1]<<8;
        color |= p[2];
#else
        color = p[2]<<16;
        color |= p[1]<<8;
        color |= p[0];
#endif

	return color;

}


#endif /* _DRAW24_HPP_ */
