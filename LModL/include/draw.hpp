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

#ifndef _DRAW_HPP_
#define _DRAW_HPP_

#include "drawcore.hpp"
#include "draw8.hpp"
#include "draw24.hpp"
#include "draw32p.hpp"

using namespace std;


void blitSurfaceClipSum8to32( const SDL_Surface* copySurface, SDL_Surface*destSurface, int x, int y, Uint32* palette );

void inline blitSurfaceClipSumPalette( SDL_Surface* copySurface, SDL_Surface* destSurface, int x, int y, Uint32* palette = NULL );

void inline blitSurfaceClipSumPalette( SDL_Surface* copySurface, SDL_Surface* destSurface, int x, int y, Uint32* palette )
{

#ifndef NDEBUG
    if( copySurface->format->BytesPerPixel!=1 ) {
        cerr << "Must apply indexed color from 8bpp surface\n";
        exit(1);
    }
#endif

    switch( destSurface->format->BytesPerPixel ) {

    case 1:
        blitSurfaceClipSum8(copySurface, destSurface, x, y);
        break;

    /// TODO /// Add these
/*
    case 2:
        break;

    case 3:
        break;
*/
    case 4:
        blitSurfaceClipSum8to32(copySurface, destSurface, x, y, palette);
        break;

    default:
        cerr << "Invalid resolution requested: " << (int) destSurface->format->BytesPerPixel << " bpp\n";
        exit(1);

    }

}

inline Uint32 getPixel( const SDL_Surface *surface, int x, int y )
{
    if( surface->format->BitsPerPixel==24 )
        return getPixel24(surface, x, y);
    else if( surface->format->BitsPerPixel==32 )
        return getPixel32(surface, x, y);
    else {
        cerr << surface->format->BitsPerPixel << " bits-per-pixel getPixel not supported" << endl;
        exit(1);
    }
}

inline void putPixel( SDL_Surface *surface, int x, int y, Uint32 color )
{
    if( surface->format->BitsPerPixel==24 )
        putPixel24(surface, x, y, color);
    else if( surface->format->BitsPerPixel==32 )
        putPixel32(surface, x, y, color);
    else {
        cerr << surface->format->BitsPerPixel << " bits-per-pixel putPixel not supported" << endl;
        exit(1);
    }
}

/*
inline void chooseBpp( SDL_Surface* surface, f8, f16, f24, f32 )
{

    choice

    if( choice==NULL ) {
        cerr << "Function not supported for " << bpp << " bits-per-pixel\n";
    }


inline void clearSurface( SDL_Surface *surface, Uint32 color = 0x00000000 )
{
    * () ( switchBpp(clearSurface8, NULL, NULL, NULL, clearSurface32) );
}
*/


#endif /* _DRAW_HPP_ */
