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

#include "draw.hpp"


void blitSurfaceClipSum8to32( const SDL_Surface *copySurface, SDL_Surface*destSurface, int x, int y, Uint32* palette )
{

    /// TODO /// consider machine code

    // Clip

    Sint32 xc, yc, xSize;
    Uint8 *cp = (Uint8*) copySurface->pixels;
    Uint8 *dp = (Uint8*) destSurface->pixels;
    Sint32 xs, ySize;

    xSize = copySurface->w;
    ySize = copySurface->h;

    // Clip lower half

    if( y<0 ) {
        yc = -y;
        cp += yc*copySurface->pitch;
    }
    else {
        yc = 0;
        dp += y*destSurface->pitch;
    }

    if( x<0 ) {
        xs = -x;
        cp += xs;
    }
    else {
        xs = 0;
        dp += x<<2;
    }

    // Clip upper half

    Sint32 temp = y+ySize-destSurface->h;
    if( temp>0 )
        yc += temp;

    temp = x+xSize-destSurface->w;
    if( temp>0 )
        xs += temp;

    Sint32 ydi = destSurface->pitch-((xSize-xs)<<2);

    // Sum memory

    for( ; yc<ySize; yc++, cp += xs, dp += ydi )

        for( xc = xs; xc<xSize; xc++, cp++, dp += 4 ) {

            Uint32 r1, r2, r3;
            r1 = palette[*cp];
            r1 &= 0xfefefefe;
            r1 >>= 1;
            r2 = *( (Uint32*) dp );
            r2 &= 0xfefefefe;
            r2 >>= 1;

            r2 += r1;
            r1 = r2;
            r2 &= 0x80808080;
            r3 = r2;
            r3 >>= 7;
            r2 -= r3;
            r1 |= r2;
            r1 &= 0x7f7f7f7f;
            r1 <<= 1;

            *( (Uint32*) dp ) = r1;
//putPixelSumClip32(destSurface, xc, yc, 0xff);
        }

}
