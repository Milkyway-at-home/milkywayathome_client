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

#include "draw8.hpp"

using namespace std;

/// TODO /// Explore auto-locking options

/// TODO /// speed test - Try using an asychronous/synchronous blit of a blank hw_surface for dirty rectangle clear

/// TODO /// Retest

void dimSurface8( SDL_Surface *surface, int dimShift )
{

    lockSurface(surface);

    if( dimShift==0 )
        return;

    int xSize = surface->w;
    int ySize = surface->h;
    Sint32 ydi = surface->pitch - xSize;
    Uint8 *dp = (Uint8*) surface->pixels;

    Uint32 dimMask = (Uint32[]) {
        0x00000000,
        0x7f7f7f7f,
        0x3f3f3f3f,
        0x1f1f1f1f,
        0x0f0f0f0f,
        0x07070707,
        0x03030303,
        0x01010101
    } [dimShift];

    for( int yc = 0; yc<ySize; yc++, dp += ydi ) {

        // 4-pixel change

        xSize -= 3;
        int xc = 0;
        for( ; xc<xSize; xc += 4, dp += 4 ) {
            Uint32 r1 = *((Uint32*) dp);
            r1 >>= dimShift;
            r1 &= dimMask;
            *((Uint32*) dp) = r1;
        }

        // 1-pixel change

        xSize += 3;
        while( xc<xSize ) {
            *dp >>= dimShift;
            xc++;
            dp++;
        }

    }

    unlockSurface(surface);

}

/// TODO /// Change to 8bpp mimic from 7bpp and apply other changes

/*
void brightAdjustSurface8( SDL_Surface *surface, float iMult )
{

    lockSurface(surface);

    Uint16 im = (int) (iMult*256.);

    if( im==256 )
        return;

    int xSize = surface->w;
    int ySize = surface->h;
    Sint32 ydi = surface->pitch - xSize;
    Uint8 *dp = (Uint8*) surface->pixels;

    for( int yc = 0; yc<ySize; yc++, dp += ydi )
        for( int xc = 0; xc<xSize; xc++, dp++ ) {
            Uint16 r1 = Uint16(*dp);
            Uint32 r2 = im*r1;
            r2 >>= 8;
            if( r2>0x7f )
                r2 = 0x7f;
            *dp = (Uint8) r2;
        }

    unlockSurface(surface);

}

void flipHorizontal8( SDL_Surface *surface )
{
    lockSurface(surface);
    /// TODO /// Optimize this and add partial transform functionality, bench against SDL-Gfx library some time
    for( int y = 0; y<surface->h; y++ )
        for( int x = 0; x<(surface->w>>1); x++ ) {
            Uint8 swap = getPixel8(surface, x, y);
            putPixel8(surface, x, y, getPixel8(surface, surface->w-x-1, y));
            putPixel8(surface, surface->w-x-1, y, swap);
        }
    unlockSurface(surface);
}

void flipVertical8( SDL_Surface *surface )
{
    lockSurface(surface);
    /// TODO /// Optimize this and add partial transform functionality, bench against SDL-Gfx library some time
    for( int y = 0; y<(surface->h>>1); y++ )
        for( int x = 0; x<surface->w; x++ ) {
            Uint8 swap = getPixel8(surface, x, y);
            putPixel8(surface, x, y, getPixel8(surface, x, surface->h-y-1));
            putPixel8(surface, x, surface->h-y-1, swap);
        }
    unlockSurface(surface);
}

void resizeCopySurfaceClip8( const SDL_Surface *copySurface, SDL_Surface *destSurface, float xPercent, float yPercent )
{
    /// TODO /// STUB
}

void resizeBlendCopySurfaceClip8( const SDL_Surface *copySurface, SDL_Surface *destSurface, int xd, int yd, int xSkip, int ySkip )
{
    /// TODO /// STUB
}

void resizeCopySurfaceClip8( const SDL_Surface *copySurface, SDL_Surface *destSurface, int xd, int yd, int xSkip, int ySkip )
{
    /// TODO /// Optimize this
    for( int ys = 0, y = 0; ys<copySurface->h; ys += ySkip, y++ )
        for( int xs = 0, x = 0; xs<copySurface->w; xs += xSkip, x++ )
            putPixelClip8(destSurface, xd+x, yd+y, getPixelClip8(copySurface, xs, ys));
}

void resizeSumSurfaceClip8( const SDL_Surface *copySurface, SDL_Surface *destSurface, int xd, int yd, int xSkip, int ySkip )
{
    /// TODO /// Optimize this
    for( int ys = 0, y = 0; ys<copySurface->h; ys += ySkip, y++ )
        for( int xs = 0, x = 0; xs<copySurface->w; xs += xSkip, x++ )
            putPixelSumClip8(destSurface, xd+x, yd+y, getPixelClip8(copySurface, xs, ys));
}
*/
