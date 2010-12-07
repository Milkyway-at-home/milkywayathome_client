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

#ifndef _DRAW8_HPP_
#define _DRAW8_HPP_

#include <iostream>
#include <cstdlib>
#include "drawcore.hpp"

using namespace std;


inline void clearSurface8( SDL_Surface *surface, Uint8 index = 0x00 );

void dimSurface8( SDL_Surface *surface, int dimShift );

inline Uint8 getPixel8( const SDL_Surface *surface, int x, int y );

inline Uint8 getPixelClip8( const SDL_Surface *surface, int x, int y );

inline void putPixel8( SDL_Surface *surface, int x, int y, Uint8 index );

inline void putPixelClip8( SDL_Surface *surface, int x, int y, Uint8 index );

inline void putPixelSum8( SDL_Surface *surface, int x, int y, Uint8 index );

inline void putPixelSumClip8( SDL_Surface *surface, int x, int y, Uint8 index );

inline void blitSurfaceClipSum8( const SDL_Surface *copySurface, SDL_Surface *destSurface, int x, int y );

inline SDL_Surface *newSurface8( int xSize, int ySize );

void brightAdjustSurface8( SDL_Surface *surface, float iMult );

inline void putPixelTrans8( SDL_Surface *surface, int x, int y, Uint8 index );

inline void putPixelTransClip8( SDL_Surface *surface, int x, int y, Uint8 index );

inline void blitSurfaceClip8( const SDL_Surface *copySurface, SDL_Surface *destSurface, int x, int y );

inline void blitSurfaceClipTrans8( const SDL_Surface *copySurface, SDL_Surface *destSurface, int x, int y );

inline void _benchTest_blitSurfaceClipSum8( const SDL_Surface *copySurface, SDL_Surface *destSurface, int x, int y );

void flipHorizontal8( SDL_Surface *surface );

void flipVertical8( SDL_Surface *surface );

inline void copySurfaceClip8( const SDL_Surface *copySurface, SDL_Surface *destSurface, int xc, int yc, int xd, int yd, int xSize, int ySize );

void resizeCopySurfaceClip8( const SDL_Surface *copySurface, SDL_Surface *destSurface, float xPercent, float yPercent );

void resizeBlendCopySurfaceClip8( const SDL_Surface *copySurface, SDL_Surface *destSurface, int xd, int yd, int xSkip, int ySkip );

void resizeCopySurfaceClip8( const SDL_Surface *copySurface, SDL_Surface *destSurface, int xd, int yd, int xSkip, int ySkip );

void resizeSumSurfaceClip8( const SDL_Surface *copySurface, SDL_Surface *destSurface, int xd, int yd, int xSkip, int ySkip );


inline void clearSurface8( SDL_Surface *surface, Uint8 index )
{
    Uint32 value = (index<<24)|(index<<16)|(index<<8)|index;
    SDL_FillRect(surface, NULL, value);
}

/// TODO /// Retest

/// TODO /// speed test - Try Uint32 for 'index'
inline Uint8 getPixel8( const SDL_Surface *surface, int x, int y )
{
#ifndef NDEBUG
    if( (unsigned int) x >= (unsigned int) surface->w || (unsigned int) y >= (unsigned int) surface->h ) {
        cerr << "Putpixel access out of range (0-" << surface->w-1 << ", 0-" << surface->h-1 << ")\n";
        exit(1);
    }
#endif
    return ((Uint8*)(surface->pixels))[ y*(surface->pitch)+x ];
}

inline Uint8 getPixelClip8( const SDL_Surface *surface, int x, int y )
{
    if( (unsigned int) x >= (unsigned int) surface->w || (unsigned int) y >= (unsigned int) surface->h )
        return 0;
    return getPixel8(surface, x, y);
}

inline void putPixel8( SDL_Surface *surface, int x, int y, Uint8 index )
{
#ifndef NDEBUG
    if( (unsigned int) x >= (unsigned int) surface->w || (unsigned int) y >= (unsigned int) surface->h ) {
        cerr << "Putpixel access out of range (0-" << surface->w-1 << ", 0-" << surface->h-1 << ")\n";
        exit(1);
    }
#endif
    ((Uint8*)(surface->pixels))[ y*(surface->pitch)+x ] = index;
}

inline void putPixelClip8( SDL_Surface *surface, int x, int y, Uint8 index )
{
    if( (unsigned int) x >= (unsigned int) surface->w || (unsigned int) y >= (unsigned int) surface->h )
        return;
    putPixel8(surface, x, y, index);
}

inline void putPixelSum8( SDL_Surface *surface, int x, int y, Uint8 index )
{
    /// TODO /// Retest
#ifndef NDEBUG
    if( (unsigned int) x >= (unsigned int) surface->w || (unsigned int) y >= (unsigned int) surface->h ) {
        cerr << "Putpixel access out of range (0-" << surface->w-1 << ", 0-" << surface->h-1 << ")\n";
        exit(1);
    }
#endif
    Uint8 *pixelPtr = ((Uint8*)(surface->pixels)) + y*(surface->pitch)+x;

    /// TODO /// speed test - try removing register
    Uint16 r1 = (*pixelPtr) + index;
    if( r1&0x0100 )
        *pixelPtr = 0xff;
    else
        *pixelPtr = r1;
}

inline void putPixelSumClip8( SDL_Surface *surface, int x, int y, Uint8 index )
{
    if( (unsigned int) x >= (unsigned int) surface->w || (unsigned int) y >= (unsigned int) surface->h )
        return;
    putPixelSum8(surface, x, y, index);
}

inline void blitSurfaceClipSum8( const SDL_Surface *copySurface, SDL_Surface *destSurface, int x, int y )
{

    /// TODO /// Consider assembly alternative
    /// TODO /// Benchmark new version

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
        dp += x;
    }

    // Clip upper half

    Sint32 temp = y+ySize-destSurface->h;
    if( temp>0 )
        yc += temp;

    temp = x+xSize-destSurface->w;
    if( temp>0 )
        xs += temp;

    Sint32 ydi = destSurface->pitch-xSize+xs;

    // Sum memory

    for( ; yc<ySize; yc++, cp += xs, dp += ydi ) {

        // 4-pixel copy

        xSize -= 3;
        for( xc = xs; xc<xSize; xc += 4, cp += 4, dp += 4 ) {

            Uint32 r1, r2, r3;
            r1 = *((Uint32*) cp);
            r1 &= 0xfefefefe;
            r1 >>= 1;
            r2 = *((Uint32*) dp);
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

            *((Uint32*) dp) = r1;

        }

        // 1-pixel copy

        xSize += 3;

        while( xc<xSize ) {

            Uint16 r1 = *cp + *dp;
            if( r1&0x0100 )
                *dp = 0xfe;
            else
                *dp = r1&0xfe;

            xc++;
            cp++;
            dp++;

        }

    }

}

inline SDL_Surface *newSurface8( int xSize, int ySize )
{
    /// TODO /// Try hw surface - speed test
    SDL_Surface *surface = SDL_CreateRGBSurface(SDL_SWSURFACE, xSize, ySize, 8, 0, 0, 0, 0);
    if( surface == NULL ) {
        cerr << "SDL_CreateRGBSurface: " << SDL_GetError() << endl;
        exit(1);
    }

    return surface;
}

/*
inline void putPixelTrans8( SDL_Surface *surface, int x, int y, Uint8 index )
{
#ifndef NDEBUG
    if( (unsigned int) x >= (unsigned int) surface->w || (unsigned int) y >= (unsigned int) surface->h ) {
        cerr << "Putpixel access out of range (0-" << surface->w-1 << ", 0-" << surface->h-1 << ")\n";
        exit(1);
    }
#endif
    Uint8 *pixelPtr = ((Uint8*)(surface->pixels)) + y*(surface->pitch)+x;
    *pixelPtr = ((*pixelPtr)*(127-index)>>7)+index;
}

inline void putPixelTransClip8( SDL_Surface *surface, int x, int y, Uint8 index )
{
    if( (unsigned int) x >= (unsigned int) surface->w || (unsigned int) y >= (unsigned int) surface->h )
        return;
    putPixelTrans8(surface, x, y, index);
}

inline void blitSurfaceClip8( const SDL_Surface *copySurface, SDL_Surface *destSurface, int x, int y )
{

    /// TODO /// update this routine to accept multiple pitches
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
        cp += yc*xSize;
    }
    else {
        yc = 0;
        dp += y*destSurface->w;
    }

    if( x<0 ) {
        xs = -x;
        cp += xs;
    }
    else {
        xs = 0;
        dp += x;
    }

    // Clip upper half

    Sint32 temp = y+ySize-destSurface->h;
    if( temp>0 )
    yc += temp;

    temp = x+xSize-destSurface->w;
    if( temp>0 )
    xs += temp;

    Sint32 ydi = destSurface->w-xSize+xs;

    // Copy memory

    for( ; yc<ySize; yc++, cp += xs, dp += ydi ) {

        // 4-pixel copy

        xSize -= 3;
        for( xc = xs; xc<xSize; xc += 4, cp += 4, dp += 4 )
            *((Uint32*) dp) = *((Uint32*) cp);

        // 1-pixel copy

        xSize += 3;
        while( xc<xSize ) {
            *dp = *cp;
            xc++;
            cp++;
            dp++;
        }

    }

}

inline void blitSurfaceClipTrans8( const SDL_Surface *copySurface, SDL_Surface *destSurface, int x, int y )
{

    /// TODO /// update this routine to accept multiple pitches
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
        cp += yc*xSize;
    }
    else {
        yc = 0;
        dp += y*destSurface->w;
    }

    if( x<0 ) {
        xs = -x;
        cp += xs;
    }
    else {
        xs = 0;
        dp += x;
    }

    // Clip upper half

    Sint32 temp = y+ySize-destSurface->h;
    if( temp>0 )
        yc += temp;

    temp = x+xSize-destSurface->w;
    if( temp>0 )
        xs += temp;

    Sint32 ydi = destSurface->w-xSize+xs;

    // Transpose memory

    for( ; yc<ySize; yc++, cp += xs, dp += ydi ) {

        xc = xs;

        while( xc<xSize ) {

            Uint16 r1 = *cp;
            Uint16 r2 = r1;
            r2 ^= 0x7f;

            Uint16 r3 = *dp;

            r3 *= r2;
            r3 >>= 7;

            r1 += r3;
            r1 &= 0x7f;

            *dp = r1;

            xc++;
            cp++;
            dp++;

        }

    }

}

inline void _benchTest_blitSurfaceClipSum8( const SDL_Surface *copySurface, SDL_Surface *destSurface, int x, int y )
{
    // This function is only for testing purposes - applications should use 'blitSurfaceClipSum7' instead
    for( int yy = 0; yy<copySurface->h; yy++)
        for( int xx = 0; xx<copySurface->w; xx++)
            putPixelSumClip8(destSurface, x+xx, y+yy, getPixel8(copySurface, xx, yy));
}

inline void copySurfaceClip8( const SDL_Surface *copySurface, SDL_Surface *destSurface, int xc, int yc, int xd, int yd, int xSize, int ySize )
{
    /// TODO /// Optimize this
    for( int yy = 0; yy<ySize; yy++ )
        for( int xx = 0; xx<xSize; xx++ )
            putPixelClip8(destSurface, xd+xx, yd+yy, getPixelClip8(copySurface, xc+xx, yc+yy));
}

*/

#endif /* _DRAW8_HPP_ */
