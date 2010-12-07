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

#ifndef _DRAW32_HPP_
#define _DRAW32_HPP_

#include <iostream>
#include <cstdlib>
#include "drawcore.hpp"

using namespace std;


void clearSurface32( SDL_Surface *surface, Uint32 color = 0x00000000 );

inline void putPixelSum32( SDL_Surface *surface, int x, int y, Uint32 color );

inline void putPixelSumClip32( SDL_Surface *surface, int x, int y, Uint32 color );

inline Uint32 getPixel32( const SDL_Surface *surface, int x, int y );

inline Uint32 getPixelClip32( const SDL_Surface *surface, int x, int y );

inline void putPixel32( SDL_Surface *surface, int x, int y, Uint32 color );

inline void putPixelClip32( SDL_Surface *surface, int x, int y, Uint32 color );

inline void blitSurfaceClipSum32( const SDL_Surface *copySurface, SDL_Surface *destSurface, int x, int y );

inline SDL_Surface * newSurface32( int xSize, int ySize, bool trans = false );

void dimSurface32( SDL_Surface *surface, int dimShift );

void brightAdjustSurface32( SDL_Surface *surface, float iMult );

inline void putPixelTrans32( SDL_Surface *surface, int x, int y, Uint8 index );

inline void putPixelTransClip32( SDL_Surface *surface, int x, int y, Uint8 index );

inline void blitSurfaceClip32( const SDL_Surface *copySurface, SDL_Surface *destSurface, int x, int y );

inline void blitSurfaceClipTrans32( const SDL_Surface *copySurface, SDL_Surface *destSurface, int x, int y );

void flipHorizontal32( SDL_Surface *surface );

void flipVertical32( SDL_Surface *surface );

inline void copySurfaceClip32( const SDL_Surface *copySurface, SDL_Surface *destSurface, int xc, int yc, int xd, int yd, int xSize, int ySize );

void resizeCopySurfaceClip32( const SDL_Surface *copySurface, SDL_Surface *destSurface, float xPercent, float yPercent );

void resizeBlendCopySurfaceClip32( const SDL_Surface *copySurface, SDL_Surface *destSurface, int xd, int yd, int xSkip, int ySkip );

void resizeCopySurfaceClip32( const SDL_Surface *copySurface, SDL_Surface *destSurface, int xd, int yd, int xSkip, int ySkip );

void resizeSumSurfaceClip32( const SDL_Surface *copySurface, SDL_Surface *destSurface, int xd, int yd, int xSkip, int ySkip );


inline void putPixelSum32( SDL_Surface *surface, int x, int y, Uint32 color )
{
#ifndef NDEBUG
    if( (unsigned int) x >= (unsigned int) surface->w || (unsigned int) y >= (unsigned int) surface->h ) {
        cerr << "Putpixel access out of range (0-" << surface->w-1 << ", 0-" << surface->h-1 << ")\n";
        exit(1);
    }
#endif

    Uint8 *pixelPtr = ((Uint8*)(surface->pixels)) + y*(surface->pitch)+(x<<2);

    Uint32 r1, r2, r3;
    r1 = color&0xfefefefe;
    r1 >>= 1;
    r2 = *( (Uint32*) pixelPtr );
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

    *( (Uint32*) pixelPtr ) = r1;

}

inline void putPixelSumClip32( SDL_Surface *surface, int x, int y, Uint32 color )
{
    if( (unsigned int) x >= (unsigned int) surface->w || (unsigned int) y >= (unsigned int) surface->h )
        return;
    putPixelSum32(surface, x, y, color);
}

inline Uint32 getPixel32( const SDL_Surface *surface, int x, int y )
{
#ifndef NDEBUG
    if( (unsigned int) x >= (unsigned int) surface->w || (unsigned int) y >= (unsigned int) surface->h ) {
        cerr << "Putpixel access out of range (0-" << surface->w-1 << ", 0-" << surface->h-1 << ")\n";
        exit(1);
    }
#endif
    return *( (Uint32*) ( (Uint8*) surface->pixels + y*(surface->pitch)+(x<<2) ) );
}

inline Uint32 getPixelClip32( const SDL_Surface *surface, int x, int y )
{
    if( (unsigned int) x >= (unsigned int) surface->w || (unsigned int) y >= (unsigned int) surface->h )
        return 0;
    return getPixel32(surface, x, y);
}

inline void putPixel32( SDL_Surface *surface, int x, int y, Uint32 color )
{
#ifndef NDEBUG
    if( (unsigned int) x >= (unsigned int) surface->w || (unsigned int) y >= (unsigned int) surface->h ) {
        cerr << "Putpixel access out of range (0-" << surface->w-1 << ", 0-" << surface->h-1 << ")\n";
        exit(1);
    }
#endif
    *( (Uint32*) ( (Uint8*)  surface->pixels + y*(surface->pitch)+(x<<2) ) ) = color;
}

inline void putPixelClip32( SDL_Surface *surface, int x, int y, Uint32 color )
{
    if( (unsigned int) x >= (unsigned int) surface->w || (unsigned int) y >= (unsigned int) surface->h )
        return;
    putPixel32(surface, x, y, color);
}

inline void blitSurfaceClip32( const SDL_Surface *copySurface, SDL_Surface *destSurface, int x, int y )
{

    SDL_Rect destRect;

    destRect.x = x;
    destRect.y = y;

    SDL_BlitSurface((SDL_Surface*)copySurface, NULL, destSurface, &destRect);

}

inline void blitSurfaceClipSum32( const SDL_Surface *copySurface, SDL_Surface *destSurface, int x, int y )
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
            r1 = *( (Uint32*) cp );
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
        }

}

inline SDL_Surface * newSurface32( int xSize, int ySize, bool trans )
{

    SDL_Surface *surface;

    Uint32 rmask, gmask, bmask, amask;

#if SDL_BYTEORDER == SDL_BIG_ENDIAN
    rmask = 0x000000ff;
    gmask = 0x0000ff00;
    bmask = 0x00ff0000;
    amask = 0xff000000;
#else
    rmask = 0xff000000;
    gmask = 0x00ff0000;
    bmask = 0x0000ff00;
    amask = 0x000000ff;
#endif

    if( !trans )
        amask = 0x00000000;

    /* Initialize image buffer */
    surface = SDL_CreateRGBSurface(SDL_SWSURFACE, xSize, ySize, 32, rmask, gmask, bmask, amask);
    if( surface == NULL ) {
        fprintf(stderr, "SDL_CreateRGBSurface: %s\n", SDL_GetError());
        exit(1);
    }

    return surface;

}

/// TODO /// complete/test routines beyond this point

/*

inline void putPixelTrans32( SDL_Surface *surface, int x, int y, Uint8 index )
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

inline void putPixelTransClip32( SDL_Surface *surface, int x, int y, Uint8 index )
{
    if( (unsigned int) x >= (unsigned int) surface->w || (unsigned int) y >= (unsigned int) surface->h )
        return;
    putPixelTrans32(surface, x, y, index);
}

inline void blitSurfaceClip32( const SDL_Surface *copySurface, SDL_Surface *destSurface, int x, int y )
{

    /// TODO /// update this segment
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

inline void blitSurfaceClipTrans32( const SDL_Surface *copySurface, SDL_Surface *destSurface, int x, int y )
{

    /// TODO /// update this segment
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

inline void copySurfaceClip32( const SDL_Surface *copySurface, SDL_Surface *destSurface, int xc, int yc, int xd, int yd, int xSize, int ySize )
{
    /// TODO /// Optimize this
    for( int yy = 0; yy<ySize; yy++ )
        for( int xx = 0; xx<xSize; xx++ )
            putPixelClip32(destSurface, xd+xx, yd+yy, getPixelClip32(copySurface, xc+xx, yc+yy));
}


*/

#endif /* _DRAW32_HPP_ */
