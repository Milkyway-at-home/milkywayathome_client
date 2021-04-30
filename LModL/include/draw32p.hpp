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

#ifndef _DRAW32P_HPP_
#define _DRAW32P_HPP_

#include <iostream>
#include <cstdlib>
#include "drawcore.hpp"
#include "draw32.hpp"



using namespace std;


void clearSurface32p( SDL_Surface *surface, Uint32 color = 0x00000000 );

inline void putPixelSum32p( SDL_Surface *surface, int x, int y, Uint32 color );

inline void putPixelSumClip32p( SDL_Surface *surface, int x, int y, Uint32 color );

inline Uint32 getPixel32p( const SDL_Surface *surface, int x, int y );

inline Uint32 getPixelClip32p( const SDL_Surface *surface, int x, int y );

inline void putPixel32p( SDL_Surface *surface, int x, int y, Uint32 color );

inline void putPixelClip32p( SDL_Surface *surface, int x, int y, Uint32 color );

inline void blitSurfaceClipSum32p( const SDL_Surface *copySurface, SDL_Surface *destSurface, int x, int y );

inline SDL_Surface * newSurface32p( int xSize, int ySize );

void dimSurface32p( SDL_Surface *surface, int dimShift );

void brightAdjustSurface32p( SDL_Surface *surface, float iMult );

inline void putPixelTrans32p( SDL_Surface *surface, int x, int y, Uint8 index );

inline void putPixelTransClip32p( SDL_Surface *surface, int x, int y, Uint8 index );

inline void blitSurfaceClip32p( const SDL_Surface *copySurface, SDL_Surface *destSurface, int x, int y );

inline void blitSurfaceClipTrans32p( const SDL_Surface *copySurface, SDL_Surface *destSurface, int x, int y );

void flipHorizontal32p( SDL_Surface *surface );

void flipVertical32p( SDL_Surface *surface );

inline void copySurfaceClip32p( const SDL_Surface *copySurface, SDL_Surface *destSurface, int xc, int yc, int xd, int yd, int xSize, int ySize );

void resizeCopySurfaceClip32p( const SDL_Surface *copySurface, SDL_Surface *destSurface, float xPercent, float yPercent );

void resizeBlendCopySurfaceClip32p( const SDL_Surface *copySurface, SDL_Surface *destSurface, int xd, int yd, int xSkip, int ySkip );

void resizeCopySurfaceClip32p( const SDL_Surface *copySurface, SDL_Surface *destSurface, int xd, int yd, int xSkip, int ySkip );

void resizeSumSurfaceClip32p( const SDL_Surface *copySurface, SDL_Surface *destSurface, int xd, int yd, int xSkip, int ySkip );


inline void putPixelSum32pf( SDL_Surface *surface, int x, int y, Uint32 color )
{

/// STUB ///

#ifndef NDEBUG
    if( (unsigned int) x >= (unsigned int) surface->w || (unsigned int) y >= (unsigned int) surface->h ) {
        cerr << "Putpixel access out of range (0-" << surface->w-1 << ", 0-" << surface->h-1 << ")\n";
        exit(1);
    }
#endif

    Uint8 *pixelPtr = ((Uint8*)(surface->pixels)) + y*(surface->pitch)+(x<<2);

    // 32-bit representation
    // rrrrrrrr rrxggggg gggggxbb bbbbbbbb

    Uint32 a, b;

    a = color;                          // reg1 = color
    b = *( (Uint32*) pixelPtr );        // reg2 = *( (Uint32*) pixelPtr )

    a += b;                             // add reg1, reg2

    // Max intensity check testing
    if( a<b )                           // jc [skip 1]
        a |= 0xffc00000;                // or reg1, 0xffc000000
    if( a&0x0020000 )                   // mov reg2, reg1 / shl  reg2, 11 / jc [skip 1]
        a |= 0x001ff800;                // or reg1, 0x001ff000
    if( a&0x00000800 )                  // shl reg2, 10 / jc [skip 1]
        a |= 0x000003ff;                // or reg1, 0x000003ff

    a &= 0xffdffbff;                    // and reg1, 0xffdffbff

    *( (Uint32*) pixelPtr ) = a;

/*

    Uint32 r1, r2, r3, r4;

    // 32-bit representation
    // 00000000 00000000 00000000 00000000
    // rrrrrrrr rrxggggg ggggggxb bbbbbbbb

r = 0xffc00000
g = 0x001ffc00
b = 0x000001ff


    r1 = color;                    // r1:   rrrrrrrr rr0ggggg gggggg0b bbbbbbbb
    r2 = *( (Uint32*) pixelPtr );  // r2:   rrrrrrrr rr0ggggg gggggg0b bbbbbbbb

    r2 += r1;                      // r2: c rrrrrrrr rrgggggg ggggggbb bbbbbbbb
    r4 = r2>r1;                    // r4:   00000000 00000000 00000000 0000000c
    r4 <<= 21;                     // r4:   00000000 0c000000 00000000 00000000
    r1 = r2;                       // r1:   rrrrrrrr rrgggggg ggggggbb bbbbbbbb
    r2 &= 0x00200200;              // r2:   00000000 00g00000 000000b0 00000000
    r3 = r2;                       // r3:   00000000 00g00000 000000b0 00000000
    r3 >>= 10;                     // r3:   00000000 00000000 00000g00 0000000b
    r3 |= r4;                      // r3:   00000000 0r000000 00000g00 0000000b
    r2 -= r3;                      // r2:   mmmmmmmm m0mmmmmm mmmmm0mm mmmmmmm0
    r1 |= r2;
    r1 &= 0x7f7f7f7f;
    r1 <<= 1;

    *( (Uint32*) pixelPtr ) = r1;
*/

}

inline void putPixelSumClip32pf( SDL_Surface *surface, int x, int y, Uint32 color )
{
    if( (unsigned int) x >= (unsigned int) surface->w || (unsigned int) y >= (unsigned int) surface->h )
        return;
    putPixelSum32p(surface, x, y, color);
}

inline void putPixelSum32p( SDL_Surface *surface, int x, int y, Uint32 color )
{

#ifndef NDEBUG
    if( (unsigned int) x >= (unsigned int) surface->w || (unsigned int) y >= (unsigned int) surface->h ) {
        cerr << "Putpixel access out of range (0-" << surface->w-1 << ", 0-" << surface->h-1 << ")\n";
        exit(1);
    }
#endif

    Uint8 *pixelPtr = ((Uint8*)(surface->pixels)) + y*(surface->pitch)+(x<<2);

    // 32-bit-p representation
    // RRRRRRRR rrxGGGGG GGGggxBB BBBBBBbb

    // 32-bit representation
    // xxrrggbb RRRRRRRR GGGGGGGG BBBBBBBB

    Uint32 a, b;

    a = color;                          // reg1 = color
    b = *( (Uint32*) pixelPtr );        // reg2 = *( (Uint32*) pixelPtr )

    a = b32Tob32p(a);
    b = b32Tob32p(b);
    a += b;                             // add reg1, reg2

    // Max intensity check testing
    if( a<b )                           // jc [skip 1]
        a |= 0xffc00000;                // or reg1, 0xffc000000
    if( a&0x00200000 )                   // mov reg2, reg1 / shl  reg2, 11 / jc [skip 1]
        a |= 0x001ff800;                // or reg1, 0x001ff000
    if( a&0x00000400 )                  // shl reg2, 10 / jc [skip 1]
        a |= 0x000003ff;                // or reg1, 0x000003ff

    a &= 0xffdffbff;                    // and reg1, 0xffdffbff

    *( (Uint32*) pixelPtr ) = b32pTob32(a);

}

inline void putPixelSumClip32p( SDL_Surface *surface, int x, int y, Uint32 color )
{
    if( (unsigned int) x >= (unsigned int) surface->w || (unsigned int) y >= (unsigned int) surface->h )
        return;
    putPixelSum32p(surface, x, y, color);
}

inline Uint32 getPixel32p( const SDL_Surface *surface, int x, int y )
{
    return getPixel32(surface, x, y);
}

inline Uint32 getPixelClip32p( const SDL_Surface *surface, int x, int y )
{
    return getPixelClip32(surface, x, y);
}

inline void putPixel32p( SDL_Surface *surface, int x, int y, Uint32 color )
{
    putPixel32(surface, x, y, color);
}

inline void putPixelClip32p( SDL_Surface *surface, int x, int y, Uint32 color )
{
    putPixelClip32(surface, x, y, color);
}

inline void blitSurfaceClip32p( const SDL_Surface *copySurface, SDL_Surface *destSurface, int x, int y )
{
    blitSurfaceClip32(copySurface, destSurface, x, y);
}

inline void blitSurfaceClipSum32p( const SDL_Surface *copySurface, SDL_Surface *destSurface, int x, int y )
{

    /// TODO /// consider machine code if cross-platform is not an issue

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

            Uint32 a, b;

            a = *( (Uint32*) cp );              // reg1 = *( (Uint32*) cp )
            b = *( (Uint32*) dp );              // reg2 = *( (Uint32*) dp )

            a += b;                             // add reg1, reg2

            // Max intensity check testing
            if( a<b )                           // jc [skip 1]
                a |= 0xffc00000;                // or reg1, 0xffc000000
            if( a&0x00200000 )                  // mov reg2, reg1 / shl  reg2, 11 / jc [skip 1]
                a |= 0x001ff800;                // or reg1, 0x001ff000
            if( a&0x00000800 )                  // shl reg2, 10 / jc [skip 1]
                a |= 0x000003ff;                // or reg1, 0x000003ff

            a &= 0xffdffbff;                    // and reg1, 0xffdffbff

            *( (Uint32*) dp ) = a;

        }

}

inline SDL_Surface * newSurface32p( int xSize, int ySize )
{

    SDL_Surface *surface;

    Uint32 rmask, gmask, bmask, amask;

    rmask = 0xff000000;
    gmask = 0x000ff000;
    bmask = 0x000000ff;
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

inline void putPixelTrans32p( SDL_Surface *surface, int x, int y, Uint8 index )
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

inline void putPixelTransClip32p( SDL_Surface *surface, int x, int y, Uint8 index )
{
    if( (unsigned int) x >= (unsigned int) surface->w || (unsigned int) y >= (unsigned int) surface->h )
        return;
    putPixelTrans32p(surface, x, y, index);
}

inline void blitSurfaceClip32p( const SDL_Surface *copySurface, SDL_Surface *destSurface, int x, int y )
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

inline void blitSurfaceClipTrans32p( const SDL_Surface *copySurface, SDL_Surface *destSurface, int x, int y )
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

inline void copySurfaceClip32p( const SDL_Surface *copySurface, SDL_Surface *destSurface, int xc, int yc, int xd, int yd, int xSize, int ySize )
{
    /// TODO /// Optimize this
    for( int yy = 0; yy<ySize; yy++ )
        for( int xx = 0; xx<xSize; xx++ )
            putPixelClip32p(destSurface, xd+xx, yd+yy, getPixelClip32p(copySurface, xc+xx, yc+yy));
}


*/

#endif /* _DRAW32P_HPP_ */
