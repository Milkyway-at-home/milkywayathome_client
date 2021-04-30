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

#ifndef _DRAWCORE_HPP_
#define _DRAWCORE_HPP_

#include <string>
#include <cstdlib>
#include "rgbconv.hpp"
#include "SDL.h"

using namespace std;

#define DISPLAY_FLAGS  (SDL_SWSURFACE)

/*
#define HUE_GRAN_SHIFT      8
#define HUE_GRAN           256
#define SATURATION_GRAN_SHIFT   8
#define SATURATION_GRAN        256
*/
/*
#define HUE_GRAN_SHIFT      6
#define HUE_GRAN           64
#define SATURATION_GRAN_SHIFT   6
#define SATURATION_GRAN        64
*/

#define HUE_GRAN_SHIFT  5
#define HUE_GRAN  32
#define SATURATION_GRAN_SHIFT  5
#define SATURATION_GRAN  32


#if SDL_BYTEORDER == SDL_LIL_ENDIAN
   #define SWAP16(X)    (X)
   #define SWAP32(X)    (X)
#else
   #define SWAP16(X)    SDL_Swap16(X)
   #define SWAP32(X)    SDL_Swap32(X)
#endif

#define  b32Tob32p(a)  ( (((a)&0x00ff0000)<<8) | (((a)&0x0000ff00)<<5) | (((a)&0x000000ff)<<2) |  (((a)&0x30000000)>>6) | (((a)&0x0c000000)>>15) | (((a)&0x03000000)>>24) )
#define  b32pTob32(a)  ( (((a)&0xff000000)>>8) | (((a)&0x001fe000)>>5) | (((a)&0x000003fc)>>2) |  (((a)&0x00c00000)<<6) | (((a)&0x00001800)<<15) | (((a)&0x00000003)<<24) )

void initDraw();

SDL_Surface* setVideo( int bpp = 32, int winXSize = -1, int winYSize = -1, bool fullScreen = false, string caption = "Application", SDL_Surface* icon = NULL );
    // Returns pointer to video surface if successful in setting up graphics display
    // Returns a null pointer if an error has ocurred

Uint32 rand32();

double randFloat();

double randdouble();


class StepTimer
{

    Uint32 alarmTime;
    Uint32 sleepTime;

public:

    inline StepTimer( double sec = 1./30. );

    inline void set();

    inline void set( double sec );

    inline void setInc();

    inline void setInc( double sec );

    inline void wait();

    inline void inc();

    inline bool getSignal();

};


/// TODO /// make this PNG's and add a screenshot check that loads filename table for directory and finds the lowest available ss number
void screenShot( string fileName );

inline void clearSurface( SDL_Surface *surface );

void restoreVideo();

inline void lockSurface( SDL_Surface* surface );

inline void unlockSurface( SDL_Surface* surface );

void setPalette8( SDL_Surface* surface, SDL_Color* palette );

extern Uint32 PALETTE_TABLE[SATURATION_GRAN][HUE_GRAN][256];

extern Uint32 PALETTE_TABLE_P[SATURATION_GRAN][HUE_GRAN][1024];

extern Uint32* GRAY_PALETTE;

void setPaletteTable();

inline Uint32* getLightnessColor32( int saturation, int hue );

extern SDL_Color PALETTE_BUFFER[256];

SDL_Color* getLightnessPalette( int saturation, int hue );


inline StepTimer::StepTimer( double sec )
{
    sleepTime = (unsigned long) (1000.*sec);
    alarmTime = SDL_GetTicks();
}

inline void StepTimer::set()
{
    alarmTime = SDL_GetTicks();
}

inline void StepTimer::set( double sec )
{
    sleepTime = (unsigned long) (1000.*sec);
    alarmTime = SDL_GetTicks();
}

inline void StepTimer::setInc()
{
    alarmTime = SDL_GetTicks();
    alarmTime += sleepTime;
}

inline void StepTimer::setInc( double sec )
{
    sleepTime = (unsigned long) (1000.*sec);
    alarmTime = SDL_GetTicks();
    alarmTime += sleepTime;
}

inline void StepTimer::wait()
{

    while( SDL_GetTicks()<alarmTime )
        ;
    alarmTime += sleepTime;
}

inline void StepTimer::inc()
{
    alarmTime += sleepTime;
}

inline bool StepTimer::getSignal()
{
    return SDL_GetTicks()>=alarmTime;
}

inline void clearSurface( SDL_Surface *surface )
{
    SDL_FillRect(surface, NULL, 0x00000000);
}

inline void lockSurface( SDL_Surface* surface )
{
    if( SDL_MUSTLOCK(surface) )
        SDL_LockSurface(surface);
}

inline void unlockSurface( SDL_Surface* surface )
{
    if( SDL_MUSTLOCK(surface) )
        SDL_UnlockSurface(surface);
}

inline Uint32* getLightnessColor32( int saturation, int hue )
{
#ifndef NDEBUG
    if( (unsigned int) (saturation>>(8-SATURATION_GRAN_SHIFT)) > (unsigned int) SATURATION_GRAN
      || (unsigned int) (hue>>(8-HUE_GRAN_SHIFT)) > (unsigned int) HUE_GRAN ) {

        cerr << "Saturation and/or hue value out of range: " << saturation << ", " << hue << endl;
        exit(1);

    }
#endif
    return PALETTE_TABLE[saturation>>(8-SATURATION_GRAN_SHIFT)][hue>>(8-HUE_GRAN_SHIFT)];
}

inline Uint32* getLightnessColor32p( int saturation, int hue )
{
#ifndef NDEBUG
    if( (unsigned int) (saturation>>(8-SATURATION_GRAN_SHIFT)) > (unsigned int) SATURATION_GRAN
      || (unsigned int) (hue>>(8-HUE_GRAN_SHIFT)) > (unsigned int) HUE_GRAN ) {

        cerr << "Saturation and/or hue value out of range: " << saturation << ", " << hue << endl;
        exit(1);

    }
#endif
    return PALETTE_TABLE_P[saturation>>(8-SATURATION_GRAN_SHIFT)][hue>>(8-HUE_GRAN_SHIFT)];
}


#endif /* _DRAWCORE_HPP_ */
