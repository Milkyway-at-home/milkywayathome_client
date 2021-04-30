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

#include "drawcore.hpp"


void initDraw()
{
    atexit(SDL_Quit);
    if( SDL_Init(SDL_INIT_VIDEO)==-1 ) {
        std::cerr << "Unable to set up SDL." << std::endl;
        exit(1);
    }
}

SDL_Surface* setVideo( int bpp, int winXSize, int winYSize, bool fullScreen, string caption, SDL_Surface* icon )
    // Returns pointer to video surface if successful in setting up graphics display
    // Returns a null pointer if an error has ocurred
{

    if( bpp!=32 && bpp!=8 ) {
        cerr << "Invalid resolution requested: " << bpp << " bpp\n";
        exit(1);
    }

    atexit(restoreVideo);

    // Choose resolution
    const SDL_VideoInfo *vInfo = SDL_GetVideoInfo();

    if( winXSize<0 || winYSize<0 ) {

        if( !fullScreen ) {
            winXSize = vInfo->current_w*3/4;
            winYSize = vInfo->current_h*3/4;
        }
        else {
            winXSize = vInfo->current_w;
            winYSize = vInfo->current_h;
        }

    }

    // Set up window icon and caption (this will not generally be seen in full screen mode) - must be done before set video mode call
    SDL_WM_SetCaption(caption.c_str(), caption.c_str());
    if( icon != NULL )
        SDL_WM_SetIcon(icon, NULL);

    // Initialize mode
    Uint32 flags = DISPLAY_FLAGS;
    if( fullScreen )
        flags |= SDL_FULLSCREEN;
    SDL_Surface *display = SDL_SetVideoMode(winXSize, winYSize, bpp, flags);

    // Turn off mouse

    if( fullScreen )
        SDL_ShowCursor(0);

    if( bpp==8 ) {
        // Set default 8bpp palette
        int saturation = 64;
        int hue = 170;
        setPalette8(display, getLightnessPalette(saturation, hue));
    }

    setPaletteTable();

    return display;

}

Uint32 rand32()
{
    Uint32 r = 0;
    for( int i = 0; i<12; i++ ){
        r <<= (i%4)+1;
        r ^= rand()/(i+1);
    }
    return r;
}

double randFloat()
{

    double rVal = (double) rand32();
    return rVal/4294967296.;

}

double randDouble()
{

    double rVal = (double) rand32();
    rVal *= 4294967296.;
    rVal += (double) rand32();
    return rVal/1.8446744073709551616e+19;

}

/// TODO /// make this PNG's and add a screenshot check that loads filename table for directory and finds the lowest available ss number
void screenShot( string fileName )
{
    SDL_Surface *display = SDL_GetVideoSurface();
    if( fileName.substr(fileName.length()-4, 4)!=".bmp" )
        fileName += ".bmp";
    SDL_SaveBMP(display, fileName.c_str());
}

void restoreVideo()
{

    // Cleanup - erase screen to avoid palette flicker
    SDL_Surface *display = SDL_GetVideoSurface();
    clearSurface(display);
    SDL_Flip(display);
    StepTimer exitTimer(1./20.);
    exitTimer.inc();
    exitTimer.wait();
    SDL_FreeSurface(display);

}

void setPalette8( SDL_Surface* surface, SDL_Color* palette )
{
    SDL_SetColors(surface, palette, 0, 0x100);
}

Uint32 PALETTE_TABLE[SATURATION_GRAN][HUE_GRAN][256];

Uint32 PALETTE_TABLE_P[SATURATION_GRAN][HUE_GRAN][1024];

Uint32* GRAY_PALETTE = PALETTE_TABLE[0][0];
Uint32* GRAY_PALETTE_P = PALETTE_TABLE_P[0][0];

void setPaletteTable()
{
    SDL_Surface *display = SDL_GetVideoSurface();
    Uint32 rmask, gmask, bmask, amask;
#if SDL_BYTEORDER == SDL_BIG_ENDIAN
    rmask = 0xff000000;
    gmask = 0x000ff000;
    bmask = 0x000000ff;
    amask = 0x00000000;
#else
    rmask = 0x000000ff;
    gmask = 0x000ff000;
    bmask = 0xff000000;
    amask = 0x00000000;
#endif
    SDL_Surface * surface32p = SDL_CreateRGBSurface(SDL_SWSURFACE, 1, 1, 32, rmask, gmask, bmask, amask);
    for( int saturation = 0; saturation<SATURATION_GRAN; saturation++ )
        for( int hue = 0; hue<HUE_GRAN; hue++ )
            for( int lightness = 0; lightness<256; lightness++ ) {
                Uint8 r, g, b;
                hslToRgb(TRIG_2PI*hue/double(HUE_GRAN), saturation/double(SATURATION_GRAN), lightness/255., r, g, b);
                Uint32 pixel = rgbToColor(display, r, g, b);
                PALETTE_TABLE[saturation][hue][lightness] = pixel;
            }
    for( int saturation = 0; saturation<SATURATION_GRAN; saturation++ )
        for( int hue = 0; hue<HUE_GRAN; hue++ )
            for( int lightness = 0; lightness<1024; lightness++ ) {
                double r, g, b;
                hslToRgb(TRIG_2PI*hue/double(HUE_GRAN), saturation/double(SATURATION_GRAN), lightness/1024., r, g, b);
#if SDL_BYTEORDER == SDL_BIG_ENDIAN
                Uint32 pixel = (((Uint32)(b*1023.))<<22) | (((Uint32)(g*1023.))<<11) | ((Uint32)(r*1023.));
#else
                Uint32 pixel = (((Uint32)(r*1023.))<<22) | (((Uint32)(g*1023.))<<11) | ((Uint32)(b*1023.));
#endif
                PALETTE_TABLE_P[saturation][hue][lightness] = b32pTob32(pixel);
            }
    SDL_FreeSurface(surface32p);
}

SDL_Color PALETTE_BUFFER[256];

/// TODO /// It may be valid to simply type cast here
SDL_Color* getLightnessPalette( int saturation, int hue )
{
    Uint8 r, g, b;
    for( int lightness = 0; lightness<256; lightness++ ) {
        hslToRgb(hue/256.*TRIG_2PI, saturation/256., lightness/255., r, g, b);
        PALETTE_BUFFER[lightness].r = r;
        PALETTE_BUFFER[lightness].g = g;
        PALETTE_BUFFER[lightness].b = b;
        PALETTE_BUFFER[lightness].unused = 0;
    }
    return PALETTE_BUFFER;
}

