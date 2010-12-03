/*****************************************************************************
 *                                                                           *
 *  Copyright (C) 2010 Shane Reilly and Rensselaer Polytechnic Institute     *
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

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdlib>

// Uncomment this line to throw an error if pixel accesses are out of bounds for non-clipping calls
//#define TEST_MODE

#include "drawhalo.hpp"

using namespace std;

class star
{

    float xv, yv, xp, yp, xc, yc, mag, grav;

    public:

    star( float x, float y, float angle, float velocity, float xCenter, float yCenter, float gravity, float magnitude )
    {
        xv = cos(angle)*velocity;
        yv = sin(angle)*velocity;
        xc = xCenter;
        yc = yCenter;
        xp = x;
        yp = y;
        mag = magnitude;
        grav = gravity;
    }

    star( int xSize, int ySize, float magnitude )
    {
        float x = xSize/2.+(rand32()%(10*xSize))/20.;
        float y = ySize/2.+(rand32()%(10*ySize))/20.;
        float angle = 2.2+(rand32()%28)/100.;
        float velocity = 2.+(rand32()%25)/5.;
        float xCenter = xSize/2.;
        float yCenter = ySize/2.;
        float gravity = 1./32.;

        xv = cos(angle)*velocity;
        yv = sin(angle)*velocity;
        xc = xCenter;
        yc = yCenter;
        xp = x;
        yp = y;
        mag = magnitude;
        grav = gravity;
    }

    void draw( SDL_Surface* surface, HaloType* sc )
    {
        sc->draw(surface, xp, yp, mag);
    }

    void fullDraw( SDL_Surface* surface, BlurCircle* bc, float radius, float intensity )
    {
        bc->draw(surface, xp, yp, radius, intensity);
    }

    void pixelDraw( SDL_Surface* surface )
    {
        putPixelClip8(surface, int(xp), int(yp), min(127, int(80.*mag)));
    }

    void step()
    {
        // Calculate change in position
        xp += xv;
        yp += yv;

        // Calculate force vector
        float xd = xc-xp;
        float yd = yc-yp;
        float d = sqrt(xd*xd+yd*yd);

        float xvp = xd/d;
        float yvp = yd/d;

        // Calculate change in velocity
        xv += grav*xvp;
        yv += grav*yvp;

    }

};


int main( int args, char **argv )
{

    // Parse parameters

    if( args<2 ) {
        cout << "Usage: ./startest stars [blur_circle_radius] [luminosity]\n";
        return 1;
    }
    int stars = atoi(argv[1]);
    float radius = 3.5;
    if( args>2 )
        radius = atof(argv[2]);
    float luminosity = .4/radius;
    if( args>3 )
        luminosity = atof(argv[3]);

    // Set up SDL

    initDraw();

    SDL_Surface* display = setVideo(8, -1, -1, true);

    SDL_Event event;

    // Set up routines

    srand(31416);
    HaloType sc(radius, 1., 0, 3);

    // Draw stars

    /// TODO take out polling here and replace with event queue ///

    star** starList = new star*[stars];
    for( int i = 0; i<stars; i++ )
    if( SDL_PollEvent(&event)!=-1 && event.type==SDL_KEYDOWN )
        exit(0);
    else
        starList[i] = new star(display->w, display->h, luminosity);

    unsigned int frameCounter = 0;
    unsigned int timeCounter = SDL_GetTicks();

    StepTimer sTimer(1./60.);
    while( SDL_PollEvent(&event)==-1 || event.type!=SDL_KEYDOWN )
    {
        StepTimer terminate(3.);
        while( sTimer.getSignal() && !terminate.getSignal()) {
            for( int i = 0; i<stars; i++ ){
                starList[i]->step();
            }
        sTimer.inc();
        }

        clearSurface(display);

        if( SDL_MUSTLOCK(display) )
            SDL_LockSurface(display);

        for( int i = 0; i<stars; i++ ) {
            starList[i]->draw(display, &sc);
            //         starList[i]->pixelDraw(display);
            starList[i]->step();
        }
        /*
        while( !(SDL_PollEvent(&event)!=-1 && event.type==SDL_KEYDOWN) )
            ;
        stepTimer s(.3);
        s.wait();
        */
        if( SDL_MUSTLOCK(display) )
            SDL_UnlockSurface(display);

        SDL_Flip(display);
        frameCounter++;
        sTimer.wait();

    }

    cout << 1000.*frameCounter/float(SDL_GetTicks()-timeCounter) << "fsp" << endl;

    screenShot("lastframe.bmp");  /// TODO /// make this keyboard-driven with multi-shot functionality
    sc._drawTest();
    screenShot("startable.bmp");

    // Cleanup - erase screen to avoid palette flicker
    clearSurface(display);
    SDL_Flip(display);
    StepTimer exitTimer(1./20.);
    exitTimer.wait();

    return 0;

}

