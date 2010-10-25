/*****************************************************************************
 *                                                                           *
 *  Copyright (C) 2010 Shane Reilly, Ben Willet, Matthew Newby, Heidi        *
 *  Newberg, Malik Magdon-Ismail, Carlos Varela, Boleslaw Szymanski, and     *
 *  Rensselaer Polytechnic Institute                                         *
 *                                                                           *
 *  This file is part of Milkway@Home.                                       *
 *                                                                           *
 *  Milkyway@Home is free software: you can redistribute it and/or modify    *
 *  it under the terms of the GNU General Public License as published by     *
 *  the Free Software Foundation, either version 3 of the License, or        *
 *  (at your option) any later version.                                      *
 *                                                                           *
 *  Milkyway@Home is distributed in the hope that it will be useful,         *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of           *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the             *
 *  GNU General Public License for more details.                             *
 *                                                                           *
 *  You should have received a copy of the GNU General Public License        *
 *  along with Milkyway@Home. If not, see <http://www.gnu.org/licenses/>.    *
 *                                                                           *
 *  Shane Reilly                                                             *
 *  reills2@cs.rpi.edu                                                       *
 *                                                                           *
 *****************************************************************************/

#include <cstdlib>

#include "drawhalo.hpp"
#include "demofile.hpp"

using namespace std;

const int STEP_TOTAL = 400;
const int MOVIE_FRAMES = 1500;

int main( int args, char **argv )
{

    // Handle arguments
    if( args<4 ){
        cout << "Usage: ./nbody n_body_file_name wedge_file_name star_brightness\n";
        return 0;
    }
    double diameter = 7.;
    double fps = 30.;
    int bpp = 32;

    string fileName = "stars_82.txt";
    if( args>1 )
        fileName = argv[1];

    // Read in wedge
    WedgeFile wf;
    int totalStars = wf.getStarTotal(fileName);
    HaloField wedge(totalStars);
    wf.readStars(fileName, wedge, .01);

    // Create display
    FieldAnimation sim(bpp, fps);

    // Read in galaxy
    ImagePlot imagePlot("eso32.bmp", 25000, 30.*1.18, .3);

    sim.add(&wedge, diameter);
    sim.add(getLastFrameNBody(), diameter);
    sim.add(imagePlot.getField(), 20.);

    sim.showCamera();
    sim.cv->setFocusPoint(sim.cv->getFocusPoint(100., 45., 30.), 0.);

    while( true )
    {

        if( sim.pollEvent() ) {
            ;
        }

    }

    return 0;

}

/*
    while( true ) {

        // Read file data and display
        nBodyFile nb(nFileName, fStar);
        starField stream(*display, cv, nb.getStarTotal());

        // Initialize stream
        int step = 0;

        while( step<STEP_TOTAL )
        {

            stream.clearStars();
            nb.readStars(stream);
            step++;

            while( true ) {

                // Display iteration
                stream.clear();
                stream.draw();
                SDL_Flip(display);
                if( cv.pollEvent(event) )
                    if( event.key.keysym.sym==SDLK_SPACE )
                        break;

            }

        }

        while( true ) {

            // Display iteration
            stream.clear();
            stream.draw();
            wedges.draw();
            SDL_Flip(display);
            if( ??.pollEvent(event) ) {

                if( event.key.keysym.sym==SDLK_BACKSPACE )
                    break;
                else if( event.key.keysym.sym==SDLK_r ) {

                    // Record consecutive screenshots
                    /// TODO /// It is possible to use partial transparency here to blend frames
                             ///   on a time-step boundary if the number of frames is not evenly
                             ///   divisible into the total number of steps
                    int frame = 0;
                    int lastFrame = 0;
                    nb.reset();
                    stream.clearStars();
                    nb.readStars(stream);
                    for( double i = 0.; i<STEP_TOTAL; i+=STEP_TOTAL/MOVIE_FRAMES ) {

                        int skipFrame = (int) i-lastFrame;
                        lastFrame = (int) i;
                        for( int j = 0; j<skipFrame; j++ ) {
                            stream.clearStars();
                            nb.readStars(stream);
                        }
                        stream.draw();
                        SDL_Flip(display);
                        string fileName = "video_frame_" + intToString(frame++, 10) + ".bmp";
                        screenShot(fileName);

                    }

                }

            }

        }

    }

    return 0;

}
*/