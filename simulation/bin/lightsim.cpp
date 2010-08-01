/*****************************************************************************
 *                                                                           *
 *  Copyright (C) 2010 Shane Reilly, Ben Willet, Matthew Newberg, Heidi      *
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

#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <ios>

//#define TEST_MODE
//#define FPS_TEST

#include "drawStar.hpp"
#include "demofile.hpp"

using namespace std;

const int STEP_TOTAL = 400;
const int MOVIE_FRAMES = 1500;

int main( int args, char **argv )
{

    // Handle arguments
    if( args<4 ){
        cout << "Usage: ./nbodysim n_body_file_name wedge_file_name star_brightness\n";
        return 0;
    }
    string nFileName = argv[1];
    string wFileName = argv[2];
    double lum = atof(argv[3]);
    double rad = 3.5;

    // Set up environment
    SDL_Surface* display = setVideo();
    SDL_Event event;
    stellarClass fStar(rad, lum, 96);

    stellarClass fStarLow(rad, lum/5., 96);

    camera cv(0.253, 0., 0., 0., 0., 0.);

    // Read in wedges
    wedgeFile wf(fStarLow);
    starField wedges(*display, cv, wf.getStarTotal(wFileName));
    wf.readStars(wFileName, wedges);

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
            if( cv.pollEvent(event) ) {

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
