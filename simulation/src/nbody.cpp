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
#include "imgplot.hpp"
#include "demofile.hpp"

const int STEP_TOTAL = 400;
const int MOVIE_FRAMES = 1500;

int main(int argc, const char* argv[])
{
    std::cout << "Initializing draw routines" << std::endl;
    initDraw();

    float diameter = 3.0f;
    float fps = 60.0f;
    int bpp = 32;
    float lum = 0.5f;

    int totalNBody = 3;

    //string fileName[totalNBody];
    std::string fileName[3];
    fileName[0] = "gd1.stoa";
    fileName[1] = "orphan.stoa";
    fileName[2] = "sgrsim.stoa";

    // Read in files
    std::cout << "Reading N-body files." << std::endl;

    bool binary = false;

    NBodyFile* nBody[totalNBody];
    for (int i = 0; i < totalNBody; i++)
        nBody[i] = new NBodyFile(fileName[i], binary);

    int totalStars = nBody[0]->getStarTotal();
    HaloField* field[totalNBody];
    for (int i = 0; i < totalNBody; i++)
        field[i] = new HaloField(totalStars);

    for (int i = 0; i < totalNBody; i++)
        nBody[i]->readStars(*field[i], lum);

    // Create display
    bool fullScreen = false;
    std::cout << "Setting up window" << std::endl;
    FieldAnimation sim(bpp, fps, fullScreen, "GD-1 + Orphan + Sgr Stream Disruption", "nbody.bmp");

    // Read in galaxy
    std::cout << "Reading/generating galaxy" << std::endl;
    ImagePlot imagePlot("eso32.bmp", 5000, 30.0 * 1.18, 0.3);

    std::cout << "Adding galaxy to field..." << std::endl;
    sim.add(imagePlot.getField(), 20.0);

    for (int i = 0; i<totalNBody; i++)
        sim.add(field[i], diameter);

    // Setup
    sim.showCamera();
    sim.showAxes();
    sim.cv->moveToPoint(Vector3d(0.0, -500.0, 0.0), 0.0, 0.0);

//    double currentTime = 0.;
//    double masterTimeStep = ;

    while (true)
    {
        if (sim.pollEvent())
        {
//          if( field.getTimeStep()>currentTimeStep)
            for (int i = 0; i<totalNBody; i++)
            {
                field[i]->clearField();
                nBody[i]->readStars(*field[i], lum);
            }
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
                    for( float i = 0.; i<STEP_TOTAL; i+=STEP_TOTAL/MOVIE_FRAMES ) {

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


