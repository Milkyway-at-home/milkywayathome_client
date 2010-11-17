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
#include "imgplot.hpp"

using namespace std;


int main( int args, char **argv )
{

    const char* fileName = argv[1];

    if( args<2 ){
//        printf("Usage: ./showstar file_name\n");
//        return 0;
        fileName = "cel2mil_xyz.txt";
    }

    cout << "Initializing draw routines\n" << flush;
    initDraw();

    double diameter = 3.;
    double fps = 30.;
    int bpp = 32;

    // Read in stars
    cout << "Reading star file\n" << flush;
    HaloField* starField = readStarFileXyz(fileName);

    // Create display
    cout << "Setting up window\n" << flush;
    FieldAnimation sim(bpp, fps, false, "Star Field Simulation", "icon.bmp");

    // Read in galaxy
//    cout << "Reading/generating galaxy\n" << flush;
//    ImagePlot imagePlot("eso32.bmp", 25000, 30.*1.18*3261.63626, .3*3261.63626);

    cout << "Adding stars to field...\n" << flush;
    sim.add(starField, diameter);
//    sim.add(imagePlot.getField(), 20.);

    cout << "Changing settings...\n" << flush;
    sim.showCamera();
//    sim.cv->moveToPoint(Vector3d(-26093.0901, 0., 0.), 0., 0.);

    cout << "Entering event loop\n" << flush;

    while( true )
    {

        if( sim.pollEvent() ) {
            ;
        }

    }

    cout << "Done\n" << flush;

    return 0;

}
