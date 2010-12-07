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
/*
    if( args<2 ){
        printf("Usage: ./mwdemo wedge_file_name\n");
        return 1;
    }
*/

    cout << "Initializing draw routines\n" << flush;

    initDraw();

    float diameter = 7.;
    float fps = 30.;
    int bpp = 32;

    string fileName = "noU_Pal5.txt";
    if( args>1 )
        fileName = argv[1];

    // Read in data
    cout << "Reading data file\n" << flush;
    HaloField* cluster;
    cluster = readStarFileRaDecGru( fileName.c_str() );

    // Create display
    cout << "Setting up window\n" << flush;
    FieldAnimation sim(bpp, fps, false, "MilkyWay@Home Screensaver Demo", "milkyway.bmp");

    // Read in galaxy
    //cout << "Reading/generating galaxy\n" << flush;
    //ImagePlot imagePlot("eso32.bmp", 25000, 30.*1.18, .3);

    cout << "Adding stars to field...\n" << flush;
    sim.add(cluster, diameter);
    //sim.add(imagePlot.getField(), 20.);

    cout << "Changing settings...\n" << flush;
    sim.showCamera();
    sim.cv->setFocusPoint(sim.cv->getFocusPoint(100., 45., 30.), 0.);

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
