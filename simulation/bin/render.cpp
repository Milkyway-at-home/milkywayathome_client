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


HaloField* getLastFrameNBody()
{
    NBodyFile nb("plummer1.bin");
    HaloField* stream = new HaloField(10000);
    nb.readStars(*stream);
    return stream;
}


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

    string fileName = "stars_82.txt";
    if( args>1 )
        fileName = argv[1];

    // Read in wedge
    cout << "Reading wedge file\n" << flush;
    WedgeFile wf;
    int totalStars = wf.getStarTotal(fileName);
    HaloField wedge(totalStars);
    wf.readStars(fileName, wedge);

    // Create display
    cout << "Setting up window\n" << flush;
    FieldAnimation sim(bpp, fps, false, "MilkyWay@Home Screensaver Demo", "milkyway.bmp");

    // Read in galaxy
    cout << "Reading/generating galaxy\n" << flush;
    ImagePlot imagePlot("eso32.bmp", 25000, 30.*1.18, .3);

    cout << "Adding stars to field...\n" << flush;
    sim.add(&wedge, diameter);
//    sim.add(getLastFrameNBody(), diameter);
//    sim.add(imagePlot.getField(), 20.);

    cout << "Changing settings...\n" << flush;
    sim.showCamera();
    sim.cv->setFocusPoint(sim.cv->getFocusPoint(100., 45., 30.), 0.);

//setFocusPoint(Vector3d{-8, 0, 0});
//moveTo(Vector3d{-8, 0, 0}, 0.);

/*
 Find earth-view focus
Vector3d sev;
Vector3d eev;
Vector3d cev;
Vector3d cevRot90 = cev;
cevRot90.rotate(TRIG_2PI/4., , )

// Find random focus
Vector3d randAngle;

float topDist = wedgeStdDev*3;
float rViewDist = topDis*2;

setFocusAngle(0., sev, 0.);
setFocusAngle(0., eev, 10.);   /// TODO /// Buffered

setFocusAngle(topDist, cEvRot90, 1.);   /// TODO /// Buffered

setFocusAngle(RandomAngle);   /// TODO /// Buffered with fade-in

*/

/////////    while( sim.pollDemo() ) ;

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
