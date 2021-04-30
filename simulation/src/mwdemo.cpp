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


HaloField* getLastFrameNBody()
{
    NBodyFile nb("plummer1.bin");
    HaloField* stream = new HaloField(10000);
    nb.readStars(*stream);
    return stream;
}

int main(int argc, const char* argv[])
{
    std::cout << "Initializing draw routines" << std::endl;

    initDraw();

    float diameter = 7.0f;
    float fps = 30.0f;
    int bpp = 32;

    std::string fileName = "stars_82.txt";
    if (argc > 1)
        fileName = argv[1];

    // Read in wedge
    std::cout << "Reading wedge file" << std::endl;
    WedgeFile wf;
    int totalStars = wf.getStarTotal(fileName);
    HaloField wedge(totalStars);
    wf.readStars(fileName, wedge);

    // Create display
    std::cout << "Setting up window" << std::endl;
    FieldAnimation sim(bpp, fps, false, "MilkyWay@Home Screensaver Demo", "milkyway.bmp");

    // Read in galaxy
    std::cout << "Reading/generating galaxy" << std::endl;
    ImagePlot imagePlot("eso32.bmp", 25000, 30.0 * 1.18, .3);

    std::cout << "Adding stars to field..." << std::endl;
    sim.add(&wedge, diameter);
    sim.add(getLastFrameNBody(), diameter);
    sim.add(imagePlot.getField(), 20.);

    std::cout << "Changing settings..." << std::endl;
    sim.showCamera();
    sim.cv->setFocusPoint(sim.cv->getFocusPoint(100.0, 45.0, 30.0), 0.0);

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

    std::cout << "Entering event loop\n" << flush;

    while (true)
    {
        if (sim.pollEvent())
        {
        }

    }

    std::cout << "Done" << std::endl;

    return 0;

}


