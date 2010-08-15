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

SDL_Color lToOStarColor( int l )
{
    SDL_Color c;
    c.r = l*9/10;
    c.g = l*9/10;
    c.b = l;
    c.unused = 0;
    return c;
}

SDL_Color* getOStarPalette( SDL_Color* palette )
{
   for( int i = 0; i<0x80; i++ )
      palette[i] = lToOStarColor((i<<1)+(i>>6));
   return palette;
}

void setOStarPalette( SDL_Surface* surface )
{
    SDL_Color pal[0x80];
    SDL_SetColors(surface, getOStarPalette(pal), 0, 0x80);
}

HaloField* getLastFrameNBody()
{
    NBodyFile nb("../sim/nbody/plummer1.bin");
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
    sim.add(getLastFrameNBody(), diameter*2);
    sim.add(imagePlot.getField(), 10.);

    sim.showCamera();
    sim.cv->setFocusPoint(sim.cv->getFocusPoint(100., 45., 30.), 0.);

//setFocusPoint(Vector3d{-8, 0, 0});
//moveTo(Vector3d{-8, 0, 0}, 0.);

/*
// Find earth-view focus
Vector3d sev;
Vector3d eev;
Vector3d cev;
Vector3d cevRot90 = cev;
cevRot90.rotate(TRIG_2PI/4., , )

// Find random focus
Vector3d randAngle;

double topDist = wedgeStdDev*3;
double rViewDist = topDis*2;

setFocusAngle(0., sev, 0.);
setFocusAngle(0., eev, 10.);   /// TODO /// Buffered

setFocusAngle(topDist, cEvRot90, 1.);   /// TODO /// Buffered

setFocusAngle(RandomAngle);   /// TODO /// Buffered with fade-in

*/

/////////    while( sim.pollDemo() ) ;

    while( true )
    {

        if( sim.pollEvent() ) {
            ;
        }

    }

    return 0;

}
