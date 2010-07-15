/****************************************************************************
 *                                                                          *
 *  Copyright (C) 2010 Shane Reilly, Ben Willet, Matthew Newberg, Heidi      *
 *  Newberg, Malik Magdon-Ismail, Carlos Varela, Boleslaw Szymanski, and     *
 *  Rensselaer Polytechnic Institute                                         *
 *                                                                          *
 *  This file is part of Milkway@Home.                                      *
 *                                                                          *
 *  Milkyway@Home is free software: you can redistribute it and/or modify   *
 *  it under the terms of the GNU General Public License as published by    *
 *  the Free Software Foundation, either version 3 of the License, or       *
 *  (at your option) any later version.                                     *
 *                                                                          *
 *  Milkyway@Home is distributed in the hope that it will be useful,        *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *  GNU General Public License for more details.                            *
 *                                                                          *
 *  You should have received a copy of the GNU General Public License       *
 *  along with Milkyway@Home. If not, see <http://www.gnu.org/licenses/>.   *
 *                                                                          *
 *  Shane Reilly                                                            *
 *  reills2@cs.rpi.edu                                                      *
 *                                                                          *
 ****************************************************************************/

#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <fstream>

//#define TEST_MODE
//#define FPS_TEST

#include "drawStar.hpp"
#include "demofile.hpp"

using namespace std


int main( int args, char **argv )
{

    if( args<3 ){
        printf("Usage: ./viewwedge file_name magnitude_mult\n");
        return 1;
    }

    wedgeFile wf(argv[1]);
    wf.read
   
   stellarClass galBulge(127.5, 2.*1.5/1.*LumMult, 96);
   star galCenter(8.5, 0., 0., 1., galBulge);

   stellarClass lrgStar(3.5, 20.*LumMult, 96);
   star sol(0., 0., 0., 1., lrgStar);

   // Calculate cluster center
   double xCtr = xs/starTotal;
   double yCtr = ys/starTotal;
   double zCtr = zs/starTotal;

//double ax = atan2((yCtr-cy), (zCtr-cz));
//double az = atan2((yCtr-cy), (xCtr-cx));
ay = 0.;/// TODO ///ay = atan2((zCtr-cz), (xCtr-cx));
ax = 0.;/// TODO ///atan2((zCtr-cz), (yCtr-cy));

//star ctrLoc(xCtr, yCtr, zCtr, 1000.);

//ax = 0.; az = 0.;

   // Calculate cluster standard deviation
   xs = 0., ys = 0., zs = 0.;
/*   for( int i = 0; i<starTotal; i++ )
   {
      double sq = (xCtr-cluster[i]->x);
      xs += sq*sq;
      sq = (xCtr-cluster[i]->y);
      ys += sq*sq;
      sq = (xCtr-cluster[i]->z);
      zs += sq*sq;
   }
   double xsd = sqrt(xs/starTotal);
   double ysd = sqrt(ys/starTotal);
   double zsd = sqrt(zs/starTotal);
*///////////////
//double temp = .0000025;  // *** Stub to determine 2d zoom

//cout << endl << xCtr << "\t" << yCtr << "\t" << zCtr << " :\t" << temp << endl;
//cout << endl << xsd << "\t" << ysd << "\t" << zsd << endl;

   // Set up routines
   SDL_Surface* display = setVideo();
   SDL_Event event;

   float sec = 1./60.;
   unsigned int tock = int(float(CLOCKS_PER_SEC)*sec);
   unsigned int tick = clock();

//int ii = 0;
   unsigned int frameCounter = 0;
   unsigned int timeCounter = SDL_GetTicks();

   bool setToggle = false;

    return 0;

}
