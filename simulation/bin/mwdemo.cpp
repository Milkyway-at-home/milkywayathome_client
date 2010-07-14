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
#include "sdlwrap.hpp"

using namespace std;

const double ANGLE_CHANGE = TRIG_2PI/750.;

int WinXSize = 1280;
int WinYSize = 800;

int XCenter = WinXSize>>1;
int YCenter = WinYSize>>1;

double Radius = 3.5;
double LumMult;
double cx = 0.;
double cy = 0.;
double cz = 0.;
double ax = 0.;
double ay = 0.;
double view = TRIG_2PI/16.;  // 22.5 degrees visible from left to right on screen
double zoom = WinXSize/tan(view);
double accInc = .0001;

double xv = 0., yv = 0., zv = 0.;

void accel( double acc )
{

   double axr = ax+TRIG_PI;
   double ayr = ay+TRIG_PI;

   double xva, yva, zva;

   xva = -sin(ayr)*acc;
   zva = cos(ayr)*acc;

   yva = sin(axr)*zva;
   zva *= cos(axr);

   xv += xva;
   yv += yva;
   zv += zva;

}

void step()
{
   cx += xv;
   cy += yv;
   cz += zv;
}

class star
{

public:

   double x, y, z, l, m;
   stellarClass* starType;

   star( double x, double y, double z, double l, stellarClass &starType )
   {
      this->x = x;
      this->y = y;
      this->z = z;
      this->l = l;
      this->starType = &starType;
   }

   void draw( SDL_Surface* surface )
   {

      double xmap, ymap, lmap;
      if( _3dTo2d(ax, ay, cx, cy, cz, x, y, z, l, zoom, xmap, ymap, lmap) )
//      if( _2dTo2d(ax, ay, cx, cy, cz, x, y, z, l, zoom, xmap, ymap, lmap) )
         starType->draw(surface, xmap+XCenter, ymap+YCenter);  /// TODO /// , lmap);
//putPixelSumClip7(surface, int(xmap+XCenter), int(ymap+YCenter), 63);  /// TODO /// , lmap);

//if((Uint32)(xmap+XCenter)<1280ul && (Uint32)(ymap+YCenter)<800ul)
//cout << ax << ", " << ay << ":" << x-cx << ", " << y-cy << endl;
 //     }
   }

};


int main( int args, char **argv )
{

   if( args<3 ){
      printf("Usage: ./viewwedge file_name magnitude_mult\n");
      return 1;
   }

   LumMult = atof(argv[2]);

   if( args>3 ){
     Radius = atof(argv[3]);
   }

   srand(0);
   stellarClass fStar(Radius, LumMult, 96);

int starTotal2 = 0;

   ifstream wedgeFile;
   wedgeFile.open(argv[1]);
   string line;
   getline(wedgeFile, line);
   int starTotal = atoi(line.c_str());
   star** cluster = new star*[starTotal];
   star** cluster2 = new star*[starTotal];
   int i = 0;
   int i2 = 0;
   double xs = 0., ys = 0., zs = 0.;
   while( !wedgeFile.eof() )
   {
      getline(wedgeFile, line);
      const char* str = line.c_str();
      char *nextPtr = NULL;
      double M = 4.2;                                             // Absolute magnitude
      if( line!="" )
      {

         // Read in star
         double ra = strtod(str, &nextPtr);                       // Right ascension
         double dec = strtod((const char*)nextPtr, &nextPtr);     // Declination
         double m = strtod((const char*)nextPtr, &nextPtr);       // G-filter apparent magnitude
//         double r = pow(10., (m-M)/5.) / 100.;               // 'r' is radial distance in kiloparsecs
         double r = m;               // 'r' is radial distance in kiloparsecs

/*bool useSet2;
// Filter data
if( m<16. || m>22.) {
starTotal--;
useSet2 = false;
}
else
{
starTotal2++;
useSet2 = true;
}
*/
         // RA/Dec to cartesian conversion
         double t0 = cos(dec*TRIG_DEG_TO_RAD)*cos(ra*TRIG_DEG_TO_RAD);
         double t1 = cos(dec*TRIG_DEG_TO_RAD)*sin(ra*TRIG_DEG_TO_RAD);
         double t2 = sin(dec*TRIG_DEG_TO_RAD);

         double mm[3][3] =
            {   { -.06699 , -.87276, -.48354 } ,
                {  .49273 , -.45035,  .74458 } ,
                { -.86760 , -.18837,  .46020 }  };

         double c0 = t0*mm[0][0] + t1*mm[0][1] + t2*mm[0][2];
         double c1 = t0*mm[1][0] + t1*mm[1][1] + t2*mm[1][2];
         double c2 = t0*mm[2][0] + t1*mm[2][1] + t2*mm[2][2];

         double x = r*c0;
         double y = r*c1;
         double z = r*c2;

         // Calculate luminosity relative to sun
         // ************** STUB get std dev of abs mag here and apply randomly
         double l = pow(10.,-M*2/5+2);

         // Track center
         xs += x;
         ys += y;
         zs += z;

         // Add star to cluster
//         if( useSet2 )
//           cluster2[i2++] = new star(x, y, z, l, fStar);
//         else
           cluster[i++] = new star(x, y, z, l, fStar);

//cout << ra << "\t" << dec << "\t" << m << "\t" << r << "\t(" << x << ", " << y << ", " << z << ") : " << d << endl;

      }

   }
   wedgeFile.close();
//   if( starTotal!=i )
//   {
//      cout << "unable to read file" << endl;
//      exit(1);
//   }

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

   while ( true ) {

      do
      {
#ifndef FPS_TEST

         if( (unsigned int)clock()>(unsigned int)tick )
         {
            tick += tock;
            step();
         }
         else
#endif
         {
      //cout << ii++ << endl;
            clearSurface(display);
            if( SDL_MUSTLOCK(display) )
               SDL_LockSurface(display);

            // Draw stars
/*            if( setToggle )
               for( int i = 0; i<starTotal2; i++ )
                  cluster2[i]->draw(display);
            else
*/               for( int i = 0; i<starTotal; i++ )
                  cluster[i]->draw(display);

   //            galCenter.draw(display);
   //         sol.draw(display);
   //ctrLoc.draw(display, fStar);

            if( SDL_MUSTLOCK(display) )
               SDL_UnlockSurface(display);
            frameCounter++;
            SDL_Flip(display);
            step();
   #ifndef FPS_TEST
            while( (unsigned int)clock()<tick )
               ;
            tick += tock;
   #endif
         }

      } while( SDL_PollEvent(&event)!=-1 && event.type!=SDL_KEYDOWN );

      switch(event.type) {

      case SDL_KEYDOWN:

         switch( event.key.keysym.sym ) {

            case SDLK_LEFT:
               ay -= ANGLE_CHANGE;
               if( ay<0. )
                  ay += TRIG_2PI;
               break;
            case SDLK_RIGHT:
               ay += ANGLE_CHANGE;
               if( ay>=TRIG_2PI )
                  ay -= TRIG_2PI;
               break;
            case SDLK_UP:
               ax -= ANGLE_CHANGE;
               if( ax<0. )
                  ax += TRIG_2PI;
               break;
            case SDLK_DOWN:
               ax += ANGLE_CHANGE;
               if( ax>=TRIG_2PI )
                  ax -= TRIG_2PI;
               break;
            case SDLK_PAGEUP:
               accel(accInc);
               break;
            case SDLK_PAGEDOWN:
               accel(-accInc);
               break;
            case SDLK_END:
               xv = 0.;
               yv = 0.;
               zv = 0.;
               break;
            case SDLK_F1:
               setToggle = false;
               break;
            case SDLK_F2:
               setToggle = true;
               break;
            case SDLK_1:
               accInc = .0001;
               break;
            case SDLK_2:
               accInc = .001;
               break;
            case SDLK_3:
               accInc = .01;
               break;
            case SDLK_4:
               accInc = .1;
               break;
            case SDLK_5:
               accInc = 1.;
               break;
            case SDLK_6:
               accInc = 10.;
               break;
            case SDLK_7:
               accInc = 100.;
               break;
            case SDLK_8:
               accInc = 1000.;
               break;
            case SDLK_9:
               accInc = 10000.;
               break;
            case SDLK_ESCAPE:
               cout << 1000.*frameCounter/float(SDL_GetTicks()-timeCounter) << "fps" << endl;
               exit(0);
               break;
            case SDLK_PRINT:
               screenShot("screenshot.bmp");  /// TODO ///  Make this save incrementally
               break;
            default:
               break;

         }

      }

   }

    return 0;

}
