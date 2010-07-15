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
#include <ios>

//#define TEST_MODE
//#define FPS_TEST

#include "drawStar.hpp"
#include "binfile.hpp"

using namespace std;

const int STEP_TOTAL = 400;
const int MOVIE_FRAMES = 1500;

class nBodyFile
{

private:

    ifstream fstrm;
    string fileName;
    int starTotal;
    bool done, binFlag;
    stellarClass* starType;

public:

    void reset()
    {
        if( fstrm.is_open() )
            fstrm.close();
        if( binFlag ) {
            fstrm.open(fileName.c_str(), ios::in|ios::binary);
            starTotal = fileGetIntBin(fstrm);
        }
        else {
            fstrm.open(fileName.c_str());
            starTotal = fileGetInt(fstrm);
        }
        done = false;
    }

    nBodyFile( string fileName, stellarClass &starType, bool binFlag = true )
    {
        this->starType = &starType;
        this->binFlag = binFlag;
        this->fileName = fileName;
        reset();
    }

    int getStarTotal() { return starTotal; }

    bool readStars( starField& stream )

        // Reads next step in 'fstrm' into stream data
        // Returns true if another step exists, false if this is the last step in the file

    {

        if( done )
            return false;

cerr << "DEBUG: Location " << (unsigned long long) fstrm.tellg() << endl;

        // Confirm number of dimensions (must be 3)
        int dimensions;
        if( binFlag )
            dimensions = fileGetIntBin(fstrm);
        else
            dimensions = fileGetInt(fstrm);

        if( dimensions!=3 ) {
            cerr << "Location " << fstrm.tellg() << ": ";
            cerr << "Error reading file - dimensions must be 3 (was " << dimensions << ").\n" ;
            exit(1);
        }

        // Get elapsed time at current step in billions of years
        double elapsedTime;
        if( binFlag )
            elapsedTime = fileGetFloatBin(fstrm);
        else
            elapsedTime = fileGetDouble(fstrm);

        // Eat garbage (unknown source) - only in 19gb file
/*      if( binFlag )
            for( int i = 0; i<starTotal; i++ )
                fileGetFloatBin(fstrm);
*/
        // Get star positions and velocity vectors at current step
        float lineArg[3];
        for( int i = 0; i<starTotal; i++ ) {

            if( binFlag )
                fileGetFloatArrayBin(fstrm, 3, lineArg);
            else
                fileGetFloatArray(fstrm, 3, lineArg);
            stream.add(lineArg[0], lineArg[1], lineArg[2], 100., *starType);

        }

        // Skip over velocities since they are not used
        for( int i = 0; i<starTotal; i++ )
            if( binFlag )
                fileGetFloatArrayBin(fstrm, 3, lineArg);
            else
                fileGetFloatArray(fstrm, 3, lineArg);

        // Check to see if there is another step
        if( fstrm.eof() ) {
            fstrm.close();
            done = true;
        }
        else {

            // Check current star number against old number - exit if they don't match
            int starVerify;
            if( binFlag )
                starVerify = fileGetIntBin(fstrm);
            else
                starVerify = fileGetInt(fstrm);
            if( starTotal != starVerify ) {
                cerr << "Location " << fstrm.tellg() << ": ";
                cerr << "Error reading file - number of stars does not match initial value (" << starVerify << " != " << starTotal << ").\n";
                exit(1);
            }

        }

        return true;

    }

};

class wedgeFile
{

private:

    ifstream fstrm;
    string fileName;
    int starTotal;
    bool binFlag;
    stellarClass* starType;

public:

    void reset()
    {
        starTotal = 0;
    }

    wedgeFile( stellarClass &starType, bool binFlag = false )
    {
        this->starType = &starType;
        this->binFlag = binFlag;
        reset();
    }

    int getStarTotal() { return starTotal; }

    int getStarTotal( string fileName )
    {

        int starTotal;

        if( binFlag ) {
            fstrm.open(fileName.c_str(), ios::in|ios::binary);
            starTotal = fileGetIntBin(fstrm);
        }
        else {
            fstrm.open(fileName.c_str());
            starTotal = fileGetInt(fstrm);
        }
        fstrm.close();

        return starTotal;

    }

    void readStars( string fileName, starField& field )

        // Reads next step in 'fstrm' into stream data
        // Returns true if another step exists, false if this is the last step in the file

    {

        int starTotal;

        if( binFlag ) {
            fstrm.open(fileName.c_str(), ios::in|ios::binary);
            starTotal = fileGetIntBin(fstrm);
        }
        else {
            fstrm.open(fileName.c_str());
            starTotal = fileGetInt(fstrm);
        }

        this->starTotal += starTotal;

        // Get star positions
        float lineArg[3];
        for( int i = 0; i<starTotal; i++ ) {

            if( binFlag )
                fileGetFloatArrayBin(fstrm, 3, lineArg);
            else
                fileGetFloatArray(fstrm, 3, lineArg);

            double l = lineArg[0];
            double b = lineArg[1];
            double r = lineArg[2];

            double x = r*cos(b*TRIG_DEG_TO_RAD)*cos(l*TRIG_DEG_TO_RAD);
            double y = r*cos(b*TRIG_DEG_TO_RAD)*sin(l*TRIG_DEG_TO_RAD);
            double z = r*sin(b*TRIG_DEG_TO_RAD);

            field.add(x, y, z, 100., *starType);

            if( fstrm.eof() ) {
                fstrm.close();
                if( starTotal!=i ) {
                    cerr << "Number of stars indicated in wedge file does not match data." << endl;
                    exit(1);
                }
            }

        }

        /// TODO /// Check to see if there is more data in the file (use a look ahead perhaps to avoid doubling the error checking)
/*      if( !fstrm.eof() ) {
            fstrm.close();
            cerr << "Number of stars indicated in wedge file does not match data." << endl;
            exit(1);
        }
*/

        fstrm.close();

    }

};

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
