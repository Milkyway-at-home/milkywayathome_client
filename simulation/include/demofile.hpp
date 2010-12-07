/*****************************************************************************
 *                                                                           *
 *  Copyright (C) 2010 Shane Reilly, Ben Willet, Matthew Newby, Heidi        *
 *  Newberg, Malik Magdon-Ismail, Carlos Varela, Boleslaw Szymanski, and     *
 *  Rensselaer Polytechnic Institute                                         *
 *                                                                           *
 *  This file is part of the MilkyWay@Home Project.                          *
 *                                                                           *
 *  This program is free software: you can redistribute it and/or modify     *
 *  it under the terms of the GNU General Public License as published by     *
 *  the Free Software Foundation, either version 3 of the License, or        *
 *  (at your option) any later version.                                      *
 *                                                                           *
 *  This program is distributed in the hope that it will be useful,          *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of           *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the             *
 *  GNU General Public License for more details.                             *
 *                                                                           *
 *  You should have received a copy of the GNU General Public License        *
 *  along with this program. If not, see <http://www.gnu.org/licenses/>.     *
 *                                                                           *
 *  Shane Reilly                                                             *
 *  reills2@cs.rpi.edu                                                       *
 *                                                                           *
 *****************************************************************************/

#ifndef _DEMOFILE_HPP_
#define _DEMOFILE_HPP_

#include <iomanip>

#include "binfile.hpp"
#include "astroconv.h"
#include "drawhalo.hpp"

using namespace std;

class NBodyFile
{

private:

    ifstream fstrm;
    string fileName;
    int starTotal;
    bool done, binFlag;
    double timeStep;

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
        timeStep = 0.;
    }

    NBodyFile( string fileName, bool binFlag = true )
    {
        this->binFlag = binFlag;
        this->fileName = fileName;
        reset();
    }

    int getStarTotal() { return starTotal; }

    double getTimeStep()
    {
        return timeStep;
    }

    bool readStars( HaloField& stream, double lum = .5 )

        // Reads next step in 'fstrm' into stream data
        // Returns true if another step exists, false if this is the last step in the file

    {

        if( done )
            return false;

        stream.clearField();

//cerr << "DEBUG: Location " << (unsigned long long) fstrm.tellg() << endl;

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
        if( binFlag )
            timeStep = fileGetFloatBin(fstrm);
        else
            timeStep = fileGetFloat(fstrm);
/*
        // Eat masses
        if( binFlag )
            for( int i = 0; i<starTotal; i++ )
                fileGetDoubleBin(fstrm);
*/
        // Get star positions and velocity vectors at current step
        float lineArg[3];
        for( int i = 0; i<starTotal; i++ ) {
            if( binFlag )
                fileGetFloatArrayBin(fstrm, 3, lineArg);
            else
                fileGetFloatArray(fstrm, 3, lineArg);
            stream.add(lineArg[0], lineArg[1], lineArg[2], lum, 120, 151);

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

struct WedgeInfo
{
    bool initialized;

    double lCenter;
    double center;
    double rCenter;
    double maxDistance;

    Vector3d crossVector;
    Vector3d depthVector;

};

class WedgeFile
{

private:

    ifstream fstrm;
    string fileName;
    int starTotal;
    bool binFlag;
    WedgeInfo wedgeInfo;

public:

    void reset()
    {
        starTotal = 0;
        wedgeInfo.initialized = false;
    }

    WedgeFile( bool binFlag = false )
    {
        this->binFlag = binFlag;
        reset();
    }

    int getStarTotal() { return starTotal; }

    int getStarTotal( string fileName )
    {

        if( binFlag ) {
            fstrm.open(fileName.c_str(), ios::in|ios::binary);
            starTotal = fileGetIntBin(fstrm);
        }
        else {
            fstrm.open(fileName.c_str());
            starTotal = fileGetInt(fstrm);
        }
        fstrm.close();
//cout << starTotal << endl;
        return starTotal;

    }

    void readStars( string fileName, HaloField& field, double lum = .5, bool removeDuplicates = false )

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
//		cout << starTotal << endl;
        this->starTotal += starTotal;

        // Get star positions
        double lineArg[3];
        double *lc = NULL;
        double *bc = NULL;

        if( removeDuplicates ) {
            lc = new double[starTotal];
            bc = new double[starTotal];
        }
//cout << starTotal << endl << flush;
        int skipTotal = 0;
        for( int i = 0; ; i++ ) {

            if( binFlag )
                fileGetDoubleArrayBin(fstrm, 3, lineArg);
            else
                fileGetDoubleArray(fstrm, 3, lineArg);

            double l = lineArg[0];
            double b = lineArg[1];
            double r = lineArg[2];

            if( removeDuplicates ) {

                lc[i] = l;
                bc[i] = b;

                bool skip = false;
                for( int j = 0; j<i; j++ )
                    if( sqrt((l-lc[j])*(l-lc[j])+(b-bc[j])*(b-bc[j]))<0.001 ) {
//cout << ": " << sqrt((l-lc[i])*(l-lc[i])+(b-bc[i])*(b-bc[i])) << endl;
                        skip = true;
                        skipTotal++;
                        break;
                    }
                if( skip )
                    continue;

            }
//cout << endl;
            double x = r*cos(b*TRIG_DEG_TO_RAD)*cos(l*TRIG_DEG_TO_RAD);
            double y = r*cos(b*TRIG_DEG_TO_RAD)*sin(l*TRIG_DEG_TO_RAD);
            double z = r*sin(b*TRIG_DEG_TO_RAD);

            // Offset x to align with galactic center
            x -= 8.;

            field.add(x, y, z, lum, 64, 45);
//cout << scientific << showpoint << setprecision(6) << l << " " << b << " " << r << endl << flush;
//cout << i << ": "<<  starTotal << endl;
            if( fstrm.eof() || i==starTotal-1 ) {
                fstrm.close();
                if( starTotal-1!=i ) {
                    cerr << "Number of stars indicated in wedge file does not match data." << endl;
                    exit(1);
                }
				break;
            }

        }


        if( removeDuplicates ) {
            delete [] lc;
            delete [] bc;
            this->starTotal -= skipTotal;
        }

        /// TODO /// Check to see if there is more data in the file (use a look ahead perhaps to avoid doubling the error checking)
/*      if( !fstrm.eof() ) {
            fstrm.close();
            cerr << "Number of stars indicated in wedge file does not match data." << endl;
            exit(1);
        }
*/

    }

    const WedgeInfo getWedgeInfo()
    {

        if( !wedgeInfo.initialized ) {

            /// TODO /// STUB
            wedgeInfo.lCenter = 0.;
            wedgeInfo.center = 0.;
            wedgeInfo.rCenter = 0.;
            wedgeInfo.maxDistance = 0.;

            wedgeInfo.crossVector = Vector3d(0., 0., 0.);
            wedgeInfo.depthVector = Vector3d(0., 0., 0.);

            wedgeInfo.initialized = true;

        }

        return (const WedgeInfo) wedgeInfo;

    }

};

inline HaloField* readStarFileXyz( const char* fileName, double lum = .5, bool binFlag = false )
{

    int starTotal;

    ifstream fstrm;

    if( binFlag ) {
        fstrm.open(fileName, ios::in|ios::binary);
        starTotal = fileGetIntBin(fstrm);
    }
    else {
        fstrm.open(fileName);
        starTotal = fileGetInt(fstrm);
    }

    HaloField* field = new HaloField(starTotal);

    // Get star positions
    float lineArg[3];

    for( int i = 0; i<starTotal; i++ ) {

        if( binFlag )
            fileGetFloatArrayBin(fstrm, 3, lineArg);
        else
            fileGetFloatArray(fstrm, 3, lineArg);

        float x = lineArg[0];
        float y = lineArg[1];
        float z = lineArg[2];

        field->add(x/*-26093.0901*/, y, z, lum, 64, 45);

        if( fstrm.eof() ) {

            if( starTotal!=i ) {
                fstrm.close();
                cerr << "Number of stars indicated in wedge file does not match data." << endl;
                exit(1);
            }

        }

    }

    fstrm.close();

    return field;

}

HaloField* readStarFileRaDecGru( const char* fileName, bool binFlag = false )
{

    int starTotal;

    ifstream fstrm;

    if( binFlag ) {
        fstrm.open(fileName, ios::in|ios::binary);
        starTotal = fileGetIntBin(fstrm);
    }
    else {
        fstrm.open(fileName);
        starTotal = fileGetInt(fstrm);
    }

    HaloField* field = new HaloField(starTotal);

    // Get star positions
    float lineArg[5];

    for( int i = 0; i<starTotal; i++ ) {

        if( binFlag )
            fileGetFloatArrayBin(fstrm, 5, lineArg);
        else
            fileGetFloatArray(fstrm, 5, lineArg);

        float Ra = lineArg[0];
        float Dec = lineArg[1];
        Uint8 g_mag = 255./lineArg[2];
        Uint8 r_mag = 255./lineArg[3];
        Uint8 u_mag = 255./lineArg[4];

        double x, y, z;
        raDecRadToCart( Ra, Dec, 1000., x, y, z );

        field->add(x-8., y, z, g_mag, r_mag, u_mag);

        if( fstrm.eof() ) {

            if( starTotal!=i ) {
                fstrm.close();
                cerr << "Number of stars indicated in wedge file does not match data." << endl;
                exit(1);
            }

        }

    }

    fstrm.close();

    return field;

}
#endif /* _DEMOFILE_HPP_ */
