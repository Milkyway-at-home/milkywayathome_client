/*****************************************************************************
 *                                                                           *
 *  Copyright (C) 2010 Shane Reilly and Rensselaer Polytechnic Institute     *
 *                                                                           *
 *  This file is part of the Light Modeling Library (LModL).                 *
 *                                                                           *
 *  This library is free software: you can redistribute it and/or modify     *
 *  it under the terms of the GNU General Public License as published by     *
 *  the Free Software Foundation, either version 3 of the License, or        *
 *  (at your option) any later version.                                      *
 *                                                                           *
 *  This library is distributed in the hope that it will be useful,          *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of           *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the             *
 *  GNU General Public License for more details.                             *
 *                                                                           *
 *  You should have received a copy of the GNU General Public License        *
 *  along with this library. If not, see <http://www.gnu.org/licenses/>.     *
 *                                                                           *
 *  Shane Reilly                                                             *
 *  reills2@cs.rpi.edu                                                       *
 *                                                                           *
 *****************************************************************************/

#include <cstdlib>

#define TEST_MODE

#include "drawhalo.hpp"

using namespace std;

class CubeField
{

    HaloField* field;
    float stepOffset;
    float lum;

public:

    CubeField( float lum = 1. )
    {

        this->lum = lum;
        int pointTotal = 9*9*9;
        field = new HaloField(pointTotal);
        stepOffset = 0;
        step();

    }

    ~CubeField()
    {
        delete field;
    }

    void step()
    {
        /// TODO /// Make these more interesting
        int c = 128;
        int h = 160;

        stepOffset += .1;
        if( stepOffset>24. )
            stepOffset -= 24.;
        float a = 1000.;
        float b = a*4.;
        field->clearField();
        for( int k = 0; k<9; k++ )
            for( int j = 0; j<9; j++ )
                for( int i = 0; i<9; i++ ) {
                    float l = (((k+j+i)*100+(int)(stepOffset*100.))%2400)/1200.;
                    l = l<1. ? l : 2.-l;
                    field->add(a*k-b, a*j-b, a*i-b, lum*l, c, h);
                }
    }

    const HaloField* getField()
    {
        return field;
    }

};

int main( int args, char **argv )
{

    initDraw();

    float lum = .5;
    float diameter = 60.;
    float fps = 60.;
    int bpp = 32;

    // Handle optional arguments
    if( args>1 )
        bpp = atoi(argv[1]);
    if( args>2 )
        lum = atof(argv[2]);
    if( args>3 )
        diameter = atof(argv[3]);

    cout << "Setting up field . . .\n" << flush;
    CubeField cube(lum);
    HaloField* cubeField = (HaloField*)cube.getField();

    // Turn on demo display
    cout << "Setting up window . . .\n" << flush;
    FieldAnimation cubeDemo(bpp, fps, false);
    cout << "Adding stars . . .\n" << flush;
    cubeDemo.add(cubeField, diameter);
    cout << "Setting parameters . . .\n" << flush;
    cubeDemo.showCamera();

    cout << "Entering event cycle . . .\n" << flush;

    while( true )
    {

        // Poll for demo event and update video display as needed
        if( cubeDemo.pollEvent() )
	{

            // Other actions to be performed, such as data modification go here
            cube.step();

        }

    }

    cout << "Done." << endl;

    return 0;

}
