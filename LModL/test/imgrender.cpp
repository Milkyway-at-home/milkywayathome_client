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

#include "drawhalo.hpp"
#include "imgplot.hpp"

using namespace std;



int main( int args, char **argv )
{

    if( args<2 ) {
        cout << "Usage: ./imgrender bmp_filename [sample_size] [blur_diameter]\n";
        exit(1);
    }

    float sample = 100000;
    float objectWidth = 30.;
    float objectThickness = .3;
    float blurDiameter = 20;
    float fps = 60.;
    int bpp = 32;

    string fileName = argv[1];

    // Handle optional arguments
    if( args>2 )
        sample = atof(argv[2]);
    if( args>3 )
        blurDiameter = atof(argv[3]);

    ImagePlot imagePlot(fileName, sample, objectWidth, objectThickness);

    FieldAnimation animate(bpp, fps);
    animate.add(imagePlot.getField(), blurDiameter);

    while( true ) {

        if( animate.pollEvent() ) {
            ;
        }

    }

    return 0;

}
