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

#include <cstdio>
#include <cstdlib>
#include "quatmath.hpp"

using namespace std;


int main( int args, char **argv )
{

    // Handle arguments
    if( args<4 ){
        printf("Usage: ./quattest x y z angle\n");
        printf("\n");
        printf("    Diplays the quaternion representing the given vector and rotation angle\n");
        return 0;
    }

    float x = atof(argv[1]);
    float y = atof(argv[2]);
    float z = atof(argv[3]);
    float a = atof(argv[4]);

    Quaternion q(Vector3d(x, y, z), a);

    cout << "x = " << q.x << endl;
    cout << "y = " << q.y << endl;
    cout << "z = " << q.z << endl << endl;
    cout << "w = " << q.w << endl;

    return 0;

}