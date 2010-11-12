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

#include "demofile.hpp"
#include "drawhalo.hpp"

using namespace std;


int main( int args, char **argv )
{
/*
    if( args<2 ){
        printf("Usage: ./mwdemo wedge_file_name\n");
        return 1;
    }
*/

    string fileName;
    if( args>1 )
        fileName = argv[1];
    else {
        cout << "Usage: rmvdup [lba text file] [output file]\n";
        exit(0);
    }

    bool earthCenter = true;

    // Read in wedge
    WedgeFile wf;

    int totalStars = wf.getStarTotal(fileName);

    HaloField wedge(totalStars);

    wf.readStars(fileName, wedge, .5, true);

    totalStars = wf.getStarTotal();
    HaloPoint **star = wedge.getPoints();
    
    cout << totalStars << endl;
    for( int i = 0; i<totalStars; i++ ) {
        if( earthCenter )
            cout << star[i]->position.x+8 << " ";
        else
            cout << star[i]->position.x << " ";
        cout << star[i]->position.y << " ";
        cout << star[i]->position.z << endl;
    }
    
    return 0;
    
}
