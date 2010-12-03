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

#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cstdlib>
#include <cmath>

using namespace std;


int main( int args, char **argv )
{

   if( args!=2 ){
      printf("Usage: trigtblm DIGITS_ACCURACY\n");
      return 1;
   }

   int digits = atoi(argv[1]);
   int gran = (int) pow(10., digits);
   int size = (int) (gran*1.25+1.);
   int offset = (int) (gran*.25);
   float dgran = gran/6.2831853071795864;
   float rdgran = 6.2831853071795864/gran;

   cout << setiosflags(ios::scientific) << setprecision(16);

   cout << "#define TRIGTBL_SIZE        " << size << endl;
   cout << "#define TRIGTBL_COS_OFFSET  " << offset << endl;
   cout << "#define TRIGTBL_GRAN        " << dgran << endl << endl;

   cout << "const float TRIG_TABLE[TRIGTBL_SIZE] = \n";
   cout << "{\n";

   for( int i = 0; i<size; i++ )
      cout << "    " << sin(float(i)*rdgran) << ",\n";

   cout << "};\n";

}
