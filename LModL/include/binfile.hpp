/*****************************************************************************
 *                                                                           *
 *  Copyright (C) 2010 Shane Reilly, Heidi Newberg, Malik Magdon-Ismail,     *
 *  Carlos Varela, Boleslaw Szymanski, and Rensselaer Polytechnic Institute  *
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

#ifndef _BINFILE_HPP_
#define _BINFILE_HPP_

#include <iostream>
#include <fstream>
#include <string>

#include "SDL.h"

using namespace std;


#if SDL_BYTEORDER == SDL_LIL_ENDIAN
#define SWAP16(X)    (X)
#define SWAP32(X)    (X)
#else
#define SWAP16(X)    SDL_Swap16(X)
#define SWAP32(X)    SDL_Swap32(X)
#endif


string intToString( Sint64 i, int base = 10 );

int isBigEndian();

string fileGetString( ifstream &fstrm );

int fileGetInt( ifstream &fstrm );

void fileGetIntArray( ifstream &fstrm, int total, int *array );

double fileGetFloat( ifstream &fstrm );

double fileGetDouble( ifstream &fstrm );

void fileGetFloatArray( ifstream &fstrm, int total, float *array );

void fileGetDoubleArray( ifstream &fstrm, int total, double *array );

void _checkFileErrorCpp();

Sint32 fileGetIntBin( ifstream &fstrm );

double fileGetDoubleBin( ifstream &fstrm );

void fileGetDoubleArrayBin( ifstream &fstrm, int total, double *array );

float fileGetFloatBin( ifstream &fstrm );

void fileGetFloatArrayBin( ifstream &fstrm, int total, float *array );


#endif /* _BINFILE_HPP_ */
