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

#include "binfile.hpp"


string intToString( Sint64 i, int base )
{

    string s, sign;

    if( i<0 ) {
        i *= -1;
        sign = "-";
    }
    do {

        int a = i%base;
        if( a>9 )
            s = char('A'+(a-10)) + s;
        else
            s = char('0'+a) + s;
        i /= base;

    } while (i);

    return sign+s;

}

int isBigEndian()
{
    int endianTest = 1;
    if ( *( (char*) &endianTest ) == 1 )
        return 0;
    else
        return 1;
}

string fileGetString( ifstream &fstrm )
{

    string line;
    do {

        if( fstrm.eof() ) {
            cout << "Error reading file.\n";
            exit(1);
        }
        getline(fstrm, line);

    } while( line=="" );

    return line;

}

int fileGetInt( ifstream &fstrm )
{

    string line = fileGetString(fstrm);
    return atoi(line.c_str());

}

void fileGetIntArray( ifstream &fstrm, int total, int *array )
{

    if( total<1 )
        return;

    string line = fileGetString(fstrm);

    const char *str = line.c_str();
    char *nextPtr = NULL;

    array[0] = strtol(str, &nextPtr, 10);

    for( int i = 1; i<total; i++ ) {

        if( *nextPtr=='\0' ) {
            cout << "Error reading file.\n";
            exit(1);
        }
        array[i] = strtol((const char*)nextPtr, &nextPtr, 10);

    }

}

double fileGetFloat( ifstream &fstrm )
{

    string line = fileGetString(fstrm);
    return atof(line.c_str());

}

double fileGetDouble( ifstream &fstrm )
{

    string line = fileGetString(fstrm);
    return atof(line.c_str());

}

void fileGetFloatArray( ifstream &fstrm, int total, float *array )
{

    if( total<1 )
        return;

    string line = fileGetString(fstrm);

    const char *str = line.c_str();
    char *nextPtr = NULL;

    array[0] = strtod(str, &nextPtr);

    for( int i = 1; i<total; i++ ) {

        if( *nextPtr=='\0' ) {
            cout << "Error reading file.\n";
            exit(1);
        }
        array[i] = strtod((const char*)nextPtr, &nextPtr);

    }

}

void fileGetDoubleArray( ifstream &fstrm, int total, double *array )
{

    if( total<1 )
        return;

    string line = fileGetString(fstrm);

    const char *str = line.c_str();
    char *nextPtr = NULL;

    array[0] = strtod(str, &nextPtr);

    for( int i = 1; i<total; i++ ) {

        if( *nextPtr=='\0' ) {
            cout << "Error reading file.\n";
            exit(1);
        }
        array[i] = strtod((const char*)nextPtr, &nextPtr);

    }

}

void _checkFileErrorCpp()
{
    if( ios_base::goodbit ) {
        cout << "Error reading file.\n";
        exit(1);
    }

}

Sint32 fileGetIntBin( ifstream &fstrm )
{

    int i;//, lose; - only problem in 19gb file
    fstrm.read((char*) &i, 4);
//    fstrm.read((char*) &lose, 4);
    _checkFileErrorCpp();
    i = SWAP32(i);
    return i;

}

double fileGetDoubleBin( ifstream &fstrm )
{

    double d;
    fstrm.read((char*) &d, 8);
    _checkFileErrorCpp();
    d = SWAP32(d);
    return d;

}

void fileGetDoubleArrayBin( ifstream &fstrm, int total, double *array )
{

#if SDL_BYTEORDER == SDL_BIG_ENDIAN
    for( int i = 0; i<total; i++ )
        array[i] = fileGetDoubleBin(fstrm);
#else
    fstrm.read((char*) array, total<<3);
    _checkFileErrorCpp();
#endif

}

float fileGetFloatBin( ifstream &fstrm )
{

    double d;
    fstrm.read((char*) &d, 4);
    _checkFileErrorCpp();
    d = SWAP16(d);
    return d;

}

void fileGetFloatArrayBin( ifstream &fstrm, int total, float *array )
{

#if SDL_BYTEORDER == SDL_BIG_ENDIAN
    for( int i = 0; i<total; i++ )
        array[i] = fileGetFloatBin(fstrm);
#else
    fstrm.read((char*) array, total<<2);
    _checkFileErrorCpp();
#endif

}

