/****************************************************************************
 *                                                                          *
 *  Copyright (C) 2010 Shane Reilly and Rensselaer Polytechnic Institute.   *
 *                                                                          *
 *  This file is part of the Light Modeling Library (LModL).                *
 *                                                                          *
 *  This library is free software: you can redistribute it and/or modify    *
 *  it under the terms of the GNU General Public License as published by    *
 *  the Free Software Foundation, either version 3 of the License, or       *
 *  (at your option) any later version.                                     *
 *                                                                          *
 *  This library is distributed in the hope that it will be useful,         *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *  GNU General Public License for more details.                            *
 *                                                                          *
 *  You should have received a copy of the GNU General Public License       *
 *  along with this library. If not, see <http://www.gnu.org/licenses/>.    *
 *                                                                          *
 *  Shane Reilly                                                            *
 *  reills2@cs.rpi.edu                                                      *
 *                                                                          *
 ****************************************************************************/

#include <iostream>
#include <string>
#include <fstream>
#include <ios>

#define TEST_MODE

#include "sdlwrap.hpp"
#include "binfile.hpp"
#include "draw8.hpp"


using namespace std;


const int PRIMITIVE_SIZE = 100;  // 100x100 pixel primitives

//   04 01 ..    01 02 ..
//   02 ** .. -> 08 ** ..
//   .. .. ..    .. .. ..
//
//   .. 01 04    .. 02 04
//   .. ** 02 -> .. ** 10
//   .. .. ..    .. .. ..
//
//   .. .. ..    .. .. ..
//   .. ** 02 -> .. ** 10
//   .. 01 04    .. 40 80
//
//   .. .. ..    .. .. ..
//   02 ** .. -> 08 ** ..
//   04 01 ..    20 40 ..

const int mapConvert[4][3] = {
    {0x02, 0x08, 0x01},
    {0x02, 0x10, 0x04},
    {0x40, 0x10, 0x80},
    {0x40, 0x08, 0x20}
};

class fontPrimitive
{

    SDL_Surface* prim[4][256];

    int translate( int cornerIndex, int masterIndex, int &matchIndex )
        // Returns number of sides matched
        // Stores translation in 'masterIndex'
    {

        int corner = cornerIndex>>3;
        int ci = cornerIndex&0x7;

        int n = 0;
        int mc;
        matchIndex = 0;

        for( int i = 0, s = 1; i<3; i++, s<<=1 )

            if( ci&s ) {
                mc = mapConvert[corner][i];
cout << "MC: " << mc << " < " << corner << " < " << i << endl;

                if( mc&~masterIndex ) {
                    matchIndex = 0;
cout << "reject\n";
                    return 0;  // Match exceeds limitations
                }

cout << "MC: " << mc << " < " << corner << " < " << i << endl;
                matchIndex |= mc;
                n++;
            }

        return n;

    }

    int findBestSolution( int masterIndex, int &match )
        // Finds first best match for glyph
        // Match may not complete a glyph and additional calls can be made to do so
    {

        int bestN = 0;
        int bestIndex = 0;
        int matchIndex = 0;
        match = 0;
        for( int i = 0; i<32; i++ ) {
cout << i << ": " << flush;
            int n = translate(i, masterIndex, matchIndex);
cout << n << endl;
            if( n>bestN ) {
                bestN = n;
                bestIndex = i;
                match = matchIndex;
cout << "M: " << match << endl;
cout << "return I: " << bestIndex << endl;
            }
        }
        return bestIndex;

    }

    void mapCorners( SDL_Surface* corner[4][32] )
    {

        for( int type = 0; type<4; type++ )
            blitSurfaceClipSum7(corner[type][0], prim[type][0], 0, 0);

        for( int i = 1; i<256; i++ ) {
cout << i << endl ; cout.flush();

            int masterIndex = i;
            while( masterIndex ) {
                int bestMatch = 0;
                int cornerIndex = findBestSolution(masterIndex, bestMatch);

                for( int type = 0; type<4; type++ )
                    blitSurfaceClipSum7(corner[type][cornerIndex], prim[type][i], 0, 0);

cout << "best match: " << bestMatch << endl;
                masterIndex &= ~bestMatch;  // Reduce problem and repeat search
cout << "master index: " << masterIndex << endl;
            }

        }
SDL_Surface* dis = SDL_GetVideoSurface();
clearSurface(dis);
for( int t = 0; t<4; t++ )
for( int i = 0; i<32; i++ )
resizeSumSurfaceClip7(corner[t][i], dis, (i%32)*22*5/2, (t*9+(i/32))*22*5/2, 2, 2);
//for( int t = 0; t<4; t++ )
//for( int i = 0; i<256; i++ )
//resizeSumSurfaceClip7(getPrimitive(t, i), dis, (i%32)*22*5/2, (t*9+(i/32))*22*5/2, 2, 2);
SDL_Flip(dis);
    }

public:

    fontPrimitive( string fileName )
    {

        // Create Primitives
        for( int t = 0; t<4; t++ )
            for( int i = 0; i<256; i++ ) {
                prim[t][i] = newSurface8(PRIMITIVE_SIZE, PRIMITIVE_SIZE);
                /// TODO /// See if this line is needed
                clearSurface(prim[t][i]);
            }

        // Create corner glyphs
        SDL_Surface* corner[4][32];
        for( int t = 0; t<4; t++ )
            for( int i = 0; i<32; i++ ) {
                corner[t][i] = newSurface8(PRIMITIVE_SIZE, PRIMITIVE_SIZE);
                /// TODO /// See if this line is needed
                clearSurface(prim[t][i]);
            }

        // Read glyph file
        SDL_Surface* surface = SDL_LoadBMP(fileName.c_str());
        if( surface==NULL ) {
            cerr << "Unable to load " << fileName << ": " << SDL_GetError();
            exit(1);
        }

        // Convert to 7bpp indexed
        for( int y = 0; y<surface->h; y++ )
            for( int x = 0; x<surface->w; x++ )
                putPixel8(surface, x, y, 0x7f-(getPixel8(surface, x, y)>>1));

        // Read up-left glyphs
        int i = 0;
        for( int t = 0, y = 0; t<4; t++, y += PRIMITIVE_SIZE )
            for( int x = 0; i<8; i++, x += PRIMITIVE_SIZE )
                copySurfaceClip8(surface, corner[t][i], x, y, 0, 0, PRIMITIVE_SIZE, PRIMITIVE_SIZE);


SDL_Surface* dis = SDL_GetVideoSurface();
clearSurface(dis);
for( int t = 0; t<4; t++ )
for( int i = 0; i<32; i++ )
resizeSumSurfaceClip7(corner[t][i], dis, (i%32)*22*5/2, (t*9+(i/32))*22*5/2, 2, 2);
//for( int t = 0; t<4; t++ )
//for( int i = 0; i<256; i++ )
//resizeSumSurfaceClip7(getPrimitive(t, i), dis, (i%32)*22*5/2, (t*9+(i/32))*22*5/2, 2, 2);
SDL_Flip(dis);
return;

        // Flip corners horizontally and apply copy
        for( int t = 0; t<4; t++ )
            for( ; i<16; i++ ) {
                copySurfaceClip8(corner[t][i-8], corner[t][i], 0, 0, 0, 0, PRIMITIVE_SIZE, PRIMITIVE_SIZE);
                flipHorizontal8(corner[t][i]);
            }

        // Flip corners vertically and apply copy
        for( int t = 0; t<4; t++ )
            for( ; i<24; i++ ) {
                copySurfaceClip8(corner[t][i-8], corner[t][i], 0, 0, 0, 0, PRIMITIVE_SIZE, PRIMITIVE_SIZE);
                flipVertical8(corner[t][i]);
            }

        // Flip corners horizontally and apply copy
        for( int t = 0; t<4; t++ )
            for( ; i<32; i++ ) {
                copySurfaceClip8(corner[t][i-8], corner[t][i], 0, 0, 0, 0, PRIMITIVE_SIZE, PRIMITIVE_SIZE);
                flipHorizontal8(corner[t][i]);
            }

        // Apply glyphs to master copy
        mapCorners(corner);

        SDL_FreeSurface(surface);
        for( int t = 0; t<4; t++ )
            for( int i = 0; i<32; i++ )
                SDL_FreeSurface(corner[t][i]);


    }

    inline SDL_Surface* getPrimitive( int type, int index )
    {
        return prim[type][index];
    }

};

class fontPrimitiveMap
{

    SDL_Surface* cSurface;

    char* cMap[0x80];
    int mapXSize, mapYSize;

    int getMapIndex( char* cMap, int xSize, int x, int y )
    {
        int mapIndex = 0;
        if( cMap[(y-1)*xSize+x-1]!='.' )
            mapIndex |= 0x01;
        if( cMap[(y-1)*xSize+x]!='.' )
            mapIndex |= 0x02;
        if( cMap[(y-1)*xSize+x+1]!='.' )
            mapIndex |= 0x04;
        if( cMap[y*xSize+x-1]!='.' )
            mapIndex |= 0x08;
        if( cMap[y*xSize+x+1]!='.' )
            mapIndex |= 0x10;
        if( cMap[(y+1)*xSize+x-1]!='.' )
            mapIndex |= 0x20;
        if( cMap[(y+1)*xSize+x]!='.' )
            mapIndex |= 0x40;
        if( cMap[(y+1)*xSize+x+1]!='.' )
            mapIndex |= 0x80;
        return mapIndex;
    }

    int getMapType( char c )
    {

        switch( c )
        {

        case '.':
            return 0;

        case 'o':
            return 1;

        case 'e':
            return 2;

        case 'x':
            return 3;

        default:
            cerr << "Invalid character in font file: " << c << endl;
            exit(1);

        }
    }

public:

    fontPrimitiveMap( string fileName )
    {

        /// TODO /// add more error checking

        // Open font primitive map
        ifstream fs;
        fs.open(fileName.c_str());
        int a[2];
        fileGetIntArray(fs, 2, a);
        mapXSize = a[0]+2;
        mapYSize = a[1]+2;
        cSurface = newSurface8((mapXSize-2)*PRIMITIVE_SIZE, (mapYSize-2)*PRIMITIVE_SIZE);

        // Set all characters to 'spaces' by default
        for( int i = 0; i<0x80; i++ ) {
            cMap[i] = new char[mapXSize*mapYSize];
            for( int y = 0; y<mapYSize; y++ )
                for( int x = 0; x<mapXSize; x++ )
                    cMap[i][y*mapXSize+x] = '.';
        }

        // Load characters
        for( int i = 0; i<0x80; i++ ) {

            for( int y = 1; y<mapYSize-1; y++ ) {
                string s = fileGetString(fs);
                for( int x = 1; x<mapXSize-1; x++ )
                    cMap[i][y*mapXSize+x] = s[x-1];
            }

        }

    }

    ~fontPrimitiveMap()
    {
        /// TODO /// STUB
    }

    void putChar( SDL_Surface* surface, double x, double y, char c, fontPrimitive &prim )
    {

        int cIndex = c&0x7f;

        int ix = (int) x;
        int iy = (int) y;

        /// TODO /// use fx, fy

        // Apply character map
        clearSurface(cSurface);
        for( int ys = 1; ys<mapYSize-1; ys++ )
            for( int xs = 1; xs<mapXSize-1; xs++ ) {
                int mt = getMapType(cMap[cIndex][ys*mapXSize+xs]);
                SDL_Surface* aPrim;
                aPrim = prim.getPrimitive(mt, getMapIndex(cMap[cIndex], mapXSize, xs, ys));
                blitSurfaceClipSum7(aPrim, cSurface, (xs-1)*PRIMITIVE_SIZE, (ys-1)*PRIMITIVE_SIZE);
            }

        resizeSumSurfaceClip7(cSurface, surface, ix, iy, 50, 50);

    }

};

int main( int args, char **argv )
{

    // Handle arguments
    if( args<4 ){
        cout << "Usage: ./cfont primitive_file fpm_file [font_size]\n";
        return 0;
    }
    string primFile = argv[1];
    string fpmFile = argv[2];
    int size = atoi(argv[3]);

    // Set up environment
    SDL_Surface* display = setVideo8();
    SDL_Event event;

    // Read in primitives
    fontPrimitive pf(primFile);
    fontPrimitiveMap pMap(fpmFile);

    StepTimer keyTimer(1./60.);

    while( true ) {

        clearSurface(display);

        if( SDL_PollEvent(&event) )

            if( event.type==SDL_KEYDOWN )
            {

                if( event.key.keysym.sym==SDLK_ESCAPE )
                    return 0;

                else if( event.key.keysym.sym==SDLK_SPACE && keyTimer.getSignal()  ) {
                    keyTimer.set();
                    for( int i = 32; i< 0x80; i++ )
                        pMap.putChar(display, (double)(i-32)*12., 0., i, pf);

                    SDL_Flip(display);
                }

                else if( event.key.keysym.sym==SDLK_BACKSPACE && keyTimer.getSignal()  ) {
                    keyTimer.set();
SDL_Surface* dis = SDL_GetVideoSurface();
clearSurface(dis);
for( int t = 0; t<4; t++ )
for( int i = 0; i<256; i++ )
resizeSumSurfaceClip7(pf.getPrimitive(t, i), dis, (i%32)*22, (t*9+(i/32))*22, 5, 5);
SDL_Flip(dis);
                }

            }

    }

    return 0;

}
