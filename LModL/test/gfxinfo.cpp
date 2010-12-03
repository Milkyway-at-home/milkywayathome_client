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
#include <fstream>
#include <string>
#include <cmath>
#include <cstdlib>

#include "SDL.h"

using namespace std;


int main( int args, char **argv )
{

   // Set up SDL routines
   if( SDL_Init(SDL_INIT_VIDEO)==-1 ) {
      cerr << "SDL_Init: " << SDL_GetError() << endl;
      SDL_Quit();
      exit(1);
   }
   atexit(SDL_Quit);

   const SDL_VideoInfo* vInfo = SDL_GetVideoInfo();

   cout << "hw_available: " << vInfo->hw_available << endl;
    cout << "wm_available: " << vInfo->wm_available << endl;
    cout << "UnusedBits1: " << vInfo->UnusedBits1 << endl;
    cout << "UnusedBits2: " << vInfo->UnusedBits2 << endl;
    cout << "blit_hw: " << vInfo->blit_hw << endl;
    cout << "blit_hw_CC: " << vInfo->blit_hw_CC << endl;
    cout << "blit_hw_A: " << vInfo->blit_hw_A << endl;
    cout << "blit_sw: " << vInfo->blit_sw << endl;
    cout << "blit_sw_CC: " << vInfo->blit_sw_CC << endl;
    cout << "blit_sw_A: " << vInfo->blit_sw_A << endl;
    cout << "blit_fill: " << vInfo->blit_fill << endl;
    cout << "UnusedBits3: " << vInfo->UnusedBits3 << endl;
    cout << "video_mem: " << vInfo->video_mem << endl;
    cout << "current_w: " << vInfo->current_w << endl;
    cout << "current_h: " << vInfo->current_h << endl << endl;

   //   SDL_Palette *palette;
    cout << "BitsPerPixel: " << (int) vInfo->vfmt->BitsPerPixel << endl;
    cout << "BytesPerPixel: " << (int) vInfo->vfmt->BytesPerPixel << endl;
    cout << "Rloss: " << (int) vInfo->vfmt->Rloss << endl;
    cout << "Gloss: " << (int) vInfo->vfmt->Gloss << endl;
    cout << "Bloss: " << (int) vInfo->vfmt->Bloss << endl;
    cout << "Aloss: " << (int) vInfo->vfmt->Aloss << endl;
    cout << "Rshift: " << (int) vInfo->vfmt->Rshift << endl;
    cout << "Gshift: " << (int) vInfo->vfmt->Gshift << endl;
    cout << "Bshift: " << (int) vInfo->vfmt->Bshift << endl;
    cout << "Ashift: " << (int) vInfo->vfmt->Ashift << endl;
    cout << "Rmask: " << (int) vInfo->vfmt->Rmask << endl;
    cout << "Gmask: " << (int) vInfo->vfmt->Gmask << endl;
    cout << "Bmask: " << (int) vInfo->vfmt->Bmask << endl;
    cout << "Amask: " << (int) vInfo->vfmt->Amask << endl;

    return 0;

}
