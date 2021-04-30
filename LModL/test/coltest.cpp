#include <iostream>
#include <cstdlib>

#include "draw8.hpp"
#include "draw32.hpp"
#include "draw32p.hpp"

using namespace std;

int main( int args, char **argv )
{
/*
cout << hex << b32pTob32(0xffdffbfe) << endl;
cout << b32Tob32p(b32pTob32(0xffdffbfe)) << endl;
cout << b32pTob32(b32Tob32p(0xffdffbfe)) << endl;
*/
    initDraw();

    SDL_Surface* display = setVideo();

    if( display==NULL ) {
        cerr << "Error setting display\n";
        exit(1);
    }

    SDL_Event event;
    StepTimer timer(.01);

    for( int t = 0;  ;t++ ) {
/*
        while(true)
            if( SDL_PollEvent(&event) && event.type==SDL_KEYDOWN && keyTimer.getSignal() )
                break;
*/

int t2;
if( t==1022 )
    t = 0;
if( t>511 )
    t2 = 1022-t;
else
    t2 = t;

clearSurface(display);
        for( int y = 0; y<512; y++ )
            for( int x = 0; x<512; x++ ) {
                int h = x>>1;
                int c = t2>>1;
                int l = y<<1;
//float r=0, g=0, b=0;
//HSL_to_RGB(h/255., c/255., l/255., r, g, b);
//putPixelSumClip32(display, x, y, (int(r*255.)<<16)|(int(g*255.)<<8)|(int(b*255.)));
//putPixelSumClip32(display, x, y, cieLchToColor32(l, c, h));
Uint32 pixel = getLightnessColor32p(c, h)[l];
                putPixelSumClip32p(display, x, y, pixel);
//cout << hex << b32pTob32(pixel) << endl;
if( pixel!=b32pTob32(b32Tob32p(pixel)) )
    cout << hex << pixel << " !=" << b32pTob32(b32Tob32p(pixel)) << endl;
//if(rand()%100==1)
  //  cout << c << ", " << h << ", " << l << ": " << hex <<(int)getLumaColor32(c, h)[l] << endl;
            }


//SDL_BlitSurface(sbuf, NULL, display, NULL);

        SDL_Flip(display);
timer.wait();
        if( SDL_PollEvent(&event) && event.type==SDL_KEYDOWN )
            break;
    }

   return 0;

}

