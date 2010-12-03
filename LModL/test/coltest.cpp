#include <iostream>
#include <cstdlib>

#include "draw8.hpp"
#include "draw32.hpp"

using namespace std;

int main( int args, char **argv )
{

    initDraw();

    SDL_Surface* display = setVideo();
    setPaletteTable();

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
    
        for( int y = 0; y<512; y++ )
            for( int x = 0; x<512; x++ ) {
                int h = x>>1;
                int c = t2>>1;
                int l = y>>1;
//float r=0, g=0, b=0;
//HSL_to_RGB(h/255., c/255., l/255., r, g, b);
//putPixelSumClip32(display, x, y, (int(r*255.)<<16)|(int(g*255.)<<8)|(int(b*255.)));
//putPixelSumClip32(display, x, y, cieLchToColor32(l, c, h));
                putPixelSumClip32(display, x, y, getLightnessColor32(c, h)[l]);
//if(rand()%100==1)
  //  cout << c << ", " << h << ", " << l << ": " << hex <<(int)getLumaColor32(c, h)[l] << endl;
            }
        SDL_Flip(display);
clearSurface(display);
timer.wait();
        if( SDL_PollEvent(&event) && event.type==SDL_KEYDOWN )
	    return 0;
    }


}

