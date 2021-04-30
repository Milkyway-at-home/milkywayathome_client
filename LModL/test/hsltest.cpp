#include <iostream>
#include <cstdlib>

#include "rgbconv.hpp"

using namespace std;

void test( Uint8 r, Uint8 g, Uint8 b, bool quiet )
{

    float x, y, l;
    float h, s;
    Uint8 ro = r;
    Uint8 go = g;
    Uint8 bo = b;
    rgbToXyl(r, g, b, x, y, l);
    if( !quiet )
        cout <<(int)r << ", "  << (int)g << ", "  << (int)b;
    xylToRgb(x, y, l, r, g, b);
    if( !quiet )
        cout << " == " <<  (int)r << ", "  << (int)g << ", "  << (int)b << endl;

    if( r!=ro || g!=go || b!=bo ) {
        cout << "Error: \n";
        cout <<(int)ro << ", "  << (int)go << ", "  << (int)bo << " != ";
        cout <<(int)r << ", "  << (int)g << ", "  << (int)b << endl;
        exit(0);
    }

    rgbToHsl(r, g, b, h, s, l);
    if( !quiet )
        cout << (int)r << ", "  << (int)g << ", "  << (int)b << " == ";
    hslToRgb(h, s, l, r, g, b);
    if( !quiet )
        cout << (int)r << ", "  << (int)g << ", "  << (int)b << "\n";

    if( r!=ro || g!=go || b!=bo ) {
        cout << "Error: \n";
        cout <<(int)ro << ", "  << (int)go << ", "  << (int)bo << " != ";
        cout <<(int)r << ", "  << (int)g << ", "  << (int)b << endl;
        exit(0);
    }


}

int main( int args, char **argv )
{

    if( args>3 ) {

        unsigned char r, g, b;

        r = atoi(argv[1]);
        g = atoi(argv[2]);
        b = atoi(argv[3]);

        test(r, g, b, false);

    }
    else {

        cout << "Testing . . .\n";
        for( int r = 0; r<255; r++ )
            for( int g = 0; g<255; g++ )
                for( int b = 0; b<255; b++ )
                    test(r, g, b, true);
        cout << "All checked\n";

    }
    return 0;

}
