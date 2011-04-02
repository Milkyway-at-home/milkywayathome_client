#include <fstream>
#include "binfile.hpp"
#include "SDL.h"


using namespace std;

int main( int args, char** argv )
{

    struct Star
    {
        Uint32 hipNumber;

        float x;
        float y;
        float z;

        Uint16 absMag ;

        Uint16 spectralClass;

    } __attribute__ ((packed));

    ifstream file;
    file.open("stars.dat", ios :: in | ios :: binary );
    if( !file ) {
        cerr << "Could not open stars.dat\n";
        exit(1);
    }

    // Remove header
    file.seekg(10);
    int totalStars = fileGetIntBin(file);
#if SDL_BYTEORDER == SDL_BIG_ENDIAN
    totalStars = SDL_Swap32(totalStars);
#endif
    Star* star = new Star[totalStars];

    for( int i = 0; i<totalStars; i++ ) {

        file.read((char*) &(star[i]), sizeof(Star));
#if SDL_BYTEORDER == SDL_BIG_ENDIAN
        star.hipNumber = SDL_Swap32(star.hipNumber);
        star.x = SDL_Swap32(star.x);
        star.y = SDL_Swap32(star.y);
        star.z = SDL_Swap32(star.z);
        star.absMag = SDL_Swap16(star.absMag);
        star.spectralClass = SDL_Swap16(star.spectralClass);
#endif
    }

    cout << showpos << showpoint << scientific << totalStars << endl;

    for( int i = 0; i<totalStars; i++ )
        cout << star[i].x << " " << star[i].y << " " << star[i].z << endl;

    return 0;

}