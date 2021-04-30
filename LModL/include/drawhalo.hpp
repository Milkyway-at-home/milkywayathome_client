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

#ifndef _DRAWHALO_HPP_
#define _DRAWHALO_HPP_

#include <cassert>

#include "draw.hpp"
#include "drawcore.hpp"
#include "3dproj.hpp"

using namespace std;


#define BLUR_GRANULARITY 0


extern const int PRINT_XSIZE;
extern const int PRINT_YSIZE;


class BlurCircle
{

    struct
    {
        int granShift;
        int granSqShift;
        int granularity;
        float granNorm;
        float granNormAdd;
        float cMult;
    } blurStats;

public:

    BlurCircle( int granShift = 0 );
        // Initialize parameters
        // 2^'granShift' represents the number of times a pixel value is sampled

    void draw( SDL_Surface *surface, float x, float y, float diameter, float intensity );
        // Surface must be locked before calling this function
        // 'x' and 'y' are the pixel coordinates of the draw
        // diameter is the diameter of the blur circle representing 3.0 standard deviations of a Gaussian distribution
        // 'hue' is the hue value of the incoming light on a scale from 0. to 1.
        // 'sat' is the saturation value of the incoming light on a scale from 0. to 1.
        // Intensity is the maximum amount of light in the blur on a scale from 0. to 1.
        // Technically the intesity should be divided by pi*radius^2 since the incoming light is distributed across blur, but
        //   this calculation is left to the calling application in order to make intensity conversions across different
        //   screen resolutions easy

};

class HaloType
{

    SDL_Surface ****pointOffset;
    float radius;
    float lumMult, lumDiv;
    float maxLum;

    int lumGranShift;
    int lumGranularity;
    int haloGranShift;
    int haloGranularity;
    float haloGranfloat;

public:

    HaloType( float diameter = 7., float maxLuminosity = 1., int lumGranShift = 7, int haloGranShift = 3 );

    ~HaloType();

    void set( float diameter = 7., float maxLuminosity = 1. );

    float getDiameter();

    float getLum();

    void draw( SDL_Surface *surface, float x, float y, float luminosity, Uint32* palette = GRAY_PALETTE );

//   void draw( SDL_Surface* surface, fix32 x, fix32 y );

   void _drawTest() const;

};

struct HaloPoint
{
    Vector3d position;
    float lightness;
    Uint32* palette;
};

class HaloField
{

private:

    int arrayTotal;
    int stackPtr, stackEndPtr;
    float lightAdd;
    float axesLimit;

    HaloPoint **field;

public:

    HaloField( int pointTotal );

    ~HaloField();

    void setLum( float lightness );

    float getLum() const;

    Vector3d getCenter() const;

    void clearField();

    void add( float x, float y, float z, float L, int C, int h );
        // Lightness (L) is a value from 0. to 1. representing the blur center-pixel L value
        // Higher L values will increase overall lightness of the blur, but all intensities are capped at 1.
        // Saturation (C) is an integer from 0 to 255
        // Hue (h) is an integer from 0 to 255

    void add( float x, float y, float z, Uint8 G, Uint8 R, Uint8 B );
        // R, G, B are integers from 0 to 255

    HaloPoint** getPoints() { return field; }

    void getNext( float x, float y, float z, float L, int C, int h );

    void set( int index, float x, float y, float z, float l, float c, float h );

    void drawAxes( SDL_Surface* surface, const Camera& cv, HaloType& lineBlur );

    void drawCamera( SDL_Surface *surface, Camera *cv, HaloType& lineBlur, HaloType& endBlur );

    void draw( SDL_Surface* surface, Camera *cv, HaloType &haloType, int skip = 1 );

};

class FieldAnimation
{

    float fps;
    bool cameraDisplay, axesDisplay;
    StepTimer *fpsTimer, *keyTimer, *fastKeyTimer;
    HaloType *lineBlur, *endBlur;
    float sizeAdjust;

public:

    SDL_Surface *display;
    Camera* cv;
    SDL_Event* event;

    struct Field
    {
        HaloType *blur;
        HaloField* data;
        float LIGHTNESSdjust;
        int skipDraw;
    };

    Field field[100];
    int totalFields;
    int fieldPtr;
    int fieldModIndex;

    FieldAnimation( int bpp = 32, float fps = 60., bool fullScreen = false, string caption = "Application", string iconFileName = "icon.bmp" );

    void resetFieldPointer();

    void add( HaloField *newField, float blurDiameter = 7. );

    void showCamera();

    void hideCamera();

    void showAxes();

    void hideAxes();

    void applyFrame( SDL_Surface* surface );

    bool pollDemo();

    bool pollEvent();

};

bool clip2d( float &x, float &y, float xc, float yc, float xMin, float xMax, float yMin, float yMax );

void drawLine( SDL_Surface *surface, HaloType *bc, float lum, float xs, float ys, float xe, float ye, Uint32* palette = GRAY_PALETTE );

void drawLine( SDL_Surface *surface, HaloType *bc, const Camera* cam, float lum, const Vector3d start, const Vector3d end, Uint32* palette = GRAY_PALETTE );



#endif /* _DRAWHALO_HPP_ */
