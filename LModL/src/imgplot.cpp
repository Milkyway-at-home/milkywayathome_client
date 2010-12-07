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

#include "imgplot.hpp"


ImagePlot::ImagePlot( string fileName, int pointTotal, float width, float thickness, float lum )
{
    field = NULL;
    readImage(fileName, pointTotal, width, thickness, lum);
}

ImagePlot::~ImagePlot()
{
    delete field;
}

void ImagePlot::readImage( string fileName, int pointTotal, float width, float thickness, float lum )
{

    delete field;
    this->pointTotal = pointTotal;
    this->lum = lum;
    field = new HaloField(pointTotal);

    SDL_Surface* image = SDL_LoadBMP(fileName.c_str());
    if( image==NULL ) {
        cout << "Unable to load file\n";
        exit(1);
    }

    float unitsPerPixel = width/((image->w+image->h)/2.);

    int maxPixel = 0;
    for( int y = 0; y<image->h; y++ )
        for( int x = 0; x<image->w; x++ )
            if( colorToL(image, getPixelClip32(image, x, y))>maxPixel )
                maxPixel += colorToL(image, getPixelClip32(image, x, y))*colorToL(image, getPixelClip32(image, x, y));
    float maxPixelDbl = sqrt((float) maxPixel);

    for( int i = 0; i<pointTotal; ) {

        int ix = rand32()%(image->w);
        int iy = rand32()%(image->h);
        int threshold = rand32()%16;
//cout << "try" << ix << ", " << iy<< "\n" << flush;
        int pixel = getPixelClip32(image, ix, iy);
//cout << "succeed\n" << flush;

        float dx = ix - image->w/2.;
        float dy = image->h/2. - iy;  // Flip orientation for graph consistency

        float x = unitsPerPixel*dx + randFloat()*unitsPerPixel;
        float y = unitsPerPixel*dy + randFloat()*unitsPerPixel;
        float z = float(colorToL(image, pixel))/maxPixelDbl*width*
                      (acos(randFloat()*2.-1.)/TRIG_PI-.5)*unitsPerPixel;
        float l = .1;

        Uint8 r, g, b;
        double h, s, l2;
        colorToRgb(image, pixel, r, g, b);

//    cout << (int) r << ", " << (int) g << ", " << (int) b << endl;
        rgbToHsl(r, g, b, h, s, l2);
if( h!=h )
    h = 0;
//cout << h << ", "<<c << ", " << l2 << endl;
        if( colorToL(image, pixel)>(255-threshold*threshold) ) {
            if( h<0. )
                h += TRIG_2PI;
            field->add(x, y, z, l, s*255., h/TRIG_2PI*255.);
            i++;
        }

    }

}

HaloField* ImagePlot::getField()
{
    return field;
}
