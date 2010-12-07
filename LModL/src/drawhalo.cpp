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

#include "drawhalo.hpp"


const int PRINT_XSIZE = 3300;
const int PRINT_YSIZE = 3300;


BlurCircle::BlurCircle( int granShift )
{
    blurStats.granShift = granShift;
    blurStats.granSqShift = 2*granShift;
    blurStats.granularity = 1<<granShift;
    blurStats.granNorm = 1./(float)blurStats.granularity;
    blurStats.granNormAdd = -.5+1./(float)(blurStats.granularity*2);
    blurStats.cMult = 1/sqrt(2*acos(-1));
}

void BlurCircle::draw( SDL_Surface *surface, float x, float y, float diameter, float intensity )
{

    /// TODO /// resolve issue with magnitude>1. not expanding size
    if( diameter<3. ) {
        intensity *= diameter;
        diameter = 3.;
    }
    float radius = diameter/2.;
    float normMult = 232.;             /// TODO /// optimize
    float amp = 5.*intensity*255.;     /// TODO /// optimize

    int ix = int(x);
    int iy = int(y);

    float fx = x - (float) ix;
    float fy = y - (float) iy;

    float xgs = .5 - fx - int(radius);
    float yg = .5 - fy - int(radius);

    int xpe = (int) (x + radius);
    int ype = (int) (y + radius);

    for( int yp = (int) (y - int(radius)); yp<=ype; yg++, yp++ ) {
        float xg = xgs;
        for( int xp = (int) (x - int(radius)); xp<=xpe; xg++, xp++ ) {
            // Integrate over a single pixel
            int mSum = 0;
            float xg2s = xg+blurStats.granNormAdd;
            float yg2 = yg+blurStats.granNormAdd;
            for( int yGran = blurStats.granularity; yGran; yGran--, yg2 += blurStats.granNorm ) {
                float xg2 = xg2s;
                for( int xGran = blurStats.granularity; xGran; xGran--, xg2 += blurStats.granNorm ) {
                    float dist = sqrt(xg2*xg2+yg2*yg2);
                    if( dist>radius )
                        continue;
                    int total = (int) ((blurStats.cMult*exp(-9.*dist*dist/(radius*radius)/2.)*normMult)*amp)>>8;
                    mSum += total;
                }
            }
            int intensity = mSum>>(blurStats.granSqShift+2);
            /*  int lossMask = (1<<(blurStats.granSqShift+2))-1;
                int loss = lossMask&mSum;
                if( rand()&lossMask<loss )
                intensity++;
            */
            //intensity += rand()%(min(6, 1+(intensity>>1)));
            putPixelSumClip8(surface, xp, yp, min(0x7f, intensity));
        }
    }

}

HaloType::HaloType( float diameter, float maxLuminosity, int lumGranShift, int haloGranShift )
{
    this->lumGranShift = lumGranShift;
    lumGranularity = 1<<lumGranShift;
    this->haloGranShift = haloGranShift;
    haloGranularity = 1<<haloGranShift;
    haloGranfloat = 1<<haloGranShift;

    pointOffset = new SDL_Surface***[lumGranularity];
    for( int l = 0; l<lumGranularity; l++ ) {
        pointOffset[l] = new SDL_Surface**[haloGranularity];
        for( int y = 0; y<haloGranularity; y++ ) {
            pointOffset[l][y] = new SDL_Surface*[haloGranularity];
            for( int x = 0; x<haloGranularity; x++ )
                pointOffset[l][y][x] = NULL;
        }
    }
    set(diameter, maxLuminosity);
}

HaloType::~HaloType()
{
    for( int l = 0; l<lumGranularity; l++ ) {
        for( int y = 0; y<haloGranularity; y++ ) {
            for( int x = 0; x<haloGranularity; x++ )
                if( pointOffset[l][y][x]!=NULL )
                    SDL_FreeSurface(pointOffset[l][y][x]);
            delete [] pointOffset[l][y];
        }
        delete [] pointOffset[l];
    }
    delete [] pointOffset;
}

void HaloType::set( float diameter, float maxLuminosity )
{

    maxLum = maxLuminosity;
    lumMult = maxLuminosity/float(lumGranularity);
    lumDiv = float(lumGranularity)/maxLuminosity;
    float radius = diameter / 2.;

    for( int l = 0; l<lumGranularity; l++ )
        for( int y = 0; y<haloGranularity; y++ )
            for( int x = 0; x<haloGranularity; x++ )
                if( pointOffset[l][y][x] != NULL )
                    SDL_FreeSurface(pointOffset[l][y][x]);

    BlurCircle bc(BLUR_GRANULARITY);

    // Round to nearest 4 pixels to speed up blitting
    int maxDiameter = int(diameter+1.99999);
    maxDiameter = ((maxDiameter+0x03)>>2)<<2;
    float step = 1./haloGranularity;

//        float lumStep = maxLuminosity/float(lumGranularity);
    float lum = lumMult;
    for( int l = 0; l<lumGranularity; l++, lum += lumMult ) {

        float yOffset = 0.+step/2.;
        for( int y = 0; y<haloGranularity; y++, yOffset+=step ) {

            float xOffset = 0.+step/2.;
            for( int x = 0; x<haloGranularity; x++, xOffset+=step ) {
                SDL_Surface *surface = newSurface8(maxDiameter, maxDiameter);
                lockSurface(surface);
                bc.draw(surface, radius+xOffset, radius+yOffset, diameter, lum);
                unlockSurface(surface);
                pointOffset[l][y][x] = surface;
            }

        }

    }

    this->radius = radius;

}

float HaloType::getDiameter()
{
    return radius*2.;
}

float HaloType::getLum()
{
    return maxLum;
}

void HaloType::draw( SDL_Surface *surface, float x, float y, float luminosity, Uint32* palette )
{
    x -= radius;
    y -= radius;
    int xf = (int) (x*haloGranfloat);
    int yf = (int) (y*haloGranfloat);
    int xi = xf>>haloGranShift;
    int yi = yf>>haloGranShift;
    xf -=  xi<<haloGranShift;
    yf -=  yi<<haloGranShift;

    unsigned int iLum = int(luminosity*lumDiv);
    iLum = min((unsigned int)lumGranularity-1, iLum);
    blitSurfaceClipSumPalette(pointOffset[iLum][yf][xf], surface, xi, yi, palette);
}

/*
void HaloType::draw( SDL_Surface* surface, fix32 x, fix32 y )
{
  blitSurfaceClipSum8(pointOffset[y.fpart(haloGranShift)][x.fpart(haloGranShift)], surface, x.ipart(), y.ipart());
}
*/

void HaloType::_drawTest() const
{
    SDL_Surface *display = SDL_GetVideoSurface();
    clearSurface(display);
    for( int y = 0; y<haloGranularity; y++ )
        for( int x = 0; x<haloGranularity; x++ )
            blitSurfaceClipSum8(pointOffset[lumGranularity-1][y][x], display, x*(pointOffset[lumGranularity-1][y][x])->w, y*(pointOffset[lumGranularity-1][y][x])->h);
}

HaloField::HaloField( int pointTotal )
{
    arrayTotal = pointTotal;
    field = new HaloPoint*[pointTotal];
    for( int i = 0; i<pointTotal; i++  )
        field[i] = new HaloPoint;
    stackPtr = 0;
    stackEndPtr = 0;
    lightAdd = 0.;
    axesLimit = 0;
}

HaloField::~HaloField()
{
    for( int i = 0; i<arrayTotal; i++  )
        delete field[i];
    delete [] field;
}

void HaloField::setLum( float lightness )
{
     lightAdd = lightness;
}

float HaloField::getLum() const
{
    return lightAdd;
}

Vector3d HaloField::getCenter() const
{
    Vector3d center(0., 0., 0.);
    for( int i = 0; i<stackEndPtr; i++ )
        center += field[i]->position;
    center /= stackEndPtr;
    return center;
}

void HaloField::clearField()
{
    stackPtr = 0;
    stackEndPtr = 0;
}

void HaloField::add( float x, float y, float z, float L, int C, int h )
{
#ifndef NDEBUG
    if( stackPtr>arrayTotal ) {
        cerr << "Halo-field stack pointer overflow" << endl;
        exit(1);
    }
#endif
    field[stackPtr]->position.x = x;
    field[stackPtr]->position.y = y;
    field[stackPtr]->position.z = z;
    field[stackPtr]->lightness = L;
//cout << "!!!" << getLightnessColor32(C, h) << " < " << L << ", " << C << ", " << h << endl << flush;
    field[stackPtr]->palette = getLightnessColor32p(C, h);
//if( (int)(field[stackPtr]->palette)==0x592040 ){
    stackPtr++;
    if( stackPtr>stackEndPtr )
        stackEndPtr = stackPtr;
    if( x>axesLimit )
        axesLimit = x;
    if( y>axesLimit )
        axesLimit = y;
    if( z>axesLimit )
        axesLimit = z;
}

void HaloField::add( float x, float y, float z, Uint8 G, Uint8 R, Uint8 B )
{
#ifndef NDEBUG
    if( stackPtr>arrayTotal ) {
        cerr << "Halo-field stack pointer overflow" << endl;
        exit(1);
    }
#endif
    field[stackPtr]->position.x = x;
    field[stackPtr]->position.y = y;
    field[stackPtr]->position.z = z;
   
    double h, s, l;
    rgbToHsl(R, G, B, h, s, l);

    field[stackPtr]->lightness = l;

//cout << "!!!" << getLightnessColor32(C, h) << " < " << L << ", " << C << ", " << h << endl << flush;
    field[stackPtr]->palette = getLightnessColor32p(s, h);
//if( (int)(field[stackPtr]->palette)==0x592040 ){
    stackPtr++;
    if( stackPtr>stackEndPtr )
        stackEndPtr = stackPtr;
    if( x>axesLimit )
        axesLimit = x;
    if( y>axesLimit )
        axesLimit = y;
    if( z>axesLimit )
        axesLimit = z;
}

void HaloField::set( int index, float x, float y, float z, float l, float c, float h )
{
    stackPtr = index;
    add(x, y, z, l, c, h);
}

void HaloField::drawAxes( SDL_Surface* surface, const Camera& cv, HaloType& lineBlur )
{
    lockSurface(surface);
    drawLine(surface, &lineBlur, &cv, 1., Vector3d(0., 0., 0.), Vector3d(axesLimit, 0., 0.));
    drawLine(surface, &lineBlur, &cv, 1., Vector3d(0., 0., 0.), Vector3d(0., axesLimit, 0.));
    drawLine(surface, &lineBlur, &cv, 1., Vector3d(0., 0., 0.), Vector3d(0., 0., axesLimit));
    unlockSurface(surface);
}

void HaloField::drawCamera( SDL_Surface *surface, Camera *cv, HaloType& lineBlur, HaloType& endBlur )
{
    lockSurface(surface);
    if( surface->h<480 )
        return;

    // Draw camera
    Vector3d eye, up, right;
    cv->getCameraDebug(eye, up, right, 50.);
    float xOffset = surface->w - 120;
    float yOffset = 120;
    drawLine(surface, &lineBlur, 1., xOffset, yOffset, xOffset+up.x, yOffset-up.z);
    drawLine(surface, &lineBlur, 1., xOffset, yOffset, xOffset+eye.x, yOffset-eye.z);
    endBlur.draw(surface, xOffset+eye.x, yOffset-eye.z, 10.);
    unlockSurface(surface);
}

void HaloField::draw( SDL_Surface* surface, Camera *cv, HaloType &haloType, int skip )
{
    lockSurface(surface);
    Vector3d map;
    int skipCount = 0;
    for( int i = 0; i<stackEndPtr; i++ ) {
        skipCount++;
        if( skipCount==skip ) {
            skipCount = 0;
            if( cv->getCameraProjection(field[i]->position, map) ) {
                cv->getDisplayOffset(map);
                haloType.draw(surface, map.x, map.y, max(0.f, (field[i]->lightness)+lightAdd), field[i]->palette);
            }
        }
    }
    unlockSurface(surface);
}

FieldAnimation::FieldAnimation( int bpp, float fps, bool fullScreen, string caption, string iconFileName )
{
    resetFieldPointer();
    event = new SDL_Event;

    this->fps = fps;
    fpsTimer = new StepTimer(1./fps);
    keyTimer = new StepTimer(1./2.);
    fastKeyTimer = new StepTimer(1./5.);

    // Set up environment
    display = setVideo(bpp, -1, -1, fullScreen, caption, NULL);
    setPaletteTable();
    cv = new Camera(display->w, display->h, fps);
    sizeAdjust = float(display->w)/1024.;
    lineBlur = new HaloType(11.*sizeAdjust, .1, 0, 3);
    endBlur = new HaloType(21.*sizeAdjust, 1., 0, 3);

    cameraDisplay = false;
    axesDisplay = false;
}

void FieldAnimation::resetFieldPointer()
{
    totalFields = 0;
    fieldPtr = 0;
    fieldModIndex = 0;
}

void FieldAnimation::add( HaloField *newField, float blurDiameter )
{

    if( 100==totalFields ) {
        cerr << "Currently only 100 fields supported\n";
        exit(1);
    }
    field[totalFields].blur = new HaloType(blurDiameter*sizeAdjust, 1., 6, 1);

    field[totalFields].data = newField;
    field[totalFields].LIGHTNESSdjust = 1.;
    field[totalFields].skipDraw = 1;
    totalFields++;

}

void FieldAnimation::showCamera()
{
    cameraDisplay = true;
}

void FieldAnimation::hideCamera()
{
    cameraDisplay = false;
}

void FieldAnimation::showAxes()
{
    axesDisplay = true;
}

void FieldAnimation::hideAxes()
{
    axesDisplay = false;
}

void FieldAnimation::applyFrame( SDL_Surface* surface )
{
    clearSurface(surface);

    if( surface->w!=display->w ) {

        cv->setSize(surface->w, surface->h);
        float adjust = (float) surface->w/display->w;
        for( int i = 0; i<totalFields; i++ ) {
            HaloType* tempBlur = new HaloType((field[i].blur->getDiameter()*adjust), field[i].blur->getLum(), 5, 1);
            field[i].data->draw(surface, cv, *tempBlur);
            delete tempBlur;
        }

        HaloType* tempBlur = NULL;
        if( cameraDisplay || axesDisplay )
            tempBlur = new HaloType((lineBlur->getDiameter()*adjust), lineBlur->getLum(), 5, 1);
        if( cameraDisplay ) {
            HaloType* tempBlur2 = new HaloType((endBlur->getDiameter()*adjust), endBlur->getLum(), 5, 1);
            field[0].data->drawCamera(surface, cv, *tempBlur, *tempBlur2);
            delete tempBlur2;
        }
        if( axesDisplay )
            field[0].data->drawAxes(surface, *cv, *tempBlur);
        if( cameraDisplay || axesDisplay )
            delete tempBlur;

        cv->setSize(display->w, display->h);

    }
    else {
        for( int i = 0; i<totalFields; i++ )
            field[i].data->draw(surface, cv, *(field[i].blur), field[i].skipDraw);
        if( cameraDisplay )
            field[0].data->drawCamera(surface, cv, *lineBlur, *endBlur);
        if( axesDisplay )
            field[0].data->drawAxes(surface, *cv, *lineBlur);
    }
}

bool FieldAnimation::pollDemo()
{

    if( fpsTimer->getSignal() ) {

        // Display iteration
        applyFrame(display);
        SDL_Flip(display);

    }

    return cv->getAutoStatus();

}


bool FieldAnimation::pollEvent()
{
    if( !fpsTimer->getSignal() )
        return false;
    pollDemo();

    // Handle events while skipping frame-draws as needed
    do {

        // Update camera position / get user input
        if( cv->pollEvent(event) ) {

            // Handle key-press (no pause)

            switch( event->key.keysym.sym ) {

            case SDLK_ESCAPE:
                exit(0);
                break;

            default:
                break;

            }

            // Handle key-press (with pause)

            if( keyTimer->getSignal() )

                switch( event->key.keysym.sym ) {

                case SDLK_c:
                    cameraDisplay = !cameraDisplay;
                    keyTimer->setInc();
                    fastKeyTimer->setInc();
                    break;

                case SDLK_a:
                    axesDisplay = !axesDisplay;
                    keyTimer->setInc();
                    fastKeyTimer->setInc();
                    break;

                case SDLK_BACKSPACE:
                    if( event->key.keysym.mod & KMOD_SHIFT )
                        cv->setFocusPoint(Vector3d(0., 0., 0.));
                    else
                        cv->setFocusPoint(Vector3d(-8., 0., 0.));

                    keyTimer->setInc();
                    fastKeyTimer->setInc();
                    break;

                case SDLK_RETURN:
                    if( event->key.keysym.mod & KMOD_SHIFT )
                        cv->moveToPoint(Vector3d(0., 0., 0.));
                    else if( event->key.keysym.mod & KMOD_ALT )
                        SDL_WM_ToggleFullScreen(display);  /// TODO /// not working yet
                    else
                        cv->moveToPoint(Vector3d(-8., 0., 0.));
                    keyTimer->setInc();
                    fastKeyTimer->setInc();
                    break;

                case SDLK_1:
                    if( fieldModIndex/10*10<totalFields )
                        fieldModIndex = fieldModIndex/10*10;
                    break;

                case SDLK_2:
                    if( fieldModIndex/10*10+1<totalFields )
                        fieldModIndex = fieldModIndex/10*10+1;
                    break;

                case SDLK_3:
                    if( fieldModIndex/10*10+2<totalFields )
                        fieldModIndex = fieldModIndex/10*10+2;
                    break;

                case SDLK_4:
                    if( fieldModIndex/10*10+3<totalFields )
                        fieldModIndex = fieldModIndex/10*10+3;
                    break;

                case SDLK_5:
                    if( fieldModIndex/10*10+4<totalFields )
                        fieldModIndex = fieldModIndex/10*10+4;
                    break;

                case SDLK_6:
                    if( fieldModIndex/10*10+5<totalFields )
                        fieldModIndex = fieldModIndex/10*10+5;
                    break;

                case SDLK_7:
                    if( fieldModIndex/10*10+6<totalFields )
                        fieldModIndex = fieldModIndex/10*10+6;
                    break;

                case SDLK_8:
                    if( fieldModIndex/10*10+7<totalFields )
                        fieldModIndex = fieldModIndex/10*10+7;
                    break;

                case SDLK_9:
                    if( fieldModIndex/10*10+8<totalFields )
                        fieldModIndex = fieldModIndex/10*10+8;
                    break;

                case SDLK_0:
                    if( fieldModIndex/10*9<totalFields )
                        fieldModIndex = fieldModIndex/10*9;
                    break;

                default:
                    break;

                }

            // Handle key-press (with short pause)

            if( fastKeyTimer->getSignal() )

                switch( event->key.keysym.sym ) {

				case SDLK_p:
				case SDLK_EQUALS:
                    if( event->key.keysym.mod & KMOD_CTRL ) {
                        if( !keyTimer->getSignal() )
                            break;
                        if( fieldModIndex<totalFields-1 )
                            fieldModIndex++;
                    }
                    else if( event->key.keysym.mod & KMOD_ALT ) {
                        if( field[fieldModIndex].skipDraw>1 )
                            field[fieldModIndex].skipDraw--;
                    }
                    else if( event->key.keysym.mod & KMOD_SHIFT ) {
                        float diameter = field[fieldModIndex].blur->getDiameter();
                        float maxLum = field[fieldModIndex].blur->getLum();
                            field[fieldModIndex].blur->set(diameter+.5, maxLum);
                    }
                    else {
                        float lum = field[fieldModIndex].data->getLum();
                        if( lum<10. )
                            field[fieldModIndex].data->setLum(lum+.01);
                    }
                    keyTimer->setInc();
                    fastKeyTimer->setInc();
                    break;

				case SDLK_m:
				case SDLK_MINUS:
                    if( event->key.keysym.mod & KMOD_CTRL ) {
                        if( !keyTimer->getSignal() )
                            break;
                        if( fieldModIndex>0 )
                            fieldModIndex--;
                    }
                    else if( event->key.keysym.mod & KMOD_ALT ) {
                        if( field[fieldModIndex].skipDraw<10 )
                            field[fieldModIndex].skipDraw++;
                    }
                    else if( event->key.keysym.mod & KMOD_SHIFT ) {
                        float diameter = field[fieldModIndex].blur->getDiameter();
                        float maxLum = field[fieldModIndex].blur->getLum();
                        if( diameter>3.5 )
                            field[fieldModIndex].blur->set(diameter-.5, maxLum);
                    }
                    else {
                        float lum = field[fieldModIndex].data->getLum();
                        field[fieldModIndex].data->setLum(lum-.01);
                    }
                    keyTimer->setInc();
                    fastKeyTimer->setInc();
                    break;

                case SDLK_PRINT:
                    if( event->key.keysym.mod & KMOD_CTRL )
                    {
                        SDL_Surface* printImage = newSurface32p(PRINT_XSIZE, PRINT_YSIZE);
//SDL_Surface* printImage = SDL_LoadBMP("hq_image.bmp"); /// TODO /// STUB
                        applyFrame(printImage);
/// TODO /// STUB for 24-bit bmp saves
SDL_Surface* img24 = SDL_CreateRGBSurface(SDL_SWSURFACE, PRINT_XSIZE, PRINT_YSIZE, 24,
#if SDL_BYTEORDER == SDL_BIG_ENDIAN
                             0xff0000,
                             0x00ff00,
                             0x0000ff,
#else
                             0x0000ff,
                             0x00ff00,
                             0xff0000,
#endif
0);

for(int y = 0; y<PRINT_YSIZE; y++)
for(int x = 0; x<PRINT_XSIZE; x++)
putPixel24(img24, x, y, getPixel32p(printImage, x, y));
//SDL_BlitSurface(printImage, NULL, img24, NULL);

                        SDL_SaveBMP(img24, "hq_image.bmp");
                        SDL_FreeSurface(printImage);
                        SDL_FreeSurface(img24);
                    }
                    else
                        screenShot("screenshot.bmp");  /// TODO ///  Make this save incrementally

                    keyTimer->setInc();
                    fastKeyTimer->setInc();
                    break;

                default:
                    break;

                }

            }

        // Report to timer that frame was applied
        fpsTimer->inc();

    } while( fpsTimer->getSignal() );

    return true;

}

bool clip2d( float &x, float &y, float xc, float yc, float xMin, float xMax, float yMin, float yMax ) {

    enum Region{cr, lt, up, ltUp, rt, ltRt, upRt, ltUpRt, dn, ltDn, upDn, ltUpDn, rtDn, ltRtDn, upRtDn, ltUpRtDn} region = cr;

    float upInt = 0., dnInt = 0., rtInt = 0., ltInt = 0.;

    if( x<xMin ) {
        region = Region(region | lt);
        float t = (xMin - x) / xc;
        ltInt = y + yc*t;
    }
    else if( x>xMax ) {
        region = Region(region | rt);
        float t = (xMax - x) / xc;
        rtInt = y + yc*t;
    }
    if( y<yMin ) {
        region = Region(region | up);
        float t = (yMin - y) / yc;
        upInt = x + xc*t;
    }
    else if( y>yMax ) {
        region = Region(region | dn);
        float t = (yMax - y) / yc;
        dnInt = x + xc*t;
    }

    switch( region )
    {

    case cr:
        // Nothing to do
        break;

    case lt:
        // Clip to left intersection
        if( ltInt<yMin || ltInt>yMax )
            return false;
        x = xMin;
        y = ltInt;
        break;

    case up:
        // Clip to up intersection
        if( upInt<xMin || upInt>xMax )
            return false;
        x = upInt;
        y = yMin;
        break;

    case ltUp:
        // Clip to either left or up intersection
        if( ltInt<yMin || ltInt>yMax ) {
            if( upInt<xMin || upInt>xMax )
                return false;
            // Use up intersection
            x = upInt;
            y = yMin;
        }
        else {
            // Use left intersection
            x = xMin;
            y = ltInt;
        }
        break;

    case rt:
        // Clip to right intersection
        if( rtInt<yMin || rtInt>yMax )
            return false;
        x = xMax;
        y = rtInt;
        break;

    case upRt:
        // Clip to either up or right intersection
        if( upInt<xMin || upInt>xMax ) {
            if( rtInt<yMin || rtInt>yMax )
                return false;
            x = xMax;
            y = rtInt;
        }
        else {
            x = upInt;
            y = yMin;
        }
        break;

    case dn:
        // Clip to down intersection
        if( dnInt<xMin || dnInt>xMax )
            return false;
        x = dnInt;
        y = yMax;
        break;

    case ltDn:
        // Clip to either left or down intersection
        if( ltInt<yMin || ltInt>yMax ) {
            if( dnInt<xMin || dnInt>xMax )
                return false;
            x = dnInt;
            y = yMax;
        }
        else {
            x = xMin;
            y = ltInt;
        }
        break;

    case rtDn:
        // Clip to either right or down intersection
        if( rtInt<yMin || rtInt>yMax ) {
            if( dnInt<xMin || dnInt>xMax )
                return false;
            x = dnInt;
            y = yMax;
        }
        else {
            x = xMax;
            y = rtInt;
        }
        break;

    default:
        assert(false);

    }

    return true;

}

void drawLine( SDL_Surface *surface, HaloType *bc, float lum, float xs, float ys, float xe, float ye, Uint32* palette )
{

    // Perform line clipping
    float xMax = (float) (surface->w-1);
    float yMax = (float) (surface->h-1);
    float xc = xe-xs;
    float yc = ye-ys;
    if( !clip2d(xs, ys, xc, yc, 0., xMax, 0., yMax) || !clip2d(xe, ye, xc, yc, 0., xMax, 0., yMax) )
        return;

    /// TODO /// Optimize with clipping
    /// TODO /// Optimize with non-redundant blur blitting
    // The blur circle passed to this routine should be about blur-diameter 2*blur_diameter times less intense than the expected line intensity

    float dist = 2. * sqrt(abs(xe-xs)*abs(xe-xs)+abs(ye-ys)*abs(ye-ys));
    int iDist = int(dist);

    float xStep = (xe-xs)/dist;
    float yStep = (ye-ys)/dist;

    float x = xs;
    float y = ys;
    for( int i = 0; i<iDist; i++, x += xStep, y += yStep )
        bc->draw(surface, x, y, lum, palette);

}

void drawLine( SDL_Surface *surface, HaloType *bc, const Camera* cam, float lum, const Vector3d start, const Vector3d end, Uint32* palette )
{

    /// TODO /// perhaps these should be 2D vectors - no apparent application for z
    /// TODO /// perhaps increasing intensity appropriately for draw repetition within the circle of confusion would create a desireable effect
    Vector3d startP, endP;
    cam->get3dTransform(start, startP);
    cam->get3dTransform(end, endP);

    if( startP.z<=0.000001 && endP.z<=0.000001 )
        return;

    // Clip line along z-axis
    if( startP.z<=0.000001 ) {

        Vector3d line = endP;
        line -= startP;
        float t = startP.z/line.z;
        startP.x -= line.x*t;
        startP.y -= line.y*t;
        endP.z = 0.000001;

    }
    else if( endP.z<=0.000001 ) {

        Vector3d line = endP;
        line -= startP;
        float t = endP.z/line.z;
        endP.x -= line.x*t;
        endP.y -= line.y*t;
        endP.z = 0.000001;

    }

    cam->get3dProjection(startP);
    cam->getDisplayOffset(startP);

    cam->get3dProjection(endP);
    cam->getDisplayOffset(endP);

    drawLine(surface, bc, lum, startP.x, startP.y, endP.x, endP.y, palette);

}
