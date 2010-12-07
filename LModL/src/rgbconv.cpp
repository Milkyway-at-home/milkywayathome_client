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

#include "rgbconv.hpp"


void rgbToXyl( Uint8 r, Uint8 g, Uint8 b, double &x, double &y, double &l )
{

    double dr = r-127.5;
    double dg = g-127.5;
    double db = b-127.5;

    // Rotate X-axis

    double cosa = sqrt(2.)/2.;
    double sina = sqrt(2.)/2.;

    double xp = dr;
    double yp = dg*cosa - db*sina;
    double lp = dg*sina + db*cosa;

    // Rotate Y-axis

    cosa = sqrt(2./3.);
    sina = -1./sqrt(3.);

    x = lp*sina + xp*cosa;
    y = yp;
    l = lp*cosa - xp*sina;

}

void xylToRgb( double x, double y, double l, Uint8 &r, Uint8 &g, Uint8 &b )
{

    // Rotate Y-axis

    double cosa = sqrt(2./3.);
    double sina = 1./sqrt(3.);

    double xp = l*sina + x*cosa;
    double yp = y;
    double lp = l*cosa - x*sina;

    // Rotate X-axis

    cosa = sqrt(2.)/2.;
    sina = -sqrt(2.)/2.;

    double yp2 = yp*cosa - lp*sina;
    lp = yp*sina + lp*cosa;


    double dr = xp+127.5 + .5;
    double dg = yp2+127.5 + .5;
    double db = lp+127.5 + .5;

    /// STUB /// This can be handled better by changing line size
    if( dr>255. )
        dr = 255.;
    if( dg>255. )
        dg = 255.;
    if( db>255. )
        db = 255.;
    if( dr<0. )
        dr = 0.;
    if( dg<0. )
        dg = 0.;
    if( db<0. )
        db = 0.;

    r = (Uint8) dr;
    g = (Uint8) dg;
    b = (Uint8) db;

}

void xylToRgb( double x, double y, double l, double &r, double &g, double &b )
{

    // Rotate Y-axis

    double cosa = sqrt(2./3.);
    double sina = 1./sqrt(3.);

    double xp = l*sina + x*cosa;
    double yp = y;
    double lp = l*cosa - x*sina;

    // Rotate X-axis

    cosa = sqrt(2.)/2.;
    sina = -sqrt(2.)/2.;

    double yp2 = yp*cosa - lp*sina;
    lp = yp*sina + lp*cosa;


    double dr = xp+127.5 + .5;
    double dg = yp2+127.5 + .5;
    double db = lp+127.5 + .5;

    /// STUB /// This can be handled better by changing line size
    if( dr>255. )
        dr = 255.;
    if( dg>255. )
        dg = 255.;
    if( db>255. )
        db = 255.;
    if( dr<0. )
        dr = 0.;
    if( dg<0. )
        dg = 0.;
    if( db<0. )
        db = 0.;

    r = dr/255.;
    g = dg/255.;
    b = db/255.;

}

void xylToHsl( double x, double y, double li, double &h, double &s, double &l )
{

    l = li/(sqrt(3.)*127.5);

    h = atan2(y, x);

    s = sqrt(x*x+y*y);
    s /= 255.*sqrt(2./3.);
    double absL = l<0. ? -l:l;
    if( absL>1./3. )
        s /= 3./2.*(1.-absL);

    l /= 2.;
    l += .5;

}

void hslToXyl( double h, double s, double li, double &x, double &y, double &l )
{

    l -= .5;
    l *= 2.;

    if( h!=h ) {
        x = 0.;
        y = 0.;
    }
    else {

        double absL = l<0. ? -l:l;
        if( absL>1./3. ) {
            s *= 3./2.*(1.-absL);
        }
        s *= 255.*sqrt(2./3.);
        x = s*cos(h);
        y = s*sin(h);

    }

    l = (li-.5)*(sqrt(3.)*255.);

}

void rgbToHsl( Uint8 r, Uint8 g, Uint8 b, double &h, double &s, double &l )
    // Returns hue, saturation and lighting
    //   h is a degree in radians
    //   s values are from 0. - 1.
    //   l values are from 0. - 1.
    //   r, g, b, are integers from 0-255
{

    double x, y;
    rgbToXyl(r, g, b, x, y, l);
    xylToHsl(x, y, l, h, s, l);
    if( r==g && r==b ) {
        h = 0./0.;
        s = 0.;
    }

}

void hslToRgb( double h, double s, double l, Uint8 &r, Uint8 &g, Uint8 &b )
    // Returns machine R, G, B values
    //   h is a degree in radians
    //   s values are from 0. - 1.
    //   l values are from 0. - 1.
    //   r, g, b, are integers from 0-255
{

    double x, y;
    hslToXyl(h, s, l, x, y, l);
    xylToRgb(x, y, l, r, g, b);

}

void hslToRgb( double h, double s, double l, double &r, double &g, double &b )
    // Returns machine R, G, B values
    //   h is a degree in radians
    //   s values are from 0. - 1.
    //   l values are from 0. - 1.
    //   r, g, b, are integers from 0-255
{

    double x, y;
    hslToXyl(h, s, l, x, y, l);
    xylToRgb(x, y, l, r, g, b);

}

