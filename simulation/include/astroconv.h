/*****************************************************************************
 *                                                                           *
 *  Copyright (C) 2010 Shane Reilly, Ben Willet, Matthew Newby, Heidi        *
 *  Newberg, Malik Magdon-Ismail, Carlos Varela, Boleslaw Szymanski, and     *
 *  Rensselaer Polytechnic Institute                                         *
 *                                                                           *
 *  This file is part of the MilkyWay@Home Project.                          *
 *                                                                           *
 *  This program is free software: you can redistribute it and/or modify     *
 *  it under the terms of the GNU General Public License as published by     *
 *  the Free Software Foundation, either version 3 of the License, or        *
 *  (at your option) any later version.                                      *
 *                                                                           *
 *  This program is distributed in the hope that it will be useful,          *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of           *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the             *
 *  GNU General Public License for more details.                             *
 *                                                                           *
 *  You should have received a copy of the GNU General Public License        *
 *  along with this program. If not, see <http://www.gnu.org/licenses/>.     *
 *                                                                           *
 *  Shane Reilly                                                             *
 *  reills2@cs.rpi.edu                                                       *
 *                                                                           *
 *****************************************************************************/

#ifndef _ASTROCONV_H_
#define _ASTROCONV_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <math.h>
#include "trigtbl.hpp"

inline void raDecRadToCart( double ra, double dec, double r, double &x, double &y, double &z )
{

    // RA/Dec to cartesian conversion
    double t0 = cos(dec*TRIG_DEG_TO_RAD)*cos(ra*TRIG_DEG_TO_RAD);
    double t1 = cos(dec*TRIG_DEG_TO_RAD)*sin(ra*TRIG_DEG_TO_RAD);
    double t2 = sin(dec*TRIG_DEG_TO_RAD);

    double mm[3][3] =
        {   { -.06699 , -.87276, -.48354 } ,
            {  .49273 , -.45035,  .74458 } ,
            { -.86760 , -.18837,  .46020 }  };

    double c0 = t0*mm[0][0] + t1*mm[0][1] + t2*mm[0][2];
    double c1 = t0*mm[1][0] + t1*mm[1][1] + t2*mm[1][2];
    double c2 = t0*mm[2][0] + t1*mm[2][1] + t2*mm[2][2];

    x = r*c0;
    y = r*c1;
    z = r*c2;

}

inline double absMagToLum( double absMagnitude )
    // Calculate luminosity relative to sun given the absolute magnitude
{
    return pow(10., 2.-.4*absMagnitude);
}

#ifdef __cplusplus
}
#endif

#endif /* _ASTROCONV_H_ */
