/*
Copyright 2008-2010 Travis Desell, Matthew Arsenault, Dave Przybylo,
Nathan Cole, Boleslaw Szymanski, Heidi Newberg, Carlos Varela, Malik
Magdon-Ismail and Rensselaer Polytechnic Institute.

This file is part of Milkway@Home.

Milkyway@Home is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Milkyway@Home is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Milkyway@Home.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef _COORDINATES_H_
#define _COORDINATES_H_

#include "separation_constants.h"
#include "separation_types.h"
#include "milkyway_cl.h"
#include "milkyway_extra.h"

//vickej2 for sgr stripes, the great circles are defined thus:
//sgr stripes run parallel to sgr longitude lines, centered on lamda=2.5*wedge number
//with mu=0 at the sgr equator and increasing in the +z direction (increasing from the equator with beta)
//and nu=0 at the center and increasing in the -y direction (inversely to lamda)
//in this manner an equatorial stripe of standard coordinate conventions is created.

/* Return eta from stripe number */
HOT CONST_F ALWAYS_INLINE OLD_GCC_EXTERNINLINE
inline real atEtaFromStripeNumber_rad(int wedge)
{
    return wedge * d2r(stripeSeparation) - d2r((real) 57.5) - (wedge > 46 ? M_PI : 0.0);
}

HOT CONST_F ALWAYS_INLINE OLD_GCC_EXTERNINLINE
inline real atEtaFromStripeNumber_deg(int wedge)
{
    return wedge * stripeSeparation - 57.5 - (wedge > 46 ? 180.0 : 0.0);
}

LB gc2lb(const int wedge, const real mu, const real nu);

#endif /* _COORDINATES_H_ */

