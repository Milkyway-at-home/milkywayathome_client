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

#include "separation_types.h"
#include "separation_constants.h"

OLD_GCC_EXTERNINLINE
inline real calcG(const real coords)
{
    return 5.0 * (mw_log10(1000.0 * coords) - 1.0) + absm;
}

mwvector stripe_normal(int wedge);

mwvector lbr2xyz(const ASTRONOMY_PARAMETERS* ap, const mwvector lbr);
mwvector xyz_mag(const ASTRONOMY_PARAMETERS* ap, mwvector point, real offset);

LB gc2lb(const int wedge, const real mu, const real nu);

#endif /* _COORDINATES_H_ */

