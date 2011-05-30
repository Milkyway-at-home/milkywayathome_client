/*
Copyright (C) 2011  Matthew Arsenault
Copyright (c) 2011 Rensselaer Polytechnic Institute.

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

#include "nbody_util.h"
#include "milkyway_math.h"

mwvector nbodyCenterOfMass(const NBodyState* st)
{
    unsigned int i;
    const Body* b;
    mwvector cm = ZERO_VECTOR;
    mwvector tmp;
    real mass = 0.0;

    for (i = 0; i < st->nbody; ++i)
    {
        b = &st->bodytab[i];

        tmp = mw_mulvs(Pos(b), Mass(b));
        mass += Mass(b);
        mw_incaddv(cm, tmp);
    }

    mw_incdivs(cm, mass);

    return cm;
}

