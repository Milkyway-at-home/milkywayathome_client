/*
 *  Copyright (c) 2008-2010 Travis Desell, Nathan Cole
 *  Copyright (c) 2008-2010 Boleslaw Szymanski, Heidi Newberg
 *  Copyright (c) 2008-2010 Carlos Varela, Malik Magdon-Ismail
 *  Copyright (c) 2008-2011 Rensselaer Polytechnic Institute
 *  Copyright (c) 2010-2011 Matthew Arsenault
 *
 *  This file is part of Milkway@Home.
 *
 *  Milkway@Home is free software: you may copy, redistribute and/or modify it
 *  under the terms of the GNU General Public License as published by the
 *  Free Software Foundation, either version 3 of the License, or (at your
 *  option) any later version.
 *
 *  This file is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_

#include "separation_types.h"

void freeStreams(Streams* streams);

void calcIntegralStepSizes(IntegralArea* i);

IntegralArea* readParameters(const char* file,
                             AstronomyParameters* ap,
                             BackgroundParameters* bgp,
                             Streams* streams);

int setParameters(AstronomyParameters* ap,
                  BackgroundParameters* bgp,
                  Streams* streams,
                  const real* parameters,
                  unsigned int numberParameters);


#endif /* _PARAMETERS_H_ */

