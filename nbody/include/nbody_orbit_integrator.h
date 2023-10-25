/* Copyright 2010 Matthew Arsenault, Travis Desell, Boleslaw
Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail and
Rensselaer Polytechnic Institute.

Copyright 2016 Siddhartha Shelton

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

#ifndef _NBODY_ORBIT_INTEGRATOR_H_
#define _NBODY_ORBIT_INTEGRATOR_H_

#include "nbody_types.h"

void nbReverseOrbit(mwvector* finalPos,
                    mwvector* finalVel,
                    const Potential* pot,
                    mwvector pos,
                    mwvector vel,
                    real tstop,
                    real dt);

void nbReverseOrbits(mwvector* finalPos,
                    mwvector* finalVel,
                    const Potential* pot,
                    mwvector* pos,
                    mwvector* vel,
                    size_t len,
                    real tstop,
                    real dt);

void nbReverseOrbit_LMC(mwvector* finalPos,
                    mwvector* finalVel,
                    mwvector* LMCfinalPos,
                    mwvector* LMCfinalVel,
                    const Potential* pot,
                    mwvector pos,
                    mwvector vel,
                    mwvector LMCposition,
                    mwvector LMCvelocity,
                    mwbool LMCDynaFric,
                    real ftime,
                    real tstop,
                    real dt,
                    real LMCmass,
                    real LMCscale,
                    real coulomb_log);

void getLMCArray(mwvector ** shiftArrayPtr, size_t * shiftSizePtr);

void getLMCPosVel(mwvector * LMCposPtr, mwvector * LMCvelPtr);

void nbPrintReverseOrbit(mwvector* finalPos,
                         mwvector* finalVel,
                         const Potential* pot,
                         mwvector pos,
                         mwvector vel,
                         real tstop,
                         real tstopforward,
                         real dt);
    
#endif /* _NBODY_ORBIT_INTEGRATOR_H_ */

