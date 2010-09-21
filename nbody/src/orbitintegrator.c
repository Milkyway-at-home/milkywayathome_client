/* Copyright 2010 Matthew Arsenault, Travis Desell, Boleslaw
Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail and
Rensselaer Polytechnic Institute.

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

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "nbody_priv.h"
#include "orbitintegrator.h"

inline void acceleration(vectorptr RESTRICT acc, const NBodyCtx* ctx, const vectorptr RESTRICT pos)
{
    /* lookup table for functions for calculating accelerations */
    static const HaloAccel haloFuncs[] = { [LogarithmicHalo] = logHaloAccel,
                                           [NFWHalo]         = nfwHaloAccel,
                                           [TriaxialHalo]    = triaxialHaloAccel };

    static const DiskAccel diskFuncs[] = { [ExponentialDisk]   = exponentialDiskAccel,
                                           [MiyamotoNagaiDisk] = miyamotoNagaiDiskAccel };

    static const SphericalAccel sphFuncs[] = { [SphericalPotential] = sphericalAccel };

    vector acctmp1 = ZERO_VECTOR;
    vector acctmp2 = ZERO_VECTOR;

    /* Use the type of potential to index into the table, and use the
     * appropriate function */

    diskFuncs[ctx->pot.disk.type](acctmp1, &ctx->pot.disk, pos);
    haloFuncs[ctx->pot.halo.type](acctmp2, &ctx->pot.halo, pos);

    ADDV(acc, acctmp1, acctmp2);

    sphFuncs[ctx->pot.sphere[0].type](acctmp1, &ctx->pot.sphere[0], pos);

    /* add the resulting vectors */
    INCADDV(acc, acctmp1);
}

/* Simple orbit integrator in user-defined potential
    Written for BOINC Nbody
    willeb 10 May 2010 */
void reverseOrbit(InitialConditions* fc, const NBodyCtx* ctx, InitialConditions* ic)
{
    vector acc, v, x;
    real t;
    int i;

    const real tstop = ctx->model.time_orbit;
    const real dt    = ctx->model.orbit_timestep;

    // Set the initial conditions
    SETV(x, ic->position);
    SETV(v, ic->velocity);
    INCNEGV(v);

    // Get the initial acceleration
    acceleration(acc, ctx, x);

    // Loop through time
    for (t = 0; t <= tstop; t += dt)
    {
        // Update the velocities and positions
        for (i = 0; i < 3; ++i)
        {
            v[i] += acc[i] * dt;
            x[i] += v[i] * dt;
        }
        // Compute the new acceleration
        acceleration(acc, ctx, x);
    }

    /* Report the final values (don't forget to reverse the velocities) */
    SETV(fc->position, x);
    SETV(fc->velocity, v);
    INCNEGV(fc->velocity);

    fc->useGalC = ic->useGalC;  /* Not actually necessary */
    fc->useRadians = ic->useRadians;

}


