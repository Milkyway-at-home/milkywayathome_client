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

#include "nbody_priv.h"
#include "orbitintegrator.h"
#include "accelerations.h"

/* Simple orbit integrator in user-defined potential
    Written for BOINC Nbody
    willeb 10 May 2010 */
static void reverseOrbit(NBodyCtx* ctx, DwarfModel* model)
{
    mwvector acc, v, x;
    real t;

    InitialConditions* ic  = &model->initialConditions;

    const real tstop = ctx->time_orbit;
    const real dt    = ctx->orbit_timestep;

    // Set the initial conditions
    x = ic->position;
    v = ic->velocity;
    mw_incnegv(v);

    // Get the initial acceleration
    acc = acceleration(ctx, x);

    // Loop through time
    for (t = 0; t <= tstop; t += dt)
    {
        // Update the velocities and positions
        mw_incaddv_s(v, acc, dt);
        mw_incaddv_s(x, v, dt);

        // Compute the new acceleration
        acc = acceleration(ctx, x);
    }

    /* Report the final values (don't forget to reverse the velocities) */
    ic->position = x;
    ic->velocity = v;
    mw_incnegv(ic->velocity);
}

/* For each of the models, if its initial conditions needs to be reverse orbited, do it. */
void reverseModelOrbits(NBodyCtx* ctx)
{
    unsigned int i;

    for (i = 0; i < ctx->modelNum; ++i)
    {
        if (ctx->models[i].initialConditions.reverseOrbit)
            reverseOrbit(ctx, &ctx->models[i]);
    }
}

