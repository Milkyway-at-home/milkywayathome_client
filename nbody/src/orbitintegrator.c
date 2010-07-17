/* Simple orbit integrator in user-defined potential
    Written for BOINC Nbody
    willeb 10 May 2010 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "nbody_priv.h"
#include "orbitintegrator.h"
#include "accelerations.h"
#include "real.h"

inline void acceleration(vectorptr restrict acc, const NBodyCtx* ctx, const vectorptr restrict pos)
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


