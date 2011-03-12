/*
Copyright (C) 2011  Matthew Arsenault

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

#include <openssl/sha.h>

#include "milkyway_util.h"
#include "nbody_priv.h"
#include "nbody_lua_types.h"
#include "milkyway_lua_marshal.h"
#include "nbody_plummer.h"

static char* showHash(const unsigned char* hash, int initializer)
{
    char* buf;
    const uint32_t* hashi = (const uint32_t*) hash;

    if (!hash)
        return NULL;

    buf = (char*) mwMalloc(2 * SHA_DIGEST_LENGTH + 1);
    buf[2 * SHA_DIGEST_LENGTH] = '\0';

    if (initializer)
        sprintf(buf, "{ 0x%x, 0x%x, 0x%x, 0x%x, 0x%x }", hashi[0], hashi[1], hashi[2], hashi[3], hashi[4]);
    else
        sprintf(buf, "%08x%08x%08x%08x%08x", hashi[0], hashi[1], hashi[2], hashi[3], hashi[4]);

    return buf;
}

static void printHash(const unsigned char* hash, int initializer)
{
    char* buf;

    buf = showHash(hash, initializer);
    puts(buf);
    free(buf);
}


/* Hash of just the bodies masses, positions and velocities */
int hashBodies(unsigned char* hash, const Body* bodies, unsigned int nbody)
{
    SHA_CTX shaCtx;
    unsigned int i;
    const Body* b;
    struct
    {
        real mass;
        mwvector pos, vel;
    } hashableBody;

    if (nbody == 0)
        return warn1("Can't hash 0 bodies\n");

    if (!SHA1_Init(&shaCtx))
        return warn1("Initializing SHA1 context failed\n");

    for (i = 0; i < nbody; ++i)
    {
        b = &bodies[i];
        hashableBody.mass = Mass(b);
        hashableBody.pos  = Pos(b);
        hashableBody.vel  = Vel(b);

        if (!SHA1_Update(&shaCtx, &hashableBody, sizeof(hashableBody)))
            return warn1("Error updating SHA1 hash for body %u\n", i);
    }

    if (!SHA1_Final(hash, &shaCtx))
        return warn1("Error finalizing SHA1 hash\n");

    return 0;
}

unsigned char* getBodyHash(const NBodyState* st, unsigned int nbody)
{
    unsigned char* bodyHash;

    bodyHash = (unsigned char*) mwCalloc(SHA_DIGEST_LENGTH, sizeof(unsigned char));
    if (hashBodies(bodyHash, st->bodytab, nbody))
    {
        free(bodyHash);
        return NULL;
    }

    return bodyHash;
}

static Body* testPlummer_1_1()
{
    uint32_t seed = 99;
    unsigned int nbody = 100;
    real mass = 20;
    mwbool ignore = FALSE;
    mwvector rShift = ZERO_VECTOR;
    mwvector vShift = ZERO_VECTOR;
    real scaleRadius = 1.0;
    dsfmt_t prng;

    Body* bodies;

    dsfmt_init_gen_rand(&prng, seed);

    bodies = (Body*) mwMalloc(sizeof(Body) * nbody);

    if (generatePlummerC(bodies, &prng, nbody,
                         mass, ignore, rShift, vShift, scaleRadius))
    {
        warn("Error generating Plummer test\n");
        free(bodies);
        return NULL;
    }

    return bodies;
}

int main(int argc, const char* argv[])
{
    Body* b;

    b = testPlummer_1_1();
    if (!b)
        return warn1("Error creating test\n");

    unsigned char hash[SHA_DIGEST_LENGTH];


    //printBodies(b, 100);
    hashBodies(hash, b, 100);
    printHash(hash, 0);


    free(b);

    return 0;
}


