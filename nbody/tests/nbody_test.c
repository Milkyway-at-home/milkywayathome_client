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

#include <openssl/evp.h>
#include <openssl/sha.h>

#include "milkyway_util.h"
#include "nbody_priv.h"
#include "nbody_lua_types.h"
#include "milkyway_lua_marshal.h"
#include "nbody_plummer.h"


/* things in NBodyCtx which influence individual steps that aren't the potential. */
typedef struct
{
    real theta;
    real eps2;
    real treeRSize;
    criterion_t criterion;

    mwbool useQuad;
    mwbool allowIncest;
} NBodyCtxRaw;

#define EMPTY_NBODYCTXRAW { 0.0, 0.0, 0.0, InvalidCriterion, FALSE, FALSE }

typedef union
{
    uint32_t mdi[5];
    unsigned char md[SHA_DIGEST_LENGTH];
} BodyHash;

#define EMPTY_BODY_HASH { .mdi = { 0x0, 0x0, 0x0, 0x0, 0x0 } }


typedef struct
{
    const char* modelName;       /* Identifier for initial distribution used */
    const Potential* potential;
    NBodyCtxRaw ctx;             /* Simulation configuration */
    uint32_t seed;               /* Seed used in for each result */
    unsigned int nbody;
    unsigned int numberSteps;    /* Steps taken for each result */
    BodyHash hash;               /* Final hash of bodies */
} NBodyResult;

#define END_NBODYRESULT { NULL, NULL, EMPTY_NBODYCTXRAW, 0, 0, 0, EMPTY_BODY_HASH }

typedef struct
{
    const NBodyResult* results;
} NBodyResultSet;

static void setNBodyCtxFromNBodyCtxRaw(NBodyCtx* ctx, const NBodyCtxRaw* ctxRaw)
{
    /* Tested, varied fields */
    ctx->theta       = ctxRaw->theta;                // 5
    ctx->eps2        = ctxRaw->eps2;                 // 5
    ctx->treeRSize   = ctxRaw->treeRSize;            // 5
    ctx->criterion   = ctxRaw->criterion;            // 4
    ctx->useQuad     = ctxRaw->useQuad;              // 2
    ctx->allowIncest = ctxRaw->allowIncest;          // 2
}

static void generateTestSet()
{
    const NBodyResultSet smallSet[] =
        {
            //{ .seed = 43, .nbody = 100, .numberSteps = 10,

        };

#if 0
    for (i = 0; i < numberSeeds; ++i)
    {
        for (j = 0; j < numberBodyTests; ++j)
        {
            for (k = 0; k < numberStepTests; ++k)
            {

            }
        }
    }
#endif

}


static void showHash(char* buf, const BodyHash* hash)
{
    sprintf(buf, "%08x%08x%08x%08x%08x", hash->mdi[0], hash->mdi[1], hash->mdi[2], hash->mdi[3], hash->mdi[4]);

}

/* Hash of just the bodies masses, positions and velocities */
static int hashBodiesRaw(EVP_MD_CTX* hashCtx, BodyHash* hash, const Body* bodies, unsigned int nbody)
{
    unsigned int i, mdLen;
    const Body* b;
    struct
    {
        mwvector pos, vel;
        real mass;
        body_t type;
        /* Padding happens */
    } hashableBody;

    if (nbody == 0)
        return warn1("Can't hash 0 bodies\n");

    if (!EVP_DigestInit_ex(hashCtx, EVP_sha1(), NULL))
        return warn1("Initializing hash digest failed\n");

    /* Prevent random garbage from getting hashed. The struct will be
     * padded and won't be the same size as 2 * sizeof(mwvector) +
     * sizeof(real) so bad things happen when hashing sizeof(hashableBody) */
    memset(&hashableBody, 0, sizeof(hashableBody));

    for (i = 0; i < nbody; ++i)
    {
        b = &bodies[i];

        hashableBody.pos  = Pos(b);
        hashableBody.vel  = Vel(b);
        hashableBody.mass = Mass(b);
        hashableBody.type = Type(b);

        if (!EVP_DigestUpdate(hashCtx, &hashableBody, sizeof(hashableBody)))
            return warn1("Error updating hash for body %u\n", i);
    }

    if (!EVP_DigestFinal_ex(hashCtx, hash->md, &mdLen))
        return warn1("Error finalizing hash\n");

    assert(mdLen == SHA_DIGEST_LENGTH);

    if (!EVP_MD_CTX_cleanup(hashCtx))
        return warn1("Error cleaning up hash context\n");

    return 0;
}

int hashBodies(BodyHash* hash, const Body* bodies, unsigned int nbody)
{
    EVP_MD_CTX hashCtx;
    int failed = 0;

    OpenSSL_add_all_digests();
    EVP_MD_CTX_init(&hashCtx);

    failed = hashBodiesRaw(&hashCtx, hash, bodies, nbody);

    if (!EVP_MD_CTX_cleanup(&hashCtx))
        failed |= warn1("Error cleaning up hash context\n");
    EVP_cleanup();

    return failed;
}

int hashSortBodies(BodyHash* hash, Body* bodies, unsigned int nbody)
{
    sortBodies(bodies, nbody);
    return hashBodies(hash, bodies, nbody);
}

BodyHash* getBodyHash(const NBodyState* st, unsigned int nbody)
{
    BodyHash* bodyHash;

    bodyHash = (BodyHash*) mwCalloc(1, sizeof(BodyHash));
    if (hashBodies(bodyHash, st->bodytab, nbody))
    {
        free(bodyHash);
        return NULL;
    }

    return bodyHash;
}

static int compareHash(const BodyHash* a, const BodyHash* b)
{
    return memcmp(a, b, sizeof(BodyHash));
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


static int luaHashOrHashSortBodies(lua_State* luaSt, int sort)
{
    NBodyState* st;
    BodyHash hash;
    char hashBuf[2 * sizeof(BodyHash) + 1] = "";

    if (lua_gettop(luaSt) != 1)
        luaL_argerror(luaSt, 1, "Expected 1 argument");

    st = checkNBodyState(luaSt, 1);

    if (sort)
        sortBodies(st->bodytab, st->nbody);

    hashBodies(&hash, st->bodytab, st->nbody);

    showHash(hashBuf, &hash);
    lua_pushstring(luaSt, hashBuf);

    return 1;
}

static int luaHashBodies(lua_State* luaSt)
{
    return luaHashOrHashSortBodies(luaSt, 0);
}

static int luaHashSortBodies(lua_State* luaSt)
{
    return luaHashOrHashSortBodies(luaSt, 1);
}

static int installHashFunctions(lua_State* luaSt)
{
    static const luaL_reg hashMethods[] =
        {
            { "hashBodies",     luaHashBodies     },
            { "hashSortBodies", luaHashSortBodies },
            { NULL, NULL }
        };

    luaL_register(luaSt, NBODYSTATE_TYPE, hashMethods);
    lua_pop(luaSt, 1);

    return 0;
}


static void testState()
{
    lua_State* luaSt;

    luaSt = nbodyLuaOpen();
    if (!luaSt)
        return;

    registerNBodyState(luaSt);
    installHashFunctions(luaSt);

    if (luaL_dofile(luaSt, "teststate.lua"))
        mw_lua_pcall_warn(luaSt, "Error evaluating script");

    lua_close(luaSt);
}

int main(int argc, const char* argv[])
{
    Body* bs;

    bs = testPlummer_1_1();
    if (!bs)
        return warn1("Error creating test\n");

    BodyHash hash;


    //printBodies(bs, 100);
    hashSortBodies(&hash, bs, 100);

    free(bs);

    testState();

    return 0;
}


