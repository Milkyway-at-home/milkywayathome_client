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

static char* showHash(const unsigned char* hash, int initializer)
{
    char* buf;
    const uint32_t* hashi = (const uint32_t*) hash;

    if (!hash)
        return NULL;

    assert(SHA_DIGEST_LENGTH / sizeof(uint32_t) == 5);

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

static int compareVectors(mwvector a, mwvector b)
{
    return X(a) < X(b) || Y(a) < Y(b) || Z(a) < Z(b);
}

/* Function for sorting bodies */
static int compareBodies(const void* _a, const void* _b)
{
    const Body* a = (const Body*) _a;
    const Body* b = (const Body*) _b;

    return Mass(a) < Mass(b) || compareVectors(Pos(a), Pos(b)) || compareVectors(Vel(a), Vel(b));
}

/* Sort the bodies. The actual order doesn't matter, it just needs to
 * be consistent when we hash. This is so when if we end up shifting
 * bodies around for the GPU, the tests will still work as
 * expected. */
static void sortBodies(Body* bodies, unsigned int nbody)
{
    qsort(bodies, nbody, sizeof(Body), compareBodies);
}


/* Hash of just the bodies masses, positions and velocities */
static int hashBodiesRaw(EVP_MD_CTX* hashCtx, unsigned char* hash, const Body* bodies, unsigned int nbody)
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

    if (!EVP_DigestFinal_ex(hashCtx, hash, &mdLen))
        return warn1("Error finalizing hash\n");

    assert(mdLen == SHA_DIGEST_LENGTH);

    if (!EVP_MD_CTX_cleanup(hashCtx))
        return warn1("Error cleaning up hash context\n");

    return 0;
}

int hashBodies(unsigned char* hash, const Body* bodies, unsigned int nbody)
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

int hashSortBodies(unsigned char* hash, Body* bodies, unsigned int nbody)
{
    sortBodies(bodies, nbody);
    return hashBodies(hash, bodies, nbody);
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

static int compareHash(const unsigned char* a, const unsigned char* b)
{
    return memcmp(a, b, sizeof(unsigned char) * SHA_DIGEST_LENGTH);
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
    Body* bs;

    bs = testPlummer_1_1();
    if (!bs)
        return warn1("Error creating test\n");

    unsigned char hash[SHA_DIGEST_LENGTH];


    //printBodies(bs, 100);
    hashSortBodies(hash, bs, 100);
    printHash(hash, 0);

    free(bs);

    return 0;
}


