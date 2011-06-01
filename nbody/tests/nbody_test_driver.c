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
#include "milkyway_lua.h"
#include "nbody_plummer.h"
#include "nbody_defaults.h"
#include "nbody_tree.h"


/* things in NBodyCtx which influence individual steps that aren't the potential. */
typedef struct
{
    real theta;
    real treeRSize;
    criterion_t criterion;

    mwbool useQuad;
    mwbool allowIncest;
} NBodyCtxTest;

#define EMPTY_NBODYCTXTEST { 0.0, 0.0, InvalidCriterion, FALSE, FALSE }

typedef struct
{
    const char* potentialName;
    const char* modelName;
    unsigned int nSteps;
    unsigned int nbody;
    uint32_t seed;
    mwbool doublePrec;
    NBodyCtxTest ctx;
} NBodyTest;

#define EMPTY_NBODYTEST { NULL, NULL, 0, 0, 0, DOUBLEPREC, EMPTY_NBODYCTXTEST }

typedef union
{
    uint32_t mdi[5];
    unsigned char md[SHA_DIGEST_LENGTH];
} MWHash;

#define EMPTY_BODY_HASH { .mdi = { 0x0, 0x0, 0x0, 0x0, 0x0 } }


static void showHash(char* buf, const MWHash* hash)
{
    sprintf(buf, "%08x%08x%08x%08x%08x", hash->mdi[0], hash->mdi[1], hash->mdi[2], hash->mdi[3], hash->mdi[4]);
}

static int hashValueFromType(lua_State* luaSt, EVP_MD_CTX* hashCtx, int type, int idx)
{
    int rc = 1;
    real n;
    int b;
    const char* str;

    switch (type)
    {
        case LUA_TNUMBER:
            n = lua_tonumber(luaSt, idx);
            rc = EVP_DigestUpdate(hashCtx, &n, sizeof(n));
            break;

        case LUA_TBOOLEAN:
            b = lua_toboolean(luaSt, idx);
            rc = EVP_DigestUpdate(hashCtx, &b, sizeof(b));
            break;

        case LUA_TSTRING:
            str = lua_tostring(luaSt, idx);
            rc = EVP_DigestUpdate(hashCtx, str, strlen(str));
            break;

        case LUA_TUSERDATA:
        case LUA_TTABLE:
        case LUA_TFUNCTION:
        case LUA_TLIGHTUSERDATA:
        case LUA_TTHREAD:
        case LUA_TNIL:
        default:
            mw_panic("Unhandled type %s (%d)\n", luaL_typename(luaSt, type), type);
    }

    if (rc == 0)
        return warn1("Error updating hash of type\n");

    return 0;
}

static int checkNBodyTestTable(lua_State* luaSt, int idx, NBodyTest* testOut)
{
    static NBodyTest test = EMPTY_NBODYTEST;
    static const char* criterionName = NULL;
    static real seedf = 0.0;
    static real nStepsf = 0.0;
    static real nbodyf = 0.0;
    static mwbool failed = FALSE;
    static const char* resultHash = NULL;
    static const char* resultName = NULL;
    static const MWNamedArg argTable[] =
        {
            { "potential",   LUA_TSTRING,  NULL, TRUE,  &test.potentialName   },
            { "model",       LUA_TSTRING,  NULL, TRUE,  &test.modelName       },
            { "nbody",       LUA_TNUMBER,  NULL, TRUE,  &nbodyf               },
            { "nSteps",      LUA_TNUMBER,  NULL, TRUE,  &nStepsf              },
            { "seed",        LUA_TNUMBER,  NULL, TRUE,  &seedf                },
            { "theta",       LUA_TNUMBER,  NULL, TRUE,  &test.ctx.theta       },
            { "treeRSize",   LUA_TNUMBER,  NULL, TRUE,  &test.ctx.treeRSize   },
            { "criterion",   LUA_TSTRING,  NULL, TRUE,  &criterionName        },
            { "useQuad",     LUA_TBOOLEAN, NULL, TRUE,  &test.ctx.useQuad     },
            { "allowIncest", LUA_TBOOLEAN, NULL, TRUE,  &test.ctx.allowIncest },

            { "doublePrec",  LUA_TBOOLEAN, NULL, FALSE, &test.doublePrec      },

            /* Unused in hash; these ones may or may not exist, just don't error if there */
            { "result",     LUA_TSTRING,   NULL,  FALSE, &resultHash          },
            { "err",        LUA_TSTRING,   NULL,  FALSE, &resultName          },
            { "failed",     LUA_TBOOLEAN,  NULL,  FALSE, &failed              },
            END_MW_NAMED_ARG
        };

    /* Oh look, we conveniently already have the table constructed
     * that we want. Just run type checking on it. */
    handleNamedArgumentTable(luaSt, argTable, idx);

    test.seed = (uint32_t) seedf;
    test.nSteps = (unsigned int) nStepsf;
    test.nbody = (unsigned int) nbodyf;

    test.ctx.criterion = readCriterion(luaSt, criterionName);
    lua_pushvalue(luaSt, lua_gettop(luaSt));

    if (testOut)
        *testOut = test;

    return 1;
}

static int hashNBodyTestCore(EVP_MD_CTX* hashCtx, MWHash* hash, const NBodyTest* t)
{
    int rc = 0;

    if (!EVP_DigestInit_ex(hashCtx, EVP_sha1(), NULL))
        return warn1("Initializing hash digest failed\n");

    rc |= !EVP_DigestUpdate(hashCtx, &t->ctx.theta,       sizeof(t->ctx.theta));
    rc |= !EVP_DigestUpdate(hashCtx, &t->ctx.treeRSize,   sizeof(t->ctx.treeRSize));
    rc |= !EVP_DigestUpdate(hashCtx, &t->ctx.criterion,   sizeof(t->ctx.criterion));
    rc |= !EVP_DigestUpdate(hashCtx, &t->ctx.useQuad,     sizeof(t->ctx.useQuad));
    rc |= !EVP_DigestUpdate(hashCtx, &t->ctx.allowIncest, sizeof(t->ctx.allowIncest));

    rc |= !EVP_DigestUpdate(hashCtx, &t->seed,            sizeof(t->seed));
    rc |= !EVP_DigestUpdate(hashCtx, &t->nSteps,          sizeof(t->nSteps));
    rc |= !EVP_DigestUpdate(hashCtx, &t->nbody,           sizeof(t->nbody));
    rc |= !EVP_DigestUpdate(hashCtx, t->modelName,        strlen(t->modelName));
    rc |= !EVP_DigestUpdate(hashCtx, t->potentialName,    strlen(t->potentialName));

    rc |= !EVP_DigestUpdate(hashCtx, &t->doublePrec,      sizeof(t->doublePrec));

    if (rc)
        return warn1("Error updating hashing for NBodyTest\n");

    if (!EVP_DigestFinal_ex(hashCtx, hash->md, NULL))
        return warn1("Error finalizing hash for NBodyTest\n");

    if (!EVP_MD_CTX_cleanup(hashCtx))
        return warn1("Error cleaning up hash context for NBodyCtxTest\n");

    return 0;
}

int hashNBodyTest(MWHash* hash, NBodyTest* test)
{
    EVP_MD_CTX hashCtx;
    int failed = 0;

    EVP_MD_CTX_init(&hashCtx);
    failed = hashNBodyTestCore(&hashCtx, hash, test);
    if (!EVP_MD_CTX_cleanup(&hashCtx))
        return warn1("Error cleaning up hash context\n");

    return failed;
}

/* Return the hash of NBodyCtxTest table to Lua */
static int hashNBodyTestTable(lua_State* luaSt)
{
    NBodyTest test;
    MWHash hash;
    char buf[SHA_DIGEST_LENGTH];

    checkNBodyTestTable(luaSt, lua_gettop(luaSt), &test);
    hashNBodyTest(&hash, &test);

    showHash(buf, &hash);
    lua_pushstring(luaSt, buf);

    return 1;
}

static int statusIsFatal(lua_State* luaSt)
{
    NBodyStatus rc = readNBodyStatus(luaSt, lua_tostring(luaSt, 1));
    lua_pushboolean(luaSt, nbodyStatusIsFatal(rc));
    return 1;
}

static void registerNBodyTestFunctions(lua_State* luaSt)
{
    lua_register(luaSt, "hashNBodyTest", hashNBodyTestTable);
    lua_register(luaSt, "statusIsFatal", statusIsFatal);
}

/* Hash of just the bodies masses, positions and velocities */
static int hashBodiesCore(EVP_MD_CTX* hashCtx, MWHash* hash, const Body* bodies, unsigned int nbody)
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

    return 0;
}

static void nbodyTestInit()
{
    OpenSSL_add_all_digests();
}

static void nbodyTestCleanup()
{
    EVP_cleanup();
}

int hashBodies(MWHash* hash, const Body* bodies, unsigned int nbody)
{
    EVP_MD_CTX hashCtx;
    int failed = 0;

    EVP_MD_CTX_init(&hashCtx);
    failed = hashBodiesCore(&hashCtx, hash, bodies, nbody);
    if (!EVP_MD_CTX_cleanup(&hashCtx))
        return warn1("Error cleaning up hash context\n");

    return failed;
}

int hashSortBodies(MWHash* hash, Body* bodies, unsigned int nbody)
{
    sortBodies(bodies, nbody);
    return hashBodies(hash, bodies, nbody);
}

MWHash* getMWHash(const NBodyState* st, unsigned int nbody)
{
    MWHash* bodyHash;

    bodyHash = (MWHash*) mwCalloc(1, sizeof(MWHash));
    if (hashBodies(bodyHash, st->bodytab, nbody))
    {
        free(bodyHash);
        return NULL;
    }

    return bodyHash;
}

static int compareHash(const MWHash* a, const MWHash* b)
{
    return memcmp(a, b, sizeof(MWHash));
}


static int luaHashOrHashSortBodies(lua_State* luaSt, int sort)
{
    NBodyState* st;
    MWHash hash;
    char hashBuf[2 * sizeof(MWHash) + 1] = "";

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


/* Create a context with everything unset, useful for testing but
 * undesirable for actual work. */
static int createTestNBodyCtx(lua_State* luaSt)
{
    return pushNBodyCtx(luaSt, &defaultNBodyCtx);
}

static void registerNBodyCtxTestMethods(lua_State* luaSt)
{
    static const luaL_reg testMethodsNBodyCtx[] =
        {
            { "createTestCtx", createTestNBodyCtx },
            { NULL, NULL }
        };

    luaL_register(luaSt, NBODYCTX_TYPE, testMethodsNBodyCtx);
}

static int runNBodyTest(const char* file, const char** args, unsigned int nArgs)
{
    int rc;
    lua_State* luaSt;

    luaSt = nbodyLuaOpen(TRUE);
    if (!luaSt)
        return 1;

    /* Register special functions used by tests but not in real things
     * for various reasons, such as avoiding depending on openssl and
     * to not include useless / and or less safe versions of
     * functions. */
    registerNBodyState(luaSt);
    installHashFunctions(luaSt);
    registerNBodyTestFunctions(luaSt);
    registerNBodyCtxTestMethods(luaSt);
    //registerFindRCrit(luaSt);

    rc = dofileWithArgs(luaSt, file, args, nArgs);
    if (rc)
        mw_lua_pcall_warn(luaSt, "Error evaluating script '%s'\n", file);

    lua_close(luaSt);

    return rc;
}

int main(int argc, const char* argv[])
{
    int rc;
    const char* testScript;

    testScript = argc > 1 ? argv[1] : NULL;

    nbodyTestInit();

    rc = runNBodyTest(testScript, &argv[2], argc - 2);

    nbodyTestCleanup();

    return rc;
}


