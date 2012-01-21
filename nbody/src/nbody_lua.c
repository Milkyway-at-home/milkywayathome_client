/*
 *  Copyright (c) 2011 Matthew Arsenault
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

#include <lua.h>
#include <lualib.h>
#include <lauxlib.h>

#include "nbody_types.h"
#include "milkyway_util.h"

#include "nbody_util.h"
#include "nbody_lua.h"
#include "nbody_lua_types.h"
#include "nbody_lua_models.h"
#include "nbody_lua_util.h"
#include "nbody_lua_misc.h"
#include "milkyway_lua.h"
#include "nbody_check_params.h"
#include "nbody_defaults.h"

static int getNBodyCtxFunc(lua_State* luaSt)
{
    return mw_lua_getglobalfunction(luaSt, "makeContext");
}

static int getNBodyPotentialFunc(lua_State* luaSt)
{
    return mw_lua_getglobalfunction(luaSt, "makePotential");
}

static int getHistogramFunc(lua_State* luaSt)
{
    return mw_lua_getglobalfunction(luaSt, "makeHistogram");
}

static int getBodiesFunc(lua_State* luaSt)
{
    return mw_lua_getglobalfunction(luaSt, "makeBodies");
}

static int bindArgSeed(lua_State* luaSt, const NBodyFlags* nbf)
{
    lua_pushinteger(luaSt, nbf->seed);
    lua_setglobal(luaSt, "argSeed");

    return 0;
}

static int bindVersionNumber(lua_State* luaSt)
{
    lua_pushinteger(luaSt, NBODY_VERSION_MAJOR);
    lua_setglobal(luaSt, "NBODY_VERSION_MAJOR");

    lua_pushinteger(luaSt, NBODY_VERSION_MINOR);
    lua_setglobal(luaSt, "NBODY_VERSION_MINOR");

    lua_pushstring(luaSt, NBODY_VERSION);
    lua_setglobal(luaSt, "NBODY_VERSION");

    return 0;
}

static int nbCheckMinVersionRequired(lua_State* luaSt)
{
    const char* version;
    int major, minor;

    lua_getglobal(luaSt, "nbodyMinVersion");
    version = luaL_optstring(luaSt, -1, "0.0");

    if (sscanf(version, "%d.%d", &major, &minor) != 2)
    {
        mw_printf("Error reading minimum required version\n");
        lua_pop(luaSt, 1);
        return 0;
    }

    lua_pop(luaSt, 1);

    if ((NBODY_VERSION_MAJOR < major) || (NBODY_VERSION_MAJOR == major && NBODY_VERSION_MINOR < minor))
    {
        mw_printf("Application version too old. Workunit requires version %d.%d, but this is %d.%d\n",
                  major, minor,
                  NBODY_VERSION_MAJOR, NBODY_VERSION_MINOR
            );
        return FALSE;
    }

    return TRUE;
}

/* Open a lua_State and load the stuff we define, but do not run anything */
lua_State* nbLuaOpen(mwbool debug)
{
    lua_State* luaSt;

    luaSt = mw_lua_newstate();
    if (!luaSt)
    {
        mw_printf("Failed to get Lua state\n");
        return NULL;
    }

    mw_lua_openlibs(luaSt, debug);

    registerNBodyTypes(luaSt);
    mwRegisterTypes(luaSt);

    registerPredefinedModelGenerators(luaSt);
    registerSphericalKinds(luaSt);
    registerDiskKinds(luaSt);
    registerHaloKinds(luaSt);

    registerModelFunctions(luaSt);
    registerUtilityFunctions(luaSt);
    nbRegisterUtilityFunctions(luaSt);

    return luaSt;
}

/* Open a lua_State, bind run information such as server arguments and
 * BOINC status, and evaluate input script. */
lua_State* nbOpenLuaStateWithScript(const NBodyFlags* nbf)
{
    char* script;
    lua_State* luaSt;
    int execFailed;

    luaSt = nbLuaOpen(nbf->debugLuaLibs);
    if (!luaSt)
        return NULL;

    bindVersionNumber(luaSt);
    bindArgSeed(luaSt, nbf);
    mwBindBOINCStatus(luaSt);

    script = mwReadFileResolved(nbf->inputFile);
    if (!script)
    {
        mwPerror("Opening Lua script '%s'", nbf->inputFile);
        lua_close(luaSt);
        return NULL;
    }

    execFailed = dostringWithArgs(luaSt, script, nbf->forwardedArgs, nbf->numForwardedArgs);
    free(script);
    if (execFailed)
    {
        mw_lua_perror(luaSt, "Error loading Lua script '%s'", nbf->inputFile);
        lua_close(luaSt);
        return NULL;
    }

    if (!nbCheckMinVersionRequired(luaSt))
    {
        lua_close(luaSt);
        return NULL;
    }

    return luaSt;
}

static int nbEvaluateContext(lua_State* luaSt, NBodyCtx* ctx)
{
    NBodyCtx* tmp;

    if (!luaSt || getNBodyCtxFunc(luaSt))
    {
        return 1;
    }

    if (lua_pcall(luaSt, 0, 1, 0))
    {
        mw_lua_perror(luaSt, "Error evaluating NBodyCtx");
        return 1;
    }

    tmp = expectNBodyCtx(luaSt, lua_gettop(luaSt));
    if (!tmp)
    {
        mw_lua_perror(luaSt, "Invalid return from makeContext()");
        return 1;
    }

    *ctx = *tmp;
    lua_pop(luaSt, 1);

    return checkNBodyCtxConstants(ctx);
}

/* Evaluate the potential function and expect a closure.  Not used on
   the first pass. Only used when using a Lua closure potential
 */
static int nbGetPotentialClosure(lua_State* luaSt)
{
    int top;
    int closure;

    if (!luaSt || getNBodyPotentialFunc(luaSt))
    {
        return LUA_NOREF;
    }

    if (lua_pcall(luaSt, 0, 1, 0))
    {
        mw_lua_perror(luaSt, "Error evaluating potential closure");
        return LUA_NOREF;
    }

    top = lua_gettop(luaSt);
    closure = mw_lua_checkluaclosure(luaSt, top);
    if (closure == LUA_NOREF)
    {
        mw_printf("Failed to get potential Lua closure\n");
    }

    return closure;
}

int nbOpenPotentialEvalStatePerThread(NBodyState* st, const NBodyFlags* nbf)
{
    int i;
    int* closures;
    lua_State** states;
    const int maxThreads = nbGetMaxThreads();

    states = mwCalloc(maxThreads, sizeof(lua_State*));
    closures = mwCalloc(maxThreads, sizeof(int));

    /* CHECKME: Is it OK to open all states in the master thread first? */
    for (i = 0; i < maxThreads; ++i)
    {
        /* FIXME: Is there a better way to copy a lua_State? */
        states[i] = nbOpenLuaStateWithScript(nbf);
        closures[i] = nbGetPotentialClosure(states[i]);
        if (closures[i] == LUA_NOREF)
        {
            free(states);
            free(closures);
            return 1;
        }
    }

    st->potEvalStates = states;
    st->potEvalClosures = closures;

    return 0;
}

/* Evaluate potential from a Lua closure. Will set the error flag on
 * the NBodyState in the event of an error.
 *
 * The closure used must of type number, number, number -> number,
 * number, number. (i.e. takes 3 numbers (x, y, z) and returns 3 numbers (a_x, a_y, a_z))
 */
void nbEvalPotentialClosure(NBodyState* st, mwvector pos, mwvector* aOut)
{
  #ifdef _OPENMP
    const int tid = omp_get_thread_num();
  #else
    const int tid = 0;
  #endif

    int top;
    mwvector a;
    static const mwvector badVector = mw_vec(REAL_MAX, REAL_MAX, REAL_MAX);
    lua_State* luaSt = st->potEvalStates[tid];

    /* Push closure */
    getLuaClosure(luaSt, &st->potEvalClosures[tid]);

    /* Push position arguments */
    lua_pushnumber(luaSt, X(pos));
    lua_pushnumber(luaSt, Y(pos));
    lua_pushnumber(luaSt, Z(pos));

    /* Call closure */
    if (lua_pcall(luaSt, 3, 3, 0))
    {
        /* Avoid spewing the same error billions of times */
        if (!st->potentialEvalError)
        {
            /* FIXME: This isn't really correct. We really need a lock. The worst that should happen*/
          #ifdef _OPENMP
            #pragma omp critical
          #endif
            {
                mw_lua_perror(luaSt, "Error evaluating potential closure");
                st->potentialEvalError = TRUE;
            }
        }

        /* Make sure we break everything */
        *aOut = badVector;
        return;
    }


    /* Retrieve acceleration. Use to* type function to avoid excessive
     * error message printing in general case. */
    top = lua_gettop(luaSt);
    if (!lua_isnumber(luaSt, top) || !lua_isnumber(luaSt, top - 1) || !lua_isnumber(luaSt, top - 2))
    {
        if (!st->potentialEvalError)
        {
          #ifdef _OPENMP
            #pragma omp critical
          #endif
            {
                mw_printf("Error in Lua potential function: "
                          "Expected number, number, number. Got %s, %s, %s\n",
                          luaL_typename(luaSt, top - 2),
                          luaL_typename(luaSt, top - 1),
                          luaL_typename(luaSt, top));
                st->potentialEvalError = TRUE;
            }
        }

        *aOut = badVector;
        return;
    }

    Z(a) = lua_tonumber(luaSt, top);
    Y(a) = lua_tonumber(luaSt, top - 1);
    X(a) = lua_tonumber(luaSt, top - 2);

    *aOut = a;
    lua_pop(luaSt, 3);
}

static int nbEvaluatePotential(lua_State* luaSt, NBodyCtx* ctx)
{
    int top;

    if (getNBodyPotentialFunc(luaSt))
    {
        return 1;
    }

    if (lua_pcall(luaSt, 0, 1, 0))
    {
        mw_lua_perror(luaSt, "Error evaluating makePotential()");
        return 1;
    }

    top = lua_gettop(luaSt);
    if (nbGetPotentialTyped(luaSt, ctx, top, "Invalid return from makePotential()"))
    {
        lua_pop(luaSt, 1);
        return 1;
    }

    lua_pop(luaSt, 1);
    return 0;
}

NBodyLikelihoodMethod nbEvaluateLikelihoodMethod(lua_State* luaSt)
{
    int top;

    static const MWEnumAssociation methodOptions[] =
        {
            { "EMD",             NBODY_EMD              },
            { "Original",        NBODY_ORIG_CHISQ       },
            { "AltOriginal",     NBODY_ORIG_ALT         },
            { "ChisqAlt",        NBODY_CHISQ_ALT        },
            { "Poisson",         NBODY_POISSON          },
            { "Kolmogorov",      NBODY_KOLMOGOROV       },
            { "KullbackLeibler", NBODY_KULLBACK_LEIBLER },
            { "Saha",            NBODY_SAHA             },
            END_MW_ENUM_ASSOCIATION
        };

    lua_getglobal(luaSt, "nbodyLikelihoodMethod");
    top = lua_gettop(luaSt);
    if (lua_isnoneornil(luaSt, top))
    {
        return DEFAULT_LIKELIHOOD_METHOD;
    }

    return (NBodyLikelihoodMethod) expectEnum(luaSt, methodOptions, top);
}

int nbEvaluateHistogramParams(lua_State* luaSt, HistogramParams* hp)
{
    HistogramParams* tmp;

    if (!luaSt || getHistogramFunc(luaSt))
    {
        return 1;
    }

    if (lua_pcall(luaSt, 0, 1, 0))
    {
        mw_lua_perror(luaSt, "Error evaluating makeHistogram()");
        return 1;
    }

    tmp = expectHistogramParams(luaSt, lua_gettop(luaSt));
    if (!tmp)
    {
        mw_lua_perror(luaSt, "Invalid return from makeHistogram()");
        return 1;
    }

    *hp = *tmp;
    return 0;
}

/* Test that the histogram params in the input from the file are OK
 * for file verification */
int nbHistogramParamsCheck(const NBodyFlags* nbf, HistogramParams* hp)
{
    lua_State* luaSt;
    int rc;
    NBodyLikelihoodMethod method;

    luaSt = nbOpenLuaStateWithScript(nbf);
    if (!luaSt)
    {
        return 1;
    }

    rc = nbEvaluateHistogramParams(luaSt, hp);
    method = nbEvaluateLikelihoodMethod(luaSt);
    lua_close(luaSt);

    return (rc || method == NBODY_INVALID_METHOD);
}

static Body* nbEvaluateBodies(lua_State* luaSt, const NBodyCtx* ctx, int* n)
{
    int level, nResults;

    level = lua_gettop(luaSt);
    if (getBodiesFunc(luaSt))
    {
        return NULL;
    }

    pushNBodyCtx(luaSt, ctx);

    if (ctx->potentialType == EXTERNAL_POTENTIAL_DEFAULT)
        pushPotential(luaSt, &ctx->pot);
    else
        lua_pushnil(luaSt);

    if (lua_pcall(luaSt, 2, LUA_MULTRET, 0))
    {
        mw_lua_perror(luaSt, "Error evaluating makeBodies()");
        return NULL;
    }

    nResults = lua_gettop(luaSt) - level;

    return readModels(luaSt, nResults, n);
}

static int nbEvaluateInitialNBodyState(lua_State* luaSt, NBodyCtx* ctx, NBodyState* st)
{
    Body* bodies;
    int nbody;

    if (nbEvaluateContext(luaSt, ctx))
        return 1;

    if (nbEvaluatePotential(luaSt, ctx))
        return 1;

    bodies = nbEvaluateBodies(luaSt, ctx, &nbody);
    if (!bodies)
        return 1;

    setInitialNBodyState(st, ctx, bodies, nbody);

    return 0;
}

int nbSetup(NBodyCtx* ctx, NBodyState* st, const NBodyFlags* nbf)
{
    int rc;
    lua_State* luaSt;

    luaSt = nbOpenLuaStateWithScript(nbf);
    if (!luaSt)
        return 1;

    rc = nbEvaluateInitialNBodyState(luaSt, ctx, st);
    lua_close(luaSt);

    return rc;
}

