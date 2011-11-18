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
#include "milkyway_lua.h"
#include "nbody_check_params.h"

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

    luaSt = lua_open();
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

    registerModelUtilityFunctions(luaSt);
    registerUtilityFunctions(luaSt);

    return luaSt;
}

/* Open a lua_State, bind run information such as server arguments and
 * BOINC status, and evaluate input script. */
lua_State* nbOpenLuaStateWithScript(const NBodyFlags* nbf)
{
    char* script;
    lua_State* luaSt;

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

    if (   dostringWithArgs(luaSt, script, nbf->forwardedArgs, nbf->numForwardedArgs)
        || !nbCheckMinVersionRequired(luaSt))
    {
        mw_lua_pcall_warn(luaSt, "Error loading Lua script '%s'", nbf->inputFile);
        lua_close(luaSt);
        luaSt = NULL;
    }

    free(script);
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
        mw_lua_pcall_warn(luaSt, "Error evaluating NBodyCtx");
        return 1;
    }

    tmp = expectNBodyCtx(luaSt, lua_gettop(luaSt));
    if (!tmp)
    {
        lua_pop(luaSt, 1);
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
        mw_lua_pcall_warn(luaSt, "Error evaluating potential closure");
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
                mw_lua_pcall_warn(luaSt, "Error evaluating potential closure");
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
    Potential* tmp;

    getNBodyPotentialFunc(luaSt);
    if (lua_pcall(luaSt, 0, 1, 0))
    {
        mw_lua_pcall_warn(luaSt, "Error evaluating Potential");
        return 1;
    }

    top = lua_gettop(luaSt);
    if (lua_isnoneornil(luaSt, top))
    {
        ctx->potentialType = EXTERNAL_POTENTIAL_NONE;
        lua_pop(luaSt, 1);
        return 0;
    }
    else if (lua_isfunction(luaSt, top))
    {
        /* Set the potential type. We will be reevaluating the script
         * later on if we're using this, so don't bother getting the
         * closure now. */
        ctx->potentialType = EXTERNAL_POTENTIAL_CUSTOM_LUA;
        lua_pop(luaSt, 1);
        return 0;
    }
    else
    {
        tmp = expectPotential(luaSt, top);
        if (!tmp)
        {
            lua_pop(luaSt, 1);
            return 1;
        }

        ctx->potentialType = EXTERNAL_POTENTIAL_DEFAULT;
        ctx->pot = *tmp;

        lua_pop(luaSt, 1);
        return checkPotentialConstants(&ctx->pot);
    }

    mw_unreachable();
    return 0;
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
        mw_lua_pcall_warn(luaSt, "Error evaluating HistogramParams");
        return 1;
    }

    tmp = expectHistogramParams(luaSt, lua_gettop(luaSt));
    if (tmp)
    {
        *hp = *tmp;
    }

    lua_pop(luaSt, 1);
    return (tmp == NULL);
}

/* Test that the histogram params in the input from the file are OK
 * for file verification */
int nbHistogramParamsCheck(const NBodyFlags* nbf, HistogramParams* hp)
{
    lua_State* luaSt;
    int rc;

    luaSt = nbOpenLuaStateWithScript(nbf);
    if (!luaSt)
    {
        return 1;
    }

    rc = nbEvaluateHistogramParams(luaSt, hp);
    lua_close(luaSt);

    return rc;
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
        mw_lua_pcall_warn(luaSt, "Error evaluating bodies");
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

