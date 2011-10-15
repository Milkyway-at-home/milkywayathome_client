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
    return mw_lua_checkglobalfunction(luaSt, "makeContext");
}

static int getNBodyPotentialFunc(lua_State* luaSt)
{
    return mw_lua_checkglobalfunction(luaSt, "makePotential");
}

static int getHistogramFunc(lua_State* luaSt)
{
    return mw_lua_checkglobalfunction(luaSt, "makeHistogram");
}

static int getBodiesFunc(lua_State* luaSt)
{
    return mw_lua_checkglobalfunction(luaSt, "makeBodies");
}

static int bindArgSeed(lua_State* luaSt, const NBodyFlags* nbf)
{
    lua_pushinteger(luaSt, nbf->seed);
    lua_setglobal(luaSt, "argSeed");

    return 0;
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

    bindArgSeed(luaSt, nbf);
    mwBindBOINCStatus(luaSt);

    script = mwReadFileResolved(nbf->inputFile);
    if (!script)
    {
        mwPerror("Opening Lua script '%s'", nbf->inputFile);
        lua_close(luaSt);
        return NULL;
    }

    if (dostringWithArgs(luaSt, script, nbf->forwardedArgs, nbf->numForwardedArgs))
    {
        mw_lua_pcall_warn(luaSt, "Error loading Lua script '%s'", nbf->inputFile);
        lua_close(luaSt);
        luaSt = NULL;
    }

    free(script);
    return luaSt;
}

static int evaluateContext(lua_State* luaSt, NBodyCtx* ctx)
{
    NBodyCtx* tmp;

    getNBodyCtxFunc(luaSt);
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

    if (!luaSt)
    {
        return LUA_NOREF;
    }

    getNBodyPotentialFunc(luaSt);
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
 * the NBodyState in the event of an error */
void nbEvalPotentialClosure(NBodyState* st, mwvector pos, mwvector* aOut)
{
  #ifdef _OPENMP
    const int tid = omp_get_thread_num();
  #else
    const int tid = 0;
  #endif

    int top;
    const mwvector* a;
    static const mwvector nanVector = mw_vec(NAN, NAN, NAN);
    lua_State* luaSt = st->potEvalStates[tid];

    /* Push closure */
    getLuaClosure(luaSt, &st->potEvalClosures[tid]);

    /* Push position argument */
    pushVector(luaSt, pos);

    /* Call closure */
    if (lua_pcall(luaSt, 1, 1, 0))
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
        *aOut = nanVector;
        return;
    }

    /* Retrieve acceleration. Use to* type function to avoid excessive
     * error message printing in general case. */
    top = lua_gettop(luaSt);
    a = toVector(luaSt, top);
    if (!a)
    {
        if (!st->potentialEvalError)
        {
          #ifdef _OPENMP
            #pragma omp critical
          #endif
            {
                mw_printf("Error in Lua potential function: ");
                mw_lua_typecheck(luaSt, top, LUA_TUSERDATA, MWVECTOR_TYPE);
                st->potentialEvalError = TRUE;
            }
        }

        lua_pop(luaSt, 1);
        *aOut = nanVector;
        return;
    }

    *aOut = *a;
    lua_pop(luaSt, 1);
}

static int evaluatePotential(lua_State* luaSt, NBodyCtx* ctx)
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

static int evaluateHistogram(lua_State* luaSt, HistogramParams* hp)
{
    HistogramParams* tmp;
    int rc = 0;

    getHistogramFunc(luaSt);
    if (lua_pcall(luaSt, 0, 1, 0))
    {
        mw_lua_pcall_warn(luaSt, "Error evaluating HistogramParams");
        return 1;
    }

    tmp = expectHistogramParams(luaSt, lua_gettop(luaSt));
    if (!tmp)
        rc = 1;
    else
        *hp = *tmp;
    lua_pop(luaSt, 1);
    return rc;
}

static Body* evaluateBodies(lua_State* luaSt, const NBodyCtx* ctx, int* n)
{
    int level, nResults;

    level = lua_gettop(luaSt);

    getBodiesFunc(luaSt);
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

static int evaluateInitialNBodyState(lua_State* luaSt, NBodyCtx* ctx, NBodyState* st, HistogramParams* hp)
{
    Body* bodies;
    int nbody;

    if (evaluateContext(luaSt, ctx))
        return 1;

    if (evaluatePotential(luaSt, ctx))
        return 1;

    if (evaluateHistogram(luaSt, hp))
        return 1;

    bodies = evaluateBodies(luaSt, ctx, &nbody);
    if (!bodies)
        return 1;

    setInitialNBodyState(st, ctx, bodies, nbody);

    return 0;
}

int nbSetup(NBodyCtx* ctx, NBodyState* st, HistogramParams* hp, const NBodyFlags* nbf)
{
    int rc;
    lua_State* luaSt;

    luaSt = nbOpenLuaStateWithScript(nbf);
    if (!luaSt)
        return 1;

    rc = evaluateInitialNBodyState(luaSt, ctx, st, hp);
    lua_close(luaSt);

    return rc;
}

