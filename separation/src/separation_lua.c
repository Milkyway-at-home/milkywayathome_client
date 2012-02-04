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

#include "milkyway_lua.h"
#include "milkyway_util.h"
#include "separation_lua.h"

#define STREAMS_NAME "streams"
#define AREAS_NAME "area"
#define CONSTANTS_NAME "constants"
#define BACKGROUND_NAME "background"

typedef struct
{
    const char* name;
    int hasDefault;
    real* value;
} SeparationConstant;

/* Default values */
static real _sun_r0 = const_sun_r0;
static real _convolve = 120.0;
static real _wedge = 0.0;

static AstronomyParameters* _ap = NULL;
static Streams* _streams = NULL;
static BackgroundParameters* _bg = NULL;

static IntegralArea* _ias = NULL;
static int _nCut = 0;

static const SeparationConstant constants[] =
{
    { "sun_r0",   TRUE,  &_sun_r0   },
    { "convolve", TRUE,  &_convolve },
    { "wedge",    FALSE, &_wedge   },
    { NULL, FALSE, NULL }
};

/* Calculate total probability calculations for checkpointing */
static uint64_t findTotalCalcProbs(const IntegralArea* cuts, int nCut)
{
    int i;
    uint64_t total_calc_probs = 0;

    for (i = 0; i < nCut; ++i)
    {
        total_calc_probs += (uint64_t) cuts[i].mu_steps * cuts[i].nu_steps * cuts[i].r_steps;
    }

    return total_calc_probs;
}


static void setAPConstants(AstronomyParameters* ap)
{
    ap->number_streams = _streams->number_streams;
    ap->number_integrals = _nCut;

    ap->convolve = (int) _convolve;
    ap->sun_r0 = _sun_r0;

    ap->wedge = (int) _wedge;

    ap->total_calc_probs = (real) findTotalCalcProbs(_ias, _nCut);
}


static int pushConstants(lua_State* luaSt)
{
    const SeparationConstant* p = constants;

    while (p->name)
    {
        if (p->hasDefault)
        {
            lua_pushnumber(luaSt, *p->value);
        }
        else
        {
            lua_pushnil(luaSt);
        }

        lua_setglobal(luaSt, p->name);
        ++p;
    }

    return 0;
}

static int tryEvaluateScript(lua_State* luaSt, const char* script, const SeparationFlags* sf)
{
    if (script[0] == '\0')
    {
        mw_printf("Parameter file '%s' is empty\n", sf->ap_file);
        return 1;
    }

    if (dostringWithArgs(luaSt, script, sf->forwardedArgs, sf->nForwardedArgs))
    {
        mw_lua_perror(luaSt, "Error loading Lua script '%s'", sf->ap_file);
        return 1;
    }

    return 0;
}

/* Open a lua_State, bind run information such as server arguments and
 * BOINC status, and evaluate input script. */
static lua_State* separationOpenLuaStateWithScript(const SeparationFlags* sf)
{
    int failed;
    char* script;
    lua_State* luaSt;

    luaSt = separationLuaOpen(FALSE);
    if (!luaSt)
        return NULL;

    script = mwReadFileResolved(sf->ap_file);
    if (!script)
    {
        mwPerror("Opening Lua script '%s'", sf->ap_file);
        lua_close(luaSt);
        return NULL;
    }

    failed = tryEvaluateScript(luaSt, script, sf);
    free(script);

    if (failed)
    {
        lua_close(luaSt);
        return NULL;
    }

    return luaSt;
}

static int readIntegralArea(lua_State* luaSt, IntegralArea* iaOut, int table)
{
    uint64_t r, mu, nu;
    static IntegralArea ia;
    static real nuStepsf, muStepsf, rStepsf;
    static const MWNamedArg iaArgTable[] =
        {
            { "nu_min",   LUA_TNUMBER, NULL, TRUE, &ia.nu_min },
            { "nu_max",   LUA_TNUMBER, NULL, TRUE, &ia.nu_max },
            { "nu_steps", LUA_TNUMBER, NULL, TRUE, &nuStepsf  },

            { "mu_min",   LUA_TNUMBER, NULL, TRUE, &ia.mu_min },
            { "mu_max",   LUA_TNUMBER, NULL, TRUE, &ia.mu_max },
            { "mu_steps", LUA_TNUMBER, NULL, TRUE, &muStepsf  },

            { "r_min",    LUA_TNUMBER, NULL, TRUE, &ia.r_min  },
            { "r_max",    LUA_TNUMBER, NULL, TRUE, &ia.r_max  },
            { "r_steps",  LUA_TNUMBER, NULL, TRUE, &rStepsf   },
            END_MW_NAMED_ARG
        };

    handleNamedArgumentTable(luaSt, iaArgTable, table);

    ia.nu_steps = (unsigned int) nuStepsf;
    ia.mu_steps = (unsigned int) muStepsf;
    ia.r_steps = (unsigned int) rStepsf;

    r = (uint64_t) ia.r_steps;
    mu = (uint64_t) ia.mu_steps;
    nu = (uint64_t) ia.nu_steps;

    if (nu == 0 || mu == 0 || r == 0)
    {
        mw_printf("Integral size { %u, %u, %u } cannot be 0\n", ia.nu_steps, ia.mu_steps, ia.r_steps);
        return 1;
    }

    if ((r > UINT64_MAX / mu) || ((r * mu) > UINT64_MAX / nu))
    {
        mw_printf("Integral size { %u, %u, %u } will overflow progress calculation\n",
                  ia.nu_steps, ia.mu_steps, ia.r_steps);
        return 1;
    }

    calcIntegralStepSizes(&ia);

    *iaOut = ia;

    return 0;
}

static int readStreamTable(lua_State* luaSt, StreamParameters* spOut, int table)
{
    static StreamParameters sp;
    static const MWNamedArg streamArgTable[] =
        {
            { "epsilon", LUA_TNUMBER, NULL, TRUE, &sp.epsilon },
            { "mu",      LUA_TNUMBER, NULL, TRUE, &sp.mu      },
            { "r",       LUA_TNUMBER, NULL, TRUE, &sp.r       },
            { "theta",   LUA_TNUMBER, NULL, TRUE, &sp.theta   },
            { "phi",     LUA_TNUMBER, NULL, TRUE, &sp.phi     },
            { "sigma",   LUA_TNUMBER, NULL, TRUE, &sp.sigma   },
            END_MW_NAMED_ARG
        };

    handleNamedArgumentTable(luaSt, streamArgTable, table);
    *spOut = sp;

    return 0;
}

static int evaluateStreams(lua_State* luaSt)
{
    int table;
    int i, n;

    lua_getglobal(luaSt, STREAMS_NAME);
    table = lua_gettop(luaSt);

    if (expectTable(luaSt, table))
        luaL_error(luaSt, "Expected '%s' to be a table", STREAMS_NAME);


    n = luaL_getn(luaSt, table);

    /* CHECKME: Is this valid? */
    if (n == 0)
    {
        lua_pop(luaSt, 1);
        return 0;
    }

    _streams->number_streams = n;
    _streams->parameters = mwMalloc(n * sizeof(StreamParameters));

    for (i = 0; i < n; ++i)
    {
        lua_rawgeti(luaSt, table, i + 1);
        readStreamTable(luaSt, &_streams->parameters[i], lua_gettop(luaSt));
        lua_pop(luaSt, 1);
    }

    lua_pop(luaSt, 1);
    return 0;
}

static int evaluateIntegralAreas(lua_State* luaSt)
{
    int i, table;

    lua_getglobal(luaSt, AREAS_NAME);

    table = lua_gettop(luaSt);
    mw_lua_checktable(luaSt, table);

    _nCut = luaL_getn(luaSt, table);

    if (_nCut == 0)
    {
        lua_pop(luaSt, 1);
        return luaL_error(luaSt, "At least one cut required");
    }

    _ias = mwMallocA(_nCut * sizeof(IntegralArea));

    for (i = 0; i < (int) _nCut; ++i)
    {
        lua_rawgeti(luaSt, table, i + 1);
        readIntegralArea(luaSt, &_ias[i], lua_gettop(luaSt));
        lua_pop(luaSt, 1);
    }

    lua_pop(luaSt, 1);
    return 0;
}

static int evaluateConstants(lua_State* luaSt)
{
    const SeparationConstant* c = constants;

    while (c->name)
    {
        lua_getglobal(luaSt, c->name);
        if (mw_lua_typecheck(luaSt, -1, LUA_TNUMBER, NULL))
        {
            return luaL_error(luaSt,
                              "constant '%s': %s",
                              c->name,
                              lua_tostring(luaSt, -1));
        }

        *c->value = lua_tonumber(luaSt, -1);
        ++c;
    }

    return 0;
}

static const BackgroundParameters defaultBG =
{
    /* .alpha   */   1.0,
    /* .r0      */   0.0,
    /* .q       */   0.0,
    /* .delta   */   1.0,
    /* .epsilon */   0.0,
    /* .a       */   0.0,
    /* .b       */   0.0,
    /* .c       */   0.0
};

static int evaluateBackground(lua_State* luaSt)
{
    int table;
    static BackgroundParameters bg = EMPTY_BACKGROUND_PARAMETERS;
    static const MWNamedArg bgArgTable[] =
        {
            { "alpha",   LUA_TNUMBER, NULL, FALSE, &bg.alpha   },
            { "r0",      LUA_TNUMBER, NULL, TRUE,  &bg.r0      },
            { "q",       LUA_TNUMBER, NULL, TRUE,  &bg.q       },
            { "delta",   LUA_TNUMBER, NULL, FALSE, &bg.delta   },

            { "epsilon", LUA_TNUMBER, NULL, FALSE, &bg.epsilon },

            { "a",       LUA_TNUMBER, NULL, FALSE, &bg.a       },
            { "b",       LUA_TNUMBER, NULL, FALSE, &bg.b       },
            { "c",       LUA_TNUMBER, NULL, FALSE, &bg.c       },
            END_MW_NAMED_ARG
        };

    lua_getglobal(luaSt, BACKGROUND_NAME);
    table = lua_gettop(luaSt);
    mw_lua_checktable(luaSt, table);

    bg = defaultBG;
    handleNamedArgumentTable(luaSt, bgArgTable, table);
    *_bg = bg;

    return 0;
}

static int evaluateGlobalName(lua_State* luaSt, lua_CFunction func, const char* name)
{
    lua_pushcfunction(luaSt, func);
    if (lua_pcall(luaSt, 0, 0, 0))
    {
        mw_lua_perror(luaSt, "Error evaluating %s", name);
        return 1;
    }

    return 0;
}

#define NUMBER_FIT_BG_PARAMETERS 2
#define NUMBER_FIT_STREAM_PARAMETERS 6

/* Set arguments using the same fit parameters as have been used */
static int luaDefaultSetBGStreamParametersFromArguments(lua_State* luaSt)
{
    int i, ptable, bgTable, streamTable, n, nStream, idx, thisStream;

    ptable = mw_lua_checktable(luaSt, lua_gettop(luaSt));
    n = luaL_getn(luaSt, ptable);  /* Number parameters */

    /* Set background parameters */
    lua_newtable(luaSt);
    bgTable = lua_gettop(luaSt);
    lua_pushvalue(luaSt, bgTable);
    lua_setglobal(luaSt, BACKGROUND_NAME);

    bgTable = lua_gettop(luaSt);

    lua_rawgeti(luaSt, ptable, 1);
    lua_setfield(luaSt, bgTable, "q");

    lua_rawgeti(luaSt, ptable, 2);
    lua_setfield(luaSt, bgTable, "r0");

    /* Not included in fit */
    lua_pushnumber(luaSt, 0.0);
    lua_setfield(luaSt, bgTable, "epsilon");

    lua_pushnumber(luaSt, 1.0);
    lua_setfield(luaSt, bgTable, "alpha");

    lua_pushnumber(luaSt, 1.0);
    lua_setfield(luaSt, bgTable, "delta");

    lua_pop(luaSt, 1);


    /* Set stream parameters */
    lua_newtable(luaSt);
    streamTable = lua_gettop(luaSt);
    lua_pushvalue(luaSt, streamTable);
    lua_setglobal(luaSt, STREAMS_NAME);

    streamTable = lua_gettop(luaSt);

    if (!mwDivisible(n - NUMBER_FIT_BG_PARAMETERS, NUMBER_FIT_STREAM_PARAMETERS))
    {
        return luaL_error(luaSt, "Parameter count (%d) inconsistent with default argument mapping\n", n);
    }

    nStream = (n - NUMBER_FIT_BG_PARAMETERS) / NUMBER_FIT_STREAM_PARAMETERS;

    for (i = 0; i < nStream; ++i)
    {
        lua_newtable(luaSt);
        thisStream = lua_gettop(luaSt);
        idx = NUMBER_FIT_BG_PARAMETERS + i * NUMBER_FIT_STREAM_PARAMETERS;

        lua_rawgeti(luaSt, ptable, idx + 1);
        lua_setfield(luaSt, thisStream, "epsilon");

        lua_rawgeti(luaSt, ptable, idx + 2);
        lua_setfield(luaSt, thisStream, "mu");

        lua_rawgeti(luaSt, ptable, idx + 3);
        lua_setfield(luaSt, thisStream, "r");

        lua_rawgeti(luaSt, ptable, idx + 4);
        lua_setfield(luaSt, thisStream, "theta");

        lua_rawgeti(luaSt, ptable, idx + 5);
        lua_setfield(luaSt, thisStream, "phi");

        lua_rawgeti(luaSt, ptable, idx + 6);
        lua_setfield(luaSt, thisStream, "sigma");

        mw_lua_assert_top_type(luaSt, LUA_TTABLE);

        lua_rawseti(luaSt, streamTable, i + 1);
    }

    lua_pop(luaSt, 1);
    return 0;
}

/* Open a lua_State and load the stuff we define, but do not run anything */
lua_State* separationLuaOpen(mwbool debug)
{
    lua_State* luaSt;

    luaSt = mw_lua_newstate();
    if (!luaSt)
    {
        mw_printf("Failed to get Lua state\n");
        return NULL;
    }

    mw_lua_openlibs(luaSt, debug);

    mwRegisterTypes(luaSt);
    pushConstants(luaSt);

    lua_register(luaSt, "defaultArgMapping", luaDefaultSetBGStreamParametersFromArguments);

    return luaSt;
}

/* It will be easier to make this less chaotic in the next release
 * when we can dump the old parameters file */
IntegralArea* setupSeparation(AstronomyParameters* ap,
                              BackgroundParameters* bg,
                              Streams* streams,
                              const SeparationFlags* sf)
{
    int rc = 0;
    lua_State* luaSt;

    luaSt = separationOpenLuaStateWithScript(sf);
    if (!luaSt)
        return NULL;

    _ap = ap;
    _bg = bg;
    _streams = streams;

    rc |= evaluateGlobalName(luaSt, evaluateConstants, CONSTANTS_NAME);
    rc |= evaluateGlobalName(luaSt, evaluateBackground, BACKGROUND_NAME);
    rc |= evaluateGlobalName(luaSt, evaluateStreams, STREAMS_NAME);
    rc |= evaluateGlobalName(luaSt, evaluateIntegralAreas, AREAS_NAME);

    lua_close(luaSt);

    if (rc)
    {
        free(_streams->parameters);
        mwFreeA(_ias);
        return NULL;
    }

    setAPConstants(ap);
    return _ias;
}


