/*
 * Copyright (c) 2011 Matthew Arsenault
 * Copyright (c) 2016-2018 Siddhartha Shelton
 * This file is part of Milkway@Home.
 *
 * Milkyway@Home is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Milkyway@Home is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Milkyway@Home.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <lua.h>
#include <lualib.h>
#include <lauxlib.h>

#include "nbody_types.h"
#include "nbody_show.h"
#include "nbody_lua_nbodyctx.h"
#include "nbody_lua_potential.h"
#include "milkyway_lua.h"
#include "milkyway_util.h"
#include "nbody_defaults.h"
#include "nbody_lua_misc.h"
#include "nbody_util.h"
#include "nbody_lua_util.h"
#include "nbody_potential.h"
#include "nbody_density.h"

static const real pi __attribute__((unused)) = 3.1415926535;

static const MWEnumAssociation criterionOptions[] =
{
    { "TreeCode",     TreeCode     },
    { "Exact",        Exact        },
    { "BH86",         BH86         },
    { "SW93",         SW93         },
    END_MW_ENUM_ASSOCIATION
};

static int getCriterionT(lua_State* luaSt, void* v)
{
    return pushEnum(luaSt, criterionOptions, *(criterion_t*) v);
}

static int setCriterionT(lua_State* luaSt, void* v)
{
    *(criterion_t*) v = checkEnum(luaSt, criterionOptions, 3);
    return 0;
}

/* eps2 is a flattened eps2_size x eps2_size table of softening lengths, and
 * eps2_index is a list of eps2_size type labels indexing into it. Both are
 * variable-length (sized by eps2_size), so unlike a fixed-length member the
 * generic offset-only getRealArray/setRealArray can't be used directly here:
 * they'd need to know the length, but only get a pointer to the member
 * itself. These wrappers recover the enclosing NBodyCtx via offsetof
 * arithmetic (the same struct these fields live in) so they can look up
 * ctx->eps2_size and use it as the length. */

static int getEps2(lua_State* luaSt, void* v)
{
    NBodyCtx* ctx = (NBodyCtx*) ((char*) v - offsetof(NBodyCtx, eps2));
    return getRealArray(luaSt, v, ctx->eps2_size * ctx->eps2_size);
}

static int setEps2(lua_State* luaSt, void* v)
{
    NBodyCtx* ctx = (NBodyCtx*) ((char*) v - offsetof(NBodyCtx, eps2));
    size_t len;

    if (!lua_istable(luaSt, 3))
    {
        return luaL_error(luaSt, "Expected table");
    }

    len = luaL_getn(luaSt, 3);
    if (len != ctx->eps2_size * ctx->eps2_size)
    {
        return luaL_error(luaSt, "eps2 table has %d entries, expected eps2_size^2 = %d",
                           (int) len, (int) (ctx->eps2_size * ctx->eps2_size));
    }

    free(*(real**) v);
    *(real**) v = (real*) mwCalloc(len, sizeof(real));
    for (size_t i = 0; i < len; ++i)
    {
        lua_rawgeti(luaSt, 3, i + 1);
        (*(real**) v)[i] = (real) lua_tonumber(luaSt, -1);
        lua_pop(luaSt, 1);
    }

    return 0;
}

static int getEps2Index(lua_State* luaSt, void* v)
{
    NBodyCtx* ctx = (NBodyCtx*) ((char*) v - offsetof(NBodyCtx, eps2_index));
    return getIntArray(luaSt, v, ctx->eps2_size);
}

static int setEps2Index(lua_State* luaSt, void* v)
{
    NBodyCtx* ctx = (NBodyCtx*) ((char*) v - offsetof(NBodyCtx, eps2_index));
    size_t len;

    if (!lua_istable(luaSt, 3))
    {
        return luaL_error(luaSt, "Expected table");
    }

    len = luaL_getn(luaSt, 3);
    if (len != ctx->eps2_size)
    {
        return luaL_error(luaSt, "eps2_index table has %d entries, expected eps2_size = %d",
                           (int) len, (int) ctx->eps2_size);
    }

    free(*(int**) v);
    *(int**) v = (int*) mwCalloc(len, sizeof(int));
    for (size_t i = 0; i < len; ++i)
    {
        lua_rawgeti(luaSt, 3, i + 1);
        (*(int**) v)[i] = (int) lua_tointeger(luaSt, -1);
        lua_pop(luaSt, 1);
    }

    return 0;
}

/* Length of the table stored under `key` in the named-argument table at
 * `table`'s stack index, without disturbing anything else on the stack.
 * Used to cross-check eps2/eps2_index against eps2_size right after they're
 * all read in by handleNamedArgumentTable(). */
static size_t luaNamedTableFieldLen(lua_State* luaSt, int table, const char* key)
{
    size_t len;
    lua_getfield(luaSt, table, key);
    len = (size_t) luaL_getn(luaSt, -1);
    lua_pop(luaSt, 1);
    return len;
}

NBodyCtx* checkNBodyCtx(lua_State* luaSt, int idx)
{
    return (NBodyCtx*) mw_checknamedudata(luaSt, idx, NBODYCTX_TYPE);
}

NBodyCtx* toNBodyCtx(lua_State* luaSt, int idx)
{
    return (NBodyCtx*) mw_tonamedudata(luaSt, idx, NBODYCTX_TYPE);
}

NBodyCtx* expectNBodyCtx(lua_State* luaSt, int idx)
{
    return (NBodyCtx*) expectType(luaSt, idx, NBODYCTX_TYPE);
}

int pushNBodyCtx(lua_State* luaSt, const NBodyCtx* p)
{
    return pushType(luaSt, NBODYCTX_TYPE, sizeof(NBodyCtx), (void*) p);
}

criterion_t readCriterion(lua_State* luaSt, const char* name)
{
    return (criterion_t) readEnum(luaSt, criterionOptions, name);
}

static int createNBodyCtx(lua_State* luaSt)
{
    static NBodyCtx ctx;
    static const char* criterionName = NULL;
    real nStepf = 0.0;

    static const MWNamedArg argTable[] =
        {  
            { "timestep",        LUA_TNUMBER, NULL,      TRUE,  &ctx.timestep,          1 },
            { "timeEvolve",      LUA_TNUMBER, NULL,      TRUE,  &ctx.timeEvolve,        1 },
            { "timeBack",        LUA_TNUMBER, NULL,      FALSE, &ctx.timeBack,          1 },
            { "theta",           LUA_TNUMBER, NULL,      FALSE, &ctx.theta,             1 },
            { "eps2",            LUA_TTABLE,  REAL_TYPE, TRUE,  &ctx.eps2,              0 }, //arrayLen 0 = variable length; sized from the actual Lua table (eps2_size^2 entries)
            { "eps2_index",      LUA_TTABLE,  INT_TYPE,  TRUE,  &ctx.eps2_index,        0 }, //arrayLen 0 = variable length; sized from the actual Lua table (eps2_size entries)
            { "eps2_size",       LUA_TNUMBER, SIZE_TYPE, TRUE,  &ctx.eps2_size,         1 },
            { "treeRSize",       LUA_TNUMBER, NULL,      FALSE, &ctx.treeRSize,         1 },
            { "sunGCDist",       LUA_TNUMBER, NULL,      FALSE, &ctx.sunGCDist,         1 },
            { "sunVelx",         LUA_TNUMBER, NULL,      FALSE, &ctx.sunVelx,           1 },
            { "sunVely",         LUA_TNUMBER, NULL,      FALSE, &ctx.sunVely,           1 },
            { "sunVelz",         LUA_TNUMBER, NULL,      FALSE, &ctx.sunVelz,           1 },
 
            { "b",               LUA_TNUMBER, NULL,      FALSE, &ctx.b,                 1 },
            { "r",               LUA_TNUMBER, NULL,      FALSE, &ctx.r,                 1 },
            { "vx",              LUA_TNUMBER, NULL,      FALSE, &ctx.vx,                1 },
            { "vy",              LUA_TNUMBER, NULL,      FALSE, &ctx.vy,                1 },
            { "vz",              LUA_TNUMBER, NULL,      FALSE, &ctx.vz,                1 },

            { "criterion",       LUA_TSTRING, NULL,      FALSE, &criterionName,         1 },
            { "SimpleOutput",    LUA_TBOOLEAN,NULL,      FALSE, &ctx.SimpleOutput,      1 },
            { "useQuad",         LUA_TBOOLEAN,NULL,      FALSE, &ctx.useQuad,           1 },
            { "allowIncest",     LUA_TBOOLEAN,NULL,      FALSE, &ctx.allowIncest,       1 },
            { "quietErrors",     LUA_TBOOLEAN,NULL,      FALSE, &ctx.quietErrors,       1 },
            { "useBestLike",     LUA_TBOOLEAN,NULL,      FALSE, &ctx.useBestLike,       1 },
            { "BestLikeStart",   LUA_TNUMBER, NULL,      FALSE, &ctx.BestLikeStart,     1 },
            { "useVelDisp",      LUA_TBOOLEAN,NULL,      FALSE, &ctx.useVelDisp,        1 },
            { "useBetaDisp",     LUA_TBOOLEAN,NULL,      FALSE, &ctx.useBetaDisp,       1 },
            { "useBetaComp",     LUA_TBOOLEAN,NULL,      FALSE, &ctx.useBetaComp,       1 },
            { "useVlos",         LUA_TBOOLEAN,NULL,      FALSE, &ctx.useVlos,           1 },
            { "useDist",         LUA_TBOOLEAN,NULL,      FALSE, &ctx.useDist,           1 },
            { "usePropMot",      LUA_TBOOLEAN,NULL,      FALSE, &ctx.usePropMot,        1 },
            { "useMomentum",     LUA_TBOOLEAN,NULL,      FALSE, &ctx.useMomentum,       1 },
            { "Ntsteps",         LUA_TNUMBER, NULL,      FALSE, &ctx.Ntsteps,           1 },
            { "Nstep_control",   LUA_TBOOLEAN,NULL,      FALSE, &ctx.Nstep_control,     1 }, 
            { "MultiOutput",     LUA_TBOOLEAN,NULL,      FALSE, &ctx.MultiOutput,       1 },
            { "InitialOutput",   LUA_TBOOLEAN,NULL,      FALSE, &ctx.InitialOutput,     1 },
            { "OutputFreq",      LUA_TNUMBER, NULL,      FALSE, &ctx.OutputFreq,        1 },
            { "BetaSigma",       LUA_TNUMBER, NULL,      TRUE,  &ctx.BetaSigma,         1 },
            { "VelSigma",        LUA_TNUMBER, NULL,      TRUE,  &ctx.VelSigma,          1 },
            { "DistSigma",       LUA_TNUMBER, NULL,      TRUE,  &ctx.DistSigma,         1 },
            { "PMSigma",         LUA_TNUMBER, NULL,      TRUE,  &ctx.PMSigma,           1 },
            { "MomentumSigma",   LUA_TNUMBER, NULL,      TRUE,  &ctx.MomentumSigma,     1 },
            { "IterMax",         LUA_TNUMBER, NULL,      FALSE, &ctx.IterMax,           1 },
            { "BetaCorrect",     LUA_TNUMBER, NULL,      TRUE,  &ctx.BetaCorrect,       1 },
            { "VelCorrect",      LUA_TNUMBER, NULL,      TRUE,  &ctx.VelCorrect,        1 },
            { "DistCorrect",     LUA_TNUMBER, NULL,      TRUE,  &ctx.DistCorrect,       1 },
            { "PMCorrect",       LUA_TNUMBER, NULL,      TRUE,  &ctx.PMCorrect,         1 },
            { "MomentumCorrect", LUA_TNUMBER, NULL,      TRUE,  &ctx.MomentumCorrect,   1 },
            { "LMC",             LUA_TBOOLEAN,NULL,      FALSE, &ctx.LMC,               1 },
	        { "LMCfunction",     LUA_TNUMBER, NULL,      FALSE, &ctx.LMCfunction,       1 },
            { "LMCmass",         LUA_TNUMBER, NULL,      FALSE, &ctx.LMCmass,           1 },
            { "LMCscale",        LUA_TNUMBER, NULL,      FALSE, &ctx.LMCscale,          1 },
	        { "LMCscale2",       LUA_TNUMBER, NULL,      FALSE, &ctx.LMCscale2,         1 },
            { "LMCDynaFric",     LUA_TBOOLEAN,NULL,      FALSE, &ctx.LMCDynaFric,       1 },
            { "coulomb_log",     LUA_TNUMBER, NULL,      FALSE, &ctx.coulomb_log,       1 },
            { "calibrationRuns", LUA_TNUMBER, UINT_TYPE, FALSE, &ctx.calibrationRuns,   1 },
            { "samplingBounds",  LUA_TTABLE,  REAL_TYPE, FALSE, &ctx.samplingBounds,    2 },
            END_MW_NAMED_ARG
        };

    criterionName = NULL;
    ctx = defaultNBodyCtx;

    if (lua_gettop(luaSt) != 1)
        return luaL_argerror(luaSt, 1, "Expected named argument table");

    handleNamedArgumentTable(luaSt, argTable, 1);

    /* eps2 must hold exactly eps2_size^2 entries (a flattened eps2_size x
     * eps2_size table of pairwise softening lengths) and eps2_index must
     * hold exactly eps2_size entries (the type label for each row/column).
     * These are read in as separate named arguments above with independent,
     * variable lengths taken from whatever tables the workunit provided, so
     * nothing else guarantees they agree with eps2_size -- catch a mismatch
     * here instead of letting nbGravity/nbGravity_Exact index out of bounds
     * with it later. */
    if (ctx.eps2_size == 0)
    {
        return luaL_argerror(luaSt, 1, "eps2_size must be at least 1");
    }

    if (luaNamedTableFieldLen(luaSt, 1, "eps2") != ctx.eps2_size * ctx.eps2_size)
    {
        return luaL_argerror(luaSt, 1, "eps2 table length does not match eps2_size^2");
    }

    if (luaNamedTableFieldLen(luaSt, 1, "eps2_index") != ctx.eps2_size)
    {
        return luaL_argerror(luaSt, 1, "eps2_index table length does not match eps2_size");
    }

    /* FIXME: Hacky handling of enum. Will result in not good error
     * messages as well as not fitting in. */
    if (criterionName) /* Not required */
    {
        ctx.criterion = readCriterion(luaSt, criterionName);
    }

    if ((ctx.criterion != Exact) && (ctx.theta < 0.0))
    {
        return luaL_argerror(luaSt, 1, "Theta argument required for criterion != 'Exact'");
    }
    else if (ctx.criterion == Exact)
    {
        /* These don't mean anything here */
        ctx.theta = 0.0;
        ctx.useQuad = FALSE;
    }

    nStepf = mw_ceil(ctx.timeEvolve / ctx.timestep);
    if (nStepf >= (real) UINT_MAX)
    {
        luaL_error(luaSt,
                   "Number of timesteps exceeds UINT_MAX: %f timesteps (%f / %f)\n",
                   nStepf,
                   ctx.timeEvolve, ctx.timestep);
    }
    
    ctx.nStep = (unsigned int) nStepf;
    
    #ifdef NBODY_DEV_OPTIONS
        if(ctx.Nstep_control)
        {
            mw_printf("BE WARNED: manually controlling time is unnatural and should be used with the utmost caution.\n");
            ctx.nStep = (int) ctx.Ntsteps;
        }
    #endif
    
    /*
    This looks really old and I don't think we have to check for WUs with version < 0.9 anymore.
    I'm going to remove it (otherwise nbody lite prints two version statements due to the nbReadMinVersion
    call here) and if it ends up breaking something, put it back later */
    // {
        // int major = 0, minor = 0;

        /* Automatically correct the timestep size so an integer
         * number of timesteps covers the evolution time.
         *
         * Only do this if we require a minimum version of 0.90 to
         * avoid not validating against currently existing workunits
         */

        // if (  !nbReadMinVersion(luaSt, &major, &minor)    /* If we fail to read version */
        //     || (major > 0 || (major == 0 && minor >= 90)) /* Version required >= 0.90 */
        //     || (major == 0 && minor == 0))                /* Min version not set */
        // {
        //     ctx.timestep = nbCorrectTimestep(ctx.timeEvolve, ctx.timestep);
        // }
        // else
        // {
        //     mw_printf("Warning: not applying timestep correction for workunit with min version %d.%d\n", major, minor);
        // }
    // }
    
    //instead, just put this here -Tom
    ctx.timestep = nbCorrectTimestep(ctx.timeEvolve, ctx.timestep);

    pushNBodyCtx(luaSt, &ctx);
    return 1;
}

static int toStringNBodyCtx(lua_State* luaSt)
{
    return toStringType(luaSt, (StructShowFunc) showNBodyCtx, (LuaTypeCheckFunc) checkNBodyCtx);
}

static int eqNBodyCtx(lua_State* luaSt)
{
    lua_pushboolean(luaSt, equalNBodyCtx(checkNBodyCtx(luaSt, 1), checkNBodyCtx(luaSt, 2)));
    return 1;
}

static int addPotential(lua_State* luaSt)
{
    NBodyCtx* ctx;

    if (lua_gettop(luaSt) != 2)
        return luaL_argerror(luaSt, 1, "Expected named 2 arguments");

    ctx = checkNBodyCtx(luaSt, 1);
    ctx->pot = *checkPotential(luaSt, 2);

    return 0;
}

static const luaL_reg metaMethodsNBodyCtx[] =
{
    { "__tostring", toStringNBodyCtx },
    { "__eq",       eqNBodyCtx       },
    { NULL, NULL }
};

static const luaL_reg methodsNBodyCtx[] =
{
    { "create", createNBodyCtx },
    { "addPotential", addPotential },
    { NULL, NULL }
};

static const Xet_reg_pre gettersNBodyCtx[] =
{
    { "timestep",        getNumber,     offsetof(NBodyCtx, timestep)       },
    { "timeEvolve",      getNumber,     offsetof(NBodyCtx, timeEvolve)     },
    { "timeBack",        getNumber,     offsetof(NBodyCtx, timeBack)       },
    { "theta",           getNumber,     offsetof(NBodyCtx, theta)          },
    { "eps2",            getEps2,       offsetof(NBodyCtx, eps2)           },
    { "eps2_index",      getEps2Index,  offsetof(NBodyCtx, eps2_index)     },
    { "eps2_size",       getSizeT,      offsetof(NBodyCtx, eps2_size)      },
    { "treeRSize",       getNumber,     offsetof(NBodyCtx, treeRSize)      },
    { "sunGCDist",       getNumber,     offsetof(NBodyCtx, sunGCDist)      },
    { "sunVelx",         getNumber,     offsetof(NBodyCtx, sunVelx)        },
    { "sunVely",         getNumber,     offsetof(NBodyCtx, sunVely)        },
    { "sunVelz",         getNumber,     offsetof(NBodyCtx, sunVelz)        },
    { "criterion",       getCriterionT, offsetof(NBodyCtx, criterion)      },
    { "SimpleOutput",    getBool,       offsetof(NBodyCtx, SimpleOutput)   },
    { "useQuad",         getBool,       offsetof(NBodyCtx, useQuad)        },
    { "allowIncest",     getBool,       offsetof(NBodyCtx, allowIncest)    },
    { "quietErrors",     getBool,       offsetof(NBodyCtx, quietErrors)    },
    { "useBestLike",     getBool,       offsetof(NBodyCtx, useBestLike)    },
    { "useVelDisp",      getBool,       offsetof(NBodyCtx, useVelDisp)     },
    { "useBetaDisp",     getBool,       offsetof(NBodyCtx, useBetaDisp)    },
    { "useBetaComp",     getBool,       offsetof(NBodyCtx, useBetaComp)    },
    { "useVlos",         getBool,       offsetof(NBodyCtx, useVlos)        },
    { "useDist",         getBool,       offsetof(NBodyCtx, useDist)        },
    { "usePropMot",      getBool,       offsetof(NBodyCtx, usePropMot)     },
    { "useMomentum",     getBool,       offsetof(NBodyCtx, useMomentum)    },
    { "BestLikeStart",   getNumber,     offsetof(NBodyCtx, BestLikeStart)  },
    { "Nstep_control",   getBool,       offsetof(NBodyCtx, Nstep_control)  },
    { "Ntsteps",         getNumber,     offsetof(NBodyCtx, Ntsteps)        },
    { "MultiOutput",     getBool,       offsetof(NBodyCtx, MultiOutput)    },
    { "OutputFreq",      getNumber,     offsetof(NBodyCtx, OutputFreq)     },
    { "InitialOutput",   getBool,       offsetof(NBodyCtx, InitialOutput)  },
    { "BetaSigma",       getNumber,     offsetof(NBodyCtx, BetaSigma)      },
    { "VelSigma",        getNumber,     offsetof(NBodyCtx, VelSigma)       },
    { "DistSigma",       getNumber,     offsetof(NBodyCtx, DistSigma)      },
    { "PMSigma",         getNumber,     offsetof(NBodyCtx, PMSigma)        },
    { "MomentumSigma",   getNumber,     offsetof(NBodyCtx, MomentumSigma)  },
    { "IterMax",         getNumber,     offsetof(NBodyCtx, IterMax)        },
    { "BetaCorrect",     getNumber,     offsetof(NBodyCtx, BetaCorrect)    },
    { "VelCorrect",      getNumber,     offsetof(NBodyCtx, VelCorrect)     },
    { "DistCorrect",     getNumber,     offsetof(NBodyCtx, DistCorrect)    },
    { "PMCorrect",       getNumber,     offsetof(NBodyCtx, PMCorrect)      },
    { "MomentumCorrect", getNumber,     offsetof(NBodyCtx, MomentumCorrect)},
    { "LMC",             getBool,       offsetof(NBodyCtx, LMC)            },
    { "LMCfunction",     getNumber,     offsetof(NBodyCtx, LMCfunction)    },
    { "LMCmass",         getNumber,     offsetof(NBodyCtx, LMCmass)        },
    { "LMCscale",        getNumber,     offsetof(NBodyCtx, LMCscale)       },
    { "LMCscale2",       getNumber,     offsetof(NBodyCtx, LMCscale2)      },
    { "LMCDynaFric",     getBool,       offsetof(NBodyCtx, LMCDynaFric)    },
    { "coulomb_log",     getNumber,     offsetof(NBodyCtx, coulomb_log)    },
    { "calibrationRuns", getNumber,     offsetof(NBodyCtx, calibrationRuns)},
    { NULL, NULL, 0 }
};

static const Xet_reg_pre settersNBodyCtx[] =
{
    { "timestep",        setNumber,     offsetof(NBodyCtx, timestep)       },
    { "timeEvolve",      setNumber,     offsetof(NBodyCtx, timeEvolve)     },
    { "timeBack",        setNumber,     offsetof(NBodyCtx, timeBack)       },
    { "theta",           setNumber,     offsetof(NBodyCtx, theta)          },
    { "eps2",            setEps2,       offsetof(NBodyCtx, eps2)           },
    { "eps2_index",      setEps2Index,  offsetof(NBodyCtx, eps2_index)     },
    { "eps2_size",       setSizeT,      offsetof(NBodyCtx, eps2_size)      },
    { "treeRSize",       setNumber,     offsetof(NBodyCtx, treeRSize)      },
    { "sunGCDist",       setNumber,     offsetof(NBodyCtx, sunGCDist)      },
    { "sunVelx",         setNumber,     offsetof(NBodyCtx, sunVelx)        },
    { "sunVely",         setNumber,     offsetof(NBodyCtx, sunVely)        },
    { "sunVelz",         setNumber,     offsetof(NBodyCtx, sunVelz)        },
    { "criterion",       setCriterionT, offsetof(NBodyCtx, criterion)      },
    { "SimpleOutput",    setBool,       offsetof(NBodyCtx, SimpleOutput)   },
    { "useQuad",         setBool,       offsetof(NBodyCtx, useQuad)        },
    { "allowIncest",     setBool,       offsetof(NBodyCtx, allowIncest)    },
    { "quietErrors",     setBool,       offsetof(NBodyCtx, quietErrors)    },
    { "useBestLike",     setBool,       offsetof(NBodyCtx, useBestLike)    },
    { "useVelDisp",      setBool,       offsetof(NBodyCtx, useVelDisp)     },
    { "useBetaDisp",     setBool,       offsetof(NBodyCtx, useBetaDisp)    },
    { "useBetaComp",     setBool,       offsetof(NBodyCtx, useBetaComp)    },
    { "useVlos",         setBool,       offsetof(NBodyCtx, useVlos)        },
    { "useDist",         setBool,       offsetof(NBodyCtx, useDist)        },
    { "usePropMot",      setBool,       offsetof(NBodyCtx, usePropMot)     },
    { "useMomentum",     setBool,       offsetof(NBodyCtx, useMomentum)    },
    { "BestLikeStart",   setNumber,     offsetof(NBodyCtx, BestLikeStart)  },
    { "Nstep_control",   setBool,       offsetof(NBodyCtx, Nstep_control)  },
    { "Ntsteps",         setNumber,     offsetof(NBodyCtx, Ntsteps)        },
    { "MultiOutput",     setBool,       offsetof(NBodyCtx, MultiOutput)    },
    { "OutputFreq",      setNumber,     offsetof(NBodyCtx, OutputFreq)     },
    { "InitialOutput",   setBool,       offsetof(NBodyCtx, InitialOutput)  },
    { "BetaSigma",       setNumber,     offsetof(NBodyCtx, BetaSigma)      },
    { "VelSigma",        setNumber,     offsetof(NBodyCtx, VelSigma)       },
    { "DistSigma",       setNumber,     offsetof(NBodyCtx, DistSigma)      },
    { "PMSigma",         setNumber,     offsetof(NBodyCtx, PMSigma)        },
    { "MomentumSigma",   setNumber,     offsetof(NBodyCtx, MomentumSigma)  },
    { "IterMax",         setNumber,     offsetof(NBodyCtx, IterMax)        },
    { "BetaCorrect",     setNumber,     offsetof(NBodyCtx, BetaCorrect)    },
    { "VelCorrect",      setNumber,     offsetof(NBodyCtx, VelCorrect)     },
    { "DistCorrect",     setNumber,     offsetof(NBodyCtx, DistCorrect)    },
    { "PMCorrect",       setNumber,     offsetof(NBodyCtx, PMCorrect)      },
    { "MomentumCorrect", setNumber,     offsetof(NBodyCtx, MomentumCorrect)},
    { "LMC",             setBool,       offsetof(NBodyCtx, LMC)            },
    { "LMCfunction",     setNumber,     offsetof(NBodyCtx, LMCfunction)    },
    { "LMCmass",         setNumber,     offsetof(NBodyCtx, LMCmass)        },
    { "LMCscale",        setNumber,     offsetof(NBodyCtx, LMCscale)       },
    { "LMCscale2",       setNumber,     offsetof(NBodyCtx, LMCscale2)      },
    { "LMCDynaFric",     setBool,       offsetof(NBodyCtx, LMCDynaFric)    },
    { "coulomb_log",     setNumber,     offsetof(NBodyCtx, coulomb_log)    },
    { "calibrationRuns", setNumber,     offsetof(NBodyCtx, calibrationRuns)},
    { NULL, NULL, 0 }
};

int registerNBodyCtx(lua_State* luaSt)
{
    return registerStruct(luaSt,
                          NBODYCTX_TYPE,
                          gettersNBodyCtx,
                          settersNBodyCtx,
                          metaMethodsNBodyCtx,
                          methodsNBodyCtx);
}
