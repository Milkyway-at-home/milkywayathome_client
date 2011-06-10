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

#include "milkyway_lua_marshal.h"
#include "milkyway_math.h"
#include "milkyway_lua_math.h"

#define DEFINE_LUA_MW_FUNC_1(cname, name)                               \
     static int lua_##cname(lua_State* luaSt)                           \
     {                                                                  \
         int nArgs;                                                     \
         real number;                                                   \
                                                                        \
         nArgs = lua_gettop(luaSt);                                     \
                                                                        \
         if (nArgs != 1)                                                \
             return luaL_argerror(luaSt, 0, "Expected 1 argument");     \
                                                                        \
         number = (real) luaL_checknumber(luaSt, 1);                    \
         lua_pop(luaSt, 1);                                             \
         lua_pushnumber(luaSt, cname(number));                          \
                                                                        \
         return 1;                                                      \
     }                                                                  \
                                                                        \
     static void register_lua_##cname(lua_State* luaSt)                 \
     {                                                                  \
         lua_register(luaSt, #name, lua_##cname);                       \
     }

#define DEFINE_LUA_MW_FUNC_2(cname, name)                               \
     static int lua_##cname(lua_State* luaSt)                           \
     {                                                                  \
         int nArgs;                                                     \
         real number1, number2;                                         \
                                                                        \
         nArgs = lua_gettop(luaSt);                                     \
                                                                        \
         if (nArgs != 2)                                                \
             return luaL_argerror(luaSt, 0, "Expected 2 arguments");    \
                                                                        \
         number1 = (real) luaL_checknumber(luaSt, 1);                   \
         number2 = (real) luaL_checknumber(luaSt, 2);                   \
         lua_pop(luaSt, 2);                                             \
         lua_pushnumber(luaSt, cname(number1, number2));                \
                                                                        \
         return 1;                                                      \
     }                                                                  \
                                                                        \
     static void register_lua_##cname(lua_State* luaSt)                 \
     {                                                                  \
         lua_register(luaSt, #name, lua_##cname);                       \
     }

DEFINE_LUA_MW_FUNC_1(sqr, sqr)
DEFINE_LUA_MW_FUNC_1(cube, cube)
DEFINE_LUA_MW_FUNC_1(d2r, d2r)
DEFINE_LUA_MW_FUNC_1(r2d, r2d)

DEFINE_LUA_MW_FUNC_1(mw_sin, sin)
DEFINE_LUA_MW_FUNC_1(mw_cos, cos)
DEFINE_LUA_MW_FUNC_1(mw_tan, tan)

DEFINE_LUA_MW_FUNC_1(mw_asin, asin)
DEFINE_LUA_MW_FUNC_1(mw_acos, atan)
DEFINE_LUA_MW_FUNC_1(mw_log, log)
DEFINE_LUA_MW_FUNC_1(mw_exp, exp)


DEFINE_LUA_MW_FUNC_1(mw_log10, log10)

DEFINE_LUA_MW_FUNC_1(mw_sinh, sinh)
DEFINE_LUA_MW_FUNC_1(mw_cosh, cosh)
DEFINE_LUA_MW_FUNC_1(mw_tanh, tanh)


DEFINE_LUA_MW_FUNC_1(mw_fabs, fabs)

DEFINE_LUA_MW_FUNC_1(mw_ceil, ceil)


DEFINE_LUA_MW_FUNC_1(mw_exp2, exp2)
DEFINE_LUA_MW_FUNC_1(mw_exp10, exp10)

DEFINE_LUA_MW_FUNC_1(mw_floor, floor)
DEFINE_LUA_MW_FUNC_2(mw_fmod, fmod)
DEFINE_LUA_MW_FUNC_2(mw_ldexp, ldexp)
DEFINE_LUA_MW_FUNC_1(mw_rsqrt, rsqrt)
DEFINE_LUA_MW_FUNC_1(mw_sqrt, sqrt)

DEFINE_LUA_MW_FUNC_1(mw_tanpi, tanpi)
DEFINE_LUA_MW_FUNC_1(mw_sinpi, sinpi)
DEFINE_LUA_MW_FUNC_1(mw_cospi, cospi)

DEFINE_LUA_MW_FUNC_2(mw_pow, pow)
DEFINE_LUA_MW_FUNC_2(mw_hypot, hypot)
DEFINE_LUA_MW_FUNC_2(mw_atan2, atan2)


#if 0
/* MSVC is missing all of these */
DEFINE_LUA_MW_FUNC_1(mw_expm1, expm1)
DEFINE_LUA_MW_FUNC_2(mw_nextafter, nextafter)
DEFINE_LUA_MW_FUNC_2(mw_remainder, remainder)
DEFINE_LUA_MW_FUNC_2(mw_fmax, fmax)
DEFINE_LUA_MW_FUNC_2(mw_fmin, fmin)
DEFINE_LUA_MW_FUNC_1(mw_erf, erf)
DEFINE_LUA_MW_FUNC_1(mw_erfc, erfc)
DEFINE_LUA_MW_FUNC_1(mw_cbrt, cbrt)
DEFINE_LUA_MW_FUNC_1(mw_log2, log2)
DEFINE_LUA_MW_FUNC_1(mw_rint, rint)

DEFINE_LUA_MW_FUNC_1(mw_tgamma, tgamma)
DEFINE_LUA_MW_FUNC_1(mw_lgamma, lgamma)
DEFINE_LUA_MW_FUNC_1(mw_logb, logb)

DEFINE_LUA_MW_FUNC_1(mw_acosh, acosh)
DEFINE_LUA_MW_FUNC_1(mw_asinh, asinh)
DEFINE_LUA_MW_FUNC_1(mw_atanh, atanh)

DEFINE_LUA_MW_FUNC_1(mw_trunc, trunc)
DEFINE_LUA_MW_FUNC_1(mw_log1p, log1p)
DEFINE_LUA_MW_FUNC_1(mw_round, round)

DEFINE_LUA_MW_FUNC_2(mw_fdim, fdim)
DEFINE_LUA_MW_FUNC_1(mw_ilogb, ilogb)

#endif /* 0 */



// these have different / annoying
//DEFINE_LUA_MW_FUNC_1(mw_fract, fract)
//DEFINE_LUA_MW_FUNC_1(mw_frexp, frexp)
//DEFINE_LUA_MW_FUNC_1(mw_modf, modf)
//DEFINE_LUA_MW_FUNC_1(mw_nan, nan)
//DEFINE_LUA_MW_FUNC_1(mw_remquo, remquo)
// #define mw_sincos sincos

static void registerConstant(lua_State* luaSt, const char* name, real val)
{
    lua_pushnumber(luaSt, val);
    lua_setglobal(luaSt, name);
}

static void registerConstants(lua_State* luaSt)
{
    registerConstant(luaSt, "pi", M_PI);
    registerConstant(luaSt, "pi_2", M_PI_2);
    registerConstant(luaSt, "pi_4", M_PI_4);
    registerConstant(luaSt, "pi_4_3", PI_4_3);
    registerConstant(luaSt, "pi_2_3", PI_2_3);
    registerConstant(luaSt, "pi_3_2", PI_3_2);
    registerConstant(luaSt, "pi2", M_2PI);
    registerConstant(luaSt, "sqrt_2pi", SQRT_2PI);
    registerConstant(luaSt, "ee", M_E);
    registerConstant(luaSt, "log2e", M_LOG2E);
    registerConstant(luaSt, "log10e", M_LOG10E);
    registerConstant(luaSt, "ln2", M_LN2);
    registerConstant(luaSt, "ln10", M_LN10);
    registerConstant(luaSt, "sqrt2", M_SQRT2);

}


void registerMilkywayMath(lua_State* luaSt)
{
    registerConstants(luaSt);

    register_lua_mw_sin(luaSt);
    register_lua_mw_cos(luaSt);
    register_lua_mw_tan(luaSt);
    register_lua_mw_sin(luaSt);
    register_lua_mw_cos(luaSt);
    register_lua_mw_tan(luaSt);
    register_lua_mw_asin(luaSt);
    register_lua_mw_acos(luaSt);
    register_lua_mw_log(luaSt);
    register_lua_mw_exp(luaSt);

    register_lua_mw_log10(luaSt);
    register_lua_mw_sinh(luaSt);
    register_lua_mw_cosh(luaSt);
    register_lua_mw_tanh(luaSt);
    register_lua_mw_fabs(luaSt);
    register_lua_mw_ceil(luaSt);
    register_lua_mw_exp2(luaSt);
    register_lua_mw_exp10(luaSt);
    register_lua_mw_floor(luaSt);
    register_lua_mw_rsqrt(luaSt);
    register_lua_mw_sqrt(luaSt);
    register_lua_mw_cospi(luaSt);
    register_lua_mw_tanpi(luaSt);

    register_lua_mw_fmod(luaSt);
    register_lua_mw_ldexp(luaSt);

    register_lua_mw_sinpi(luaSt);
    register_lua_mw_pow(luaSt);
    register_lua_mw_hypot(luaSt);
    register_lua_mw_atan2(luaSt);

    register_lua_sqr(luaSt);
    register_lua_cube(luaSt);
    register_lua_r2d(luaSt);
    register_lua_d2r(luaSt);

#if 0
    /* The missing from MSVC functions */
    register_lua_mw_log1p(luaSt);
    register_lua_mw_expm1(luaSt);
    register_lua_mw_logb(luaSt);
    register_lua_mw_log2(luaSt);
    register_lua_mw_acosh(luaSt);
    register_lua_mw_asinh(luaSt);
    register_lua_mw_atanh(luaSt);
    register_lua_mw_cbrt(luaSt);
    register_lua_mw_erfc(luaSt);
    register_lua_mw_erf(luaSt);
    register_lua_mw_tgamma(luaSt);
    register_lua_mw_lgamma(luaSt);
    register_lua_mw_rint(luaSt);
    register_lua_mw_round(luaSt);
    register_lua_mw_trunc(luaSt);
    register_lua_mw_fdim(luaSt);
    register_lua_mw_fmax(luaSt);
    register_lua_mw_fmin(luaSt);
    register_lua_mw_remainder(luaSt);
    register_lua_mw_nextafter(luaSt);
    register_lua_mw_ilogb(luaSt);
#endif /* 0 */
}
