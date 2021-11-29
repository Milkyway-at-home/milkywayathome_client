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

#include "milkyway_math.h"
#include "milkyway_lua_marshal.h"
#include "milkyway_lua_math.h"

real checkReal(lua_State* luaSt, int idx)
{
    real* num = mw_checknamedudata(luaSt, idx, REAL_TYPE);
    return *num;
}

real toReal(lua_State* luaSt, int idx)
{
    real* num = mw_tonamedudata(luaSt, idx, REAL_TYPE);
    return *num;
}

real expectReal(lua_State* luaSt, int idx)
{
    real* num = expectType(luaSt, idx, REAL_TYPE);
    return *num;
}

int pushReal(lua_State* luaSt, real realIn)
{
    real* realNew;

    realNew = (real*) lua_newuserdata(luaSt, sizeof(real));
    *realNew = realIn;

    luaL_getmetatable(luaSt, REAL_TYPE);
    lua_setmetatable(luaSt, -2);

    return 1;
}

int getReal(lua_State* luaSt, void* r)
{
    pushReal(luaSt, *(real*) r);
    return 1;
}

int setReal(lua_State* luaSt, void* r)
{
    *(real*)r = checkReal(luaSt, 3);
    return 0;
}


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

DEFINE_LUA_MW_FUNC_1(sqr_0, sqr)
DEFINE_LUA_MW_FUNC_1(cube_0, cube)
DEFINE_LUA_MW_FUNC_1(d2r_0, d2r)
DEFINE_LUA_MW_FUNC_1(r2d_0, r2d)

DEFINE_LUA_MW_FUNC_1(mw_sin_0, sin)
DEFINE_LUA_MW_FUNC_1(mw_cos_0, cos)
DEFINE_LUA_MW_FUNC_1(mw_tan_0, tan)

DEFINE_LUA_MW_FUNC_1(mw_asin_0, asin)
DEFINE_LUA_MW_FUNC_1(mw_acos_0, atan)
DEFINE_LUA_MW_FUNC_1(mw_log_0, log)
DEFINE_LUA_MW_FUNC_1(mw_exp_0, exp)


DEFINE_LUA_MW_FUNC_1(mw_log10_0, log10)

DEFINE_LUA_MW_FUNC_1(mw_sinh_0, sinh)
DEFINE_LUA_MW_FUNC_1(mw_cosh_0, cosh)
DEFINE_LUA_MW_FUNC_1(mw_tanh_0, tanh)


DEFINE_LUA_MW_FUNC_1(mw_fabs_0, fabs)

DEFINE_LUA_MW_FUNC_1(mw_ceil_0, ceil)


DEFINE_LUA_MW_FUNC_1(mw_exp2_0, exp2)
DEFINE_LUA_MW_FUNC_1(mw_exp10_0, exp10)

DEFINE_LUA_MW_FUNC_1(mw_floor_0, floor)
DEFINE_LUA_MW_FUNC_2(mw_fmod_0, fmod)
DEFINE_LUA_MW_FUNC_2(mw_ldexp_0, ldexp)
DEFINE_LUA_MW_FUNC_1(mw_rsqrt_0, rsqrt)
DEFINE_LUA_MW_FUNC_1(mw_sqrt_0, sqrt)

DEFINE_LUA_MW_FUNC_1(mw_tanpi_0, tanpi)
DEFINE_LUA_MW_FUNC_1(mw_sinpi_0, sinpi)
DEFINE_LUA_MW_FUNC_1(mw_cospi_0, cospi)

DEFINE_LUA_MW_FUNC_2(mw_pow_0, pow)
DEFINE_LUA_MW_FUNC_2(mw_hypot_0, hypot)
DEFINE_LUA_MW_FUNC_2(mw_atan2_0, atan2)


#if 0
/* MSVC is missing all of these */
DEFINE_LUA_MW_FUNC_1(mw_expm1_0, expm1)
DEFINE_LUA_MW_FUNC_2(mw_nextafter_0, nextafter)
DEFINE_LUA_MW_FUNC_2(mw_remainder_0, remainder)
DEFINE_LUA_MW_FUNC_2(mw_fmax_0, fmax)
DEFINE_LUA_MW_FUNC_2(mw_fmin_0, fmin)
DEFINE_LUA_MW_FUNC_1(mw_erf_0, erf)
DEFINE_LUA_MW_FUNC_1(mw_erfc_0, erfc)
DEFINE_LUA_MW_FUNC_1(mw_cbrt_0, cbrt)
DEFINE_LUA_MW_FUNC_1(mw_log2_0, log2)
DEFINE_LUA_MW_FUNC_1(mw_rint_0, rint)

DEFINE_LUA_MW_FUNC_1(mw_tgamma_0, tgamma)
DEFINE_LUA_MW_FUNC_1(mw_lgamma_0, lgamma)
DEFINE_LUA_MW_FUNC_1(mw_logb_0, logb)

DEFINE_LUA_MW_FUNC_1(mw_acosh_0, acosh)
DEFINE_LUA_MW_FUNC_1(mw_asinh_0, asinh)
DEFINE_LUA_MW_FUNC_1(mw_atanh_0, atanh)

DEFINE_LUA_MW_FUNC_1(mw_trunc_0, trunc)
DEFINE_LUA_MW_FUNC_1(mw_log1p_0, log1p)
DEFINE_LUA_MW_FUNC_1(mw_round_0, round)

DEFINE_LUA_MW_FUNC_2(mw_fdim_0, fdim)
DEFINE_LUA_MW_FUNC_1(mw_ilogb_0, ilogb)

#endif /* 0 */



// these have different / annoying
//DEFINE_LUA_MW_FUNC_1(mw_fract_0, fract)
//DEFINE_LUA_MW_FUNC_1(mw_frexp_0, frexp)
//DEFINE_LUA_MW_FUNC_1(mw_modf_0, modf)
//DEFINE_LUA_MW_FUNC_1(mw_nan_0, nan)
//DEFINE_LUA_MW_FUNC_1(mw_remquo_0, remquo)
// #define mw_sincos_0 sincos

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

    register_lua_mw_sin_0(luaSt);
    register_lua_mw_cos_0(luaSt);
    register_lua_mw_tan_0(luaSt);
    register_lua_mw_asin_0(luaSt);
    register_lua_mw_acos_0(luaSt);
    register_lua_mw_log_0(luaSt);
    register_lua_mw_exp_0(luaSt);

    register_lua_mw_log10_0(luaSt);
    register_lua_mw_sinh_0(luaSt);
    register_lua_mw_cosh_0(luaSt);
    register_lua_mw_tanh_0(luaSt);
    register_lua_mw_fabs_0(luaSt);
    register_lua_mw_ceil_0(luaSt);
    register_lua_mw_exp2_0(luaSt);
    register_lua_mw_exp10_0(luaSt);
    register_lua_mw_floor_0(luaSt);
    register_lua_mw_rsqrt_0(luaSt);
    register_lua_mw_sqrt_0(luaSt);
    register_lua_mw_cospi_0(luaSt);
    register_lua_mw_tanpi_0(luaSt);

    register_lua_mw_fmod_0(luaSt);
    register_lua_mw_ldexp_0(luaSt);

    register_lua_mw_sinpi_0(luaSt);
    register_lua_mw_pow_0(luaSt);
    register_lua_mw_hypot_0(luaSt);
    register_lua_mw_atan2_0(luaSt);

    register_lua_sqr_0(luaSt);
    register_lua_cube_0(luaSt);
    register_lua_r2d_0(luaSt);
    register_lua_d2r_0(luaSt);

#if 0
    /* The missing from MSVC functions */
    register_lua_mw_log1p_0(luaSt);
    register_lua_mw_expm1_0(luaSt);
    register_lua_mw_logb_0(luaSt);
    register_lua_mw_log2_0(luaSt);
    register_lua_mw_acosh_0(luaSt);
    register_lua_mw_asinh_0(luaSt);
    register_lua_mw_atanh_0(luaSt);
    register_lua_mw_cbrt_0(luaSt);
    register_lua_mw_erfc_0(luaSt);
    register_lua_mw_erf_0(luaSt);
    register_lua_mw_tgamma_0(luaSt);
    register_lua_mw_lgamma_0(luaSt);
    register_lua_mw_rint_0(luaSt);
    register_lua_mw_round_0(luaSt);
    register_lua_mw_trunc_0(luaSt);
    register_lua_mw_fdim_0(luaSt);
    register_lua_mw_fmax_0(luaSt);
    register_lua_mw_fmin_0(luaSt);
    register_lua_mw_remainder_0(luaSt);
    register_lua_mw_nextafter_0(luaSt);
    register_lua_mw_ilogb_0(luaSt);
#endif /* 0 */
}
