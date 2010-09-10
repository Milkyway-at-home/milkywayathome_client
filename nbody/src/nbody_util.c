/* Copyright 2010 Matthew Arsenault, Travis Desell, Boleslaw
Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail and
Rensselaer Polytechnic Institute.

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

#include "nbody_priv.h"
#include "milkyway_util.h"

void cartesianToLbr_rad(const NBodyCtx* ctx, vectorptr restrict lbR, const vectorptr restrict r)
{
    const real xp = X(r) + ctx->sunGCDist;

    L(lbR) = mw_atan2(Y(r), xp);
    B(lbR) = mw_atan2( Z(r), mw_sqrt( sqr(xp) + sqr(Y(r)) ) );
    R(lbR) = mw_sqrt(sqr(xp) + sqr(Y(r)) + sqr(Z(r)));

    if (L(lbR) < 0.0)
        L(lbR) += 2 * M_PI;
}

void cartesianToLbr(const NBodyCtx* ctx, vectorptr restrict lbR, const vectorptr restrict r)
{
    cartesianToLbr_rad(ctx, lbR, r);
    L(lbR) = r2d(L(lbR));
    B(lbR) = r2d(B(lbR));
}

inline static void _lbrToCartesian(vectorptr cart, const real l, const real b, const real r, const real sun)
{
    X(cart) = r * mw_cos(l) * mw_cos(b) - sun;
    Y(cart) = r * mw_sin(l) * mw_cos(b);
    Z(cart) = r * mw_sin(b);
}

void lbrToCartesian_rad(const NBodyCtx* ctx, vectorptr cart, const vectorptr lbr)
{
    _lbrToCartesian(cart, L(lbr), B(lbr), R(lbr), ctx->sunGCDist);
}

void lbrToCartesian(const NBodyCtx* ctx, vectorptr cart, const vectorptr lbr)
{
    _lbrToCartesian(cart, d2r(L(lbr)), d2r(B(lbr)), R(lbr), ctx->sunGCDist);
}

