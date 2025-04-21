/*
This code tests the functions in nbody_coordinates.c to transform 
between galactocentric Cartesian coordinates (pos and vel) and 
proper motions in RA and Dec
*/

#include "nbody_coordinates.h"
#include "nbody_types.h"
#include "nbody_histogram.h"
#include "nbody_mass.h"
#include "milkyway_util.h"
#include "nbody_defaults.h"

struct {
    mwvector v;
    mwvector r;
} TestObj;

TestObj* createTestObj(real x, real y, real z, real vx, real vy, real vz) {
    TestObj* to;
    mwvector* vto = to->v;
    mwvector* rto = to->r;

    X(vto) = vx;
    Y(vto) = vy;
    Z(vto) = vz;

    X(rto) = x;
    Y(rto) = y;
    Z(rto) = z;

    return to;
}

int testObject(TestObj* to, real compareMuDec, real compareMuRA) {
    ctx = defaultNBodyCtx

    real epsilon = 0.5;
    real mudec, mura;

    mudec = nbVXVYVZtomuDec(a->r,a->v,ctx->sunVelx,ctx->sunVely,ctx->sunVelz,ctx->sunGCDist,ctx->NGPdec,ctx->lNCP);
    mura = nbVXVYVZtomuRA(a->r,a->v,ctx->sunVelx,ctx->sunVely,ctx->sunVelz,ctx->sunGCDist,ctx->NGPdec,ctx->lNCP);
    mw_printf("mura: %.15f, mudec: %.15f", , -mudec, -mura);

    if ((mw_abs(mudec - compareMuDec) >= epsilon) || (mw_abs(mura - compareMuRA) >= epsilon)) {
        return 1;
    }

    return 0;
}

int main() {
    TestObj* a = createTestObj(-10.0,0.0,0.0,10.3,229.2,6.9); //matching sun velocity
    TestObj* b = createTestObj(3.3,2.8,-3.1,10.3,229.2,6.9); //matching sun velocity
    TestObj* c = createTestObj(-10.0,0.0,0.0,40.1,229.2,6.9); //moving only in line of sight relative to sun
    TestObj* gc = createTestObj(0.0,0.0,0.0,0.0,0.0,0.0); //Galactic center

    mw_printf("Proper motion of test objects. Should be (0.0,0.0), (0.0,0.0), (0.0,0.0), and (-2.7,-5.6) mura,mudec")
    if (testObject(a,0.0,0.0) == 1) {
        mw_printf("Proper Motion Test Failed: An object matching solar velocity should have 0 proper motion");
        return 1;
    } else if (testObject(b,0.0,0.0) == 1) {
        mw_printf("Proper Motion Test Failed: An object matching solar velocity should have 0 proper motion");
        return 1;
    } else if (testObject(c,0.0,0.0) == 1) {
        mw_printf("Proper Motion Test Failed: An object moving along line of sight should have 0 proper motion");
        return 1;
    } else if (testObject(gc,-5.6,-2.7) == 1) {
        mw_printf("Proper Motion Test Failed: The Galactic center should have -2.7 mura, -5.6 mudec");
        return 1;
    }

    return 0;
}
