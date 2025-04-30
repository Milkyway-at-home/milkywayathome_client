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

mwvector* createTestObj(real x, real y, real z, real vx, real vy, real vz) {
    mwvector* vr = mwCalloc(2, sizeof(mwvector));

    if (vr == NULL) {
        printf("Allocation Failed");
        exit(0);
    }

    X(vr[1]) = vx;
    Y(vr[1]) = vy;
    Z(vr[1]) = vz;

    X(vr[0]) = x;
    Y(vr[0]) = y;
    Z(vr[0]) = z;

    return vr;
}

int testObject(mwvector* to, real compareMuDec, real compareMuRA) {
    NBodyCtx ctx = defaultNBodyCtx;

    real epsilon = 0.1;
    real mudec, mura;

    mudec = nbVXVYVZtomuDec(to[0],to[1],ctx.sunVelx,ctx.sunVely,ctx.sunVelz,ctx.sunGCDist,ctx.NGPdec,ctx.lNCP);
    mura = nbVXVYVZtomuRA(to[0],to[1],ctx.sunVelx,ctx.sunVely,ctx.sunVelz,ctx.sunGCDist,ctx.NGPdec,ctx.lNCP);
    mw_printf("mura: %.15f, mudec: %.15f\n", -mudec, -mura);

    if ((mw_abs(mudec - compareMuDec) >= epsilon) || (mw_abs(mura - compareMuRA) >= epsilon)) {
        return 1;
    }

    return 0;
}

int main() {
    int failed = 0;
    mwvector* a = createTestObj(-10.0,0.0,0.0,10.3,229.2,6.9); //matching sun velocity
    mwvector* b = createTestObj(3.3,2.8,-3.1,10.3,229.2,6.9); //matching sun velocity
    mwvector* c = createTestObj(-10.0,0.0,0.0,40.1,229.2,6.9); //moving only in line of sight relative to sun
    mwvector* gc = createTestObj(0.0,0.0,0.0,0.0,0.0,0.0); //Galactic center

    mw_printf("Rotation matrix test. Should rotate a length 1 vector into a length 1 vector.");
    mwvector test_vec;
    X(test_vec) = 1;
    Y(test_vec) = 0;
    Z(test_vec) = 0;
    test_vec = cartesianalign(test_vec, d2r(27.4), d2r(192), d2r(123));
    real epsilon = 0.001;
    if (mw_abs(mw_length(test_vec)-1) >= epsilon ) {
        mw_printf("Rotation matrix not unitary\n");
        failed = 1;
    }

    mw_printf("Proper motion of test objects. Should be (0.0,0.0), (0.0,0.0), (0.0,0.0), and (-2.7,-5.6) mura,mudec\n");
    if (testObject(a,0.0,0.0) == 1) {
        mw_printf("Proper Motion Test Failed: An object matching solar velocity should have 0 proper motion\n");
        failed = 1;
    } if (testObject(b,0.0,0.0) == 1) {
        mw_printf("Proper Motion Test Failed: An object matching solar velocity should have 0 proper motion\n");
        failed = 1;
    } if (testObject(c,0.0,0.0) == 1) {
        mw_printf("Proper Motion Test Failed: An object moving along line of sight should have 0 proper motion\n");
        failed = 1;
    } if (testObject(gc,-5.34,-2.54) == 1) {
        mw_printf("Proper Motion Test Failed: The Galactic center should have about -2.6 mura, -5.4 mudec\n");
        failed = 1;
    }

    free(a);
    free(b);
    free(c);
    free(gc);

    return failed;
}
