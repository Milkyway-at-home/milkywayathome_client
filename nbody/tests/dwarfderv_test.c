/* This program tests the derivative calculations of a 2-component dwarf galaxy.*/

#include "milkyway_util.h"
#include "nbody_mixeddwarf.h"
#include "dSFMT.h"
#include <stdio.h>

#define autodiff_thresh (0.0001)

static inline real_0 Qfunc(real_0 v_s, real_0 v_esc)
{
    real_0 ratio = v_s/v_esc;
    return v_s * mw_pow_0(1.0 - sqr_0(ratio), -7/2) * (24 - 78*sqr_0(ratio) + 80*fourth_0(ratio) - 27*sixth_0(ratio));
}

static inline real_0 Pfunc(real_0 v_s, real_0 v_esc)
{
    real_0 ratio = v_s/v_esc;
    real_0 part1 = (24 - 378*sqr_0(ratio) + 712*fourth_0(ratio) - 349*sixth_0(ratio)) * (24 - 78*sqr_0(ratio) + 80*fourth_0(ratio) - 27*sixth_0(ratio));
    real_0 part2 = 24 * sqr_0(ratio) * mw_pow_0(1.0 - sqr_0(ratio), 7/2) * (324 - 710*sqr_0(ratio) + 402*fourth_0(ratio) - 27*sixth_0(ratio));
    return v_s * mw_pow_0(1.0 - sqr_0(ratio), -8) * (part1 + part2);
}

static inline int checkRadDerv(mwvector* pos, real_0 a)
{
    int failed = 0;
    real radius = mw_hypot(&X(pos), &Y(pos));
    radius = mw_hypot(&radius, &Z(pos));

    real_0 expect_derv = radius.value / a;
    real_0 actual_derv = radius.gradient[BARYON_RADIUS_POS];
    real_0 check = mw_abs_0(1.0 - expect_derv/actual_derv);
    if(check > autodiff_thresh)
    {
        failed += 1;
        printf("\t Radial Gradient test failed, |1 + %1f / %1f| = %1f\n", expect_derv, actual_derv, check);
    }

    real_0 expect_hess = 0.0;
    int hess_indx = BARYON_RADIUS_POS*(BARYON_RADIUS_POS+1)/2 + BARYON_RADIUS_POS;
    real_0 actual_hess = radius.hessian[hess_indx];
    check = mw_abs_0(actual_hess - expect_hess);
    if(check > autodiff_thresh)
    {
        failed += 1;
        printf("\t Radial Hessian test failed, |%1f - %1f| = %1f\n", expect_hess, actual_hess, check);
    }

    return failed;
}

static inline int checkVelDerv(mwvector* vel, mwvector* pos, real_0 a, real_0 M)
{
    int failed = 0;
    real radius = mw_hypot(&X(pos), &Y(pos));
    radius = mw_hypot(&radius, &Z(pos));
    real_0 r_s = radius.value;

    real velocity = mw_hypot(&X(pos), &Y(pos));
    velocity = mw_hypot(&velocity, &Z(pos));
    real_0 v_s = velocity.value;

    real_0 v_esc = mw_sqrt_0(2*M / mw_hypot_0(r_s, a));

    real_0 expect_derv_a = -r_s / 24.0 / (r_s*r_s + a*a) * Qfunc(v_s, v_esc);
    real_0 actual_derv_a = velocity.gradient[BARYON_RADIUS_POS];
    real_0 check = mw_abs_0(1.0 - expect_derv_a/actual_derv_a);
    if(check > autodiff_thresh)
    {
        failed += 1;
        printf("\t Velocity Gradient (a) test failed, |1 + %1f / %1f| = %1f\n", expect_derv_a, actual_derv_a, check);
    }

    real_0 expect_derv_M = Qfunc(v_s, v_esc)/ 48.0 / M;
    real_0 actual_derv_M = velocity.gradient[BARYON_MASS_POS];
    check = mw_abs_0(1.0 - expect_derv_M/actual_derv_M);
    if(check > autodiff_thresh)
    {
        failed += 1;
        printf("\t Velocity Gradient (M) test failed, |1 + %1f / %1f| = %1f\n", expect_derv_M, actual_derv_M, check);
    }

    real_0 expect_hess_a2 = r_s / 24.0 / (r_s*r_s + a*a) * (Qfunc(v_s, v_esc)/a + r_s / 24.0 / (r_s*r_s + a*a) * Pfunc(v_s, v_esc));
    int hess_indx = BARYON_RADIUS_POS*(BARYON_RADIUS_POS+1)/2 + BARYON_RADIUS_POS;
    real_0 actual_hess_a2 = velocity.hessian[hess_indx];
    check = mw_abs_0(1.0 - expect_hess_a2/actual_hess_a2);
    if(check > autodiff_thresh)
    {
        failed += 1;
        printf("\t Velocity Hessian (a2) test failed, |1 + %1f / %1f| = %1f\n", expect_hess_a2, actual_hess_a2, check);
    }

    real_0 expect_hess_aM = -r_s / 1152.0 / M /(r_s*r_s + a*a) * Pfunc(v_s, v_esc);
    hess_indx = BARYON_RADIUS_POS*(BARYON_RADIUS_POS+1)/2 + BARYON_MASS_POS;
    real_0 actual_hess_aM = velocity.hessian[hess_indx];
    check = mw_abs_0(1.0 - expect_hess_aM/actual_hess_aM);
    if(check > autodiff_thresh)
    {
        failed += 1;
        printf("\t Velocity Hessian (aM) test failed, |1 + %1f / %1f| = %1f\n", expect_hess_aM, actual_hess_aM, check);
    }

    real_0 expect_hess_M2 = (Pfunc(v_s, v_esc)/48.0 - Qfunc(v_s, v_esc)) / 48.0 / M / M;
    hess_indx = BARYON_MASS_POS*(BARYON_MASS_POS+1)/2 + BARYON_MASS_POS;
    real_0 actual_hess_M2 = velocity.hessian[hess_indx];
    check = mw_abs_0(1.0 - expect_hess_M2/actual_hess_M2);
    if(check > autodiff_thresh)
    {
        failed += 1;
        printf("\t Velocity Hessian (M2) test failed, |1 + %1f / %1f| = %1f\n", expect_hess_M2, actual_hess_M2, check);
    }

    return failed;
}

//make a two component plummer-plummer dwarf and check it
int testDwarfDerivatives()
{
	int failed = 0;
	
	unsigned int numBodies = 500;
        unsigned int halfBodies = numBodies/2;
	mwvector* positions    = mwCalloc(numBodies, sizeof(mwvector));
	mwvector* velocities   = mwCalloc(numBodies, sizeof(mwvector));
	real* masses           = mwCalloc(numBodies, sizeof(real));
	
	//we want the dwarf to be at the origin
	mwvector rshift, vshift;
	rshift.x = ZERO_REAL;
	rshift.y = ZERO_REAL;
	rshift.z = ZERO_REAL;
	vshift.x = ZERO_REAL;
	vshift.y = ZERO_REAL;
	vshift.z = ZERO_REAL;

        Dwarf* comp1       = mwMalloc(sizeof(Dwarf));
	comp1->type        = Plummer;
	comp1->mass        = 12.0;
	comp1->scaleLength = 0.2;

	Dwarf* comp2       = mwMalloc(sizeof(Dwarf));
	comp2->type        = comp1->type;
	comp2->mass        = 0.0001;
	comp2->scaleLength = 20.0;

	dsfmt_t prng;
	dsfmt_init_gen_rand(&prng, 1234); //initialize the random variable

	//Actually generate the dwarf bodies by calling a special version of the actual generation function from nbody_mixeddwarf.c
	nbGenerateMixedDwarfCore_TESTVER(positions, velocities, masses, &prng, numBodies, comp1, comp2, &rshift, &vshift);
	
        for(int i = 0; i < numBodies; i++)
        {
            if(i < halfBodies)
            {
                failed += checkRadDerv(&positions[i], comp1->scaleLength);
                failed += checkVelDerv(&velocities[i], &positions[i], comp1->scaleLength, comp1->mass);
            }
        }

	free(positions);
	free(velocities);
	free(masses);
	free(comp1);
	free(comp2);

	return failed;
}

int main()
{

    int failed = 0;

    failed += testDwarfDerivatives();

    if(failed == 0)
    {
        printf("All tests successful!\n");
    }

    return failed;
}
