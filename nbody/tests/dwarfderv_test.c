/* This program tests the derivative calculations of a 2-component dwarf galaxy.*/

#include "milkyway_util.h"
#include "nbody_mixeddwarf.h"
#include "dSFMT.h"
#include <stdio.h>

#define autodiff_thresh (0.0)

static inline int checkRadDerv(mwvector* pos, real_0 a)
{
    printf("POS = [%.15f, %.15f, %.15f]\n", X(pos).value, Y(pos).value, Z(pos).value);
    int failed = 0;
    real radius = mw_hypot(&X(pos), &Y(pos));
    radius = mw_hypot(&radius, &Z(pos));

    real_0 expect_derv = radius.value / a;
    real_0 actual_derv = radius.gradient[BARYON_RADIUS_POS];
    real_0 check = mw_abs_0(1.0 - expect_derv/actual_derv);
    if(check > autodiff_thresh)
    {
        failed += 1;
        printf("    Radial Gradient test failed, |1 - %.15f / %.15f| = %.15f\n", expect_derv, actual_derv, check);
    }

    real_0 expect_hess = 0.0;
    int hess_indx = BARYON_RADIUS_POS*(BARYON_RADIUS_POS+1)/2 + BARYON_RADIUS_POS;
    real_0 actual_hess = radius.hessian[hess_indx];
    check = mw_abs_0(actual_hess - expect_hess);
    if(check > autodiff_thresh)
    {
        failed += 1;
        printf("    Radial Hessian test failed, |%.15f - %.15f| = %.15f\n", expect_hess, actual_hess, check);
    }

    return failed;
}

static inline int checkVelDerv(mwvector* vel, mwvector* pos, real_0 a, real_0 M)
{
    printf("VEL = [%.15f, %.15f, %.15f]\n", X(vel).value, Y(vel).value, Z(vel).value);
    int failed = 0;
    real radius = mw_hypot(&X(pos), &Y(pos));
    radius = mw_hypot(&radius, &Z(pos));
    real_0 r_s = radius.value;

    real velocity = mw_hypot(&X(vel), &Y(vel));
    velocity = mw_hypot(&velocity, &Z(vel));
    real_0 v_s = velocity.value;

    real_0 v_esc = mw_sqrt_0(2.0*M / mw_hypot_0(r_s, a));

    real_0 expect_derv_a = -v_s / 2.0 / a;
    real_0 actual_derv_a = velocity.gradient[BARYON_RADIUS_POS];
    real_0 check = mw_abs_0(1.0 - expect_derv_a/actual_derv_a);
    if(check > autodiff_thresh)
    {
        failed += 1;
        printf("    Velocity Gradient (a) test failed, |1 - %.15f / %.15f| = %.15f\n", expect_derv_a, actual_derv_a, check);
    }

    real_0 expect_derv_M = v_s / 2.0 / M;
    real_0 actual_derv_M = velocity.gradient[BARYON_MASS_POS];
    check = mw_abs_0(1.0 - expect_derv_M/actual_derv_M);
    if(check > autodiff_thresh)
    {
        failed += 1;
        printf("    Velocity Gradient (M) test failed, |1 - %.15f / %.15f| = %.15f\n", expect_derv_M, actual_derv_M, check);
    }

    real_0 expect_hess_a2 = 0.75 * v_s / a / a;
    int hess_indx = BARYON_RADIUS_POS*(BARYON_RADIUS_POS+1)/2 + BARYON_RADIUS_POS;
    real_0 actual_hess_a2 = velocity.hessian[hess_indx];
    check = mw_abs_0(1.0 - expect_hess_a2/actual_hess_a2);
    if(check > autodiff_thresh)
    {
        failed += 1;
        printf("    Velocity Hessian (a2) test failed, |1 - %.15f / %.15f| = %.15f\n", expect_hess_a2, actual_hess_a2, check);
    }

    int eff_i;
    int eff_j;
    if(BARYON_RADIUS_POS<BARYON_MASS_POS)
    {
        eff_i = BARYON_MASS_POS;
        eff_j = BARYON_RADIUS_POS;
    }
    else
    {
        eff_i = BARYON_RADIUS_POS;
        eff_j = BARYON_MASS_POS;
    }
    real_0 expect_hess_aM = -0.25 * v_s / a / M;
    hess_indx = eff_i*(eff_i+1)/2 + eff_j;
    real_0 actual_hess_aM = velocity.hessian[hess_indx];
    check = mw_abs_0(1.0 - expect_hess_aM/actual_hess_aM);
    if(check > autodiff_thresh)
    {
        failed += 1;
        printf("    Velocity Hessian (aM) test failed, |1 - %.15f / %.15f| = %.15f\n", expect_hess_aM, actual_hess_aM, check);
    }

    real_0 expect_hess_M2 = -0.25 * v_s / M / M;
    hess_indx = BARYON_MASS_POS*(BARYON_MASS_POS+1)/2 + BARYON_MASS_POS;
    real_0 actual_hess_M2 = velocity.hessian[hess_indx];
    check = mw_abs_0(1.0 - expect_hess_M2/actual_hess_M2);
    if(check > autodiff_thresh)
    {
        failed += 1;
        printf("    Velocity Hessian (M2) test failed, |1 - %.15f / %.15f| = %.15f\n", expect_hess_M2, actual_hess_M2, check);
    }

    return failed;
}

//make a two component plummer-plummer dwarf and check it
int testDwarfDerivatives()
{
	int failed = 0;
	
	unsigned int numBodies = 4;
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
	comp2->mass        = REAL_EPSILON*10.0;
	comp2->scaleLength = 100.0;

	dsfmt_t prng;
	dsfmt_init_gen_rand(&prng, 1234); //initialize the random variable

	//Actually generate the dwarf bodies by calling a special version of the actual generation function from nbody_mixeddwarf.c
	nbGenerateMixedDwarfCore_TESTVER(positions, velocities, masses, &prng, numBodies, comp1, comp2, &rshift, &vshift, FALSE);
        printf("Dwarf Generated!\n");
	
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
