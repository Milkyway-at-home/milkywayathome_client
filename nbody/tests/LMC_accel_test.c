/* This program tests LMC acceleration fnctions */

#include <math.h>
#include "nbody_potential.h"

/* Computes the acceleration of LMC functions */
int test_LMC_accel()
{
    int failed = 0;

    //test parameters
    mwvector x = mw_vec(10, 10, 10);
    mwvector LMCx = mw_vec(0, 0, 0);
    real LMCmass = 202439.64970382853;

    //call LMC functions
    A1 = LMCAcceleration(1, x, LMCx, LMCmass, 8.32, 1.0);//Plummer LMC
    A2 = LMCAcceleration(2, x, LMCx, LMCmass, 5.45, 1.0);//Henrquist LMC
    A3 = LMCAcceleration(3, mw_vec(5, 5, 5), LMCx, 2336886.4425964956, 39.8, 16.6);//before cutoff 
    A4 = LMCAcceleration(3, x, LMCx, 2336886.4425964956, 39.8, 16.6);//after cutoff

    //tolerance value
    delta = 0.1

    if (fabs(A1.x - (-285.34)) < delta)
	{
	    printf("Plummer LMC test failed: accel=%f\n",A1);
            failed = 1;
	}
    if (fabs(A2.x - (-225.42)) < delta)
	{
	    printf("Hernquist LMC test failed: accel=%f\n",A2);
            failed = 1;
	}
    if (fabs(A3.x - (-574.52)) < delta)
        {
	    printf("cutoff Hernquist (before) test failed: accel=%f\n",A3)
            failed = 1;
	}
    if (fabs(A4.x - (-1319.8)) < delta)
	{
	    printf("cutoff Hernquist (after) test failed: accel=%f\n",A4)
            failed = 1;
	}

    printf("\tPassed LMC accel test\n");
    fflush(stdout);

    return failed;

}


