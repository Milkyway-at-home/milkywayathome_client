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
    mwvector x_cutoff = mw_vec(5, 5, 5);

    //declare acceleration vectors, tolerance, functions
    mwvector A1, A2, A3, A4;
    real delta;
    int F1 = 1;
    int F2 = 2;
    int F3 = 3;

    //call LMC functions
    A1 = LMCAcceleration(F1, x, LMCx, LMCmass, 8.32, 1.0);//Plummer LMC
    A2 = LMCAcceleration(F2, x, LMCx, LMCmass, 5.45, 1.0);//Henrquist LMC
    A3 = LMCAcceleration(F3, x_cutoff, LMCx, 2336886.4425964956, 39.8, 16.6);//before cutoff 
    A4 = LMCAcceleration(F3, x, LMCx, 2336886.4425964956, 39.8, 16.6);//after cutoff

    //tolerance value
    delta = 0.1;

    if (fabs(A1.x - (-285.34)) > delta)
	{
	    printf("Plummer LMC test failed: accel=%f\n",A1.x);
            failed = 1;
	}
    if (fabs(A2.x - (-225.42)) > delta)
	{
	    printf("Hernquist LMC test failed: accel=%f\n",A2.x);
            failed = 1;
	}
    if (fabs(A3.x - (-574.52)) > delta)
        {
	    printf("cutoff Hernquist (before) test failed: accel=%f\n",A3.x);
            failed = 1;
	}
    if (fabs(A4.x - (-392.69)) > delta)
	{
	    printf("cutoff Hernquist (after) test failed: accel=%f\n",A4.x);
            failed = 1;
	}

    return failed;

}

/* Main function to run the test */
int main()
{   
    int failed = test_LMC_accel();

    if(failed == 0)
    {
        printf("All tests successful!\n");
    }

    return failed;
}


