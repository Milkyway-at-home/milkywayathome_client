/*****************************************************\
* Test for the the plummer momentum offset module     *
* - see nbody_plummer_momentum.c for details          *
*                                                     *
* Aurora Chen - Apr. 2025 - Written for MilkyWay@Home *
\*****************************************************/

#include "nbody_priv.h"
#include "nbody_plummer_momentum.h"

int main() 
{
    const real errLimit = 0.00001;
    real a, M0, velScale, velExp;

    int test_fails = 0;

    // SMC
    mwvector pos = mw_vec(15.9462, -37.5198, -43.6455);
    mwvector vel = mw_vec(21.99, -201.36, 171.25);
    
    M0 = 29241.283;
    a = 2.9;
    real velScaleExp[7] = {0.0, 0.0, 1.00102375871, 1.00161587857, 1.00211796381, 1.00254456226, 1.00291232756};

    for (int k = 2; k <= 6; k++)
    {
        velScale = velocityAdj (a, pos, vel, M0, (real)k);
        velExp = velScaleExp[k];
        if (mw_abs(velScale - velExp) > errLimit)
        {
            mw_printf("SMC: k = %.2f\n", (real)k);
            mw_printf("Velocity offset factor: %.9fx\n", velScale);
            mw_printf("|Error| = %.9f > %.9f\n", mw_abs(velScale - velExp), errLimit);
            test_fails += 1;
        }
    }

    // Sagittarius Dwarf
    mwvector pos2 = mw_vec(15.9872, 2.3527, -6.1609);
    mwvector vel2 = mw_vec(223.97, -5.34, 185.78);

    M0 = 1799.464;
    a = 1.53;
    real velScaleExp2[7] = {0.0, 0.0, 1.00339852947, 1.0053641653, 1.00703091692, 1.00844707818, 1.00966793345};

    for (int k = 2; k <= 6; k++)
    {
        velScale = velocityAdj (a, pos2, vel2, M0, (real)k);
        velExp = velScaleExp2[k];
        if (mw_abs(velScale - velExp) > errLimit)
        {
            mw_printf("Sgr: k = %.2f\n", (real)k);
            mw_printf("Velocity offset factor: %.9fx\n", velScale);
            mw_printf("|Error| = %.9f > %.9f\n", mw_abs(velScale - velExp), errLimit);
            test_fails += 1;
        }
    }
    return test_fails;
}
