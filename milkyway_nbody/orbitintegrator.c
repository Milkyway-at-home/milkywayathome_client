/* Simple orbit integrator in user-defined potential
    Written for BOINC Nbody
    willeb 10 May 2010 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "defs.h"
#include "code.h"

//#define X 0
//#define Y 1
//#define Z 2

void integrate()
{
    vector acc, v, x;
    real time;
    int i;

    // Set the initial conditions
    x[0] = ps.Xinit;
    x[1] = ps.Yinit;
    x[2] = ps.Zinit;
    v[0] = -ps.Xinit;
    v[1] = -ps.Yinit;
    v[2] = -ps.Zinit;

    //printf("%f %f %f %f %f %f\n", ps.Xinit, ps.Yinit, ps.Zinit, ps.Xinit, ps.Yinit, ps.Zinit);

    // Get the initial acceleration
    acceleration(x, acc);

    // Loop through time
    for (time = 0; time <= ps.orbittstop; time += ps.dtorbit)
    {
        // Update the velocities and positions
        for (i = 0; i < 3; ++i)
        {
            v[i] += acc[i] * ps.dtorbit;
            x[i] += v[i] * ps.dtorbit;
            //printf(" %f", x[i]);
        }

        /*
        for (i = 0; i <= 2; i++)
        {
            //printf(" %f", v[i]);
        }
        */
        //printf("orbittime = %f\n", time);

        // Compute the new acceleration
        acceleration(x, acc);
    }

    // Report the final values (don't forget to reverse the velocities)
    ps.XC = x[0];
    ps.YC = x[1];
    ps.ZC = x[2];
    ps.XC = -v[0];
    ps.YC = -v[1];
    ps.ZC = -v[2];

}

void acceleration(real* pos, real* acc)
{
    // The external potential

    real miya_ascal = 6.5;
    real miya_bscal = 0.26;
    real miya_mass = 4.45865888E5;
    real plu_rc = 0.7;
    real bulge_mass = 1.52954402E5;
    real vhalo = 73;
    real haloq = 1.0;
    real halod = 12.0;

    real apar, qpar, spar, ppar, lpar, rpar, rcyl;
    real vhalosqr2, rppPar, spar15;

    rcyl = sqrt(pos[0] * pos[0] + pos[1] * pos[1]);
    qpar = sqrt(pos[2] * pos[2] + miya_bscal * miya_bscal);
    apar = miya_ascal + qpar;
    spar = (pos[0] * pos[0]) + (pos[1] * pos[1]) + ((miya_ascal + qpar) * (miya_ascal + qpar));
    ppar = sqrt ((pos[0] * pos[0]) + (pos[1] * pos[1]) + (pos[2] * pos[2])) + plu_rc;
    rpar = ppar - plu_rc;
    lpar = (rcyl * rcyl) + ((pos[2] / haloq) * (pos[2] / haloq)) + (halod * halod);


    vhalosqr2 = 2.0 * vhalo * vhalo;
    rppPar = rpar * ppar * ppar;
    spar15 = pow(spar, 1.5);


    /* This function returns the acceleration vector */
    acc[0] = - ( ( (vhalosqr2 * pos[0]) / lpar )
                 + ((bulge_mass * pos[0]) / rppPar )
                 + ((miya_mass * pos[0]) / spar15 ) );

    acc[1] = - ( ( (vhalosqr2 * pos[1]) / lpar )
                 + ((bulge_mass * pos[1]) / rppPar )
                 + ((miya_mass * pos[1]) / spar15 ) );

    acc[2] = - ( ( (vhalosqr2 * pos[2]) / (haloq * haloq * lpar) )
                 + ((bulge_mass * pos[2]) / rppPar)
                 + ((miya_mass * pos[2] * apar) / (qpar * spar15) ) );
}

