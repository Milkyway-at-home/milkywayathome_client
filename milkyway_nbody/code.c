/* ************************************************************************** */
/* CODE.C: hierarchical N-body code. */
/* */
/* Copyright (c) 1993 by Joshua E. Barnes, Honolulu, HI. */
/* It's free because it's yours. */
/* ************************************************************************** */

#define extern /* nada */

#include "code.h"
#include <time.h>
#include <stdlib.h>

/*  * Default values for input parameters.
 */

const char* defv[] =
{

    /* file names for input/output */
    "in=",          /* Input file with initial conditions */
    "out=",         /* Output file of N-body frames */

    /* params to control N-body integration */
    "dtime=0.03125",        /* Integration time-step */
    "eps=0.025",        /* Potential softening parameter */
    "theta=1.0",        /* Cell subdivision tolerence */
    "usequad=false",        /* If true, use quad moments */
    "rsize=4.0",        /* Size of initial root cell */
    "options=",         /* Various control options */
    "tstop=2.0",        /* Time to stop integration */
    "dtout=0.25",       /* Data output interval */

    /* params used if no input specified to make a Plummer model */
    "nbody=1024",       /* Number of particles for test run */
    "seed=123",         /* Random number seed for test run */

    NULL,
};

char* headline = "Hierarchical N-body Code";   /* default id for run */

static void startrun(void);          /* initialize system state */
static void stepsystem(void);            /* advance by one time-step */
static void testdata(void);          /* generate test data */
static void pickshell(vector, real);     /* pick point on shell */

/*  * MAIN: toplevel routine for hierarchical N-body code.
 */

int main(int argc, char* argv[])
{

    printf("Initializing parameters...");
    static real dtnbody;
    static real kmax;
    float chisqans = 0.0;

    PluMass = atof(argv[1]);    //mass
    r0 = atof(argv[2]);     //radius
    orbittstop = atof(argv[3]); //how far back the orbit goes
    tstop = atof(argv[4]);      //evolution time
    nbody = atoi(argv[5]);      //number of bodies

    // Constants
    options = "";
    outfile = "out";
    rsize = 4;

    // Calculate starting galactic coordinates
    // The l,b,r and vx, vy, vz are hard coded for the Orphan project
    // In the future, these will change
    lstart = 218.0;
    bstart = 53.5;

    // From the vhalo = 73 model result from Newberg et al 2009
    Rstart = 28.6;
    VXinit = -156;
    VYinit = 79;
    VZinit = 107;

    Xinit = Rstart * cos(lstart * 3.14159 / 180.0) * cos(bstart * 3.14159 / 180.0) - 8.0; // 8.0 is sun-gc distance (TODO: make par)
    Yinit = Rstart * sin(lstart * 3.14159 / 180.0) * cos(bstart * 3.14159 / 180.0);
    Zinit = Rstart * sin(bstart * 3.14159 / 180.0);

    eps = r0 / (10 * sqrt((real)nbody));
    dtnbody = (1 / 10.0) * (1 / 10.0) * sqrt(((4 / 3) * 3.14159 * r0 * r0 * r0) / (PluMass));
    //dtnbody = pow(2.718,log(0.5)*kmax);
    freq = 1.0 / dtnbody;
    freqout = freq;

    //printf("eps = %f dtnbody = %f\n", eps, dtnbody);
    printf("done\n");

    // Calculate the reverse orbit (use dt = dtnbody/2 to make sure we get enough orbit precision, this puts the results in XC, YC, ZC, VXC, VYC, VZC, which is then used by the testdata routine to do the shift, be aware: changing the mass changes the orbit results, but this is OK)
    printf("Calculating reverse orbit...");
    dtorbit = dtnbody / 2.0;
    integrate();
    printf("done\n");

    printf("Beginning run...\n");
    startrun();                 /* set params, input data */
    initoutput();               /* begin system output */
    while (tnow < tstop - 1.0 / (1024 * freq)) /* while not past tstop */
        stepsystem();               /* advance N-body system */

    // Get the likelihood
    chisqans = chisq();
    printf("Run finished. chisq = %f\n", chisqans);

    stopoutput();               /* finish up output */

    return 0;
}

/*  * STARTRUN: startup hierarchical N-body code.
 */

static void startrun(void)
{
    if (nbody < 1)              /* check input value */
        error("startrun: nbody = %d is absurd\n", nbody);

    //srand48((long) time(NULL));   /* set random generator */
    srand48((long) 0.0);    /* set random generator */
    testdata();             /* make test model */

    nstep = 0;                  /* start counting steps */
    tout = tnow;                /* schedule first output */
}

/*  * TESTDATA: generate Plummer model initial conditions for test runs,
 * scaled to units such that M = -4E = G = 1 (Henon, Hegge, etc).
 * See Aarseth, SJ, Henon, M, & Wielen, R (1974) Astr & Ap, 37, 183.
 */

#define MFRAC  0.999                /* cut off 1-MFRAC of mass */

static void testdata(void)
{
    printf("Initializing plummer model...");
    real rsc, vsc, r, v, x, y;
    vector cmr, cmv;
    vector rshift, vshift, scaledrshift, scaledvshift;

    // The coordinates to shift the plummer sphere by
    (rshift)[0] = XC;
    (rshift)[1] = YC;
    (rshift)[2] = ZC;
    (vshift)[0] = VXC;
    (vshift)[1] = VYC;
    (vshift)[2] = VZC;

    printf("Shifting plummer sphere to r = (%f, %f, %f) v = (%f, %f, %f)...\n", rshift[0], rshift[1], rshift[2], vshift[0], vshift[1], vshift[2]);

    bodyptr p;

    headline = "Hierarchical code: Plummer model";
    /* supply default headline */
    tnow = 0.0;                 /* reset elapsed model time */
    bodytab = (bodyptr) allocate(nbody * sizeof(body));
    /* alloc space for bodies */
    rsc = r0;               /* set length scale factor */
    vsc = rsqrt(PluMass / rsc);         /* and recip. speed scale */
    CLRV(cmr);                  /* init cm pos, vel */
    CLRV(cmv);
    CLRV(scaledrshift);
    CLRV(scaledvshift);
    MULVS(scaledrshift, rshift, rsc);   /* Multiply shift by scale factor */
    MULVS(scaledvshift, vshift, vsc);   /* Multiply shift by scale factor */

    for (p = bodytab; p < bodytab + nbody; p++) /* loop over particles */
    {
        Type(p) = BODY;             /* tag as a body */
        Mass(p) = PluMass / (real)nbody;            /* set masses equal */
        r = 1 / rsqrt(rpow(xrandom(0.0, MFRAC), /* pick r in struct units */
                           -2.0 / 3.0) - 1);
        pickshell(Pos(p), rsc * r);     /* pick scaled position */
        ADDV(Pos(p), Pos(p), rshift);       /* move the position */
        ADDV(cmr, cmr, Pos(p));         /* add to running sum */
        do                      /* select from fn g(x) */
        {
            x = xrandom(0.0, 1.0);      /* for x in range 0:1 */
            y = xrandom(0.0, 0.1);      /* max of g(x) is 0.092 */
        }
        while (y > x * x * rpow(1 - x * x, 3.5)); /* using von Neumann tech */
        v = rsqrt(2.0) * x / rpow(1 + r * r, 0.25); /* find v in struct units */
        pickshell(Vel(p), vsc * v);     /* pick scaled velocity */
        ADDV(Vel(p), Vel(p), vshift);       /* move the velocity */
        ADDV(cmv, cmv, Vel(p));         /* add to running sum */
    }
    DIVVS(cmr, cmr, (real) nbody);      /* normalize cm coords */
    DIVVS(cmv, cmv, (real) nbody);

    printf("done\n");
}

/*  * PICKSHELL: pick a random point on a sphere of specified radius.
 */

static void pickshell(vector vec, real rad)
{
    int k;
    real rsq, rsc;

    do                      /* pick point in NDIM-space */
    {
        for (k = 0; k < NDIM; k++)      /* loop over dimensions */
            vec[k] = xrandom(-1.0, 1.0);        /* pick from unit cube */
        DOTVP(rsq, vec, vec);           /* compute radius squared */
    }
    while (rsq > 1.0);                      /* reject if outside sphere */
    rsc = rad / rsqrt(rsq);         /* compute scaling factor */
    MULVS(vec, vec, rsc);           /* rescale to radius given */
}

/*  * STEPSYSTEM: advance N-body system one time-step.
 */

static void stepsystem(void)
{
    bodyptr p;
    real dt;
    vector dvel, dpos;
    if (nstep == 0)                 /* about to take 1st step? */
    {
        printf("Building tree...Starting Nbody simulation...\n");
        maketree(bodytab, nbody);       /* build tree structure */
        nfcalc = n2bcalc = nbccalc = 0;     /* zero counters */
        for (p = bodytab; p < bodytab + nbody; p++)
        {
            /* loop over all bodies */
            hackgrav(p, Mass(p) > 0.0);     /* get force on each */
            nfcalc++;               /* count force calcs */
            n2bcalc += n2bterm;         /* and 2-body terms */
            nbccalc += nbcterm;         /* and body-cell terms */
        }
        output();               /* do initial output */
    }
    dt = 1.0 / freq;                /* set basic time-step */
    for (p = bodytab; p < bodytab + nbody; p++) /* loop over all bodies */
    {
        MULVS(dvel, Acc(p), 0.5 * dt);      /* get velocity increment */
        ADDV(Vel(p), Vel(p), dvel);     /* advance v by 1/2 step */
        MULVS(dpos, Vel(p), dt);        /* get positon increment */
        ADDV(Pos(p), Pos(p), dpos);     /* advance r by 1 step */
    }
    maketree(bodytab, nbody);           /* build tree structure */
    nfcalc = n2bcalc = nbccalc = 0;     /* zero counters */
    for (p = bodytab; p < bodytab + nbody; p++) /* loop over bodies */
    {
        hackgrav(p, Mass(p) > 0.0);     /* get force on each */
        nfcalc++;               /* count force calcs */
        n2bcalc += n2bterm;         /* and 2-body terms */
        nbccalc += nbcterm;         /* and body-cell terms */
    }
    for (p = bodytab; p < bodytab + nbody; p++) /* loop over all bodies */
    {
        MULVS(dvel, Acc(p), 0.5 * dt);          /* get velocity increment */
        ADDV(Vel(p), Vel(p), dvel);             /* advance v by 1/2 step */
    }
    nstep++;                    /* count another time step */
    tnow = tnow + dt;               /* finally, advance time */
    output();                   /* do major or minor output */
}
