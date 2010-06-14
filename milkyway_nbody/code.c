/* ************************************************************************** */
/* CODE.C: hierarchical N-body code. */
/* */
/* Copyright (c) 1993 by Joshua E. Barnes, Honolulu, HI. */
/* It's free because it's yours. */
/* ************************************************************************** */

#include "code.h"
#include "defs.h"
#include <time.h>
#include <stdlib.h>
#include <string.h>

#include "json_params.h"

NBodyParams ps = { 0, };
Tree t        = EMPTY_TREE;
NBodyCtx ctx  = EMPTY_CTX;
NBodyState st = EMPTY_STATE;
InitialConditions ic = EMPTY_INITIAL_CONDITIONS;


char* headline = "Hierarchical N-body Code";   /* default id for run */

static void startrun(NBodyCtx*, NBodyState*);           /* initialize system state */
static void stepsystem(void);         /* advance by one time-step */
static void testdata(void);           /* generate test data */
static void pickshell(vector, real);  /* pick point on shell */

/* MAIN: toplevel routine for hierarchical N-body code. */

int main(int argc, char* argv[])
{
    float chisqans = 0.0;

    initNBody(&ctx, &ic, argc, (const char**) argv);
    printf("Reached end of init\n");
    exit(EXIT_SUCCESS);


    // Calculate the reverse orbit
    printf("Calculating reverse orbit...");
    integrate();
    printf("done\n");

    printf("Beginning run...\n");
    startrun(&ctx, &st);                 /* set params, input data */

    while (st.tnow < ctx.tstop - 1.0 / (1024 * ctx.freq)) /* while not past ctx.tstop */
        stepsystem();               /* advance N-body system */

    // Get the likelihood
    chisqans = chisq();
    printf("Run finished. chisq = %f\n", chisqans);

    nbody_ctx_destroy(&ctx);               /* finish up output */

    return 0;
}

/* STARTRUN: startup hierarchical N-body code. */

static void startrun(NBodyCtx* ctx, NBodyState* st)
{
    if (ctx->model.nbody < 1)              /* check input value */
        error("startrun: ctx.model.nbody = %d is absurd\n", ctx->model.nbody);

    //srand48((long) time(NULL));   /* set random generator */
    srand48((long) 0.0);    /* set random generator */
    testdata();             /* make test model */

    st->nstep = 0;                  /* start counting steps */
    st->tout = st->tnow;            /* schedule first output */
}

/* TESTDATA: generate Plummer model initial conditions for test runs,
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
    (rshift)[0] = ps.XC;
    (rshift)[1] = ps.YC;
    (rshift)[2] = ps.ZC;
    (vshift)[0] = ps.XC;
    (vshift)[1] = ps.YC;
    (vshift)[2] = ps.ZC;

    printf("Shifting plummer sphere to r = (%f, %f, %f) v = (%f, %f, %f)...\n",
           rshift[0],
           rshift[1],
           rshift[2],
           vshift[0],
           vshift[1],
           vshift[2]);

    bodyptr p;

    /* TODO: This doesn't belong here*/
    ctx.headline = strdup("Hierarchical code: Plummer model");
    /* supply default ctx.headline */
    st.tnow = 0.0;                 /* reset elapsed model time */
    st.bodytab = (bodyptr) allocate(ctx.model.nbody * sizeof(body));
    /* alloc space for bodies */
    rsc = ps.r0;               /* set length scale factor */
    vsc = rsqrt(ps.PluMass / rsc);         /* and recip. speed scale */
    CLRV(cmr);                  /* init cm pos, vel */
    CLRV(cmv);
    CLRV(scaledrshift);
    CLRV(scaledvshift);
    MULVS(scaledrshift, rshift, rsc);   /* Multiply shift by scale factor */
    MULVS(scaledvshift, vshift, vsc);   /* Multiply shift by scale factor */

    for (p = st.bodytab; p < st.bodytab + ctx.model.nbody; p++) /* loop over particles */
    {
        Type(p) = BODY;             /* tag as a body */
        Mass(p) = ps.PluMass / (real)ctx.model.nbody;            /* set masses equal */
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
    DIVVS(cmr, cmr, (real) ctx.model.nbody);      /* normalize cm coords */
    DIVVS(cmv, cmv, (real) ctx.model.nbody);

    printf("done\n");
}

/* PICKSHELL: pick a random point on a sphere of specified radius. */

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

/* STEPSYSTEM: advance N-body system one time-step. */

static void stepsystem()
{
    bodyptr p;
    real dt;
    vector dvel, dpos;
    int nfcalc;          /* count force calculations */
    int n2bcalc;         /* count body-body interactions */
    int nbccalc;         /* count body-cell interactions */


    if (st.nstep == 0)                 /* about to take 1st step? */
    {
        printf("Building tree...Starting Nbody simulation...\n");
        maketree(st.bodytab, ctx.model.nbody);       /* build tree structure */
        nfcalc = n2bcalc = nbccalc = 0;     /* zero counters */
        for (p = st.bodytab; p < st.bodytab + ctx.model.nbody; p++)
        {
            /* loop over all bodies */
            hackgrav(p, Mass(p) > 0.0);     /* get force on each */
            nfcalc++;               /* count force calcs */
            n2bcalc += st.n2bterm;         /* and 2-body terms */
            nbccalc += st.nbcterm;         /* and body-cell terms */
        }
        output();               /* do initial output */
    }
    dt = 1.0 / ctx.freq;                /* set basic time-step */
    for (p = st.bodytab; p < st.bodytab + ctx.model.nbody; p++) /* loop over all bodies */
    {
        MULVS(dvel, Acc(p), 0.5 * dt);      /* get velocity increment */
        ADDV(Vel(p), Vel(p), dvel);     /* advance v by 1/2 step */
        MULVS(dpos, Vel(p), dt);        /* get positon increment */
        ADDV(Pos(p), Pos(p), dpos);     /* advance r by 1 step */
    }

    maketree(st.bodytab, ctx.model.nbody);           /* build tree structure */
    nfcalc = n2bcalc = nbccalc = 0;     /* zero counters */
    for (p = st.bodytab; p < st.bodytab + ctx.model.nbody; p++) /* loop over bodies */
    {
        hackgrav(p, Mass(p) > 0.0);     /* get force on each */
        nfcalc++;               /* count force calcs */
        n2bcalc += st.n2bterm;         /* and 2-body terms */
        nbccalc += st.nbcterm;         /* and body-cell terms */
    }
    for (p = st.bodytab; p < st.bodytab + ctx.model.nbody; p++) /* loop over all bodies */
    {
        MULVS(dvel, Acc(p), 0.5 * dt);          /* get velocity increment */
        ADDV(Vel(p), Vel(p), dvel);             /* advance v by 1/2 step */
    }
    st.nstep++;                 /* count another time step */
    st.tnow = st.tnow + dt;     /* finally, advance time */
    output();                   /* do major or minor output */
}

