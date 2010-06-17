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

/* MAIN: toplevel routine for hierarchical N-body code. */

int main(int argc, char* argv[])
{
    float chisqans = 0.0;

    initNBody(&ctx, &ic, argc, (const char**) argv);

    printf("Reached end of init\n");

    /* FIXME: This is wrong */
    ps.PluMass    = ctx.model.mass;
    ps.r0         = ctx.model.scale_radius;
    ps.Xinit      = ic.position[0];
    ps.Yinit      = ic.position[1];
    ps.Zinit      = ic.position[2];
    ps.sunGCDist  = ic.sunGCDist;
    ps.VXinit     = ic.velocity[0];
    ps.VYinit     = ic.velocity[1];
    ps.VZinit     = ic.velocity[2];
    ps.orbittstop = ctx.model.time_orbit;
    ps.dtorbit    = ctx.model.timestep / 2.0;

    ctx.tstop = ctx.model.time_dwarf;

    ctx.criterion = NEWCRITERION;

    t.rsize = 4.0;

    initoutput(&ctx);

    printContext(&ctx);
    printInitialConditions(&ic);

    printf("arstarstarst: %g %g %g\n", ps.Xinit, ps.Yinit, ps.Zinit);
    printf("varst: %g %g %g\n", ps.VXinit, ps.VYinit, ps.VZinit);


    // Calculate the reverse orbit
    printf("Calculating reverse orbit...");
    integrate();
    printf("done\n");

    printf("Beginning run...\n");
    startrun(&ctx, &st);                 /* set params, input data */

    printf("st.tnow = %g, ctx.tstop = %g, ctx.freq = %g\n",
           st.tnow, ctx.tstop, ctx.freq);

	printf("eps = %f dtnbody = %f\n", ctx.model.eps, ctx.model.timestep);

    printf("tnow init = %g, tstop = %g\n", st.tnow, ctx.tstop);

    printf("tstop = %g, while < %g or %g\n", ctx.tstop,
           ctx.tstop - 1.0 / (1024.0 * ctx.freq),
           ctx.tstop - 1.0 / (1024 * ctx.freq));




    while (st.tnow < ctx.tstop - 1.0 / (1024.0 * ctx.freq)) /* while not past ctx.tstop */
        stepsystem();               /* advance N-body system */

    printf("nstep final = %d\n", st.nstep);
    printf("Step system done\n");
    // Get the likelihood
    //chisqans = chisq();
    //printf("Run finished. chisq = %f\n", chisqans);

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
    generatePlummer();             /* make test model */

    st->nstep = 0;                  /* start counting steps */
    st->tout = st->tnow;            /* schedule first output */
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
        printf("Tree made\n");
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

