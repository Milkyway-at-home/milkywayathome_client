/* ************************************************************************** */
/* nbody_simulation.c: hierarchical N-body code. */
/* Copyright (c) 1993 by Joshua E. Barnes, Honolulu, HI. */
/* It's free because it's yours. */
/* ************************************************************************** */

#include "nbody.h"
#include "util.h"
#include "json_params.h"

/* startrun: startup hierarchical N-body code. */
inline static void startrun(const NBodyCtx* ctx, const InitialConditions* ic, NBodyState* st)
{
    if (ctx->model.nbody < 1)              /* check input value */
        error("startrun: ctx.model.nbody = %d is absurd\n", ctx->model.nbody);

    /* FIXME: Make seed be long anyway */
    srand48((long) ctx->seed);       /* set random generator */
    generatePlummer(ctx, ic, st);    /* make test model */

    st->nstep = 0;                  /* start counting steps */
    st->tout  = st->tnow;           /* schedule first output */
    st->tree.rsize = ctx->tree_rsize;

}

/* stepsystem: advance N-body system one time-step. */
inline static void stepsystem(const NBodyCtx* ctx, NBodyState* st)
{
    bodyptr p;
    vector dvel, dpos;
    int nfcalc;          /* count force calculations */
    int n2bcalc;         /* count body-body interactions */
    int nbccalc;         /* count body-cell interactions */

    const bodyptr bodytab = st->bodytab;
    const bodyptr endp    = bodytab + ctx->model.nbody;
    const real dt         = 1.0 / ctx->freq;         /* set basic time-step */

    if (st->nstep == 0)                 /* about to take 1st step? */
    {
        printf("Building tree...Starting Nbody simulation...\n");
        maketree(ctx, st, bodytab, ctx->model.nbody);       /* build tree structure */
        printf("Tree made\n");
        nfcalc = n2bcalc = nbccalc = 0;     /* zero counters */
        for (p = bodytab; p < endp; p++)
        {
            /* loop over all bodies */
            hackgrav(ctx, st, p, Mass(p) > 0.0);     /* get force on each */
            ++nfcalc;                       /* count force calcs */
            n2bcalc += st->n2bterm;          /* and 2-body terms */
            nbccalc += st->nbcterm;          /* and body-cell terms */
        }
        output(ctx, st);                /* do initial output */
    }

    for (p = bodytab; p < endp; p++)    /* loop over all bodies */
    {
        MULVS(dvel, Acc(p), 0.5 * dt);  /* get velocity increment */
        ADDV(Vel(p), Vel(p), dvel);     /* advance v by 1/2 step */
        MULVS(dpos, Vel(p), dt);        /* get positon increment */
        ADDV(Pos(p), Pos(p), dpos);     /* advance r by 1 step */
    }

    maketree(ctx, st, bodytab, ctx->model.nbody);     /* build tree structure */
    nfcalc = n2bcalc = nbccalc = 0;          /* zero counters */
    for (p = st->bodytab; p < endp; p++)     /* loop over bodies */
    {
        hackgrav(ctx, st, p, Mass(p) > 0.0);   /* get force on each */
        ++nfcalc;                     /* count force calcs */
        n2bcalc += st->n2bterm;           /* and 2-body terms */
        nbccalc += st->nbcterm;           /* and body-cell terms */
    }

    for (p = bodytab; p < endp; p++) /* loop over all bodies */
    {
        MULVS(dvel, Acc(p), 0.5 * dt);          /* get velocity increment */
        ADDV(Vel(p), Vel(p), dvel);             /* advance v by 1/2 step */
    }

    st->nstep++;           /* count another time step */
    st->tnow += dt;        /* finally, advance time */
    output(ctx, st);             /* do major or minor output */
}


void runSystem(const NBodyCtx* ctx, const InitialConditions* ic, NBodyState* st)
{
    const real tstop = ctx->model.time_dwarf - 1.0 / (1024.0 * ctx->freq);

    startrun(ctx, ic, st);

    while (st->tnow < tstop)
        stepsystem(ctx, st);               /* advance N-body system */
}

