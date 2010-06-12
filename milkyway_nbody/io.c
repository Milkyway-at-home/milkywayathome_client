/* ************************************************************************** */
/* IO.C: I/O routines for export version of hierarchical N-body code. */
/* Public routines: inputdata(), inictx.toutput(), stopoutput(), output(). */
/* */
/* Copyright (c) 1993 by Joshua E. Barnes, Honolulu, HI. */
/* It's free because it's yours. */
/* ************************************************************************** */

#include "code.h"

static void diagnostics(void);
static void in_int(FILE*, int*);
static void in_real(FILE*, real*);
static void in_vector(FILE*, vector);
static void out_int(FILE*, int);
static void out_real(FILE*, real);
static void out_vector(FILE*, vector);
static void out_2vectors(FILE*, vector, vector);
static void printvec(char*, vector);

/* INPUTDATA: read initial conditions from input file. */

#if 0
void inputdata(void)
{
    FILE* instr;
    static char headbuf[128];
    int ndim;
    bodyptr p;

    instr = fopen(ctx.infile, "r");         /* open input FILE* */
    if (instr == NULL)
        error("inputdata: cannot find file %s\n", ctx.infile);
    sprintf(headbuf, "Hierarchical code: input file %s", ctx.infile);
    ctx.headline = headbuf;
    in_int(instr, &ctx.nbody);
    if (ctx.nbody < 1)
        error("inputdata: ctx.nbody = %d is absurd\n", ctx.nbody);
    in_int(instr, &ndim);
    if (ndim != NDIM)
        error("inputdata: ndim = %d is absurd\n", ndim);
    in_real(instr, &st.tnow);
    st.bodytab = (bodyptr) allocate(ctx.nbody * sizeof(body));
    for (p = st.bodytab; p < st.bodytab + ctx.nbody; p++) /* loop over new bodies */
        Type(p) = BODY;             /* init body type */
    for (p = st.bodytab; p < st.bodytab + ctx.nbody; p++)
        in_real(instr, &Mass(p));
    for (p = st.bodytab; p < st.bodytab + ctx.nbody; p++)
        in_vector(instr, Pos(p));
    for (p = st.bodytab; p < st.bodytab + ctx.nbody; p++)
        in_vector(instr, Vel(p));
    fclose(instr);              /* close input FILE* */
}
#endif

/* INITOUTPUT: initialize output routines. */

void initoutput(NBodyCtx* ctx)
{
    if (ctx->outfilename)                       /* output file specified? */
    {
        ctx->outfile = fopen(ctx->outfilename, "w");           /* setup output FILE* */
        if (ctx->outfile == NULL)
            error("initoutput: cannot open file %s\n", ctx->outfilename);
    }
    else
        ctx->outfile = NULL;              /* turn off data output */

}

/*  * Counters and accumulators for output routines.
 */

static real mtot;                /* total mass of N-body system */
static real etot[3];             /* binding, kinetic, potential energy */
static matrix keten;     /* kinetic energy tensor */
static matrix peten;     /* potential energy tensor */
static vector cmphase[2];    /* center of mass coordinates */
static vector amvec;     /* angular momentum vector */

/*  * OUTPUT: compute diagnostics and output data.
 */

void output(void)
{
    bodyptr p;
    vector lbR;
    diagnostics();              /* compute std diagnostics */
    //printf("st.tnow = %f\n", st.tnow);
    if (ctx.tstop - st.tnow < 0.01 / ctx.freq)
    {
        printf("st.tnow = %f\n", st.tnow);
        for (p = st.bodytab; p < st.bodytab + ctx.nbody; p++)
        {
            (lbR)[2] = sqrt(Pos(p)[0] * Pos(p)[0] + Pos(p)[1] * Pos(p)[1] + Pos(p)[2] * Pos(p)[2]);
            (lbR)[1] = r2d(atan2(Pos(p)[2], sqrt((Pos(p)[0]) * (Pos(p)[0]) + Pos(p)[1] * Pos(p)[1])));
            (lbR)[0] = r2d(atan2(Pos(p)[1], Pos(p)[0]));

            if ((lbR)[0] < 0)
            {
                (lbR)[0] += 360.0;
            }

            out_2vectors(ctx.outfile, lbR, Vel(p));
        }
        printf("\tParticle data written to file %s\n\n", ctx.outfilename);
        fflush(ctx.outfile);             /* drain output buffer */
        st.tout += 1 / ctx.freqout;     /* schedule next data out */
    }
}

/*  * DIAGNOSTICS: compute various dynamical diagnostics.
 */

static void diagnostics(void)
{
    register bodyptr p;
    real velsq;
    vector tmpv;
    matrix tmpt;

    mtot = 0.0;                 /* zero total mass */
    etot[1] = etot[2] = 0.0;            /* zero total KE and PE */
    CLRM(keten);                /* zero ke tensor */
    CLRM(peten);                /* zero pe tensor */
    CLRV(cmphase[0]);               /* zero c. of m. position */
    CLRV(cmphase[1]);               /* zero c. of m. velocity */
    CLRV(amvec);                /* zero am vector */
    for (p = st.bodytab; p < st.bodytab + ctx.nbody; p++) /* loop over all particles */
    {
        mtot += Mass(p);                        /* sum particle masses */
        DOTVP(velsq, Vel(p), Vel(p));       /* square vel vector */
        etot[1] += 0.5 * Mass(p) * velsq;   /* sum current KE */
        etot[2] += 0.5 * Mass(p) * Phi(p);  /* and current PE */
        MULVS(tmpv, Vel(p), 0.5 * Mass(p)); /* sum 0.5 m v_i v_j */
        OUTVP(tmpt, tmpv, Vel(p));
        ADDM(keten, keten, tmpt);
        MULVS(tmpv, Pos(p), Mass(p));       /* sum m r_i a_j */
        OUTVP(tmpt, tmpv, Acc(p));
        ADDM(peten, peten, tmpt);
        MULVS(tmpv, Pos(p), Mass(p));       /* sum cm position */
        ADDV(cmphase[0], cmphase[0], tmpv);
        MULVS(tmpv, Vel(p), Mass(p));       /* sum cm momentum */
        ADDV(cmphase[1], cmphase[1], tmpv);
        CROSSVP(tmpv, Pos(p), Vel(p));      /* sum angular momentum */
        MULVS(tmpv, tmpv, Mass(p));
        ADDV(amvec, amvec, tmpv);
    }
    etot[0] = etot[1] + etot[2];                /* sum KE and PE */
    DIVVS(cmphase[0], cmphase[0], mtot);        /* normalize cm coords */
    DIVVS(cmphase[1], cmphase[1], mtot);
}

/*  * Low-level input and output operations.
 */

static void in_int(FILE* str, int* iptr)
{
    if (fscanf(str, "%d", iptr) != 1)
        error("in_int: input conversion error\n");
}

static void in_real(FILE* str, real* rptr)
{
    double tmp;

    if (fscanf(str, "%lf", &tmp) != 1)
        error("in_real: input conversion error\n");
    *rptr = tmp;
}

static void in_vector(FILE* str, vector vec)
{
    double tmpx, tmpy, tmpz;

    if (fscanf(str, "%lf%lf%lf", &tmpx, &tmpy, &tmpz) != 3)
        error("in_vector: input conversion error\n");
    vec[0] = tmpx;
    vec[1] = tmpy;
    vec[2] = tmpz;
}

static void out_int(FILE* str, int ival)
{
    fprintf(str, "  %d\n", ival);
}

static void out_real(FILE* str, real rval)
{
    fprintf(str, " %21.14E\n", rval);
}

static void out_vector(FILE* str, vector vec)
{
    fprintf(str, " %21.14E %21.14E %21.14E\n", vec[0], vec[1], vec[2]);
}

static void out_2vectors(FILE* str, vector vec1, vector vec2)
{
    fprintf(str, " %21.14E %21.14E %21.14E %21.14E %21.14E %21.14E\n", vec1[0], vec1[1], vec1[2], vec2[0], vec2[1], vec2[2]);
}

static void printvec(char* name, vector vec)
{
    printf("          %10s%10.4f%10.4f%10.4f\n",
           name, vec[0], vec[1], vec[2]);
}

void printContext(NBodyCtx* ctx)
{

    printf("ctx = { \n"
           "  nbody              = %d\n"
           "  outfilename        = %s\n"
           "  headline           = %s\n"
           "  criterion          = %d\n"
           "  usequad            = %d\n"
           "  allowIncest        = %d\n"
           "  seed               = %d\n"
           "  accuracy parameter = %g\n"
           "};\n",
           ctx->nbody,
           ctx->outfilename,
           ctx->headline,
           ctx->criterion,
           ctx->usequad,
           ctx->allowIncest,
           ctx->seed,
           ctx->theta
        );
}

