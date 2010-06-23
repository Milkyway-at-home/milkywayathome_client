/* ************************************************************************** */
/* IO.C: I/O routines for export version of hierarchical N-body code. */
/* Public routines: inputdata(), inictx.toutput(), stopoutput(), output(). */
/* */
/* Copyright (c) 1993 by Joshua E. Barnes, Honolulu, HI. */
/* It's free because it's yours. */
/* ************************************************************************** */

#include "nbody_priv.h"

void initoutput(NBodyCtx* ctx)
{
    if (ctx->outfilename)                       /* output file specified? */
    {
        ctx->outfile = fopen(ctx->outfilename, "w");           /* setup output FILE* */
        if (ctx->outfile == NULL)
            error("initoutput: cannot open file %s\n", ctx->outfilename);
    }
    else
        ctx->outfile = stdout;
}

/* Low-level input and output operations. */

static void out_2vectors(FILE* str, vector vec1, vector vec2)
{
    fprintf(str, " %21.14E %21.14E %21.14E %21.14E %21.14E %21.14E\n", vec1[0], vec1[1], vec1[2], vec2[0], vec2[1], vec2[2]);
}

/* DIAGNOSTICS: compute various dynamical diagnostics.  */

static void diagnostics(const NBodyCtx* ctx, const NBodyState* st)
{
    register bodyptr p;
    real velsq;
    vector tmpv;
    matrix tmpt;
    real mtot = 0.0;                                     /* total mass of N-body system */
    real etot[3] = { 0.0, 0.0, 0.0 };                    /* binding, kinetic, potential energy */
    matrix keten = ZERO_MATRIX;                          /* kinetic energy tensor */
    matrix peten = ZERO_MATRIX;                          /* potential energy tensor */
    vector cmphase[2] = { ZERO_VECTOR, ZERO_VECTOR };    /* center of mass coordinates */
    vector amvec = ZERO_VECTOR;                         /* angular momentum vector */

    const bodyptr endp = st->bodytab + ctx->model.nbody;
    for (p = st->bodytab; p < endp; p++) /* loop over all particles */
    {
        mtot += Mass(p);                    /* sum particle masses */
        SQRV(velsq, Vel(p));                /* square vel vector */
        etot[1] += 0.5 * Mass(p) * velsq;   /* sum current KE */
        etot[2] += 0.5 * Mass(p) * Phi(p);  /* and current PE */
        MULVS(tmpv, Vel(p), 0.5 * Mass(p)); /* sum 0.5 m v_i v_j */
        OUTVP(tmpt, tmpv, Vel(p));
        ADDM(keten, keten, tmpt);
        MULVS(tmpv, Pos(p), Mass(p));       /* sum m r_i a_j */
        OUTVP(tmpt, tmpv, Acc(p));
        ADDM(peten, peten, tmpt);
        MULVS(tmpv, Pos(p), Mass(p));       /* sum cm position */
        INCADDV(cmphase[0], tmpv);
        MULVS(tmpv, Vel(p), Mass(p));       /* sum cm momentum */
        INCADDV(cmphase[1], tmpv);
        CROSSVP(tmpv, Pos(p), Vel(p));      /* sum angular momentum */
        INCMULVS(tmpv, Mass(p));
        INCADDV(amvec, tmpv);
    }
    etot[0] = etot[1] + etot[2];                /* sum KE and PE */
    INCDIVVS(cmphase[0], mtot);                 /* normalize cm coords */
    INCDIVVS(cmphase[1], mtot);
}

/* OUTPUT: compute diagnostics and output data. */
void output(const NBodyCtx* ctx, NBodyState* st)
{
    bodyptr p;
    vector lbR;
    const bodyptr endp = st->bodytab + ctx->model.nbody;

    diagnostics(ctx, st);              /* compute std diagnostics */
    if (ctx->model.time_dwarf - st->tnow < 0.01 / ctx->freq)
    {
        printf("st.tnow = %f\n", st->tnow);
        for (p = st->bodytab; p < endp; p++)
        {
            lbR[2] = rsqrt(Pos(p)[0] * Pos(p)[0] + Pos(p)[1] * Pos(p)[1] + Pos(p)[2] * Pos(p)[2]);
            lbR[1] = r2d(ratan2(Pos(p)[2], rsqrt((Pos(p)[0]) * (Pos(p)[0]) + Pos(p)[1] * Pos(p)[1])));
            lbR[0] = r2d(ratan2(Pos(p)[1], Pos(p)[0]));

            if (lbR[0] < 0)
                lbR[0] += 360.0;

            out_2vectors(ctx->outfile, lbR, Vel(p));
        }

        if (ctx->outfile != stdout)
            printf("\tParticle data written to file %s\n\n", ctx->outfilename);
        fflush(ctx->outfile);             /* drain output buffer */
        st->tout += 1.0 / ctx->freqout;     /* schedule next data out */
    }
}


/* A bunch of boilerplate for debug printing */

const char* showBool(const bool x)
{
    switch (x)
    {
        case FALSE:
            return "false";
        case TRUE:
            return "true";
        default:
            return "invalid boolean (but true)";
    }
}

const char* showCriterionT(const criterion_t x)
{
    static const char* table[] =
        {
            [EXACT]        = "exact",
            [NEWCRITERION] = "new criterion",
            [BH86]         = "bh86",
            [SW93]         = "sw93"
        };

    if (x > SW93)
        return "invalid criterion_t";
    else
        return table[x];

}

const char* showSphericalT(const spherical_t x)
{
    static const char* table[] =
        {
            [SphericalPotential]        = "SphericalPotential",
        };

    if (x > SphericalPotential)
        return "invalid spherical_t";
    else
        return table[x];
}

const char* showDiskT(const disk_t x)
{
    static const char* table[] =
        {
            [MiyamotoNagaiDisk] = "MiyamotoNagaiDisk",
            [ExponentialDisk]   = "ExponentialDisk"
        };

    if (x > ExponentialDisk)
        return "invalid disk_t";
    else
        return table[x];
}

const char* showHaloT(const halo_t x)
{
    static const char* table[] =
        {
            [LogarithmicHalo] = "LogarithmicHalo",
            [NFWHalo]         = "NFWHalo",
            [TriaxialHalo]    = "TriaxialHalo"
        };

    if (x > TriaxialHalo)
        return "invalid halo_t";
    else
        return table[x];
}

const char* showDwarfModelT(const dwarf_model_t x)
{
    static const char* table[] =
        {
            [DwarfModelPlummer] = "DwarfModelPlummer",
            [DwarfModelKing]    = "DwarfModelKing",
            [DwarfModelDehnen]  = "DwarfModelDehnen"
        };

    if (x > DwarfModelDehnen)
        return "invalid dwarf_model_t";
    else
        return table[x];
}

char* showSpherical(const Spherical* s)
{
    char* buf;

    if (0 > asprintf(&buf,
                     "{\n"
                     "      type  = %s\n"
                     "      mass  = %g\n"
                     "      scale = %g\n"
                     "    };\n",
                     showSphericalT(s->type),
                     s->mass,
                     s->scale))
    {
        fail("asprintf() failed\n");
    }

    return buf;
}

char* showHalo(const Halo* h)
{
    char* buf;

    if (0 > asprintf(&buf,
                     "{ \n"
                     "      type         = %s\n"
                     "      vhalo        = %g\n"
                     "      scale_length = %g\n"
                     "      flattenX     = %g\n"
                     "      flattenY     = %g\n"
                     "      flattenZ     = %g\n"
                     "      triaxAngle   = %g\n"
                     "    };\n",
                     showHaloT(h->type),
                     h->vhalo,
                     h->scale_length,
                     h->flattenX,
                     h->flattenY,
                     h->flattenZ,
                     h->triaxAngle))
    {
        fail("asprintf() failed\n");
    }

    return buf;
}

char* showDisk(const Disk* d)
{
    char* buf;

    if (0 > asprintf(&buf,
                     "{ \n"
                     "      type         = %s\n"
                     "      mass         = %g\n"
                     "      scale_length = %g\n"
                     "      scale_height = %g\n"
                     "    };\n",
                     showDiskT(d->type),
                     d->mass,
                     d->scale_length,
                     d->scale_height))
    {
        fail("asprintf() failed\n");
    }

    return buf;
}

/* For debugging. Need to make this go away for release since it uses
 * GNU extensions */
char* showPotential(const Potential* p)
{
    int rc;
    char* buf;
    char* sphBuf;
    char* diskBuf;
    char* haloBuf;

    sphBuf  = showSpherical(&p->sphere[0]);
    diskBuf = showDisk(&p->disk);
    haloBuf = showHalo(&p->halo);

    rc = asprintf(&buf,
                  "{\n"
                  "    sphere = %s\n"
                  "    disk = %s\n"
                  "    halo = %s\n"
                  "    rings  = { unused pointer %p }\n"
                  "  };\n",
                  sphBuf,
                  diskBuf,
                  haloBuf,
                  p->rings);

    if (rc < 0)
        fail("asprintf() failed\n");

    free(sphBuf);
    free(diskBuf);
    free(haloBuf);

    return buf;
}

char* showDwarfModel(const DwarfModel* d)
{
    char* buf;

    if (0 > asprintf(&buf,
                     "{ \n"
                     "      type           = %s\n"
                     "      nbody          = %d\n"
                     "      mass           = %g\n"
                     "      scale_radius   = %g\n"
                     "      timestep       = %g\n"
                     "      orbit_timestep = %g\n"
                     "      time_dwarf     = %g\n"
                     "      time_orbit     = %g\n"
                     "      eps            = %g\n"
                     "    };\n",
                     showDwarfModelT(d->type),
                     d->nbody,
                     d->mass,
                     d->scale_radius,
                     d->timestep,
                     d->orbit_timestep,
                     d->time_orbit,
                     d->time_dwarf,
                     d->eps))
    {
        fail("asprintf() failed\n");
    }

    return buf;
}

char* showInitialConditions(const InitialConditions* ic)
{
    char* buf;
    if (0 > asprintf(&buf,
                     "initial-conditions = { \n"
                     "  useGalC    = %s\n"
                     "  useRadians = %s\n"
                     "  sunGCDist  = %g\n"
                     "  position   = { %g, %g, %g }\n"
                     "  velocity   = { %g, %g, %g }\n"
                     "};\n",
                     showBool(ic->useGalC),
                     showBool(ic->useRadians),
                     ic->sunGCDist,
                     ic->position[0],
                     ic->position[1],
                     ic->position[2],
                     ic->velocity[0],
                     ic->velocity[1],
                     ic->velocity[2]))
    {
        fail("asprintf() failed\n");
    }

    return buf;
}

char* showContext(const NBodyCtx* ctx)
{
    char* buf;
    char* potBuf;
    char* modelBuf;

    potBuf   = showPotential(&ctx->pot);
    modelBuf = showDwarfModel(&ctx->model);

    if (0 > asprintf(&buf,
                     "ctx = { \n"
                     "  pot = %s\n"
                     "  model = %s\n"
                     "  headline    = %s\n"
                     "  outfilename = %s\n"
                     "  outfile     = %p\n"
                     "  criterion   = %s\n"
                     "  usequad     = %s\n"
                     "  allowIncest = %s\n"
                     "  seed        = %d\n"
                     "  tree_rsize  = %g\n"
                     "  theta       = %g\n"
                     "  freq        = %g\n"
                     "  freqout     = %g\n"
                     "};\n",
                     potBuf,
                     modelBuf,
                     ctx->headline,
                     ctx->outfilename,
                     ctx->outfile,
                     showCriterionT(ctx->criterion),
                     showBool(ctx->usequad),
                     showBool(ctx->allowIncest),
                     ctx->seed,
                     ctx->tree_rsize,
                     ctx->theta,
                     ctx->freq,
                     ctx->freqout))
    {
        fail("asprintf() failed\n");
    }

    free(potBuf);
    free(modelBuf);

    return buf;
}

void printContext(const NBodyCtx* ctx)
{
    char* buf = showContext(ctx);
    puts(buf);
    free(buf);
}

void printInitialConditions(const InitialConditions* ic)
{
    char* buf = showInitialConditions(ic);
    puts(buf);
    free(buf);
}

char* showVector(const vector v)
{
    char* buf;

    if (asprintf(&buf, "{ %g, %g, %g }", v[0], v[1], v[2]) < 0)
        fail ("asprintf() failed\n");

    return buf;

}

