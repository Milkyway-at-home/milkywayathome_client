/* Copyright 2010 Matthew Arsenault, Travis Desell, Dave Przybylo,
Nathan Cole, Boleslaw Szymanski, Heidi Newberg, Carlos Varela, Malik
Magdon-Ismail and Rensselaer Polytechnic Institute.

This file is part of Milkway@Home.

Milkyway@Home is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Milkyway@Home is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Milkyway@Home.  If not, see <http://www.gnu.org/licenses/>.
*/


#include "nbody_types.h"
#include "milkyway_util.h"
#include "nbody_show.h"

/* A bunch of boilerplate for debug printing */

const char* showBool(mwbool x)
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

const char* showCriterionT(criterion_t x)
{
    switch (x)
    {
        case Exact:
            return "Exact";
        case NewCriterion:
            return "NewCriterion";
        case BH86:
            return "BH86";
        case SW93:
            return "SW93";
        case InvalidCriterion:
            return "InvalidCriterion";
        default:
            return "Bad criterion_t";
    }
}

const char* showSphericalT(spherical_t x)
{
    switch (x)
    {
        case SphericalPotential:
            return "SphericalPotential";
        case InvalidSpherical:
            return "InvalidSpherical";
        default:
            return "Bad spherical_t";
    }
}

const char* showDiskT(disk_t x)
{
    switch (x)
    {
        case MiyamotoNagaiDisk:
            return "MiyamotoNagaiDisk";
        case ExponentialDisk:
            return "ExponentialDisk";
        case InvalidDisk:
            return "InvalidDisk";
        default:
            return "Bad disk_t";
    }
}

const char* showHaloT(halo_t x)
{
    switch (x)
    {
        case LogarithmicHalo:
            return "LogarithmicHalo";
        case NFWHalo:
            return "NFWHalo";
        case TriaxialHalo:
            return "TriaxialHalo";
        case InvalidHalo:
            return "InvalidHalo";
        default:
            return "Bad halo_t";
    }
}

const char* showNBodyStatus(NBodyStatus x)
{
    switch (x)
    {
        case NBODY_SUCCESS:
            return "NBODY_SUCCESS";
        case NBODY_ERROR:
            return "NBODY_ERROR";
        case NBODY_ASSERTION_FAILURE:
            return "NBODY_ASSERTION_FAILURE";
        case NBODY_TREE_STRUCTURE_ERROR:
            return "NBODY_TREE_STRUCTURE_ERROR";
        case NBODY_IO_ERROR:
            return "NBODY_IO_ERROR";
        case NBODY_CHECKPOINT_ERROR:
            return "NBODY_CHECKPOINT_ERROR";
        case NBODY_CL_ERROR:
            return "NBODY_CL_ERROR";
        case NBODY_CAPABILITY_ERROR:
            return "NBODY_CAPABILITY_ERROR";
        case NBODY_CONSISTENCY_ERROR:
            return "NBODY_CONSISTENCY_ERROR";
        case NBODY_UNIMPLEMENTED:
            return "NBODY_UNIMPLEMENTED";
        case NBODY_UNSUPPORTED:
            return "NBODY_UNSUPPORTED";
        case NBODY_USER_ERROR:
            return "NBODY_USER_ERROR";
        case NBODY_PARAM_FILE_ERROR:
            return "NBODY_PARAM_FILE_ERROR";
        case NBODY_LUA_POTENTIAL_ERROR:
            return "NBODY_LUA_POTENTIAL_ERROR";
        case NBODY_LIKELIHOOD_ERROR:
            return "NBODY_LIKELIHOOD_ERROR";
        case NBODY_MAX_DEPTH_ERROR:
            return "NBODY_MAX_DEPTH_ERROR";
        case NBODY_CELL_OVERFLOW_ERROR:
            return "NBODY_CELL_OVERFLOW_ERROR";
        case NBODY_RESERVED_ERROR_1:
            return "NBODY_RESERVED_ERROR_1";
        case NBODY_RESERVED_ERROR_2:
            return "NBODY_RESERVED_ERROR_2";
        case NBODY_RESERVED_ERROR_3:
            return "NBODY_RESERVED_ERROR_3";
        case NBODY_RESERVED_ERROR_4:
            return "NBODY_RESERVED_ERROR_4";
        case NBODY_RESERVED_ERROR_5:
            return "NBODY_RESERVED_ERROR_5";
        case NBODY_RESERVED_ERROR_6:
            return "NBODY_RESERVED_ERROR_6";

        case NBODY_TREE_INCEST_NONFATAL:
            return "NBODY_TREE_INCEST_NONFATAL";
        case NBODY_TREE_INCEST_FATAL:
            return "NBODY_TREE_INCEST_FATAL";
        case NBODY_RESERVED_WARNING_1:
            return "NBODY_RESERVED_WARNING_1";
        case NBODY_RESERVED_WARNING_2:
            return "NBODY_RESERVED_WARNING_2";
        case NBODY_RESERVED_WARNING_3:
            return "NBODY_RESERVED_WARNING_3";
        case NBODY_RESERVED_WARNING_4:
            return "NBODY_RESERVED_WARNING_4";
        case NBODY_RESERVED_WARNING_5:
            return "NBODY_RESERVED_WARNING_5";
        case NBODY_RESERVED_WARNING_6:
            return "NBODY_RESERVED_WARNING_6";
        case NBODY_RESERVED_WARNING_7:
            return "NBODY_RESERVED_WARNING_7";

        default:
            return "Unknown NBodyStatus";
    }
}

const char* showNBodyKernelError(NBodyKernelError x)
{
    if ((int) x > 0)
    {
        return "Exceeded maximum depth";
    }

    switch (x)
    {
        case NBODY_KERNEL_OK:
            return "NBODY_KERNEL_OK";
        case NBODY_KERNEL_CELL_OVERFLOW:
            return "NBODY_KERNEL_CELL_OVERFLOW";
        case NBODY_KERNEL_TREE_INCEST:
            return "NBODY_KERNEL_CELL_TREE_INCEST";
        case NBODY_KERNEL_TREE_STRUCTURE_ERROR:
            return "NBODY_KERNEL_TREE_STRUCTURE_ERROR";
        case NBODY_KERNEL_ERROR_OTHER:
            return "NBODY_KERNEL_ERROR_OTHER";
        default:
            return "Unknown NBodyKernelError";
    }
}

const char* showExternalPotentialType(ExternalPotentialType x)
{
    switch (x)
    {
        case EXTERNAL_POTENTIAL_DEFAULT:
            return "Milkyway@Home N-body potential";
        case EXTERNAL_POTENTIAL_NONE:
            return "None";
        case EXTERNAL_POTENTIAL_CUSTOM_LUA:
            return "Lua";
        default:
            return "Unknown ExternalPotentialType";
    }
}

char* showSpherical(const Spherical* s)
{
    char* buf;

    if (!s)
        return NULL;

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
        mw_fail("asprintf() failed\n");
    }

    return buf;
}

char* showHalo(const Halo* h)
{
    char* buf;

    if (!h)
        return NULL;

    if (0 > asprintf(&buf,
                     "{ \n"
                     "      type         = %s\n"
                     "      vhalo        = %g\n"
                     "      scaleLength  = %g\n"
                     "      flattenX     = %g\n"
                     "      flattenY     = %g\n"
                     "      flattenZ     = %g\n"
                     "      c1           = %g\n"
                     "      c2           = %g\n"
                     "      c3           = %g\n"
                     "      triaxAngle   = %g\n"
                     "    };\n",
                     showHaloT(h->type),
                     h->vhalo,
                     h->scaleLength,
                     h->flattenX,
                     h->flattenY,
                     h->flattenZ,
                     h->c1,
                     h->c2,
                     h->c3,
                     h->triaxAngle))
    {
        mw_fail("asprintf() failed\n");
    }

    return buf;
}

char* showDisk(const Disk* d)
{
    char* buf;

    if (!d)
        return NULL;

    if (0 > asprintf(&buf,
                     "{ \n"
                     "      type        = %s\n"
                     "      mass        = %g\n"
                     "      scaleLength = %g\n"
                     "      scaleHeight = %g\n"
                     "    };\n",
                     showDiskT(d->type),
                     d->mass,
                     d->scaleLength,
                     d->scaleHeight))
    {
        mw_fail("asprintf() failed\n");
    }

    return buf;
}

char* showCell(const NBodyCell* c)
{
    char* buf;
    char* posBuf;

    if (!c)
        return NULL;


    posBuf = showVector(Pos(c));

    if (0 > asprintf(&buf,
                     "NBodyCell = {\n"
                     "  cellnode = {\n"
                     "    pos  = %s\n"
                     "    next = %p\n"
                     "    mass = %f\n"
                     "    type = %d\n"
                     "  }\n"
                     "  rcrit2   = %f\n"
                     "  more     = %p\n"
                     "  stuff    = {\n"
                     "    .quad = {\n"
                     "      .xx = %f, .xy = %f, .xz = %f,\n"
                     "      .yy = %f, .yz = %f,\n"
                     "      .zz = %f\n"
                     "    },\n"
                     "\n"
                     "    .subp = {\n"
                     "      %p, %p, %p, %p,\n"
                     "      %p, %p, %p, %p\n"
                     "    }\n"
                     "  }\n"
                     "}\n",
                     posBuf,
                     (void*) Next(c),
                     Mass(c),
                     Type(c),

                     Rcrit2(c),
                     (void*) More(c),

                     Quad(c).xx, Quad(c).xy, Quad(c).xz,
                     Quad(c).yy, Quad(c).yz,
                     Quad(c).zz,

                     (void*) Subp(c)[0], (void*) Subp(c)[1], (void*) Subp(c)[2], (void*) Subp(c)[3],
                     (void*) Subp(c)[5], (void*) Subp(c)[5], (void*) Subp(c)[6], (void*) Subp(c)[7]
            ))
    {
        mw_fail("asprintf() failed\n");
    }

    free(posBuf);

    return buf;
}

void printCell(const NBodyCell* c)
{
    char* buf = showCell(c);
    puts(buf);
    free(buf);
}

char* showBody(const Body* p)
{
    char* buf;
    char* vel;
    char* pos;

    if (!p)
        return NULL;

    vel = showVector(Vel(p));
    pos = showVector(Pos(p));

    if (0 > asprintf(&buf,
                     "body { \n"
                     "      mass     = %g\n"
                     "      position = %s\n"
                     "      velocity = %s\n"
                     "      ignore   = %s\n"
                     "    };\n",
                     Mass(p),
                     pos,
                     vel,
                     showBool(ignoreBody(p))))

    {
        mw_fail("asprintf() failed\n");
    }

    free(vel);
    free(pos);

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

    if (!p)
        return NULL;

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
        mw_fail("asprintf() failed\n");

    free(sphBuf);
    free(diskBuf);
    free(haloBuf);

    return buf;
}

/* Most efficient function ever */
char* showNBodyCtx(const NBodyCtx* ctx)
{
    char* buf;
    char* potBuf;

    if (!ctx)
        return NULL;

    potBuf = showPotential(&ctx->pot);
    if (!potBuf)
        return NULL;

    if (0 > asprintf(&buf,
                     "ctx = { \n"
                     "  eps2            = %f\n"
                     "  theta           = %f\n"
                     "  timestep        = %f\n"
                     "  timeEvolve      = %f\n"
                     "  treeRSize       = %f\n"
                     "  sunGCDist       = %f\n"
                     "  criterion       = %s\n"
                     "  useQuad         = %s\n"
                     "  allowIncest     = %s\n"
                     "  checkpointT     = %d\n"
                     "  nStep           = %u\n"
                     "  potentialType   = %s\n"
                     "  pot = %s\n"
                     "};\n",
                     ctx->eps2,
                     ctx->theta,
                     ctx->timestep,
                     ctx->timeEvolve,
                     ctx->treeRSize,
                     ctx->sunGCDist,
                     showCriterionT(ctx->criterion),
                     showBool(ctx->useQuad),
                     showBool(ctx->allowIncest),
                     (int) ctx->checkpointT,
                     ctx->nStep,
                     showExternalPotentialType(ctx->potentialType),
                     potBuf
            ))
    {
        mw_fail("asprintf() failed\n");
    }

    free(potBuf);

    return buf;
}

void printNBodyCtx(const NBodyCtx* ctx)
{
    char* buf = showNBodyCtx(ctx);
    puts(buf);
    free(buf);
}

char* showHistogramParams(const HistogramParams* hp)
{
    char* buf;
    if (0 > asprintf(&buf,
                     "histogram-params = { \n"
                     "  phi      = %g\n"
                     "  theta    = %g\n"
                     "  psi      = %g\n"
                     "  startRaw = %g\n"
                     "  endRaw   = %g\n"
                     "  binSize  = %g\n"
                     "  center   = %g\n"
                     "};\n",
                     hp->phi,
                     hp->theta,
                     hp->psi,
                     hp->startRaw,
                     hp->endRaw,
                     hp->binSize,
                     hp->center))

    {
        mw_fail("asprintf() failed\n");
    }

    return buf;
}

void printHistogramParams(const HistogramParams* hp)
{
    char* buf = showHistogramParams(hp);
    puts(buf);
    free(buf);
}

void printBody(const Body* p)
{
    char* buf = showBody(p);
    puts(buf);
    free(buf);
}

void printHalo(const Halo* h)
{
    char* buf = showHalo(h);
    puts(buf);
    free(buf);
}

void printDisk(const Disk* d)
{
    char* buf = showDisk(d);
    puts(buf);
    free(buf);
}

void printPotential(const Potential* p)
{
    char* buf = showPotential(p);
    puts(buf);
    free(buf);
}

void printBodies(const Body* bs, int n)
{
    int i;

    for (i = 0; i < n; ++i)
        printBody(&bs[i]);
}

char* showNBodyTree(const NBodyTree* t)
{
    char* buf;

    if (0 > asprintf(&buf,
                     "  Tree %p = {\n"
                     "    root     = %p\n"
                     "    rsize    = %g\n"
                     "    cellUsed = %u\n"
                     "    maxDepth = %u\n"
                     "  };\n",
                     t,
                     t->root,
                     t->rsize,
                     t->cellUsed,
                     t->maxDepth))
    {
        mw_fail("asprintf() failed\n");
    }

    return buf;
}

void printNBodyTree(const NBodyTree* t)
{
    char* buf = showNBodyTree(t);
    puts(buf);
    free(buf);
}

char* showNBodyState(const NBodyState* st)
{
    char* buf;
    char* treeBuf;

    treeBuf = showNBodyTree(&st->tree);

    if (0 > asprintf(&buf,
                     "NBodyState %p = {\n"
                     "  tree           = %s\n"
                     "  freeCell       = %p\n"
                     "  lastCheckpoint = %d\n"
                     "  step           = %u\n"
                     "  nbody          = %u\n"
                     "  bodytab        = %p\n"
                     "  acctab         = %p\n"
                     "  treeIncest     = %s\n"
                     "};\n",
                     st,
                     treeBuf,
                     st->freeCell,
                     (int) st->lastCheckpoint,
                     st->step,
                     st->nbody,
                     st->bodytab,
                     st->acctab,
                     showBool(st->treeIncest)
            ))
    {
        mw_fail("asprintf() failed\n");
    }

    free(treeBuf);

    return buf;

}

void printNBodyState(const NBodyState* st)
{
    char* buf = showNBodyState(st);
    puts(buf);
    free(buf);
}

const char* showNBodyLikelihoodMethod(NBodyLikelihoodMethod x)
{
    switch (x)
    {
        case NBODY_INVALID_METHOD:
            return "InvalidMethod";
        case NBODY_EMD:
            return "EMD";
        case NBODY_ORIG_CHISQ:
            return "Original";
        case NBODY_ORIG_ALT:
            return "AltOriginal";
        case NBODY_CHISQ_ALT:
            return "ChisqAlt";
        case NBODY_POISSON:
            return "Poisson";
        case NBODY_KOLMOGOROV:
            return "Kolmogorov";
        case NBODY_KULLBACK_LEIBLER:
            return "KullbackLeibler";
        case NBODY_SAHA:
            return "Saha";
        default:
            return "Invalid NBodyLikelihoodMethod";
    }
}

