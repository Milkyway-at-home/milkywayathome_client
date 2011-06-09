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
#include "mw_asprintf.h"

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
        case NBODY_TREE_INCEST_NONFATAL:
            return "NBODY_TREE_INCEST_NONFATAL";
        case NBODY_SUCCESS:
            return "NBODY_SUCCESS";
        case NBODY_ERROR:
            return "NBODY_ERROR";
        case NBODY_TREE_STRUCTURE_ERROR:
            return "NBODY_TREE_STRUCTURE_ERROR";
        case NBODY_TREE_INCEST_FATAL:
            return "NBODY_TREE_INCEST_FATAL";
        case NBODY_IO_ERROR:
            return "NBODY_IO_ERROR";
        case NBODY_CHECKPOINT_ERROR:
            return "NBODY_CHECKPOINT_ERROR";
        default:
            return "Invalid NBodyStatus";
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
        fail("asprintf() failed\n");
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
        fail("asprintf() failed\n");
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
        fail("asprintf() failed\n");
    }

    return buf;
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
                     "      mass   = %g\n"
                     "      pos    = %s\n"
                     "      vel    = %s\n"
                     "      ignore = %s\n"
                     "    };\n",
                     Mass(p),
                     pos,
                     vel,
                     showBool(ignoreBody(p))))

    {
        fail("asprintf() failed\n");
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
        fail("asprintf() failed\n");

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

    if (0 > asprintf(&buf,
                     "ctx = { \n"
                     "  pot = %s\n"
                     "  timeEvolve      = %g\n"
                     "  timestep        = %g\n"
                     "  sunGCDist       = %g\n"
                     "  criterion       = %s\n"
                     "  useQuad         = %s\n"
                     "  allowIncest     = %s\n"
                     "  treeRSize       = %g\n"
                     "  theta           = %g\n"
                     "  eps2            = %g\n"
                     "  checkpointT     = %u\n"
                     "  freqOut         = %u\n"
                     "};\n",
                     potBuf,
                     ctx->timeEvolve,
                     ctx->timestep,
                     ctx->sunGCDist,
                     showCriterionT(ctx->criterion),
                     showBool(ctx->useQuad),
                     showBool(ctx->allowIncest),
                     ctx->treeRSize,
                     ctx->theta,
                     ctx->eps2,
                     (int) ctx->checkpointT,
                     ctx->freqOut))
    {
        fail("asprintf() failed\n");
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
        fail("asprintf() failed\n");
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

void printBodies(const Body* bs, unsigned int n)
{
    unsigned int i;

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
                     "    cellused = %u\n"
                     "    maxlevel = %u\n"
                     "  };\n",
                     t,
                     t->root,
                     t->rsize,
                     t->cellused,
                     t->maxlevel))
    {
        fail("asprintf() failed\n");
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
                     "  freecell       = %p\n"
                     "  lastCheckpoint = %d\n"
                     "  tnow           = %.15g\n"
                     "  nbody          = %u\n"
                     "  bodytab        = %p\n"
                     "  acctab         = %p\n"
                     "  treeIncest     = %s\n"
                     "  outFile        = %p\n"
                     "};\n",
                     st,
                     treeBuf,
                     st->freecell,
                     (int) st->lastCheckpoint,
                     st->tnow,
                     st->nbody,
                     st->bodytab,
                     st->acctab,
                     showBool(st->treeIncest),
                     st->outFile))
    {
        fail("asprintf() failed\n");
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

