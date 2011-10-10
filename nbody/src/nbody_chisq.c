/* Copyright 2010 Ben Willett, Matthew Arsenault, Boleslaw Szymanski,
Heidi Newberg, Carlos Varela, Malik Magdon-Ismail and Rensselaer
Polytechnic Institute.

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "nbody_priv.h"
#include "nbody_chisq.h"
#include "milkyway_util.h"

/* Calculate chisq from data read from histData and the histogram
 * generated from the simulation, histogram, with maxIdx bins. */
static real nbCalcChisq(const HistData* histData,
                        const unsigned int* histogram,
                        const unsigned int maxIdx,
                        const real totalNum)
{
    unsigned int i;
    real tmp, chisqval = 0.0;

    for (i = 0; i < maxIdx; ++i)
    {
        if (!histData[i].useBin)  /* Skip bins with missing data */
            continue;

        tmp = (histData[i].count - ((real) histogram[i] / totalNum)) / histData[i].err;
        chisqval += sqr(tmp);
    }

    // MAXIMUM likelihood, multiply by -1
    return -chisqval;
}

static void nbPrintHistogramHeader(FILE* f,
                                   const NBodyCtx* ctx,
                                   const HistogramParams* hp,
                                   int nbody,
                                   uint32_t seed,
                                   real chisq)
{
    char tBuf[256];
    const Potential* p = &ctx->pot;

    mwLocalTimeFull(tBuf, sizeof(tBuf));

    fprintf(f,
            "#\n"
            "# Generated %s\n"
            "#\n"
            "# Likelihood = %.15f\n"
            "#\n"
            "# (phi,   theta,  psi) = (%f, %f, %f)\n"
            "# (start, center, end) = (%f, %f, %f)\n"
            "# Bin size = %f\n"
            "#\n"
            "#\n",
            tBuf,
            chisq,
            hp->phi, hp->theta, hp->psi,
            hp->startRaw, hp->center, hp->endRaw,
            hp->binSize);


    fprintf(f,
            "# Nbody = %d\n"
            "# Evolve time = %f\n"
            "# Timestep = %f\n"
            "# Sun GC Dist = %f\n"
            "# Criterion = %s\n"
            "# Theta = %f\n"
            "# Quadrupole Moments = %s\n"
            "# Eps = %f\n"
            "# Seed = %u\n"
            "#\n",
            nbody,
            ctx->timeEvolve,
            ctx->timestep,
            ctx->sunGCDist,
            showCriterionT(ctx->criterion),
            ctx->theta,
            showBool(ctx->useQuad),
            mw_sqrt(ctx->eps2),
            seed
        );


    fprintf(f,
            "#\n"
            "# Potential: (%s)\n"
            "#\n",
            showExternalPotentialType(ctx->potentialType)
        );

    if (ctx->potentialType != EXTERNAL_POTENTIAL_DEFAULT)
    {
        return;
    }


    switch (p->disk.type)
    {
        case MiyamotoNagaiDisk:

            fprintf(f,
                    "# Disk: MiaymotoNagai\n"
                    "#   mass = %f\n"
                    "#   a = %f\n"
                    "#   b = %f\n"
                    "#\n",
                    p->disk.mass,
                    p->disk.scaleLength,
                    p->disk.scaleHeight);
            break;

        case ExponentialDisk:
            fprintf(f,
                    "# Disk: Exponential\n"
                    "#   mass = %f\n"
                    "#   b = %f\n"
                    "#\n",
                    p->disk.mass,
                    p->disk.scaleLength);
            break;

        case InvalidDisk:
        default:
            fprintf(f,
                    "# Disk: ???\n"
                    "#\n");
    }


    switch (p->halo.type)
    {
        case LogarithmicHalo:
            fprintf(f,
                    "# Halo: Logarithmic\n"
                    "#   vhalo = %f\n"
                    "#   d = %f\n"
                    "#\n",
                    p->halo.vhalo,
                    p->halo.scaleLength);
            break;

        case NFWHalo:
            fprintf(f,
                    "# Halo: NFW\n"
                    "#   vhalo = %f\n"
                    "#   q = %f\n"
                    "#\n",
                    p->halo.vhalo,
                    p->halo.scaleLength);
            break;

        case TriaxialHalo:
            fprintf(f,
                    "# Halo: Triaxial\n"
                    "#   vhalo = %f\n"
                    "#   rhalo = %f\n"
                    "#   qz = %f\n"
                    "#   q1 = %f\n"
                    "#   q2 = %f\n"
                    "#   phi = %f\n"
                    "#\n",
                    p->halo.vhalo,
                    p->halo.scaleLength,
                    p->halo.flattenZ,
                    p->halo.flattenX,
                    p->halo.flattenY,
                    p->halo.triaxAngle);
            break;

        case InvalidHalo:
        default:
            fprintf(f,
                    "# Halo: ???\n"
                    "#\n");
    }

    fprintf(f,
            "#\n"
            "# UseBin  Lambda  Probability  Error\n"
            "#\n"
            "\n"
        );
}

static void nbPrintHistogram(FILE* f,
                             const HistogramParams* hp,
                             const HistData* histData,
                             const unsigned int* histogram,
                             unsigned int maxIdx,
                             real start,
                             real totalNum,
                             int binUseIsValid)  /* If */
{
    unsigned int i;

    mw_boinc_print(f, "<histogram>\n");
    for (i = 0; i < maxIdx; ++i)
    {
        fprintf(f,
                "%d %12.10f %12.10f %12.10f\n",  /* Report center of the bins */
                binUseIsValid ? histData[i].useBin : TRUE,
                ((real) i  + 0.5) * hp->binSize + start,
                ((real) histogram[i]) / totalNum,
                histogram[i] == 0 ? inv(totalNum) : mw_sqrt(histogram[i]) / totalNum);
    }

    mw_boinc_print(f, "</histogram>\n");
}

static void nbWriteHistogram(const NBodyCtx* ctx,
                             NBodyState* st,
                             const NBodyFlags* nbf,
                             const HistogramParams* hp,
                             const HistData* histData,      /* Read histogram data */
                             const unsigned int* histogram, /* Binned simulation data */
                             unsigned int maxIdx,           /* number of bins */
                             real start,                    /* Calculated low point of bin range */
                             real chisq,
                             real totalNum)                 /* Total number in range */
{
    FILE* f = DEFAULT_OUTPUT_FILE;
    int binUseIsValid;

    if (nbf->histoutFileName && strcmp(nbf->histoutFileName, ""))  /* If file specified, try to open it */
    {
        f = mwOpenResolved(nbf->histoutFileName, "w+");
        if (f == NULL)
        {
            mwPerror("Writing histogram '%s' Using output file instead", nbf->histoutFileName);
            f = DEFAULT_OUTPUT_FILE;
        }
    }

    /* If we have a nan chisq, it means that either there was no
       histogram read to generate a chisq or the bins were mismatched
       in some way, so we can't tell which bins should be ignored
     */
    binUseIsValid = !isnan(chisq);
    nbPrintHistogramHeader(f, ctx, hp, st->nbody, nbf->seed, chisq);
    nbPrintHistogram(f, hp, histData, histogram, maxIdx, start, totalNum, binUseIsValid);

    if (f != DEFAULT_OUTPUT_FILE)
        fclose(f);
}

/* We create a raw histogram from the simulation.
   To compare it, we need to remove the unused bins
 */
static unsigned int nbCorrectTotalNumberInHistogram(const unsigned int* histogram,
                                                    unsigned int nBin,
                                                    unsigned int totalNum,
                                                    const HistData* histData)

{
    unsigned int i;

    for (i = 0; i < nBin; ++i)
    {
        if (!histData[i].useBin)
        {
            totalNum -= histogram[i];
        }
    }

    return totalNum;
}

/*
Takes a treecode position, converts it to (l,b), then to (lambda,
beta), and then constructs a histogram of the density in lambda.

Then calculates the cross correlation between the model histogram and
the data histogram A maximum correlation means the best fit */

/* Bin the bodies from the simulation into maxIdx bins.
   Returns null on failure
 */
static unsigned int* nbCreateHistogram(const NBodyCtx* ctx,       /* Simulation context */
                                       const NBodyState* st,      /* Final state of the simulation */
                                       const unsigned int maxIdx, /* Total number of bins */
                                       const real start,          /* Calculated start point of bin range */
                                       const HistogramParams* hp,
                                       unsigned int* totalNumOut) /* Out: Number of particles in range */
{
    real lambda;
    unsigned int idx;
    unsigned int totalNum = 0;
    Body* p;
    unsigned int* histogram;
    NBHistTrig histTrig;
    const Body* endp = st->bodytab + st->nbody;

    nbGetHistTrig(&histTrig, hp);

    histogram = (unsigned int*) mwCalloc(maxIdx, sizeof(unsigned int));

    for (p = st->bodytab; p < endp; ++p)
    {
        /* Only include bodies in models we aren't ignoring */
        if (ignoreBody(p))
            continue;

        lambda = nbXYZToLambda(&histTrig, Pos(p), ctx->sunGCDist);
        idx = (unsigned int) mw_floor((lambda - start) / hp->binSize);
        if (idx < maxIdx)
        {
            ++histogram[idx];
            ++totalNum;
        }
    }

    *totalNumOut = totalNum;
    return histogram;
}


/* The chisq is calculated by reading a histogram file of normalized data.
   Returns null on failure.
 */
static HistData* nbReadHistData(const char* histogram, unsigned int* maxIdxOut)
{
    FILE* f;
    int rc = 0;
    size_t fsize;
    HistData* histData;
    unsigned int fileCount = 0;
    unsigned int lineNum = 0;
    mwbool error = FALSE;
    mwbool readParams = FALSE;
    char lineBuf[1024];

    f = mwOpenResolved(histogram, "r");
    if (f == NULL)
    {
        mw_printf("Opening histogram file '%s'\n", histogram);
        return NULL;
    }

    fsize = mwCountLinesInFile(f);
    if (fsize == 0)
    {
        mw_printf("Histogram line count = 0\n");
        return NULL;
    }

    histData = (HistData*) mwCalloc(sizeof(HistData), fsize);

    while (fgets(lineBuf, (int) sizeof(lineBuf), f))
    {
        ++lineNum;

        if (strlen(lineBuf) + 1 >= sizeof(lineBuf))
        {
            mw_printf("Error reading histogram line %d (Line buffer too small): %s", lineNum, lineBuf);
            error = TRUE;
            break;
        }

        /* Skip comments and blank lines */
        if (lineBuf[0] == '#' || lineBuf[0] == '\n')
            continue;

        if (!readParams)  /* One line is allowed for information on the histogram */
        {
            real phi, theta, psi;

            rc = sscanf(lineBuf,
                        DOUBLEPREC ? " phi = %lf , theta = %lf , psi = %lf \n" : " phi = %f , theta = %f , psi = %f \n",
                        &phi, &theta, &psi);
            if (rc == 3)
            {
                readParams = TRUE;
                continue;
            }
        }

        rc = sscanf(lineBuf,
                    DOUBLEPREC ? "%d %lf %lf %lf \n" : "%d %f %f %f \n",
                    &histData[fileCount].useBin,
                    &histData[fileCount].lambda,
                    &histData[fileCount].count,
                    &histData[fileCount].err);
        if (rc != 4)
        {
            mw_printf("Error reading histogram line %d: %s", lineNum, lineBuf);
            error = TRUE;
            break;
        }

        ++fileCount;
    }

    fclose(f);

    if (error)
    {
        free(histData);
        return NULL;
    }

    if (maxIdxOut)
    {
        *maxIdxOut = fileCount;
    }

    return histData;
}

/* Calculate the likelihood from the final state of the simulation */
real nbChisq(const NBodyCtx* ctx, NBodyState* st, const NBodyFlags* nbf, const HistogramParams* hp)
{
    real chisqval = NAN;
    unsigned int totalNum = 0;
    unsigned int* histogram = NULL;
    HistData* histData = NULL;
    unsigned int dataMaxIdx = 0;

    /* Calculate the bounds of the bin range, making sure to use a
     * fixed bin size which spans the entire range, and is symmetric
     * around 0 */
    const real rawCount = (hp->endRaw - hp->startRaw) / hp->binSize;
    const unsigned int maxIdx = (unsigned int) ceil(rawCount);

    const real start = mw_ceil(hp->center - hp->binSize * (real) maxIdx / 2.0);

    histogram = nbCreateHistogram(ctx, st, maxIdx, start, hp, &totalNum);
    if (!histogram)
    {
        mw_printf("Failed to create histogram\n");
        return NAN;
    }

    if (nbf->histogramFileName) /* If we have an input histogram to match */
    {
        histData = nbReadHistData(nbf->histogramFileName, &dataMaxIdx);
        if (!histData)
        {
            mw_printf("Failed to read histogram\n");
        }
        else
        {
            if (dataMaxIdx != maxIdx)
            {
                mw_printf("Number of bins does not match those in histogram file. "
                          "Expected %u, got %u\n",
                          maxIdx,
                          dataMaxIdx);
            }
            else
            {
                real effTotalNum;

                if (totalNum != 0)
                {
                    effTotalNum = (real) nbCorrectTotalNumberInHistogram(histogram, maxIdx, totalNum, histData);
                    chisqval = nbCalcChisq(histData, histogram, maxIdx, effTotalNum);
                }
                else
                {
                    chisqval = -INFINITY;
                }
            }
        }
    }

    /* We want to write something whether or not the likelihood can be
     * calculated (i.e. given a histogram) */
    if (nbf->printHistogram)
    {
        nbWriteHistogram(ctx, st, nbf, hp, histData, histogram, maxIdx, start, chisqval, (real) totalNum);
    }

    free(histogram);
    free(histData);

    return chisqval;
}

