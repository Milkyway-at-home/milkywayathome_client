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
static real calcChisq(const HistData* histData,
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

static inline void printHistogram(FILE* f,
                                  const HistogramParams* hp,
                                  const HistData* histData,
                                  const unsigned int* histogram,
                                  const unsigned int maxIdx,
                                  const real start,
                                  const real totalNum)
{
    unsigned int i;

    mw_boinc_print(f, "<histogram>\n");
    for (i = 0; i < maxIdx; ++i)
    {
        fprintf(f, "%d %2.10f %2.10f %2.10f\n",  /* Report center of the bins */
                histData[i].useBin,
                ((real) i  + 0.5) * hp->binSize + start,
                ((real) histogram[i]) / totalNum,
                histogram[i] == 0 ? inv(totalNum) : mw_sqrt(histogram[i]) / totalNum);
    }

    mw_boinc_print(f, "</histogram>\n");
}

static void writeHistogram(NBodyState* st,
                           const NBodyFlags* nbf,
                           const HistogramParams* hp,
                           const HistData* histData,      /* Read histogram data */
                           const unsigned int* histogram, /* Binned simulation data */
                           const unsigned int maxIdx,     /* number of bins */
                           const real start,              /* Calculated low point of bin range */
                           const real totalNum)           /* Total number in range */
{
    FILE* f = st->outFile;

    if (nbf->histoutFileName && strcmp(nbf->histoutFileName, ""))  /* If file specified, try to open it */
    {
        f = mwOpenResolved(nbf->histoutFileName, "w");
        if (f == NULL)
        {
            perror("Writing histout. Using output file instead");
            f = st->outFile;
        }
    }

    printHistogram(f, hp, histData, histogram, maxIdx, start, totalNum);

    if (f != st->outFile)
        fclose(f);
}

/*
Takes a treecode position, converts it to (l,b), then to (lambda,
beta), and then constructs a histogram of the density in lambda.

Then calculates the cross correlation between the model histogram and
the data histogram A maximum correlation means the best fit */

/* Bin the bodies from the simulation into maxIdx bins.
   Returns null on failure
 */
static unsigned int* createHistogram(const NBodyCtx* ctx,       /* Simulation context */
                                     const NBodyState* st,      /* Final state of the simulation */
                                     const unsigned int maxIdx, /* Total number of bins */
                                     const real start,          /* Calculated start point of bin range */
                                     const HistogramParams* hp,
                                     const HistData* histData,  /* Data histogram; which bins to skip */
                                     unsigned int* totalNumOut) /* Out: Number of particles in range */
{
    real lambda;
    real bcos, bsin, lsin, lcos;
    mwvector lbr;
    unsigned int idx;
    unsigned int totalNum = 0;
    Body* p;
    unsigned int* histogram;

    real rphi = d2r(hp->phi);
    real rpsi = d2r(hp->psi);
    real rth  = d2r(hp->theta);

    const real cosphi = mw_cos(rphi);
    const real sinphi = mw_sin(rphi);
    const real sinpsi = mw_sin(rpsi);
    const real cospsi = mw_cos(rpsi);
    const real costh  = mw_cos(rth);
    const real sinth  = mw_sin(rth);

    const Body* endp = st->bodytab + st->nbody;
    histogram = (unsigned int*) mwCalloc(maxIdx, sizeof(unsigned int));

    for (p = st->bodytab; p < endp; ++p)
    {
        /* Only include bodies in models we aren't ignoring */
        if (ignoreBody(p))
            continue;

        // Convert to (l,b) (involves convert x to Sun-centered)
        // Leave in radians to make rotation easier
        lbr = cartesianToLbr_rad(Pos(p), ctx->sunGCDist);

        // Convert to (lambda, beta) (involves a rotation using the
        // Newberg et al (2009) rotation matrices)
        bcos = mw_cos(B(lbr));
        bsin = mw_sin(B(lbr));
        lsin = mw_sin(L(lbr));
        lcos = mw_cos(L(lbr));

        lambda = r2d(mw_atan2(
                         - (sinpsi * cosphi + costh * sinphi * cospsi) * bcos * lcos
                         + (-sinpsi * sinphi + costh * cosphi * cospsi) * bcos * lsin
                         + cospsi * sinth * bsin,

                           (cospsi * cosphi - costh * sinphi * sinpsi) * bcos * lcos
                         + (cospsi * sinphi + costh * cosphi * sinpsi) * bcos * lsin
                         + sinpsi * sinth * bsin ));

        idx = (unsigned int) mw_floor((lambda - start) / hp->binSize);
        if (idx < maxIdx && histData[idx].useBin)
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
static HistData* readHistData(const char* histogram, const unsigned int maxIdx)
{
    FILE* f;
    int rc = 0;
    size_t fsize;
    HistData* histData;
    unsigned int fileCount = 0;

    f = mwOpenResolved(histogram, "r");
    if (f == NULL)
    {
        perror("Opening histogram");
        return NULL;
    }

    fseek(f, 0L, SEEK_END);
    /* Make sure it's big enough, avoid fun with integer division */
    fsize = (size_t) ceil((real) (ftell(f) + 1) / 3);
    fseek(f, 0L, SEEK_SET);

    histData = (HistData*) mwCalloc(sizeof(HistData), fsize);

    while ( (rc = fscanf(f,
                       #if DOUBLEPREC
                         "%d %lf %lf %lf\n",
                       #else
                         "%d %f %f %f\n",
                       #endif
                         &histData[fileCount].useBin,
                         &histData[fileCount].lambda,
                         &histData[fileCount].count,
                         &histData[fileCount].err))
            && rc != EOF)
    {
        ++fileCount;
    }

    fclose(f);

    if (fileCount != maxIdx)
    {
        warn("Number of bins does not match those in histogram file. "
             "Expected %u, got %u\n",
             maxIdx,
             fileCount);
        free(histData);
        return NULL;
    }

    return histData;
}

/* Calculate the likelihood from the final state of the simulation */
real nbodyChisq(const NBodyCtx* ctx, NBodyState* st, const NBodyFlags* nbf, const HistogramParams* hp)
{
    real chisqval;
    unsigned int totalNum = 0;
    unsigned int* histogram;
    HistData* histData;

    /* Calculate the bounds of the bin range, making sure to use a
     * fixed bin size which spans the entire range, and is symmetric
     * around 0 */
    const real rawCount = (hp->endRaw - hp->startRaw) / hp->binSize;
    const unsigned int maxIdx = (unsigned int) ceil(rawCount);

    const real start = mw_ceil(hp->center - hp->binSize * (real) maxIdx / 2.0);

    histData = readHistData(nbf->histogramFileName, maxIdx);
    if (!histData)
    {
        warn("Failed to read histogram\n");
        return NAN;
    }

    histogram = createHistogram(ctx, st, maxIdx, start, hp, histData, &totalNum);
    if (!histogram)
    {
        warn("Failed to create histogram\n");
        return NAN;
    }

    if (nbf->printHistogram)
        writeHistogram(st, nbf, hp, histData, histogram, maxIdx, start, (real) totalNum);

    if (totalNum != 0)
        chisqval = calcChisq(histData, histogram, maxIdx, (real) totalNum);
    else
        chisqval = -INFINITY;

    free(histData);
    free(histogram);

    return chisqval;
}

