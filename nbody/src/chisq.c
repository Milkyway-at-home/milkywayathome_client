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
#include "nbody_util.h"
#include "chisq.h"

/* FIXME: sort of magic numbers */
#define phi d2r(128.79)
#define theta d2r(54.39)
#define psi d2r(90.70)
#define startRaw ((real) -50)
#define endRaw ((real) 50)
#define binsize ((real) 2.9411764705882355)
#define center ((real) 0.0)

typedef struct
{
    int useBin;
    real lambda;
    real err;
    real count;
} HistData;


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

inline static void printHistogram(FILE* f,
                                  const HistData* histData,
                                  const unsigned int* histogram,
                                  const unsigned int maxIdx,
                                  const real start,
                                  const real totalNum)
{
    unsigned int i;

  #if BOINC_APPLICATION
    fprintf(f, "<histogram>\n");
  #endif

    for (i = 0; i < maxIdx; ++i)
    {
        fprintf(f, "%d %2.10f %2.10f %2.10f\n",  /* Report center of the bins */
                histData[i].useBin,
                ((real) i  + 0.5) * binsize + start,
                ((real) histogram[i]) / totalNum,
                histogram[i] == 0 ? inv(totalNum) : sqrt(histogram[i]) / totalNum);
    }

  #if BOINC_APPLICATION
    fprintf(f, "</histogram>\n");
  #endif

}

static void writeHistogram(const char* histout,           /* Filename to write histogram to */
                           const HistData* histData,      /* Read histogram data */
                           const unsigned int* histogram, /* Binned simulation data */
                           const unsigned int maxIdx,     /* number of bins */
                           const real start,              /* Calculated low point of bin range */
                           const real totalNum)           /* Total number in range */
{
    FILE* f = DEFAULT_OUTPUT_FILE;

    if (histout && strcmp(histout, ""))  /* If file specified, try to open it */
    {
        f = nbodyOpenResolved(histout, "w");
        if (f == NULL)
        {
            perror("Writing histout. Using stderr instead");
            f = DEFAULT_OUTPUT_FILE;
        }
    }

    printHistogram(f, histData, histogram, maxIdx, start, totalNum);

    if (f != DEFAULT_OUTPUT_FILE)
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
                                     unsigned int* totalNumOut) /* Out: Number of particles in range */
{
    real lambda;
    real bcos, bsin, lsin, lcos;
    vector lbr;
    unsigned int idx;
    unsigned int totalNum = 0;
    bodyptr p;
    const real cosphi = rcos(phi);
    const real sinphi = rsin(phi);
    const real sinpsi = rsin(psi);
    const real cospsi = rcos(psi);
    const real costh  = rcos(theta);
    const real sinth  = rsin(theta);

    unsigned int* histogram = callocSafe(maxIdx, sizeof(unsigned int));

    const bodyptr endp = st->bodytab + ctx->model.nbody;

    for (p = st->bodytab; p < endp; ++p)
    {
        // Convert to (l,b) (involves convert x to Sun-centered)
        // Leave in radians to make rotation easier
        cartesianToLbr_rad(ctx, lbr, Pos(p));

        // Convert to (lambda, beta) (involves a rotation using the
        // Newberg et al (2009) rotation matrices)
        bcos = rcos(B(lbr));
        bsin = rsin(B(lbr));
        lsin = rsin(L(lbr));
        lcos = rcos(L(lbr));

        lambda = r2d(atan2(
                         - (sinpsi * cosphi + costh * sinphi * cospsi) * bcos * lcos
                         + (-sinpsi * sinphi + costh * cosphi * cospsi) * bcos * lsin
                         + cospsi * sinth * bsin,

                           (cospsi * cosphi - costh * sinphi * sinpsi) * bcos * lcos
                         + (cospsi * sinphi + costh * cosphi * sinpsi) * bcos * lsin
                         + sinpsi * sinth * bsin ));

        idx = (unsigned int) rfloor((lambda - start) / binsize);
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
static HistData* readHistData(const char* histogram, const unsigned int maxIdx)
{
    FILE* f;
    int rc = 0;
    size_t fsize;
    HistData* histData;
    unsigned int fileCount = 0;

    f = nbodyOpenResolved(histogram, "r");
    if (f == NULL)
    {
        perror("Opening histogram");
        return NULL;
    }

    fseek(f, 0L, SEEK_END);
    /* Make sure it's big enough, avoid fun with integer division */
    fsize = ceil((real) (ftell(f) + 1) / 3);
    fseek(f, 0L, SEEK_SET);

    histData = callocSafe(sizeof(HistData), fsize);

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
real nbodyChisq(const NBodyCtx* ctx, const NBodyState* st)
{
    real chisqval;
    unsigned int totalNum = 0;
    unsigned int* histogram;
    HistData* histData;

    /* Calculate the bounds of the bin range, making sure to use a
     * fixed bin size which spans the entire range, and is symmetric
     * around 0 */
    const real rawCount = (endRaw - startRaw) / binsize;
    const unsigned int maxIdx = (unsigned int) ceil(rawCount);

    const real start = rceil(center - binsize * (real) maxIdx / 2.0);

    histogram = createHistogram(ctx, st, maxIdx, start, &totalNum);
    histData  = readHistData(ctx->histogram, maxIdx);

    if (histogram && histData)
    {
        if (ctx->outputHistogram)
            writeHistogram(ctx->histout, histData, histogram, maxIdx, start, (real) totalNum);

        if (totalNum != 0)
            chisqval = calcChisq(histData, histogram, maxIdx, (real) totalNum);
        else
            chisqval = -INFINITY;
    }
    else
        chisqval = NAN;

    free(histData);
    free(histogram);

    return chisqval;
}

