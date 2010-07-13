// A BOINC compatible fitness calculation routine for the Orphan FILE* density profile project
// B. Willett Feb. 25, 2010
// Adapted for use with B&H treecode May 11, 2010


// Takes a treecode position, converts it to (l,b), then to (lambda, beta), and then constructs a histogram of the density in lambda.

// Then calculates the cross correlation between the model histogram and the data histogram
// A maximum correlation means the best fit

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "nbody_priv.h"
#include "chisq.h"

/* FIXME: Magic numbers */
#define phi d2r(128.79)
#define theta d2r(54.39)
#define psi d2r(90.70)
#define startRaw ((real) -50.0)
#define endRaw ((real) 50.0)
#define binsize ((real) 3.0)
#define center ((real) 0.0)

typedef struct
{
    real lambda;
    real err;
    real count;
} HistData;

real chisq(const NBodyCtx* ctx, NBodyState* st)
{
    real chisqval = 0.0;
    unsigned int totalNum = 0;
    unsigned int i, idx;

    const int nbody = ctx->model.nbody;
    const bodyptr endp = st->bodytab + nbody;
    bodyptr p;
    FILE* f;

    const real rawCount = (endRaw - startRaw) / binsize;

    const unsigned int maxIdx = (unsigned int) ceil(rawCount);
    const real maxIdxf = (real) maxIdx;

    const real start = rceil(center - binsize * maxIdxf / 2.0);
    const real end   = rfloor(center + binsize * maxIdxf / 2.0);

    printf("start = %g  end = %g\n", start, end);

    unsigned int* histogram = callocSafe(maxIdx, sizeof(unsigned int));

    real bcos, bsin, lsin, lcos;
    vector lbr;
    real lambda;
    real tmp;
    const real cosphi = rcos(phi);
    const real sinphi = rsin(phi);
    const real sinpsi = rsin(psi);
    const real cospsi = rcos(psi);
    const real costh  = rcos(theta);
    const real sinth  = rsin(theta);

    printf("Binning simulation results.\n");
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

        /* CHECKME: Do we need fma here? */
        lambda = r2d(atan2(
                         - (sinpsi * cosphi + costh * sinphi * cospsi) * bcos * lcos
                         + (-sinpsi * sinphi + costh * cosphi * cospsi) * bcos * lsin
                         + cospsi * sinth * bsin,

                           (cospsi * cosphi - costh * sinphi * sinpsi) * bcos * lcos
                         + (cospsi * sinphi + costh * cosphi * sinpsi) * bcos * lsin
                         + sinpsi * sinth * bsin ));

        //printf("L = %g\n", lambda);

        idx = (unsigned int) rfloor((lambda - start) / binsize);
        if (idx < maxIdx)
        {
            ++histogram[idx];
            ++totalNum;
        }
    }

    printf("Total num in range: %u\n", totalNum);

    f = ctx->histout ? nbody_fopen(ctx->histout, "w") : stdout;
    if (f == NULL)
    {
        perror("Writing histout");
        return NAN;
    }

    const real totalNumf = (real) totalNum;
    // Print out the histogram
    for (i = 0; i < maxIdx; ++i)
    {
        fprintf(f, "%2.10f %2.10f %2.10f\n",  /* Report center of the bins */
                fma(((real) i  + 0.5), binsize, start),
                ((real) histogram[i]) / totalNumf,
                sqrt(histogram[i]) / totalNumf
            );
    }

    if (f != stdout)
        fclose(f);

    /* Calculate the chisq value by reading a file.
       (the real data histogram should already be normalized) */
    f = nbody_fopen(ctx->histogram, "r");
    if (f == NULL)
    {
        perror("Opening histogram");
        return NAN;
    }

    size_t fsize;
    unsigned int filecount = 0;

    fseek(f, 0L, SEEK_END);
    fsize = ceil((real) (ftell(f) + 1) / 3);  /* Make sure it's big enough, avoid fun with integer division */
    fseek(f, 0L, SEEK_SET);

    HistData* histData = callocSafe(sizeof(HistData), fsize);

    int rc = 0;
    while ( (rc = fscanf(f,
                         #ifdef DOUBLEPREC
                           "%lf %lf %lf\n",
                         #else
                         "%f %f %f\n",
                         #endif
                         &histData[filecount].lambda,
                         &histData[filecount].count,
                         &histData[filecount].err))
            && rc != EOF)
    {
        ++filecount;
    }

    if (f != stdout)
        fclose(f);

    if (filecount != maxIdx)
    {
        warn("Number of bins does not match those in histogram file. "
             "Expected %u, got %u\n",
             maxIdx,
             filecount);
        free(histData);
        free(histogram);
        return NAN;
    }

    // Calculate the chisq
    for (i = 0; i < maxIdx; ++i)
    {
        tmp = (histData[i].count - ((real) histogram[i] / totalNumf)) / histData[i].err;
        chisqval += sqr(tmp);
    }

    free(histData);
    free(histogram);

    // MAXIMUM likelihood, multiply by -1
    return -chisqval;
}

