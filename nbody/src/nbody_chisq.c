/*
 *  Copyright (c) 2010-2011 Rensselaer Polytechnic Institute
 *  Copyright (c) 2010-2011 Ben Willett
 *  Copyright (c) 2010-2011 Matthew Arsenault
 *
 *  This file is part of Milkway@Home.
 *
 *  Milkway@Home is free software: you may copy, redistribute and/or modify it
 *  under the terms of the GNU General Public License as published by the
 *  Free Software Foundation, either version 3 of the License, or (at your
 *  option) any later version.
 *
 *  This file is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "nbody_priv.h"
#include "nbody_chisq.h"
#include "milkyway_util.h"


typedef enum
{
    NBODY_LIKELIHOOD,
    NBODY_ALT_LIKELIHOOD
} NBodyLikelihoodMethod;


/* From the range of a histogram, find the number of bins */
static unsigned int nbHistogramNBin(const HistogramParams* hp)
{
    double rawCount = (hp->endRaw - hp->startRaw) / hp->binSize;
    return (unsigned int) ceil(rawCount);   /* Total number of bins */
}

/* Find the corrected starting point for the histogram */
static real nbHistogramStart(const HistogramParams* hp)
{
    unsigned int nBin = nbHistogramNBin(hp);
    return mw_ceil(hp->center - hp->binSize * (real) nBin / 2.0);
}


/* We create a raw histogram from the simulation.
   To compare it, we need to remove the unused bins totals.

   Return the total number in the generated histogram with the counts
   subtracted from the bins corresponding to ignored bins in the data
   histogram
 */
static unsigned int nbCorrectTotalNumberInHistogram(const NBodyHistogram* histogram, /* Generated histogram */
                                                    const NBodyHistogram* data)      /* Data histogram */
{
    unsigned int i;
    unsigned int nBin = data->nBin;
    unsigned int totalNum = histogram->totalNum;

    assert(histogram->hasRawCounts);
    assert(data->nBin == histogram->nBin);

    for (i = 0; i < nBin; ++i)
    {
        if (!data->data[i].useBin)
        {
            totalNum -= histogram->data[i].rawCount;
        }
    }

    return totalNum;
}

/* Calculate chisq from read data histogarm and the generated histogram */
static real nbCalcChisq(const NBodyHistogram* data,        /* Data histogram */
                        const NBodyHistogram* histogram,   /* Generated histogram */
                        NBodyLikelihoodMethod method)
{
    unsigned int i;
    real tmp;
    real effTotalNum;
    real chiSq = 0.0;
    real n;
    real err;
    unsigned int nBin = data->nBin;


    assert(histogram->hasRawCounts);
    assert(nBin == histogram->nBin);

    if (histogram->totalNum == 0)
    {
        return -INFINITY;
    }

    effTotalNum = (real) nbCorrectTotalNumberInHistogram(histogram, data);
    for (i = 0; i < nBin; ++i)
    {
        if (data->data[i].useBin)  /* Skip bins with missing data */
        {
            n = (real) histogram->data[i].rawCount;
            err = data->data[i].err;

            switch (method)
            {
                case NBODY_LIKELIHOOD:
                    tmp = (data->data[i].count - (n / effTotalNum)) / err;
                    chiSq += sqr(tmp);
                    break;

                case NBODY_ALT_LIKELIHOOD:
                    /* effective error = sqrt( (data error)^2 + (sim count error)^2 ) */
                    err = sqrt(sqr(err) + n / effTotalNum);
                    tmp = (data->data[i].count - (n / effTotalNum)) / err;
                    chiSq += sqr(tmp);
                    break;

                default:
                    mw_fail("Invalid likelihood method\n");
            }
        }
    }

    return -chiSq;     /* MAXIMUM likelihood, multiply by -1 */
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

static void nbPrintHistogram(FILE* f, const NBodyHistogram* histogram)
{
    unsigned int i;
    const HistData* data;
    unsigned int nBin = histogram->nBin;

    mw_boinc_print(f, "<histogram>\n");
    for (i = 0; i < nBin; ++i)
    {
        data = &histogram->data[i];
        fprintf(f,
                "%d %12.10f %12.10f %12.10f\n",
                data->useBin,
                data->lambda,
                data->count,
                data->err);
    }

    mw_boinc_print(f, "</histogram>\n");
}

static void nbWriteHistogram(const NBodyCtx* ctx,
                             NBodyState* st,
                             const NBodyFlags* nbf,
                             NBodyHistogram* histogram,
                             real chisq)
{
    FILE* f = DEFAULT_OUTPUT_FILE;

    if (nbf->histoutFileName && strcmp(nbf->histoutFileName, ""))  /* If file specified, try to open it */
    {
        f = mwOpenResolved(nbf->histoutFileName, "w+");
        if (f == NULL)
        {
            mwPerror("Writing histogram '%s' Using output file instead", nbf->histoutFileName);
            f = DEFAULT_OUTPUT_FILE;
        }
    }

    nbPrintHistogramHeader(f, ctx, &histogram->params, st->nbody, nbf->seed, chisq);
    nbPrintHistogram(f, histogram);

    if (f != DEFAULT_OUTPUT_FILE)
        fclose(f);
}


/* Get normalized histogram counts and errors */
static void nbNormalizeHistogram(NBodyHistogram* histogram)
{
    unsigned int i;
    real count;

    unsigned int nBin = histogram->nBin;
    const HistogramParams* hp = &histogram->params;
    real totalNum = (real) histogram->totalNum;
    HistData* histData = histogram->data;
    real start = nbHistogramStart(&histogram->params);


    for (i = 0; i < nBin; ++i)
    {
        count = (real) histData[i].rawCount;
        histData[i].lambda = ((real) i  + 0.5) * hp->binSize + start;  /* Report center of the bins */
        histData[i].count = count / totalNum;
        histData[i].err = (histData[i].rawCount == 0) ? inv(totalNum) : mw_sqrt(count) / totalNum;
    }
}

/*
Takes a treecode position, converts it to (l,b), then to (lambda,
beta), and then constructs a histogram of the density in lambda.

Then calculates the cross correlation between the model histogram and
the data histogram A maximum correlation means the best fit */

/* Bin the bodies from the simulation into maxIdx bins.
   Returns null on failure
 */
static NBodyHistogram* nbCreateHistogram(const NBodyCtx* ctx,        /* Simulation context */
                                         const NBodyState* st,       /* Final state of the simulation */
                                         const HistogramParams* hp)  /* Range of histogram to create */
{
    real lambda;
    unsigned int i;
    unsigned int idx;
    unsigned int totalNum = 0;
    Body* p;
    NBodyHistogram* histogram;
    HistData* histData;
    NBHistTrig histTrig;
    const Body* endp = st->bodytab + st->nbody;

    real start = nbHistogramStart(hp);
    unsigned int nBin = nbHistogramNBin(hp);


    nbGetHistTrig(&histTrig, hp);

    histogram = mwCalloc(sizeof(NBodyHistogram) + nBin * sizeof(HistData), sizeof(char));
    histogram->nBin = nBin;
    histogram->hasRawCounts = TRUE;
    histogram->params = *hp;
    histData = histogram->data;

    /* It does not make sense to ignore bins in a generated histogram */
    for (i = 0; i < nBin; ++i)
    {
        histData[i].useBin = TRUE;
    }


    for (p = st->bodytab; p < endp; ++p)
    {
        /* Only include bodies in models we aren't ignoring */
        if (!ignoreBody(p))
        {
            lambda = nbXYZToLambda(&histTrig, Pos(p), ctx->sunGCDist);
            idx = (unsigned int) mw_floor((lambda - start) / hp->binSize);
            if (idx < nBin)
            {
                histData[idx].rawCount++;
                ++totalNum;
            }
        }
    }

    histogram->totalNum = totalNum; /* Total particles in range */

    nbNormalizeHistogram(histogram);

    return histogram;
}


/* The chisq is calculated by reading a histogram file of normalized data.
   Returns null on failure.
 */
static NBodyHistogram* nbReadHistogram(const char* histogramFile)
{
    FILE* f;
    int rc = 0;
    size_t fsize = 0;
    NBodyHistogram* histogram = NULL;
    HistData* histData = NULL;
    unsigned int fileCount = 0;
    unsigned int lineNum = 0;
    mwbool error = FALSE;
    mwbool readParams = FALSE;
    char lineBuf[1024];

    f = mwOpenResolved(histogramFile, "r");
    if (f == NULL)
    {
        mw_printf("Opening histogram file '%s'\n", histogramFile);
        return NULL;
    }

    fsize = mwCountLinesInFile(f);
    if (fsize == 0)
    {
        mw_printf("Histogram line count = 0\n");
        return NULL;
    }

    histogram = (NBodyHistogram*) mwCalloc(sizeof(NBodyHistogram) + fsize * sizeof(HistData), sizeof(char));
    histogram->hasRawCounts = FALSE;     /* Do we want to include these? */
    histData = histogram->data;


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
        free(histogram);
        return NULL;
    }


    histogram->nBin = fileCount;
    return histogram;
}


/* Calculate the likelihood from the final state of the simulation */
real nbChisq(const NBodyCtx* ctx, NBodyState* st, const NBodyFlags* nbf)
{
    real chiSq = NAN;
    real altChiSq = NAN;
    NBodyHistogram* data = NULL;
    NBodyHistogram* histogram = NULL;
    lua_State* luaSt = NULL;
    HistogramParams hp;


    luaSt = nbOpenLuaStateWithScript(nbf);
    if (!luaSt)
    {
        return NAN;
    }

    if (nbEvaluateHistogramParams(luaSt, &hp))
    {
        lua_close(luaSt);
        return NAN;
    }

    lua_close(luaSt);



    /* Calculate the bounds of the bin range, making sure to use a
     * fixed bin size which spans the entire range, and is symmetric
     * around 0 */

    histogram = nbCreateHistogram(ctx, st, &hp);
    if (!histogram)
    {
        mw_printf("Failed to create histogram\n");
        return NAN;
    }

    if (nbf->histogramFileName) /* If we have an input histogram to match */
    {
        data = nbReadHistogram(nbf->histogramFileName);
        if (!data)
        {
            mw_printf("Failed to read histogram\n");
        }
        else
        {
            if (data->nBin != histogram->nBin)
            {
                mw_printf("Number of bins does not match those in histogram file. "
                          "Expected %u, got %u\n",
                          histogram->nBin,
                          data->nBin);
            }
            else
            {
                chiSq = nbCalcChisq(data, histogram, NBODY_LIKELIHOOD);
                altChiSq = nbCalcChisq(data, histogram, NBODY_ALT_LIKELIHOOD);
            }
        }
    }

    mw_printf("<alt_likelihood>%.15f</alt_likelihood>\n", altChiSq);

    /* We want to write something whether or not the likelihood can be
     * calculated (i.e. given a histogram) */
    if (nbf->printHistogram)
    {
        nbWriteHistogram(ctx, st, nbf, histogram, chiSq);
    }

    free(data);
    free(histogram);

    return chiSq;
}

