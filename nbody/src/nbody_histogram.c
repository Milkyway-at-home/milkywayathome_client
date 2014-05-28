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

#include "nbody_config.h"

#include "nbody_histogram.h"
#include "nbody_chisq.h"
#include "nbody_emd.h"
#include "nbody_mass.h"
#include "nbody_defaults.h"
#include "milkyway_util.h"
#include "nbody_coordinates.h"
#include "nbody_show.h"

/* From the range of a histogram, find the number of bins */
static unsigned int nbHistogramNBin(const HistogramParams* hp)
{
    real rawCount = (hp->endRaw - hp->startRaw) / hp->binSize;
    return (unsigned int) mw_ceil(rawCount);   /* Total number of bins */
}

/* Find the corrected starting point for the histogram */
static double nbHistogramStart(const HistogramParams* hp)
{
    unsigned int nBin = nbHistogramNBin(hp);
    return (double) mw_ceil(hp->center - hp->binSize * (real) nBin / 2.0);
}

double nbNormalizedHistogramError(unsigned int n, double total)
{
    return (n == 0) ? inv(total) : sqrt((double) n) / total;
}

double nbCorrectRenormalizedInHistogram(const NBodyHistogram* histogram, const NBodyHistogram* data)
{
    unsigned int i;
    unsigned int nBin = data->nBin;
    double total = 0.0;

    for (i = 0; i < nBin; ++i)
    {
        if (data->data[i].useBin)
        {
            total += histogram->data[i].count;
        }
    }

    return total;
}

/* We create a raw histogram from the simulation.
   To compare it, we need to remove the unused bins totals.

   Return the total number in the generated histogram with the counts
   subtracted from the bins corresponding to ignored bins in the data
   histogram
 */
unsigned int nbCorrectTotalNumberInHistogram(const NBodyHistogram* histogram, /* Generated histogram */
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

static void nbPrintHistogramHeader(FILE* f,
                                   const NBodyCtx* ctx,
                                   const HistogramParams* hp,
                                   int nbody)
{
    char tBuf[256];
    const Potential* p = &ctx->pot;

    mwLocalTimeFull(tBuf, sizeof(tBuf));

    fprintf(f,
            "#\n"
            "# Generated %s\n"
            "#\n"
            "# (phi,   theta,  psi) = (%f, %f, %f)\n"
            "# (start, center, end) = (%f, %f, %f)\n"
            "# Bin size = %f\n"
            "#\n"
            "#\n",
            tBuf,
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
            "#\n",
            nbody,
            ctx->timeEvolve,
            ctx->timestep,
            ctx->sunGCDist,
            showCriterionT(ctx->criterion),
            ctx->theta,
            showBool(ctx->useQuad),
            mw_sqrt(ctx->eps2)
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

/* Print the histogram without a header. */
void nbPrintHistogram(FILE* f, const NBodyHistogram* histogram)
{
    unsigned int i;
    const HistData* data;
    unsigned int nBin = histogram->nBin;

    mw_boinc_print(f, "<histogram>\n");
    fprintf(f, "n = %u\n", histogram->totalNum);
    fprintf(f, "massPerParticle = %12.10f\n", histogram->massPerParticle);
    fprintf(f, "totalSimulated = %u\n", histogram->totalSimulated);

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

/* Write histogram to given file name along with descriptive header */
void nbWriteHistogram(const char* histoutFileName,
                      const NBodyCtx* ctx,
                      const NBodyState* st,
                      const NBodyHistogram* histogram)
{
    FILE* f = DEFAULT_OUTPUT_FILE;

    if (histoutFileName && strcmp(histoutFileName, ""))  /* If file specified, try to open it */
    {
        f = mwOpenResolved(histoutFileName, "w+");
        if (f == NULL)
        {
            mwPerror("Error opening histogram '%s'. Using default output instead.", histoutFileName);
            f = DEFAULT_OUTPUT_FILE;
        }
    }

    nbPrintHistogramHeader(f, ctx, &histogram->params, st->nbody);
    nbPrintHistogram(f, histogram);

    if (f != DEFAULT_OUTPUT_FILE)
        fclose(f);
}

/* Get normalized histogram counts and errors */
static void nbNormalizeHistogram(NBodyHistogram* histogram)
{
    unsigned int i;
    double count;

    unsigned int nBin = histogram->nBin;
    const HistogramParams* hp = &histogram->params;
    double binSize = (double) hp->binSize;
    double totalNum = (double) histogram->totalNum;
    HistData* histData = histogram->data;
    double start = nbHistogramStart(&histogram->params);

    for (i = 0; i < nBin; ++i)
    {
        count = (double) histData[i].rawCount;
        histData[i].lambda = ((double) i + 0.5) * binSize + start;  /* Report center of the bins */
        histData[i].count = count / totalNum;
        histData[i].err = nbNormalizedHistogramError(histData[i].rawCount, totalNum);
    }
}

/* Find the center of mass of the normalized histogram. If
 * useBinIndex, use the position in the histogram rather than the
 * position in lambda  */
static double nbHistogramCenterOfMass(const NBodyHistogram* hist, int useBinIndex)
{
    unsigned int i;
    unsigned int n = hist->nBin;
    const HistData* data = hist->data;
    double cm = 0.0;

    if (useBinIndex)
    {
        for (i = 0; i < n; ++i)
        {
            cm += (double) i * data[i].count;
        }
    }
    else
    {
        for (i = 0; i < n; ++i)
        {
            cm += data[i].lambda * data[i].count;
        }
    }

    /* cm /= (total mass = 1.0) */

    return cm;
}


/*
Takes a treecode position, converts it to (l,b), then to (lambda,
beta), and then constructs a histogram of the density in lambda.

Then calculates the cross correlation between the model histogram and
the data histogram A maximum correlation means the best fit */

/* Bin the bodies from the simulation into maxIdx bins.
   Returns null on failure
 */
NBodyHistogram* nbCreateHistogram(const NBodyCtx* ctx,        /* Simulation context */
                                  const NBodyState* st,       /* Final state of the simulation */
                                  const HistogramParams* hp)  /* Range of histogram to create */
{
    double lambda;
    mwvector lambdaBetaR;
    unsigned int i;
    unsigned int idx;
    unsigned int totalNum = 0;
    Body* p;
    NBodyHistogram* histogram;
    HistData* histData;
    NBHistTrig histTrig;
    const Body* endp = st->bodytab + st->nbody;
    double binSize = (double) hp->binSize;

    /* Calculate the bounds of the bin range, making sure to use a
     * fixed bin size which spans the entire range, and is symmetric
     * around 0 */

    double start = nbHistogramStart(hp);
    unsigned int nBin = nbHistogramNBin(hp);


    nbGetHistTrig(&histTrig, hp);

    histogram = mwCalloc(sizeof(NBodyHistogram) + nBin * sizeof(HistData), sizeof(char));
    histogram->nBin = nBin;
    histogram->hasRawCounts = TRUE;
    histogram->params = *hp;
    histogram->totalSimulated = (unsigned int) st->nbody;
    histogram->massPerParticle = (double) st->bodytab->bodynode.mass;
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

	    // Use vector instead
            lambdaBetaR = nbXYZToLambdaBeta(&histTrig, Pos(p), ctx->sunGCDist);
            lambda = L(lambdaBetaR);

	    idx = (unsigned int) floor((lambda - start) / binSize);
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


/* Read in a histogram from a file for calculating a likelihood value.
 */
NBodyHistogram* nbReadHistogram(const char* histogramFile)
{
    FILE* f;
    int rc = 0;
    size_t fsize = 0;
    NBodyHistogram* histogram = NULL;
    HistData* histData = NULL;
    unsigned int fileCount = 0;
    unsigned int lineNum = 0;
    mwbool error = FALSE;
    mwbool readParams = FALSE; /* Read some other optional histogram params */
    mwbool readNGen = FALSE;  /* Read the scale for the histogram (particles in data bin) */
    mwbool readTotalSim = FALSE; /*Read the total number of particles simulated for the histogram */
    mwbool readMass = FALSE; /*Read the mass per particle for the histogram*/
    unsigned int nGen = 0;    /* Number of particles read from the histogram */
    unsigned int totalSim = 0;	/*Total number of simulated particles read from the histogram */
    double mass = 0;			/*mass per particle read from the histogram */
    char lineBuf[1024];

    f = mwOpenResolved(histogramFile, "r");

    if (f == NULL)
    {
        mw_printf("Error opening histogram file '%s'\n", histogramFile);
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
            double phi, theta, psi;

            rc = sscanf(lineBuf,
                        " phi = %lf , theta = %lf , psi = %lf \n",
                        &phi, &theta, &psi);
            if (rc == 3)
            {
                readParams = TRUE;
                continue;
            }
        }

        if (!readNGen)
        {
            rc = sscanf(lineBuf, " n = %u \n", &nGen);
            if (rc == 1)
            {
                readNGen = TRUE;
                continue;
            }
        }
        if (!readMass)
        {
            rc = sscanf(lineBuf, " massPerParticle = %lf \n", &mass);
            if (rc == 1)
            {
                readMass = TRUE;
                continue;
            }
        }
        if (!readTotalSim)
        {
            rc = sscanf(lineBuf, " totalSimulated = %u \n", &totalSim);
            if (rc == 1)
            {
                readTotalSim = TRUE;
                continue;
            }
        }

        rc = sscanf(lineBuf,
                    "%d %lf %lf %lf \n",
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
    histogram->totalNum = nGen;
    histogram->totalSimulated = totalSim;
    histogram->massPerParticle = mass;
    return histogram;
}
