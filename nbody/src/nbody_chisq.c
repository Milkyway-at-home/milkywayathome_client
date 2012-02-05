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

#include "nbody_priv.h"
#include "nbody_chisq.h"
#include "milkyway_util.h"
#include "nbody_emd.h"


/* From the range of a histogram, find the number of bins */
static unsigned int nbHistogramNBin(const HistogramParams* hp)
{
    double rawCount = (hp->endRaw - hp->startRaw) / hp->binSize;
    return (unsigned int) ceil(rawCount);   /* Total number of bins */
}

/* Find the corrected starting point for the histogram */
static double nbHistogramStart(const HistogramParams* hp)
{
    unsigned int nBin = nbHistogramNBin(hp);
    return mw_ceil(hp->center - hp->binSize * (double) nBin / 2.0);
}

static inline double nbNormalizedHistogramError(unsigned int n, double total)
{
    return (n == 0) ? inv(total) : sqrt((double) n) / total;
}


static double nbCorrectRenormalizedInHistogram(const NBodyHistogram* histogram, const NBodyHistogram* data)
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


static double nbSahaTerm(double m, double s)
{
    /*
      W = product_{i = 1 .. B }

                  (m_i + s_i)!
                 --------------
                   m_i! s_i!

       Using Stirling's approximation with an additional term
       ln(n!) ~= ln(2 *pi * n) / 2 + n * ln(n) - n
       ln(n!) ~= ln(2 pi) / 2 + (n + 1/2) * ln(n) - n

       This reduces to

       ln(W) = sum_{i = 1 .. B}


     */

    /* Stirling's approximation fails miserably at n = 0

       If s_i == 0 or m_i == 0
       the sum term reduces to 0.

       Suppose m_i == 0
       ln(0!) == 0

       ln(0!s_i!) - ln(0!) - ln(s_i) == ln(s_i!) - ln(s_i!) == 0

       If both == 0
       ln(0!0!) - ln(0!) - ln(0!) == 0
     */

    if (fabs(m) < 1.0e-10 || fabs(s) < 1.0e-10)
    {
        return 0.0;
    }
    else
    {
        return -0.5 * log(M_2PI) + (m + s + 0.5) * log(m + s) - (m + 0.5) * log(m) - (s + 0.5) * log(s);
    }
}

static double nbPoissonTerm(double f, double y)
{
    /*
      Fitting a data set y(y1, y2, .. yn)  to a model function f(f1, f2, .. fn)
      sum = i .. N bins,
      2 sum(f_i - y_i) - sum(i != 1, y_i != 0) y_i * ln(f_i / y_i))
     */

    if (fabs(f) < 1.0e-10 || fabs(y) < 1.0e-10)
    {
        return 2.0 * f;
    }
    else
    {
        return 2.0 * ((f - y) - y * log(f / y));
    }
}

static double nbKullbackLeiblerTerm(double h, double k)
{
    /* Symmetrized version. (Jeffrey divergence?) */
    double m;

    if (fabs(h) < 1.0e-10 || fabs(k) < 1.0e-10)
    {
        return 0.0;
    }
    else
    {
        m = (h + k) / 2.0;
        return h * log(h / m) + k * log(k / m);
    }


#if 0
    double p = h;
    double q = k;
    /* Non-symmetrized version */
    if (fabs(p) < 1.0e-10 || fabs(q) < 1.0e-10)
    {
        /* Not sure this is really correct for q == 0 */
        return 0.0;
    }
    else
    {
        return p * (nbLog2(p) / nbLog2(q));
    }
#endif
}

static double nbChisqAlt(double p, double q)
{
    return 0.5 * sqr(p - q) / (p + q);
}

/* Calculate chisq from read data histogarm and the generated histogram */
double nbCalcChisq(const NBodyHistogram* data,        /* Data histogram */
                   const NBodyHistogram* histogram,   /* Generated histogram */
                   NBodyLikelihoodMethod method)
{
    unsigned int i;
    double tmp;
    double effTotalNum;
    double chiSq = 0.0;
    double n;
    double err;
    double simErr;
    double scale = 1.0;
    unsigned int nBin = data->nBin;

    assert(nBin == histogram->nBin);

    if (!histogram->hasRawCounts)
    {
        mw_printf("FIXME: other likelihoods need raw count on generated histogram\n");
        return NAN;
    }

    if (method == NBODY_SAHA)
    {
        /* We need to have the total number to scale to the correct
         * numbers for Saha likelihood */
        scale = (double) data->totalNum;
        if (data->totalNum == 0 || histogram->totalNum == 0)
        {
            mw_printf("Histogram scales required for Saha likelihood but missing\n");
            return NAN;
        }
    }

    if (histogram->totalNum == 0)
    {
        return INFINITY;
    }

    effTotalNum = (double) nbCorrectTotalNumberInHistogram(histogram, data);

    for (i = 0; i < nBin; ++i)
    {
        if (data->data[i].useBin)  /* Skip bins with missing data */
        {
            n = (double) histogram->data[i].rawCount;
            err = data->data[i].err;

            switch (method)
            {
                case NBODY_ORIG_CHISQ:
                    tmp = (data->data[i].count - (n / effTotalNum)) / err;
                    chiSq += sqr(tmp);
                    break;

                case NBODY_ORIG_ALT:
                    /* We already have errors from the simulation, but
                     * we need to correct the errors in case there
                     * were any bins we are skipping for matching to
                     * the data */
                    simErr = nbNormalizedHistogramError(histogram->data[i].rawCount, effTotalNum);

                    /* effective error = sqrt( (data error)^2 + (sim count error)^2 ) */
                    err = sqrt(sqr(err) + sqr(simErr));
                    tmp = (data->data[i].count - (n / effTotalNum)) / err;
                    chiSq += sqr(tmp);
                    break;

                case NBODY_CHISQ_ALT:
                    chiSq += nbChisqAlt(data->data[i].count, n / effTotalNum);
                    break;

                case NBODY_POISSON:
                    /* Poisson one */
                    chiSq += nbPoissonTerm(data->data[i].count, n / effTotalNum);
                    break;

                case NBODY_KOLMOGOROV:
                    chiSq = mw_fmax(chiSq, fabs(data->data[i].count - (n / effTotalNum)));
                    break;

                case NBODY_KULLBACK_LEIBLER:
                    /* "Relative entropy" */
                    chiSq += nbKullbackLeiblerTerm(data->data[i].count, n / effTotalNum);
                    break;

                case NBODY_SAHA:
                    /* This will actually find ln(W). W is an unreasonably large number. */
                    chiSq += nbSahaTerm(n, scale * data->data[i].count);
                    break;

                case NBODY_INVALID_METHOD:
                case NBODY_EMD:
                default:
                    mw_fail("Invalid likelihood method\n");
            }
        }
    }

    return chiSq;
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
    double totalNum = (double) histogram->totalNum;
    HistData* histData = histogram->data;
    double start = nbHistogramStart(&histogram->params);


    for (i = 0; i < nBin; ++i)
    {
        count = (double) histData[i].rawCount;
        histData[i].lambda = ((double) i  + 0.5) * hp->binSize + start;  /* Report center of the bins */
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
    unsigned int i;
    unsigned int idx;
    unsigned int totalNum = 0;
    Body* p;
    NBodyHistogram* histogram;
    HistData* histData;
    NBHistTrig histTrig;
    const Body* endp = st->bodytab + st->nbody;

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
    unsigned int nGen = 0;    /* Number of particles read from the */
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
    return histogram;
}

static double nbWorstCaseEMD(const NBodyHistogram* hist)
{
    return fabs(hist->data[0].lambda - hist->data[hist->nBin - 1].lambda);
}

double nbMatchEMD(const NBodyHistogram* data, const NBodyHistogram* histogram)
{
    unsigned int i;
    unsigned int n = data->nBin;
    unsigned int effTotalNum = 0;
    double emd;
    WeightPos* hist;
    WeightPos* dat;
    double renormalize = 0.0;

    if (data->nBin != histogram->nBin)
    {
        /* FIXME?: We could have mismatched histogram sizes, but I'm not sure what to do with ignored bins and renormalization */
        return NAN;
    }

    if (histogram->hasRawCounts)
    {
        effTotalNum = nbCorrectTotalNumberInHistogram(histogram, data);
        mw_printf("Histogram has %u bodies with %u in accepted bins\n",
                  histogram->totalNum,
                  effTotalNum
            );

        if (effTotalNum == 0)
        {
            /* If the histogram is totally empty, it is worse than the worst case */
            return INFINITY;
        }
    }
    else
    {
        renormalize = nbCorrectRenormalizedInHistogram(histogram, data);
    }

    /* We need a correctly normalized histogram with the missing bins filtered out */
    hist = mwCalloc(n, sizeof(WeightPos));
    dat = mwCalloc(n, sizeof(WeightPos));

    for (i = 0; i < n; ++i)
    {
        if (data->data[i].useBin)
        {
            dat[i].weight = data->data[i].count;
            if (histogram->hasRawCounts)
            {
                double correctedCount = (double) histogram->data[i].rawCount / (double) effTotalNum;
                hist[i].weight = correctedCount;
            }
            else
            {
                hist[i].weight = histogram->data[i].count / renormalize;
            }
        }
        /* Otherwise weight is 0.0 */

        hist[i].pos = histogram->data[i].lambda;
        dat[i].pos = data->data[i].lambda;
    }

    emd = emdCalc((const float*) dat, (const float*) hist, n, n, NULL);

    free(hist);
    free(dat);

    return emd;
}

/*
  Load information necessary to calculate the likelihood from input lua script.

  Load parameters where in the sky the histogram goes and likelihood type.

  Return TRUE on failure.
*/
int nbGetLikelihoodInfo(const NBodyFlags* nbf, HistogramParams* hp, NBodyLikelihoodMethod* method)
{
    lua_State* luaSt = NULL;

    luaSt = nbOpenLuaStateWithScript(nbf);
    if (!luaSt)
    {
        return TRUE;
    }

    if (nbEvaluateHistogramParams(luaSt, hp))
    {
        lua_close(luaSt);
        return TRUE;
    }

    *method = nbEvaluateLikelihoodMethod(luaSt);

    lua_close(luaSt);
    return FALSE;
}

/* Calculate the likelihood from the final state of the simulation */
double nbSystemChisq(const NBodyState* st,
                     const NBodyHistogram* data,
                     const NBodyHistogram* histogram,
                     NBodyLikelihoodMethod method)
{
    if (data->nBin != histogram->nBin)
    {
        mw_printf("Number of bins does not match those in histogram file. "
                  "Expected %u, got %u\n",
                  histogram->nBin,
                  data->nBin);
        return NAN;
    }

    if (method == NBODY_EMD)
    {
        /* We could have far crazier variations in the distance in cases
         * where the number of particles resulting in the bins is very
         * small, such as when a few particles on the edge are thrown out
         * and happen to end up in the binning range.
         *
         * Make sure that at least 1% of the total particles are being
         * counted to hopefully smooth away potential issues.
         *
         * If there are truly no particles in useful bins, the EMD will
         * return infinity. Having more than 0 particles should be better
         * than infinity, so use something a bit worse than the case where
         * 100% is located in opposite bins.
         */
        if (histogram->totalNum < 0.01 * (double) st->nbody)
        {
            double worstEMD;

            mw_printf("Number of particles in bins is very small compared to total. "
                      "(%u << %u). Skipping distance calculation\n",
                      histogram->totalNum,
                      st->nbody
                );
            worstEMD = nbWorstCaseEMD(histogram);
            return 2.0 * worstEMD;
        }

        return nbMatchEMD(data, histogram);
    }
    else
    {
        return nbCalcChisq(data, histogram, method);
    }
}

double nbMatchHistogramFiles(const char* datHist, const char* matchHist)
{
    NBodyHistogram* dat;
    NBodyHistogram* match;
    double emd = NAN;

    dat = nbReadHistogram(datHist);
    match = nbReadHistogram(matchHist);

    if (dat && match)
    {
        emd = nbMatchEMD(dat, match);
    }

    free(dat);
    free(match);

    return emd;
}

