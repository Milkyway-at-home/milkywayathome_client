/*
 *  Copyright (c) 2010-2011 Rensselaer Polytechnic Institute
 *  Copyright (c) 2010-2011 Ben Willett
 *  Copyright (c) 2010-2011 Matthew Arsenault
 *  Copyright (c) 2016-2018 Siddhartha Shelton
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
#include "nbody_potential_types.h"
#include "milkyway_util.h"
#include "nbody_coordinates.h"
#include "nbody_show.h"

/*Calculates the center of two numbers */
static real nbHistogramCenter(real start, real end)
{
    return (start + end)/2;
}

/* From the range of a histogram, find the bin size in Lambda */
static real nbHistogramLambdaBinSize(const HistogramParams* hp)
{
    real binSize = (hp->lambdaEnd - hp->lambdaStart) / (real) hp->lambdaBins;
    return binSize;   /* Size of bins */
}

/* From the range of a histogram, find the bin size in Beta */
static real nbHistogramBetaBinSize(const HistogramParams* hp)
{
    real binSize = (hp->betaEnd - hp->betaStart) / (real) hp->betaBins;
    return binSize;
}

real nbNormalizedHistogramError(unsigned int n, real total)
{
    return (n == 0) ? inv(total) : mw_sqrt((real) n) / total;
}


real nbCorrectRenormalizedInHistogram(const NBodyHistogram* histogram, const NBodyHistogram* data)
{
    unsigned int i;
    unsigned int nBin = data->lambdaBins * data->betaBins;
    real total = 0.0;

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
    unsigned int nBin = data->lambdaBins * data->betaBins;
    unsigned int totalNum = histogram->totalNum;

    assert(histogram->hasRawCounts);
    assert(histogram->lambdaBins == data->lambdaBins);
    assert(histogram->betaBins == data->betaBins);

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
                                   int nbody,
                                   real bestLikelihood_time,
                                   real bestLikelihood)
{
    char tBuf[256];
    const Potential* p = &ctx->pot;

    if (bestLikelihood_time == 0.0)
    {
        bestLikelihood_time = ctx->timeEvolve;
    }

    mwLocalTimeFull(tBuf, sizeof(tBuf));

    fprintf(f,
            "#\n"
            "# Generated %s\n"
            "#\n"
            "# (phi,   theta,  psi) = (%f, %f, %f)\n"
            "# (lambdaStart, lambdaCenter, lambdaEnd) = (%f, %f, %f)\n"
            "# (betaStart, betaCenter, betaEnd) = (%f, %f, %f)\n"
            "# Lambda Bin size = %f\n"
            "# Beta Bin size = %f\n"
            "#\n"
            "#\n",
            tBuf,
            hp->phi, hp->theta, hp->psi,
            hp->lambdaStart, nbHistogramCenter(hp->lambdaStart, hp->lambdaEnd), hp->lambdaEnd,
            hp->betaStart, nbHistogramCenter(hp->betaStart, hp->betaEnd), hp->betaEnd,
            nbHistogramLambdaBinSize(hp),
            nbHistogramBetaBinSize(hp));
    
    fprintf(f,
            "# Nbody = %d\n"
            "# Evolve backward time = %f\n"
            "# Evolve forward time = %f\n"
            "# Best Likelihood = %f\n"
            "# Timestep = %f\n"
            "# Sun GC Dist = %f\n"
            "# Criterion = %s\n"
            "# Theta = %f\n"
            "# Quadrupole Moments = %s\n"
            "# Eps = %f\n"
            "#\n",
            nbody,
            ctx->timeEvolve,
            bestLikelihood_time,
            -bestLikelihood,
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

    switch (p->sphere[1].type)
    {
        case HernquistSpherical:

            fprintf(f,
                    "# Spherical: Hernquist\n"
                    "#   mass = %f\n"
                    "#   a = %f\n"
                    "#\n",
                    p->sphere[0].mass,
                    p->sphere[0].scale);
            break;

        case PlummerSpherical:

            fprintf(f,
                    "# Spherical: Plummer\n"
                    "#   mass = %f\n"
                    "#   a = %f\n"
                    "#\n",
                    p->sphere[0].mass,
                    p->sphere[0].scale);
            break;

        case NoSpherical:

            fprintf(f,
                    "# Spherical: None\n");
            break;

        case InvalidSpherical:
        default:
            fprintf(f,
                    "# Spherical: ???\n"
                    "#\n");
    }

    switch (p->disk.type)
    {
        case MiyamotoNagaiDisk:

            fprintf(f,
                    "# Primary Disk: MiaymotoNagai\n"
                    "#   mass = %f\n"
                    "#   a = %f\n"
                    "#   b = %f\n"
                    "#\n",
                    p->disk.mass,
                    p->disk.scaleLength,
                    p->disk.scaleHeight);
            break;

        case DoubleExponentialDisk:

            fprintf(f,
                    "# Primary Disk: DoubleExponential\n"
                    "#   mass = %f\n"
                    "#   Rd = %f\n"
                    "#   zd = %f\n"
                    "#\n",
                    p->disk.mass,
                    p->disk.scaleLength,
                    p->disk.scaleHeight);
            break;

        case Sech2ExponentialDisk:

            fprintf(f,
                    "# Primary Disk: Sech2Exponential\n"
                    "#   mass = %f\n"
                    "#   Rd = %f\n"
                    "#   zd = %f\n"
                    "#\n",
                    p->disk.mass,
                    p->disk.scaleLength,
                    p->disk.scaleHeight);
            break;

        case FreemanDisk:
            fprintf(f,
                    "# Primary Disk: Freeman\n"
                    "#   mass = %f\n"
                    "#   b = %f\n"
                    "#\n",
                    p->disk.mass,
                    p->disk.scaleLength);
            break;

        case NoDisk:
            fprintf(f,
                    "# Primary Disk: None\n");
            break;

        case InvalidDisk:
        default:
            fprintf(f,
                    "# Primary Disk: ???\n"
                    "#\n");
    }

    switch (p->disk2.type)
    {
        case MiyamotoNagaiDisk:

            fprintf(f,
                    "# Seconday Disk: MiaymotoNagai\n"
                    "#   mass = %f\n"
                    "#   a = %f\n"
                    "#   b = %f\n"
                    "#\n",
                    p->disk.mass,
                    p->disk.scaleLength,
                    p->disk.scaleHeight);
            break;

        case DoubleExponentialDisk:

            fprintf(f,
                    "# Secondary Disk: DoubleExponential\n"
                    "#   mass = %f\n"
                    "#   Rd = %f\n"
                    "#   zd = %f\n"
                    "#\n",
                    p->disk.mass,
                    p->disk.scaleLength,
                    p->disk.scaleHeight);
            break;

        case Sech2ExponentialDisk:

            fprintf(f,
                    "# Secondary Disk: Sech2Exponential\n"
                    "#   mass = %f\n"
                    "#   Rd = %f\n"
                    "#   zd = %f\n"
                    "#\n",
                    p->disk.mass,
                    p->disk.scaleLength,
                    p->disk.scaleHeight);
            break;

        case FreemanDisk:
            fprintf(f,
                    "# Secondary Disk: Freeman\n"
                    "#   mass = %f\n"
                    "#   b = %f\n"
                    "#\n",
                    p->disk.mass,
                    p->disk.scaleLength);
            break;

        case OrbitingPointMassBar:
            fprintf(f,
                    "# Secondary Disk: OrbitingPointMassBar\n"
                    "#   mass = %f\n"
                    "#   b = %f\n"
                    "#   pattern speed = %f\n"
                    "#   start angle = %f\n"
                    "#\n",
                    p->disk2.mass,
                    p->disk2.scaleLength,
                    p->disk2.patternSpeed,
                    p->disk2.startAngle);
            break;

        case NoDisk:
            fprintf(f,
                    "# Secondary Disk: None\n"
                    "#\n");
            break;

        case InvalidDisk:
        default:
            fprintf(f,
                    "# Secondary Disk: ???\n"
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
        case CausticHalo:
            fprintf(f,
                    "# Halo: Caustic\n"
                    "#\n");
            break;

        case AllenSantillanHalo:
            fprintf(f,
                    "# Halo: Allen-Santillan\n"
                    "#   mass = %f\n"
                    "#   a = %f\n"
                    "#   gamma = %f\n"
                    "#   lambda = %f\n"
                    "#\n",
                    p->halo.mass,
                    p->halo.scaleLength,
                    p->halo.gamma,
                    p->halo.lambda);
            break;

        case WilkinsonEvansHalo:
            fprintf(f,
                    "# Halo: Wilkinson-Evans\n"
                    "#   mass = %f\n"
                    "#   a = %f\n"
                    "#\n",
                    p->halo.mass,
                    p->halo.scaleLength);
            break;

        case NFWMassHalo:
            fprintf(f,
                    "# Halo: NFW\n"
                    "#   mass = %f\n"
                    "#   a = %f\n"
                    "#\n",
                    p->halo.mass,
                    p->halo.scaleLength);
            break;

        case PlummerHalo:
            fprintf(f,
                    "# Halo: Plummer\n"
                    "#   mass = %f\n"
                    "#   a = %f\n"
                    "#\n",
                    p->halo.mass,
                    p->halo.scaleLength);
            break;

        case HernquistHalo:
            fprintf(f,
                    "# Halo: Hernquist\n"
                    "#   mass = %f\n"
                    "#   a = %f\n"
                    "#\n",
                    p->halo.mass,
                    p->halo.scaleLength);
            break;

        case NinkovicHalo:
            fprintf(f,
                    "# Halo: Ninkovic\n"
                    "#   rho0 = %f\n"
                    "#   a = %f\n"
                    "#   rl = %f\n"
                    "#\n",
                    p->halo.rho0,
                    p->halo.scaleLength,
                    p->halo.lambda);
            break;

        case NoHalo:
            fprintf(f,
                    "# Halo: None\n");
            break;

        case InvalidHalo:
        default:
            fprintf(f,
                    "# Halo: ???\n"
                    "#\n");
    }

    fprintf(f,
            "#\n"
            "#Column Headers:\n"
            "# UseBin,  Lambda,  Beta,  Normalized Counts, Count Error, "
            "Beta Dispersion,  Beta Dispersion Error,"
            "LOS Velocity Dispersion, Velocity Dispersion Error\n"
            "#\n"
            "\n"
        );
}

/* Print the histogram without a header. */
void nbPrintHistogram(FILE* f, const NBodyHistogram* histogram)
{
    unsigned int i;
    const HistData* data;
    unsigned int nBin = histogram->lambdaBins * histogram->betaBins;

    mw_boinc_print(f, "<histogram>\n");
    fprintf(f, "n = %u\n", histogram->totalNum);
    fprintf(f, "massPerParticle = %12.15f\n", histogram->massPerParticle);
    fprintf(f, "totalSimulated = %u\n", histogram->totalSimulated);
    fprintf(f, "lambdaBins = %u\n", histogram->lambdaBins);
    fprintf(f, "betaBins = %u\n", histogram->betaBins);

    
    for (i = 0; i < nBin; ++i)
    {
        data = &histogram->data[i];
        fprintf(f,
                "%d %12.15f %12.15f %12.15f %12.15f %12.15f %12.15f %12.15f %12.15f\n",
                data->useBin,
                data->lambda,
                data->beta,
                data->count,
                data->err,
                data->beta_disp,
                data->beta_disperr,
                data->vdisp,
                data->vdisperr);

    /* Print blank lines for plotting histograms in gnuplot pm3d */
        if(i % histogram->betaBins == histogram->betaBins-1)
        {
            fprintf(f, "\n");
        }

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

    nbPrintHistogramHeader(f, ctx, &histogram->params, st->nbody, st->bestLikelihood_time, st->bestLikelihood);
    nbPrintHistogram(f, histogram);

    if (f != DEFAULT_OUTPUT_FILE)
        fclose(f);
}

/* Get normalized histogram counts and errors */
static void nbNormalizeHistogram(NBodyHistogram* histogram)
{
    unsigned int i;
    unsigned int j;
    unsigned int Histindex;
    real count;

    const HistogramParams* hp = &histogram->params;

    unsigned int lambdaBins = histogram->lambdaBins;
    unsigned int betaBins = histogram->betaBins;
    real lambdaSize = nbHistogramLambdaBinSize(hp);
    real betaSize = nbHistogramBetaBinSize(hp);
    real lambdaStart = hp->lambdaStart;
    real betaStart = hp->betaStart;

    real totalNum = (real) histogram->totalNum;
    HistData* histData = histogram->data;


    for (i = 0; i < lambdaBins; ++i)
    {
        for(j = 0; j < betaBins; ++j)
        {
            Histindex = i * betaBins + j;
            count = (real) histData[Histindex].rawCount;
            
            /* Report center of the bins */
            histData[Histindex].lambda = ((real) i + 0.5) * lambdaSize + lambdaStart;
            histData[Histindex].beta   = ((real) j + 0.5) * betaSize + betaStart;
            histData[Histindex].count  = count / totalNum;
            histData[Histindex].err    = nbNormalizedHistogramError(histData[i].rawCount, totalNum);
        }
    }
}


/*
Takes a treecode position, converts it to (l,b), then to (lambda,
beta), and then constructs a histogram of the density in lambda and beta.

Then calculates the cross correlation between the model histogram and
the data histogram A maximum correlation means the best fit */

/* Returns null on failure */
NBodyHistogram* nbCreateHistogram(const NBodyCtx* ctx,        /* Simulation context */
                                  const NBodyState* st,       /* Final state of the simulation */
                                  const HistogramParams* hp)  /* Range of histogram to create */
{
    real lambda;
    real beta;
    mwvector lambdaBetaR;
    unsigned int lambdaIndex;
    unsigned int betaIndex;
    unsigned int Histindex;
    unsigned int totalNum = 0;
    Body* p;
    NBodyHistogram* histogram;
    HistData* histData;
    NBHistTrig histTrig;
    const Body* endp = st->bodytab + st->nbody;
    real lambdaSize = nbHistogramLambdaBinSize(hp);
    real betaSize = nbHistogramBetaBinSize(hp);
    /* Calculate the bounds of the bin range, making sure to use a
     * fixed bin size which spans the entire range, and is symmetric
     * around 0 */

    real lambdaStart = hp->lambdaStart;
    real betaStart = hp->betaStart;
    unsigned int lambdaBins = hp->lambdaBins;
    unsigned int betaBins = hp->betaBins;
    unsigned int IterMax = ctx->IterMax;
    /*unsigned int IterMax = 6;*/	/*Default value for IterMax*/
    unsigned int nBin = lambdaBins * betaBins;
    unsigned int body_count = 0;
    unsigned int ub_counter = 0;
    
    real Nbodies = st->nbody;
    mwbool islight = FALSE;//is it light matter?
    
    
    nbGetHistTrig(&histTrig, hp);
    histogram = mwCalloc(sizeof(NBodyHistogram) + nBin * sizeof(HistData), sizeof(char));
    histogram->lambdaBins = lambdaBins;
    histogram->betaBins = betaBins;
    histogram->hasRawCounts = TRUE;
    histogram->params = *hp;
    
    real v_line_of_sight;
    
    for (int i = 0; i < Nbodies; i++)
    {
        const Body* b = &st->bodytab[i];
        if(Type(b) == BODY(islight))
        {
            histogram->massPerParticle = Mass(b);
            body_count++;
        }
    }

    real * use_velbody  = mwCalloc(body_count, sizeof(real));
    real * use_betabody  = mwCalloc(body_count, sizeof(real));
    real * vlos      = mwCalloc(body_count, sizeof(real));
    real * betas     = mwCalloc(body_count, sizeof(real));
    
    histogram->totalSimulated = (unsigned int) body_count;
    histData = histogram->data;
    
    /* It does not make sense to ignore bins in a generated histogram */
    for (Histindex = 0; Histindex < nBin; ++Histindex)
    {
        histData[Histindex].rawCount = 0;
        histData[Histindex].v_sum    = 0.0;
        histData[Histindex].vsq_sum  = 0.0;
        histData[Histindex].vdisp    = 0.0;
        histData[Histindex].vdisperr = 0.0;
        
        histData[Histindex].beta_sum    = 0.0;
        histData[Histindex].betasq_sum  = 0.0;
        histData[Histindex].beta_disp    = 0.0;
        histData[Histindex].beta_disperr = 0.0;
        
        histData[Histindex].outliersBetaRemoved = 0.0;
        histData[Histindex].outliersVelRemoved = 0.0;
        histData[Histindex].useBin = TRUE;
    }


    for (p = st->bodytab; p < endp; ++p)
    {
        /* Only include bodies in models we aren't ignoring (like dark matter) */
        if (!ignoreBody(p))
        {
            
            /* Get the position in lbr coorinates */
            lambdaBetaR = nbXYZToLambdaBeta(&histTrig, Pos(p), ctx->sunGCDist);
            lambda = L(lambdaBetaR);
            beta = B(lambdaBetaR);
            
            use_betabody[ub_counter] = DEFAULT_NOT_USE;//defaulted to not use body
            use_velbody[ub_counter] = DEFAULT_NOT_USE;//defaulted to not use body
            
            vlos[ub_counter]     = DEFAULT_NOT_USE;//default vlos
            betas[ub_counter]    = DEFAULT_NOT_USE;

            /* Find the indices */
            lambdaIndex = (unsigned int) mw_floor((lambda - lambdaStart) / lambdaSize);
            betaIndex = (unsigned int) mw_floor((beta - betaStart) / betaSize);

            /* Check if the position is within the bounds of the histogram */
            if (lambdaIndex < lambdaBins && betaIndex < betaBins)   
            {   
                Histindex = lambdaIndex * betaBins + betaIndex;
                use_betabody[ub_counter] = Histindex;//if body is in hist, mark which hist bin
                use_velbody[ub_counter] = Histindex;//if body is in hist, mark which hist bin
                
                histData[Histindex].rawCount++;
                ++totalNum;
                
                
                v_line_of_sight = calc_vLOS(Vel(p), Pos(p), ctx->sunGCDist);//calc the heliocentric line of sight vel
                vlos[ub_counter] = v_line_of_sight;//store the vlos's so as to not have to recalc
                betas[ub_counter] = beta;
                /* each of these are components of the vel disp */
                histData[Histindex].v_sum += v_line_of_sight;
                histData[Histindex].vsq_sum += sqr(v_line_of_sight);
                
                /* each of these are components of the beta disp */
                histData[Histindex].beta_sum += beta;
                histData[Histindex].betasq_sum += sqr(beta);
            }
            ub_counter++;
        }
    }
    histogram->totalNum = totalNum; /* Total particles in range */

    nbCalcVelDisp(histogram, TRUE, ctx->VelCorrect);
    nbCalcBetaDisp(histogram, TRUE, ctx->BetaCorrect);
    /* this converges somewhere between 3 and 6 iterations */
    for(int i = 0; i < IterMax; i++)
    {
        nbRemoveBetaOutliers(st, histogram, use_betabody, betas, ctx->BetaSigma);
        nbCalcBetaDisp(histogram, FALSE, ctx->BetaCorrect);
        
        nbRemoveVelOutliers(st, histogram, use_velbody, vlos, ctx->VelSigma);
        nbCalcVelDisp(histogram, FALSE, ctx->VelCorrect);
        
    }
    
    nbNormalizeHistogram(histogram);
    
    free(use_velbody);
    free(use_betabody);
    free(vlos);
    free(betas);
    
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
    mwbool readLambdaBins = FALSE; /* Read the number of bins in the lambda direction */
    mwbool readBetaBins = FALSE; /* Read the number of bins the beta direction */
    mwbool readOpeningTag = FALSE; /* Read the <histogram> tag */
    mwbool readClosingTag = FALSE; /* Read the </histogram> tag */
    unsigned int nGen = 0;    /* Number of particles read from the histogram */
    unsigned int totalSim = 0;  /*Total number of simulated particles read from the histogram */
    unsigned int lambdaBins = 0; /* Number of bins in lambda direction */
    unsigned int betaBins = 0; /* Number of bins in beta direction */
    real mass = 0;            /*mass per particle read from the histogram */
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

        /* Skip <histogram> tags */
        if(!readOpeningTag)
        {
            if(strcmp("<histogram>\n", lineBuf) == 0)
            {
                readOpeningTag = TRUE;
                continue;
            }
        }

        if(!readClosingTag)
        {
            if(strcmp("</histogram>\n", lineBuf) == 0)
            {
                readClosingTag = TRUE;
                continue;
            }
        }

        if (!readParams)  /* One line is allowed for information on the histogram */
        {
            real phi, theta, psi;

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
        
        if (!readLambdaBins)
        {
            rc = sscanf(lineBuf, " lambdaBins = %u \n", &lambdaBins);
            if(rc == 1)
            {
                readLambdaBins = TRUE;
                continue;
            }
        }

        if (!readBetaBins)
        {
            rc = sscanf(lineBuf, " betaBins = %u \n", &betaBins);
            if(rc == 1)
            {
                readBetaBins = TRUE;
                continue;
            }
        }

        rc = sscanf(lineBuf,
                    "%d %lf %lf %lf %lf %lf %lf %lf %lf\n",
                    &histData[fileCount].useBin,
                    &histData[fileCount].lambda,
                    &histData[fileCount].beta,
                    &histData[fileCount].count,
                    &histData[fileCount].err,
                    &histData[fileCount].beta_disp,
                    &histData[fileCount].beta_disperr,
                    &histData[fileCount].vdisp,
                    &histData[fileCount].vdisperr);
        
        
        /* new standard for histograms is being enforced. Two extra columns for vel and beta dispersion 
         * and their errors. If not using them, can input zeros in the Columns
         */
        if (rc != 9)
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

    histogram->lambdaBins = lambdaBins;
    histogram->betaBins = betaBins;
    histogram->totalNum = nGen;
    histogram->totalSimulated = totalSim;
    histogram->massPerParticle = mass;
    return histogram;
}
