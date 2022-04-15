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
static real_0 nbHistogramCenter(real_0 start, real_0 end)
{
    return (start + end)/2;
}

/* From the range of a histogram, find the bin size in Lambda */
static real_0 nbHistogramLambdaBinSize(const HistogramParams* hp)
{
    real_0 binSize = (hp->lambdaEnd - hp->lambdaStart) / (real_0) hp->lambdaBins;
    return binSize;   /* Size of bins */
}

/* From the range of a histogram, find the bin size in Beta */
static real_0 nbHistogramBetaBinSize(const HistogramParams* hp)
{
    real_0 binSize = (hp->betaEnd - hp->betaStart) / (real_0) hp->betaBins;
    return binSize;
}

real nbNormalizedHistogramError(real* n, real* total)
{
    real tmp1, tmp2;
    real norm_n = mw_div(n, total);
    tmp1 = mw_mul_s(&norm_n,-2.0);
    tmp1 = mw_add_s(&tmp1, 1.0);
    tmp1 = mw_mul(&tmp1, n);
    tmp2 = sqr(&norm_n);
    tmp2 = mw_mul(&tmp2, total);
    tmp1 = mw_add(&tmp1, &tmp2);
    tmp1 = mw_sqrt(&tmp1);
    return (showRealValue(n) == 0) ? inv(total) : mw_div(&tmp1, total);
}

real nbCorrectRenormalizedInHistogram(const NBodyHistogram* histogram, const NBodyHistogram* data)
{
    unsigned int i;
    unsigned int nBin = data->lambdaBins * data->betaBins;
    real total = ZERO_REAL;

    for (i = 0; i < nBin; ++i)
    {
        if (data->data[i].useBin)
        {
            total = mw_add(&total, &histogram->data[i].variable);
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
real nbCorrectTotalNumberInHistogram(const NBodyHistogram* histogram, /* Generated histogram */
                                                    const NBodyHistogram* data)      /* Data histogram */
{
    unsigned int i;
    unsigned int nBin = data->lambdaBins * data->betaBins;
    real totalNum = histogram->totalNum;

    assert(histogram->hasRawCounts);
    assert(histogram->lambdaBins == data->lambdaBins);
    assert(histogram->betaBins == data->betaBins);

    for (i = 0; i < nBin; ++i)
    {
        if (!data->data[i].useBin)
        {
            totalNum = mw_sub(&totalNum, &histogram->data[i].rawCount);
        }
    }

    return totalNum;
}

static void nbPrintHistogramHeader(FILE* f,
                                   const NBodyCtx* ctx,
                                   const HistogramParams* hp,
                                   NBodyState* st)
{
    int nbody = st->nbody;
    real_0 bestLikelihood_time = st->bestLikelihood_time;
    real bestLikelihood = st->bestLikelihood;
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
    
    //calculate how far off the bar time was for different situations
    
    //bar on calibration on
    real_0 barTimeError = bestLikelihood_time - st->previousForwardTime;
    if(ctx->pot.disk2.type == _NO_DISK){//no bar
        barTimeError = 0;
    }else if(ctx->calibrationRuns == 0){//no calibration but bar on
        barTimeError = bestLikelihood_time - ctx->timeEvolve;
    }
    real_0 barAngleError = barTimeError * ctx->pot.disk2.patternSpeed;
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
            "# Bar Time Error = %f\n"
            "# Bar Angle Error = %f rad\n"
            "#\n",
            nbody,
            ctx->timeBack,
            bestLikelihood_time,
            showRealValue(&bestLikelihood),
            ctx->timestep,
            ctx->sunGCDist,
            showCriterionT(ctx->criterion),
            ctx->theta,
            showBool(ctx->useQuad),
            mw_sqrt_0(ctx->eps2),
            barTimeError,
            barAngleError

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

        case OrbitingBar:
            fprintf(f,
                    "# Secondary Disk: OrbitingBar\n"
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
            "# UseBin,  Lambda,  Beta, Normalized Counts, Count Error, "
            "Beta Dispersion,  Beta Dispersion Error, "
            "LOS Velocity Dispersion, Velocity Dispersion Error, "
            "LOS Velocity, LOS Velocity Error, Beta Average, Beta Average Error, "
            "Distance, Distance Error\n"
            "#\n"
            "\n"
        );
}

/* Print the histogram without a header. */
void nbPrintHistogram(FILE* f, const MainStruct* all)
{
    real_0 output[12]; // for outputting the data

    unsigned int nBin;
    nBin = all->histograms[0]->lambdaBins * all->histograms[0]->betaBins;
    real_0 missingFrac = mw_abs_0(showRealValue(&all->histograms[0]->totalNum) - mw_round_0(showRealValue(&all->histograms[0]->totalNum)));

    // outputs these numbers from one of the histograms
    // at this point, these numbers should be the same for all histograms anyway
    mw_boinc_print(f, "<histogram>\n");
    if(missingFrac < REAL_EPSILON)
    {
        fprintf(f, "n = %u\n", (int) showRealValue(&all->histograms[0]->totalNum));
    }
    else
    {
        fprintf(f, "n = %12.15f\n", showRealValue(&all->histograms[0]->totalNum));
    }
    fprintf(f, "massPerParticle = %12.15f\n", showRealValue(&all->histograms[0]->massPerParticle));
    fprintf(f, "totalSimulated = %u\n", all->histograms[0]->totalSimulated);
    fprintf(f, "lambdaBins = %u\n", all->histograms[0]->lambdaBins);
    fprintf(f, "betaBins = %u\n", all->histograms[0]->betaBins);
    if(all->usage[3]) fprintf(f, "!\n");
    
    
    for (unsigned int i = 0; i < nBin; ++i)
    {
        // should only output if the histogram is being used
        const HistData storedData = all->histograms[0]->data[i];   // for generic output parameters
        int k = 0;
        for(unsigned int j = 0; j < 6; j++)
        {
            if(all->usage[j])
            {
                output[k]   = showRealValue(&all->histograms[j]->data[i].variable);
                output[k+1] = showRealValue(&all->histograms[j]->data[i].err);
            }
            k+=2;
        }

        if(all->usage[3]) // if this histogram has been used, all histograms were calculated
        {
            fprintf(f,
                    "%d %12.15f %12.15f %12.15f %12.15f %12.15f %12.15f %12.15f %12.15f %12.15f %12.15f %12.15f %12.15f %12.15f %12.15f\n",
                    storedData.useBin,
                    showRealValue(&storedData.lambda),
                    showRealValue(&storedData.beta),
                    output[0],      // normalized counts
                    output[1],      // normalized counts error
                    output[2],      // beta disp
                    output[3],      // beta disp error
                    output[4],      // vlos disp
                    output[5],      // vlos disp error
                    output[6],      // vlos average
                    output[7],      // vlos avg error
                    output[8],      // beta average
                    output[9],      // beta average error
                    output[10],     // distance
                    output[11]);    // distance error
        }
        else // only fitting original histograms
        {
            fprintf(f,
                    "%d %12.15f %12.15f %12.15f %12.15f %12.15f %12.15f %12.15f %12.15f\n",
                    storedData.useBin,
                    showRealValue(&storedData.lambda),
                    showRealValue(&storedData.beta),
                    output[0],      // normalized counts
                    output[1],      // normalized counts error
                    output[2],      // beta disp
                    output[3],      // beta disp error
                    output[4],      // vlos disp
                    output[5]);     // vlos disp error  
        }  

    /* Print blank lines for plotting histograms in gnuplot pm3d */
        if(i % all->histograms[0]->betaBins == (all->histograms[0]->betaBins)-1)
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
                      const MainStruct* all)
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

    nbPrintHistogramHeader(f, ctx, &all->histograms[0]->params, st);
    nbPrintHistogram(f, all);

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
    real_0 lambdaSize = nbHistogramLambdaBinSize(hp);
    real_0 betaSize = nbHistogramBetaBinSize(hp);
    real_0 lambdaStart = hp->lambdaStart;
    real_0 betaStart = hp->betaStart;

    real totalNum = histogram->totalNum;
    HistData* histData = histogram->data;


    for (i = 0; i < lambdaBins; ++i)
    {
        for(j = 0; j < betaBins; ++j)
        {
            Histindex = i * betaBins + j;
            count = (real) histData[Histindex].rawCount;
            
            /* Report center of the bins */
            histData[Histindex].lambda = mw_real_const(((real_0) i + 0.5) * lambdaSize + lambdaStart);
            histData[Histindex].beta   = mw_real_const(((real_0) j + 0.5) * betaSize + betaStart);
            histData[Histindex].variable  = mw_div(&count, &totalNum);
            histData[Histindex].err    = nbNormalizedHistogramError(&histData[i].rawCount, &totalNum);
        }
    }
}


/*
Takes a treecode position, converts it to (l,b), then to (lambda,
beta), and then constructs a histogram of the density in lambda and beta.

Then calculates the cross correlation between the model histogram and
the data histogram A maximum correlation means the best fit */

/* Returns null on failure */
MainStruct* nbCreateHistogram(const NBodyCtx* ctx,        /* Simulation context */
                                  const NBodyState* st,       /* Final state of the simulation */
                                  const HistogramParams* hp)  /* Range of histogram to create */
{
    real tmp;
    real location;
    real lambda;
    real beta;
    real v_line_of_sight;
    real bodyBinFrac;
    real fullBodyFrac;
    mwvector lambdaBetaR;
    unsigned int lambdaIndex;
    unsigned int betaIndex;
    unsigned int Histindex;
    unsigned int i, j;
    real totalNum = ZERO_REAL;
    HistData* histData;
    Body* p;
    NBHistTrig histTrig;
    const Body* endp = st->bodytab + st->nbody;
    real_0 lambdaSize = nbHistogramLambdaBinSize(hp);
    real_0 betaSize = nbHistogramBetaBinSize(hp);
    /* Calculate the bounds of the bin range, making sure to use a
     * fixed bin size which spans the entire range, and is symmetric
     * around 0 */

    real_0 lambdaStart = hp->lambdaStart;
    real_0 betaStart = hp->betaStart;
    unsigned int lambdaBins = hp->lambdaBins;
    unsigned int betaBins = hp->betaBins;
    unsigned int IterMax = ctx->IterMax;
    /*unsigned int IterMax = 6;*/   /*Default value for IterMax*/
    unsigned int nBin = lambdaBins * betaBins;
    unsigned int body_count = 0;
    
    MainStruct* all = mwCalloc(6*(sizeof(NBodyHistogram) + nBin * sizeof(HistData)), sizeof(char)); 

    real_0 Nbodies = st->nbody;
    mwbool islight = FALSE;//is it light matter?
    
    nbGetHistTrig(&histTrig, hp, ctx->leftHanded);

    NBodyHistogram* hist0 = mwCalloc(sizeof(NBodyHistogram) + nBin * sizeof(HistData), sizeof(char));
    NBodyHistogram* hist1 = mwCalloc(sizeof(NBodyHistogram) + nBin * sizeof(HistData), sizeof(char));
    NBodyHistogram* hist2 = mwCalloc(sizeof(NBodyHistogram) + nBin * sizeof(HistData), sizeof(char));
    NBodyHistogram* hist3 = mwCalloc(sizeof(NBodyHistogram) + nBin * sizeof(HistData), sizeof(char));
    NBodyHistogram* hist4 = mwCalloc(sizeof(NBodyHistogram) + nBin * sizeof(HistData), sizeof(char));
    NBodyHistogram* hist5 = mwCalloc(sizeof(NBodyHistogram) + nBin * sizeof(HistData), sizeof(char));

    // this only being for hist0 should be fine as long as the first histogram
    // (normalized counts) is always the one this info is accessed from
    /*  NEEDS TO BE SET FOR ALL OF THEM! */
    hist0->lambdaBins = lambdaBins;
    hist0->betaBins = betaBins;
    hist0->hasRawCounts = TRUE;
    hist0->params = *hp;
    mwbool lightMassUnlogged = TRUE;

    for (i = 0; i < Nbodies; i++)
    {
        const Body* b = &st->bodytab[i];
        if(Type(b) == BODY(islight))
        {
            if(lightMassUnlogged)
            {
                hist0->massPerParticle = Mass(b);
                hist1->massPerParticle = Mass(b);
                hist2->massPerParticle = Mass(b);
                hist3->massPerParticle = Mass(b);
                hist4->massPerParticle = Mass(b);
                hist5->massPerParticle = Mass(b);
                lightMassUnlogged = FALSE;
            }
            body_count++;
        }
    }

    hist0->totalSimulated = (unsigned int) body_count;

    // create the histograms for each variable only if they're wanted/needed
    // otherwise mark them as unused (false)

    /// always use normalized counts
    all->usage[0] = TRUE;
    all->histograms[0] = hist0; 

    // always use beta and vlos dispersion (for backwards compatibility)
    all->usage[1] = TRUE;
    all->histograms[1] = hist1;
    hist1->lambdaBins = lambdaBins;
    hist1->betaBins = betaBins;
    hist1->hasRawCounts = TRUE;
    hist1->params = *hp;
    hist1->totalSimulated = (unsigned int) body_count;

    all->usage[2] = TRUE;
    all->histograms[2] = hist2;
    hist2->lambdaBins = lambdaBins;
    hist2->betaBins = betaBins;
    hist2->hasRawCounts = TRUE;
    hist2->params = *hp;
    hist2->totalSimulated = (unsigned int) body_count;

    if(st->useVlos || st->useBetaComp || st-> useDist)
    {
        all->usage[3] = TRUE;
        all->histograms[3] = hist3;
        hist3->lambdaBins = lambdaBins;
        hist3->betaBins = betaBins;
        hist3->hasRawCounts = TRUE;
        hist3->params = *hp;
        hist3->totalSimulated = (unsigned int) body_count;

        all->usage[4] = TRUE;
        all->histograms[4] = hist4;
        hist4->lambdaBins = lambdaBins;
        hist4->betaBins = betaBins;
        hist4->hasRawCounts = TRUE;
        hist4->params = *hp;
        hist4->totalSimulated = (unsigned int) body_count;

        all->usage[5] = TRUE;
        all->histograms[5] = hist5;
        hist5->lambdaBins = lambdaBins;
        hist5->betaBins = betaBins;
        hist5->hasRawCounts = TRUE;
        hist5->params = *hp;
        hist5->totalSimulated = (unsigned int) body_count;
    }
    else
    {
        all->usage[3] = FALSE;
        free(hist3);
        all->usage[4] = FALSE;
        free(hist4);
        all->usage[5] = FALSE;
        free(hist5);
    }

    real_0 * use_velbody   = mwCalloc(body_count, sizeof(real_0));
    real_0 * use_betabody  = mwCalloc(body_count, sizeof(real_0));
    real_0 * use_distbody  = mwCalloc(body_count, sizeof(real_0));

    real * vlos      = mwCalloc(body_count, sizeof(real));       
    real * betas     = mwCalloc(body_count, sizeof(real));
    real * distances = mwCalloc(body_count, sizeof(real));

    real * bodyFrac  = mwCalloc(body_count*nBin, sizeof(real));
    for(i = 0; i < body_count*nBin; i++) bodyFrac[i] = ZERO_REAL;
       
    
    /* It does not make sense to ignore bins in a generated histogram */
    for (unsigned int i = 0; i < 6; ++i)
    {
        if(!all->usage[i]) continue;
        histData = all->histograms[i]->data;

        for (Histindex = 0; Histindex < nBin; ++Histindex)
        {
            histData[Histindex].rawCount    = ZERO_REAL;
            histData[Histindex].variable    = ZERO_REAL;
            histData[Histindex].sum         = ZERO_REAL;
            histData[Histindex].sq_sum      = ZERO_REAL;
            histData[Histindex].err         = ZERO_REAL;
            histData[Histindex].outliersRemoved = ZERO_REAL;
            histData[Histindex].useBin = TRUE;
        }
    }

 /* This code takes a single body and extrapolates it into a probability distribution.
    This code is only really needed for when we are running with AUTODIFF since the
    normal code does not propagate derivative information. If useContBins is FALSE,
    the original code runs instead. */
    if (ctx->useContBins)
    {
      #ifdef _OPENMP
        #pragma omp parallel for private(p, lambdaBetaR, lambda, beta, lambdaIndex, betaIndex, v_line_of_sight, location, i, j, bodyBinFrac, fullBodyFrac, tmp) shared(totalNum, use_betabody, use_velbody, use_distbody, vlos, betas, distances) schedule(dynamic)
      #endif
        for (p = st->bodytab; p < endp; ++p)
        {
            /* Only include bodies in models we aren't ignoring (like dark matter) */
            if (!ignoreBody(p))
            {

                /* Get the position in lbr coorinates */
                //mw_printf("P POS = [%.15f,%.15f,%.15f]\n", showRealValue(&X(&Pos(p))), showRealValue(&Y(&Pos(p))), showRealValue(&Z(&Pos(p))));
                lambdaBetaR = nbXYZToLambdaBeta(&histTrig, &Pos(p), ctx->sunGCDist);
                lambda = lambdaBetaR.x;
                beta = lambdaBetaR.y;
            
                use_betabody[p - st->bodytab] = DEFAULT_NOT_USE;//defaulted to not use body
                use_velbody[p - st->bodytab] = DEFAULT_NOT_USE;//defaulted to not use body
                use_distbody[p - st->bodytab] = DEFAULT_NOT_USE;
            
                vlos[p - st->bodytab]      = mw_real_const(DEFAULT_NOT_USE);//default vlos
                betas[p - st->bodytab]     = mw_real_const(DEFAULT_NOT_USE);
                distances[p - st->bodytab] = mw_real_const(DEFAULT_NOT_USE);

                /* Find the indices */
                lambdaIndex = (unsigned int) mw_floor_0((showRealValue(&lambda) - lambdaStart) / lambdaSize);
                betaIndex = (unsigned int) mw_floor_0((showRealValue(&beta) - betaStart) / betaSize);

                /* Calculate vlos and distance */
                v_line_of_sight = calc_vLOS(&Vel(p), &Pos(p), ctx->sunGCDist);
                location = calc_distance(&Pos(p), ctx->sunGCDist);

                /* Log the Disp quantities */
                vlos[p - st->bodytab] = v_line_of_sight;  
                betas[p - st->bodytab] = beta;
                distances[p - st->bodytab] = location;

                /* For each bin, calculate a FB5 distribution for the body and determine what fraction goes in each bin */
                fullBodyFrac = ZERO_REAL;
                for(i = 0; i < nBin; i++)
                {
                    unsigned int fracIndex = (p - st->bodytab)*nBin + i;
                    int test_beta_index = i % betaBins;
                    int test_lambda_index = (i - test_beta_index)/betaBins;
                    /* Only perform the body fraction calculation for nearby bins to save time */
                    if((mw_abs_0((real_0) test_lambda_index - (real_0) lambdaIndex) <= ctx->bleedInRange) && (mw_abs_0((real_0) test_beta_index - (real_0) betaIndex) <= ctx->bleedInRange))
                    {
                        bodyBinFrac = getBodyBinFrac(ctx, hp, p, i);  //Returns log of body fraction
                        bodyBinFrac = mw_exp(&bodyBinFrac);
                        //mw_printf("bodyBinFrac = %.15f\n", showRealValue(&bodyBinFrac));
                    }
                    else
                    {
                        bodyBinFrac = ZERO_REAL;
                    }
                    bodyFrac[fracIndex] = bodyBinFrac;
                    for(j = 0; j < 6; j++)
                        if(all->usage[j]) all->histograms[j]->data[i].rawCount = mw_add(&all->histograms[j]->data[i].rawCount, &bodyBinFrac);
                    fullBodyFrac = mw_add(&fullBodyFrac, &bodyBinFrac);

                    /* Sums and Square Sums for chi^2 calcs must be weighted by the body fraction in the bin */
                    if(all->usage[1])
                    {
                        /* each of these are components of the beta disp */
                        tmp = mw_mul(&bodyBinFrac, &beta);
                        all->histograms[1]->data[i].sum = mw_add(&all->histograms[1]->data[i].sum, &tmp);
                        tmp = sqr(&beta);
                        tmp = mw_mul(&bodyBinFrac, &tmp);
                        all->histograms[1]->data[i].sq_sum = mw_add(&all->histograms[1]->data[i].sq_sum, &tmp);
                    }

                    if(all->usage[2])
                    {
                        /* each of these are components of the vel disp */
                        tmp = mw_mul(&bodyBinFrac, &v_line_of_sight);
                        all->histograms[2]->data[i].sum = mw_add(&all->histograms[2]->data[i].sum, &tmp);
                        tmp = sqr(&v_line_of_sight);
                        tmp = mw_mul(&bodyBinFrac, &tmp);
                        all->histograms[2]->data[i].sq_sum = mw_add(&all->histograms[2]->data[i].sq_sum, &tmp);
                    }

                    if(all->usage[3])
                    {
                        /* each of these are components of the vel disp, which is used for vel avg */
                        tmp = mw_mul(&bodyBinFrac, &v_line_of_sight);
                        all->histograms[3]->data[i].sum = mw_add(&all->histograms[3]->data[i].sum, &tmp);
                        tmp = sqr(&v_line_of_sight);
                        tmp = mw_mul(&bodyBinFrac, &tmp);
                        all->histograms[3]->data[i].sq_sum = mw_add(&all->histograms[3]->data[i].sq_sum, &tmp);
                    }

                    if(all->usage[4])
                    {
                        /* each of these are components of the beta disp, which is used for beta avg */
                        tmp = mw_mul(&bodyBinFrac, &beta);
                        all->histograms[4]->data[i].sum = mw_add(&all->histograms[4]->data[i].sum, &tmp);
                        tmp = sqr(&beta);
                        tmp = mw_mul(&bodyBinFrac, &tmp);
                        all->histograms[4]->data[i].sq_sum = mw_add(&all->histograms[4]->data[i].sq_sum, &tmp);
                    }

                    if(all->usage[5])
                    {
                        /* average distance */
                        tmp = mw_mul(&bodyBinFrac, &location);
                        all->histograms[5]->data[i].sum = mw_add(&all->histograms[5]->data[i].sum, &tmp);
                        tmp = sqr(&location);
                        tmp = mw_mul(&bodyBinFrac, &tmp);
                        all->histograms[5]->data[i].sq_sum = mw_add(&all->histograms[5]->data[i].sq_sum, &tmp);
                    }
                }
                totalNum = mw_add(&totalNum, &fullBodyFrac);
            }
        }
    }
    else //ORIGINAL DISCRETE BINNING CODE
    {
        for (p = st->bodytab; p < endp; ++p)
        {
            /* Only include bodies in models we aren't ignoring (like dark matter) */
            if (!ignoreBody(p))
            {

                /* Get the position in lbr coorinates */
                //mw_printf("P POS = [%.15f,%.15f,%.15f]\n", showRealValue(&X(&Pos(p))), showRealValue(&Y(&Pos(p))), showRealValue(&Z(&Pos(p))));
                lambdaBetaR = nbXYZToLambdaBeta(&histTrig, &Pos(p), ctx->sunGCDist);
                lambda = lambdaBetaR.x;
                beta = lambdaBetaR.y;
            
                use_betabody[p - st->bodytab] = DEFAULT_NOT_USE;//defaulted to not use body
                use_velbody[p - st->bodytab] = DEFAULT_NOT_USE;//defaulted to not use body
                use_distbody[p - st->bodytab] = DEFAULT_NOT_USE;
            
                vlos[p - st->bodytab]      = mw_real_const(DEFAULT_NOT_USE);//default vlos
                betas[p - st->bodytab]     = mw_real_const(DEFAULT_NOT_USE);
                distances[p - st->bodytab] = mw_real_const(DEFAULT_NOT_USE);

                /* Find the indices */
                lambdaIndex = (unsigned int) mw_floor_0((showRealValue(&lambda) - lambdaStart) / lambdaSize);
                betaIndex = (unsigned int) mw_floor_0((showRealValue(&beta) - betaStart) / betaSize);

                /* Check if the position is within the bounds of the histogram */
                if (lambdaIndex < lambdaBins && betaIndex < betaBins)   
                {
                    Histindex = lambdaIndex * betaBins + betaIndex;
                    use_betabody[p - st->bodytab] = Histindex;//if body is in hist, mark which hist bin
                    use_velbody[p - st->bodytab] = Histindex;
                    use_distbody[p - st->bodytab] = Histindex;

                    for(i = 0; i < 6; i++)
                        if(all->usage[i]) all->histograms[i]->data[Histindex].rawCount = mw_add_s(&all->histograms[i]->data[Histindex].rawCount, 1.0);

                    totalNum = mw_add_s(&totalNum, 1.0);
                
                    v_line_of_sight = calc_vLOS(&Vel(p), &Pos(p), ctx->sunGCDist);//calc the heliocentric line of sight vel
                    location = calc_distance(&Pos(p), ctx->sunGCDist);

                    vlos[p - st->bodytab] = v_line_of_sight;//store the vlos's so as to not have to recalc  
                    betas[p - st->bodytab] = beta;
                    distances[p - st->bodytab] = location;

                    if(all->usage[1])
                    {
                        /* each of these are components of the beta disp */
                        all->histograms[1]->data[Histindex].sum = mw_add(&all->histograms[1]->data[Histindex].sum, &beta);
                        tmp = sqr(&beta);
                        all->histograms[1]->data[Histindex].sq_sum = mw_add(&all->histograms[1]->data[Histindex].sq_sum, &tmp);
                    }
                    if(all->usage[2])
                    {
                        /* each of these are components of the vel disp */
                        all->histograms[2]->data[Histindex].sum = mw_add(&all->histograms[2]->data[Histindex].sum, &v_line_of_sight);
                        tmp = sqr(&v_line_of_sight);
                        all->histograms[2]->data[Histindex].sq_sum = mw_add(&all->histograms[2]->data[Histindex].sq_sum, &tmp);
                    }
                    if(all->usage[3])
                    {
                        /* each of these are components of the vel disp, which is used for vel avg */
                        all->histograms[3]->data[Histindex].sum = mw_add(&all->histograms[3]->data[Histindex].sum, &v_line_of_sight);
                        tmp = sqr(&v_line_of_sight);
                        all->histograms[3]->data[Histindex].sq_sum = mw_add(&all->histograms[3]->data[Histindex].sq_sum, &tmp);
                    }
                    if(all->usage[4])
                    {
                        /* each of these are components of the beta disp, which is used for beta avg */
                        all->histograms[4]->data[Histindex].sum = mw_add(&all->histograms[4]->data[Histindex].sum, &beta);
                        tmp = sqr(&beta);
                        all->histograms[4]->data[Histindex].sq_sum = mw_add(&all->histograms[4]->data[Histindex].sq_sum, &tmp);
                    }
                    if(all->usage[5])
                    {
                        /* average distance */
                        all->histograms[5]->data[Histindex].sum = mw_add(&all->histograms[5]->data[Histindex].sum, &location);
                        tmp = sqr(&location);
                        all->histograms[5]->data[Histindex].sq_sum = mw_add(&all->histograms[5]->data[Histindex].sq_sum, &tmp);
                    }
                }
            }
        }
    }

    for(i = 0; i < 6; i++)
        if(all->usage[i]) all->histograms[i]->totalNum = totalNum; /* Total particles in range */

    if(all->usage[1])    // if using beta disp
        nbCalcDisp(all->histograms[1], TRUE, ctx->BetaCorrect);
    if(all->usage[2])    // if using vel disp
        nbCalcDisp(all->histograms[2], TRUE, ctx->VelCorrect);
    if(all->usage[3])    // if using vlos average
        nbCalcDisp(all->histograms[3], TRUE, ctx->VelCorrect);
    if(all->usage[4])    // if using beta average
        nbCalcDisp(all->histograms[4], TRUE, ctx->BetaCorrect);
    if(all->usage[5])    // if using distance average
        nbCalcDisp(all->histograms[5], TRUE, ctx->DistCorrect);

    /* these converge somewhere between 3 and 6 iterations */
    if(all->usage[1])
    {
        for(i = 0; i < IterMax; i++)
        {
            nbRemoveOutliers(st, all->histograms[1], use_betabody, betas, ctx->BetaSigma, nBin, ctx->useContBins, bodyFrac);
            nbCalcDisp(all->histograms[1], FALSE, ctx->BetaCorrect);
        }
    }
    if(all->usage[2])
    {
        for(i = 0; i < IterMax; i++)
        {
            nbRemoveOutliers(st, all->histograms[2], use_velbody, vlos, ctx->VelSigma, nBin, ctx->useContBins, bodyFrac);
            nbCalcDisp(all->histograms[2], FALSE, ctx->VelCorrect);
        }
    }

    // calculation of average velocity and average beta values
    // dispersions are already calculated and in histogram - this is used to calculate error
    if(all->usage[3]) // vlos average
    {
        for(i = 0; i < IterMax; i++)
        {
            nbRemoveOutliers(st, all->histograms[3], use_velbody, vlos, ctx->VelSigma, nBin, ctx->useContBins, bodyFrac);
            nbCalcDisp(all->histograms[3], FALSE, ctx->VelCorrect);
        }
        for (i = 0; i < nBin; ++i)
        {
            real vdenom = mw_sub(&all->histograms[3]->data[i].rawCount, &all->histograms[3]->data[i].outliersRemoved);
            if(showRealValue(&vdenom) > 10) // no data for the bin
            {
                // calculates error first because the dispersion is stored as the variable at the moment
                // dispersion is used for error calc, then variable is overwritten as the average vlos (as it should be)
                tmp = mw_sqrt(&vdenom);
                all->histograms[3]->data[i].err = mw_div(&all->histograms[3]->data[i].variable, &tmp);
                all->histograms[3]->data[i].variable = mw_div(&all->histograms[3]->data[i].sum, &vdenom);
            }
            else
            {
                all->histograms[3]->data[i].err = mw_real_const(-1);
                all->histograms[3]->data[i].variable = mw_real_const(0);
            }
        }
    }
    if(all->usage[4]) // beta average
    {
        for(i = 0; i < IterMax; i++)
        {
            nbRemoveOutliers(st, all->histograms[4], use_betabody, betas, ctx->BetaSigma, nBin, ctx->useContBins, bodyFrac);
            nbCalcDisp(all->histograms[4], FALSE, ctx->BetaCorrect);
        }
        for (i = 0; i < nBin; ++i)
        {
            real bdenom = mw_sub(&all->histograms[4]->data[i].rawCount, &all->histograms[4]->data[i].outliersRemoved);
            if(showRealValue(&bdenom) > 10) // no data for the bin
            {
                // calculates error first because the dispersion is stored as the variable at the moment
                // dispersion is used for error calc, then variable is overwritten as the average beta (as it should be)
                tmp = mw_sqrt(&bdenom);
                all->histograms[4]->data[i].err = mw_div(&all->histograms[4]->data[i].variable, &tmp);
                all->histograms[4]->data[i].variable = mw_div(&all->histograms[4]->data[i].sum, &bdenom);
            }
            else
            {
                all->histograms[4]->data[i].err = mw_real_const(-1);
                all->histograms[4]->data[i].variable = mw_real_const(0);
            }
        }
    }
    if(all->usage[5]) //distance calculation
    {
        for(i = 0; i < IterMax; ++i)
        {
            nbRemoveOutliers(st, all->histograms[5], use_distbody, distances, ctx->DistSigma, nBin, ctx->useContBins, bodyFrac);
            nbCalcDisp(all->histograms[5], FALSE, ctx->DistCorrect);
        }
        for (i = 0; i < nBin; ++i)
        {
            real ddenom = mw_sub(&all->histograms[5]->data[i].rawCount, &all->histograms[5]->data[i].outliersRemoved);
            if(showRealValue(&ddenom) > 10)
            {
                // calculates error first because the dispersion is stored as the variable at the moment
                // dispersion is used for error calc, then variable is overwritten as the average distance (as it should be)
                tmp = mw_sqrt(&ddenom);
                all->histograms[5]->data[i].err = mw_div(&all->histograms[5]->data[i].variable, &tmp);
                all->histograms[5]->data[i].variable = mw_div(&all->histograms[5]->data[i].sum, &ddenom);
            }
            else
            {
                all->histograms[5]->data[i].err = mw_real_const(-1);
                all->histograms[5]->data[i].variable = mw_real_const(0);
            }
        }
    }
    
    nbNormalizeHistogram(all->histograms[0]); // sets up normalized counts histogram

    //freeing up mwcallocs (which is essentially a malloc)
    free(use_velbody);
    free(use_betabody);
    free(use_distbody);
    free(vlos);
    free(betas);
    free(distances);
    free(bodyFrac);
    
    return all;
}


/* Read in a histogram from a file for calculating a likelihood value.
 */
MainStruct* nbReadHistogram(const char* histogramFile)
{
    FILE* f;
    int rc = 0;
    size_t fsize = 0;
    MainStruct* all = NULL;
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
    mwbool buildHist = FALSE; /* only want to build the histogrma once */
    real_0 nGen = 0.0;    /* Number of particles read from the histogram */
    unsigned int totalSim = 0;  /*Total number of simulated particles read from the histogram */
    unsigned int lambdaBins = 0; /* Number of bins in lambda direction */
    unsigned int betaBins = 0; /* Number of bins in beta direction */
    mwbool used = FALSE;  /* indicates whether or not to expect extra histogram parameters */
    real_0 mass = 0.0;            /*mass per particle read from the histogram */
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
        // if using the new histogram parameters
        // there will be a line beginning with a ! to indicate
        // the extra information to be read in
        if (lineBuf[0] == '!')
        {
            used = TRUE;
            continue;
        }

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
            real_0 phi, theta, psi;

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
            rc = sscanf(lineBuf, " n = %lf \n", &nGen);
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

        if(readBetaBins && !buildHist) // only build the histogram once
        {
            unsigned int nBin = lambdaBins * betaBins;
            all = mwCalloc(6*(sizeof(NBodyHistogram) + nBin * sizeof(HistData)), sizeof(char));

            NBodyHistogram* hist0 = mwCalloc(sizeof(NBodyHistogram) + nBin * sizeof(HistData), sizeof(char));
            all->usage[0] = TRUE;
            all->histograms[0] = hist0;

            NBodyHistogram* hist1 = mwCalloc(sizeof(NBodyHistogram) + nBin * sizeof(HistData), sizeof(char));
            all->usage[1] = TRUE;
            all->histograms[1] = hist1;

            NBodyHistogram* hist2 = mwCalloc(sizeof(NBodyHistogram) + nBin * sizeof(HistData), sizeof(char));
            all->usage[2] = TRUE;
            all->histograms[2] = hist2;

            if(used)
            {
                NBodyHistogram* hist3 = mwCalloc(sizeof(NBodyHistogram) + nBin * sizeof(HistData), sizeof(char));
                all->usage[3] = TRUE;
                all->histograms[3] = hist3;

                NBodyHistogram* hist4 = mwCalloc(sizeof(NBodyHistogram) + nBin * sizeof(HistData), sizeof(char));
                all->usage[4] = TRUE;
                all->histograms[4] = hist4;

                NBodyHistogram* hist5 = mwCalloc(sizeof(NBodyHistogram) + nBin * sizeof(HistData), sizeof(char));
                all->usage[5] = TRUE;
                all->histograms[5] = hist5;
            }
            else
            {
                all->usage[3] = FALSE;
                all->usage[4] = FALSE;
                all->usage[5] = FALSE;
            }

            buildHist = TRUE;
        }

        unsigned int useBin = 0;
        double lambda = 0;
        double beta = 0;
        double variable[6];
        double errors[6];

        if(used) // new histogram output, there are more parameters to read in
        {
            rc = sscanf(lineBuf,
                        "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
                        &useBin,
                        &lambda,
                        &beta,
                        &variable[0],
                        &errors[0],
                        &variable[1],
                        &errors[1],
                        &variable[2],
                        &errors[2],
                        &variable[3],
                        &errors[3],
                        &variable[4],
                        &errors[4],
                        &variable[5],
                        &errors[5]);


            // only save the numbers that are used
            for(int i = 0; i < 6; i++)
            {
                if(all->usage[i])
                {
                    all->histograms[i]->data[fileCount].useBin = useBin;
                    all->histograms[i]->data[fileCount].lambda = mw_real_const(lambda);
                    all->histograms[i]->data[fileCount].beta = mw_real_const(beta);
                    all->histograms[i]->data[fileCount].variable = mw_real_const(variable[i]);
                    all->histograms[i]->data[fileCount].err = mw_real_const(errors[i]);
                }
            }
            
            if (rc != 15)
            {
                mw_printf("Error reading histogram line %d: %s", lineNum, lineBuf);
                error = TRUE;
                break;
            }
        }

        else // default histogram output (no new parameters)
        {
            rc = sscanf(lineBuf,
                        "%d %lf %lf %lf %lf %lf %lf %lf %lf\n",
                        &useBin,
                        &lambda,
                        &beta,
                        &variable[0],
                        &errors[0],
                        &variable[1],
                        &errors[1],
                        &variable[2],
                        &errors[2]);


            // only save the numbers that are used
            for(int i = 0; i < 3; i++)
            {
                if(all->usage[i])
                {
                    all->histograms[i]->data[fileCount].useBin = useBin;
                    all->histograms[i]->data[fileCount].lambda = mw_real_const(lambda);
                    all->histograms[i]->data[fileCount].beta = mw_real_const(beta);
                    all->histograms[i]->data[fileCount].variable = mw_real_const(variable[i]);
                    all->histograms[i]->data[fileCount].err = mw_real_const(errors[i]);
                }
            }
            
            if (rc != 9)
            {
                mw_printf("Error reading histogram line %d: %s", lineNum, lineBuf);
                error = TRUE;
                break;
            }
        }

        ++fileCount;
    }

    fclose(f);

    // checking to make sure there's something in the histogram
    unsigned int num = 3;
    if(used) num = 6; // more histogram parameters for the new output

    unsigned int used_hist = 0;
    for(unsigned int i = 0; i < num; i++)
        if(all->usage[i]) used_hist++;

    if (error || used_hist == 0)
    {
        for(int i = 0; i < num; i++)
            free(all->histograms[i]);
        free(all);
        return NULL;
    }

    for(int i = 0; i < num; i++)
    {
        if(all->usage[i])
        {
            all->histograms[i]->lambdaBins = lambdaBins;
            all->histograms[i]->betaBins = betaBins;
            all->histograms[i]->totalNum = mw_real_const(nGen);
            all->histograms[i]->totalSimulated = totalSim;
            all->histograms[i]->massPerParticle = mw_real_const(mass);
        }
    }
    
    return all;
}
/*
mwvector getHistogramCenter(const NBodyHistogram* hist, HistogramParams* hp, NBodyState* st, NBodyCtx* ctx){
    int i;
    real bestLambda, highestCount = 0, tmpCount;
    mwvector center, pos;
    center.x = 0;
    center.y = 0;
    center.z = 0;

    //get lambda from histogram
    for(i = 0; i < hist->lambdaBins; i++){
        tmpCount = hist->data[i].count;
        if(tmpCount > highestCount){
            highestCount = tmpCount;
            bestLambda = hist->data[i].lambda;
        }
    }

    real lambdaBinSize = getAngleDiffDegrees(hp->lambdaStart, hp->lambdaEnd)/lambdaBins;
    mwvector lambdaBeta;
    NBHistTrig histTrig;
    real maxLambda, minLambda, betaCount, betas = 0;
    nbGetHistTrig(&histTrig, hp);
    maxLambda = bestLambda + lambdaBinSize/2;
    minLambda = bestLambda - lambdaBinSize/2;
    if(maxLambda >= 360)
        maxLambda -= 360;
    if(minLambda < 0)
        minLambda += 360;
    //get beta from bodies (mean beta in the lambda bin)
    for(i = 0; i < st->nbody; i++){
        pos = st->bodytab[i].bodynode.pos;
        lambdaBeta = nbXYZToLambdaBeta(&histTrig, pos, ctx->sunGCDist);
        if(angleIsBetween(minLambda, maxLambda, L(lambdaBeta))){
            betas += B(lambdaBeta);
            betaCount += 1;
        }
    }

    center = nb
    
    return center;
}*/
