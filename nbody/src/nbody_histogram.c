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
            total += histogram->data[i].variable;
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
                                   NBodyState* st)
{
    int nbody = st->nbody;
    real bestLikelihood_time = st->bestLikelihood_time;
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
    real barTimeError = bestLikelihood_time - st->previousForwardTime;
    if(ctx->pot.disk2.type == _NO_DISK){//no bar
        barTimeError = 0;
    }else if(ctx->calibrationRuns == 0){//no calibration but bar on
        barTimeError = bestLikelihood_time - ctx->timeEvolve;
    }
    real barAngleError = barTimeError * ctx->pot.disk2.patternSpeed;
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
            "# Bar/LMC Calibration Time Error = %f\n"
            "# Bar Angle Error = %f rad\n"
            "#\n",
            nbody,
            ctx->timeBack,
            bestLikelihood_time,
            bestLikelihood,
            ctx->timestep,
            ctx->sunGCDist,
            showCriterionT(ctx->criterion),
            ctx->theta,
            showBool(ctx->useQuad),
            mw_sqrt(ctx->eps2),
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
    real output[12]; // for outputting the data

    unsigned int nBin;
    nBin = all->histograms[0]->lambdaBins * all->histograms[0]->betaBins;

    // outputs these numbers from one of the histograms
    // at this point, these numbers should be the same for all histograms anyway
    mw_boinc_print(f, "<histogram>\n");
    fprintf(f, "n = %u\n", all->histograms[0]->totalNum);
    fprintf(f, "massPerParticle = %12.15f\n", all->histograms[0]->massPerParticle);
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
                output[k] = all->histograms[j]->data[i].variable;
                output[k+1] = all->histograms[j]->data[i].err;
            }
            k+=2;
        }

        if(all->usage[3]) // if this histogram has been used, all histograms were calculated
        {
            fprintf(f,
                    "%d %12.15f %12.15f %12.15f %12.15f %12.15f %12.15f %12.15f %12.15f %12.15f %12.15f %12.15f %12.15f %12.15f %12.15f\n",
                    storedData.useBin,
                    storedData.lambda,
                    storedData.beta,
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
                    storedData.lambda,
                    storedData.beta,
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
            histData[Histindex].variable  = count / totalNum;
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
MainStruct* nbCreateHistogram(const NBodyCtx* ctx,        /* Simulation context */
                                  const NBodyState* st,       /* Final state of the simulation */
                                  const HistogramParams* hp)  /* Range of histogram to create */
{
    real location;
    real lambda;
    real beta;
    real v_line_of_sight;
    mwvector lambdaBetaR;
    unsigned int lambdaIndex;
    unsigned int betaIndex;
    unsigned int Histindex;
    unsigned int totalNum = 0;
    HistData* histData;
    Body* p;
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
    /*unsigned int IterMax = 6;*/   /*Default value for IterMax*/
    unsigned int nBin = lambdaBins * betaBins;
    unsigned int body_count = 0;
    unsigned int ub_counter = 0;
    
    MainStruct* all = mwCalloc(6*(sizeof(NBodyHistogram) + nBin * sizeof(HistData)), sizeof(char)); 

    real Nbodies = st->nbody;
    mwbool islight = FALSE;//is it light matter?
    
    nbGetHistTrig(&histTrig, hp);

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

    for (int i = 0; i < Nbodies; i++)
    {
        const Body* b = &st->bodytab[i];
        if(Type(b) == BODY(islight))
        {
            hist0->massPerParticle = Mass(b);
            hist1->massPerParticle = Mass(b);
            hist2->massPerParticle = Mass(b);
            hist3->massPerParticle = Mass(b);
            hist4->massPerParticle = Mass(b);
            hist5->massPerParticle = Mass(b);
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
        all->usage[4] = TRUE;
        all->usage[5] = TRUE;
        all->histograms[3] = hist3;
        hist3->lambdaBins = lambdaBins;
        hist3->betaBins = betaBins;
        hist3->hasRawCounts = TRUE;
        hist3->params = *hp;
        hist3->totalSimulated = (unsigned int) body_count;

        all->histograms[4] = hist4;
        hist4->lambdaBins = lambdaBins;
        hist4->betaBins = betaBins;
        hist4->hasRawCounts = TRUE;
        hist4->params = *hp;
        hist4->totalSimulated = (unsigned int) body_count;

        all->usage[3] = TRUE;
        all->usage[4] = TRUE;
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

    real * use_velbody   = mwCalloc(body_count, sizeof(real));
    real * use_betabody  = mwCalloc(body_count, sizeof(real));
    real * use_distbody  = mwCalloc(body_count, sizeof(real));

    real * vlos      = mwCalloc(body_count, sizeof(real));       
    real * betas     = mwCalloc(body_count, sizeof(real));
    real * distances = mwCalloc(body_count, sizeof(real));   
       
    
    /* It does not make sense to ignore bins in a generated histogram */
    for (unsigned int i = 0; i < 6; ++i)
    {
        if(!all->usage[i]) continue;
        histData = all->histograms[i]->data;

        for (Histindex = 0; Histindex < nBin; ++Histindex)
        {
            histData[Histindex].rawCount    = 0;
            histData[Histindex].variable    = 0.0;
            histData[Histindex].sum         = 0.0;
            histData[Histindex].sq_sum      = 0.0;
            histData[Histindex].err         = 0.0;
            histData[Histindex].outliersRemoved = 0.0;
            histData[Histindex].useBin = TRUE;
        }
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
            use_distbody[ub_counter] = DEFAULT_NOT_USE;
            
            vlos[ub_counter]     = DEFAULT_NOT_USE;//default vlos
            betas[ub_counter]    = DEFAULT_NOT_USE;
            distances[ub_counter] = DEFAULT_NOT_USE;

            /* Find the indices */
            lambdaIndex = (unsigned int) mw_floor((lambda - lambdaStart) / lambdaSize);
            betaIndex = (unsigned int) mw_floor((beta - betaStart) / betaSize);

            /* Check if the position is within the bounds of the histogram */
            if (lambdaIndex < lambdaBins && betaIndex < betaBins)   
            {   
                Histindex = lambdaIndex * betaBins + betaIndex;
                use_betabody[ub_counter] = Histindex;//if body is in hist, mark which hist bin
                use_velbody[ub_counter] = Histindex;
                use_distbody[ub_counter] = Histindex;
                
                for(int i = 0; i < 6; i++)
                    if(all->usage[i]) all->histograms[i]->data[Histindex].rawCount++;

                ++totalNum;
                
                v_line_of_sight = calc_vLOS(Vel(p), Pos(p), ctx->sunGCDist);//calc the heliocentric line of sight vel
                location = calc_distance(Pos(p), ctx->sunGCDist);

                vlos[ub_counter] = v_line_of_sight;//store the vlos's so as to not have to recalc  
                betas[ub_counter] = beta;
                distances[ub_counter] = location;

                if(all->usage[1])
                {
                    /* each of these are components of the beta disp */
                    all->histograms[1]->data[Histindex].sum += beta;
                    all->histograms[1]->data[Histindex].sq_sum += sqr(beta);
                }
                if(all->usage[2])
                {
                    /* each of these are components of the vel disp */
                    all->histograms[2]->data[Histindex].sum += v_line_of_sight;
                    all->histograms[2]->data[Histindex].sq_sum += sqr(v_line_of_sight);
                }
                if(all->usage[3])
                {
                    /* each of these are components of the vel disp, which is used for vel avg */
                    all->histograms[3]->data[Histindex].sum += v_line_of_sight;
                    all->histograms[3]->data[Histindex].sq_sum += sqr(v_line_of_sight);
                }
                if(all->usage[4])
                {
                    /* each of these are components of the beta disp, which is used for beta avg */
                    all->histograms[4]->data[Histindex].sum += beta;
                    all->histograms[4]->data[Histindex].sq_sum += sqr(beta);
                }
                if(all->usage[5])
                {
                    /* average distance */
                    all->histograms[5]->data[Histindex].sum += location;
                    all->histograms[5]->data[Histindex].sq_sum += sqr(location);
                }
            
            }
            ub_counter++;
        }
    }
   
    for(int i = 0; i < 6; i++)
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
        nbCalcDisp(all->histograms[5], TRUE, ctx->DistCorrect); //using beta correct for now, will need to add a dist correct

    /* these converge somewhere between 3 and 6 iterations */
    if(all->usage[1])
    {
        for(unsigned int i = 0; i < IterMax; i++)
        {
            nbRemoveOutliers(st, all->histograms[1], use_betabody, betas, ctx->BetaSigma, ctx->sunGCDist);
            nbCalcDisp(all->histograms[1], FALSE, ctx->BetaCorrect);
        }
    }
    if(all->usage[2])
    {
        for(unsigned int i = 0; i < IterMax; i++)
        {
            nbRemoveOutliers(st, all->histograms[2], use_velbody, vlos, ctx->VelSigma, ctx->sunGCDist);
            nbCalcDisp(all->histograms[2], FALSE, ctx->VelCorrect);
        }
    }

    // calculation of average velocity and average beta values
    // dispersions are already calculated and in histogram - this is used to calculate error
    if(all->usage[3]) // vlos average
    {
        for(unsigned int i = 0; i < IterMax; i++)
        {
            nbRemoveOutliers(st, all->histograms[3], use_velbody, vlos, ctx->VelSigma, ctx->sunGCDist);
            nbCalcDisp(all->histograms[3], FALSE, ctx->VelCorrect);
        }
        for (unsigned int i = 0; i < nBin; ++i)
        {
            int vdenom = all->histograms[3]->data[i].rawCount - all->histograms[3]->data[i].outliersRemoved;
            if(vdenom != 0) // no data for the bin
            {
                // calculates error first because the dispersion is stored as the variable at the moment
                // dispersion is used for error calc, then variable is overwritten as the average vlos (as it should be)
                all->histograms[3]->data[i].err = all->histograms[3]->data[i].variable / sqrt(vdenom);
                all->histograms[3]->data[i].variable = all->histograms[3]->data[i].sum / vdenom;
            }
            else
            {
                all->histograms[3]->data[i].err = 0;
                all->histograms[3]->data[i].variable = 0;
            }
        }
    }
    if(all->usage[4]) // beta average
    {
        for(unsigned int i = 0; i < IterMax; i++)
        {
            nbRemoveOutliers(st, all->histograms[4], use_betabody, betas, ctx->BetaSigma, ctx->sunGCDist);
            nbCalcDisp(all->histograms[4], FALSE, ctx->BetaCorrect);
        }
        for (unsigned int i = 0; i < nBin; ++i)
        {
            int bdenom = all->histograms[4]->data[i].rawCount - all->histograms[4]->data[i].outliersRemoved;
            if(bdenom != 0) // no data for the bin
            {
                // calculates error first because the dispersion is stored as the variable at the moment
                // dispersion is used for error calc, then variable is overwritten as the average beta (as it should be)
                all->histograms[4]->data[i].err = all->histograms[4]->data[i].variable / sqrt(bdenom);
                all->histograms[4]->data[i].variable = all->histograms[4]->data[i].sum / bdenom;
            }
            else
            {
                all->histograms[4]->data[i].err = 0;
                all->histograms[4]->data[i].variable = 0;
            }
        }
    }
    if(all->usage[5]) //distance calculation
    {
        for(unsigned int i = 0; i < IterMax; ++i)
        {
            nbRemoveOutliers(st, all->histograms[5], use_distbody, distances, ctx->DistSigma, ctx->sunGCDist);
            nbCalcDisp(all->histograms[5], FALSE, ctx->DistCorrect);
        }
        for (unsigned int i = 0; i < nBin; ++i)
        {
            int ddenom = all->histograms[5]->data[i].rawCount - all->histograms[5]->data[i].outliersRemoved;
            if(ddenom != 0)
            {
                all->histograms[5]->data[i].err = all->histograms[5]->data[i].variable / sqrt(ddenom);
                all->histograms[5]->data[i].variable  = all->histograms[5]->data[i].sum / ddenom;
            }
            else
            {
                all->histograms[5]->data[i].err = 0;
                all->histograms[5]->data[i].variable = 0;
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
    unsigned int nGen = 0;    /* Number of particles read from the histogram */
    unsigned int totalSim = 0;  /*Total number of simulated particles read from the histogram */
    unsigned int lambdaBins = 0; /* Number of bins in lambda direction */
    unsigned int betaBins = 0; /* Number of bins in beta direction */
    mwbool used = FALSE;  /* indicates whether or not to expect extra histogram parameters */
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
                    all->histograms[i]->data[fileCount].lambda = lambda;
                    all->histograms[i]->data[fileCount].beta = beta;
                    all->histograms[i]->data[fileCount].variable = variable[i];
                    all->histograms[i]->data[fileCount].err = errors[i];
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
                    all->histograms[i]->data[fileCount].lambda = lambda;
                    all->histograms[i]->data[fileCount].beta = beta;
                    all->histograms[i]->data[fileCount].variable = variable[i];
                    all->histograms[i]->data[fileCount].err = errors[i];
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
            all->histograms[i]->totalNum = nGen;
            all->histograms[i]->totalSimulated = totalSim;
            all->histograms[i]->massPerParticle = mass;
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
