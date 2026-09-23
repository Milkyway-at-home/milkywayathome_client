/*
 * Copyright (c) 2012 Rensselaer Polytechnic Institute
 * Copyright (c) 2016-2018 Siddhartha Shelton
 * 
 * This file is part of Milkway@Home.
 *
 * Milkyway@Home is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Milkyway@Home is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Milkyway@Home.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "nbody_mass.h"
#include "nbody_defaults.h"
#include "milkyway_math.h"
#include "nbody_types.h"
#include "nbody_coordinates.h"


// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // 
// There functions are involved in calculating a binomial distribution
/*In order to decrease the size of the numbers
 * computed all these functions are
 * calculated in log space*/
static real factorial(int n)
{
     int counter;
     real result = 0.0;

     for (counter = n; counter >= 1; counter--)
       {
          result += mw_log((real) counter);
       }

     return result;
}


static real choose(int n, int c)
{
    unsigned int i;
    real result = 0.0;
    
    /* This for loop calulates log(n!/(n-c)!) */
    for (i = n - c + 1; i <= (unsigned int) n; ++i)
    {
        result += mw_log(i);
    }
    result -= factorial(c);
    return result;
}

real probability_match(int n, real ktmp, real pobs)
{
    real result = 0.0;

    /*
     * Previously, this function took in k as an int. Bad move.
     * This function was called twice, one of which sent a real valued k: (int) k1 and (real) k2
     * That real k2 was then converted to int. Could result in converted (int) k1 != (int) k2 when k1 = k2. 
     * Special result was poor likelihood for some histograms when check against themselves!
     * General results: unknown. But probably not good. (most likely caused different machines to report
     * different likelihood values).
     * 
     */
    int k = (int) mw_round(ktmp);    //patch. See above. 
    //The previous calculation does not return the right values.  Furthermore, we need a zeroed metric.                                                                                              
    result =  (real) choose(n, k);
    result += k * mw_log(pobs); 
    result += (n - k) * mw_log(1.0 - pobs);
    
    
    return mw_exp(result);
}
// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // 
// IMPLEMENTATION OF GAMMA FUNCTIONS. COMPLETE AND INCOMPLETE
// These functions are adapted from the Numerical Recipes in C 2nd ed except for gammln.
/* Returns ln(Gamma(z)) via the Lanczos approximation (NR 3rd ed). */
real gammln(const real z) 
{
    //Alogrithm for the calculation of the Lanczos Approx of the complete Gamma function 
    //as implemented in Numerical Recipes 3rd ed, 2007.
    real g = 4.7421875; //g parameter for the gamma function
    real x, tmp, y, A_g;
    
    //these are the cn's
    static const real coeff[14] = {57.1562356658629235,-59.5979603554754912,
                                14.1360979747417471,-0.491913816097620199,.339946499848118887e-4,
                                .465236289270485756e-4,-.983744753048795646e-4,.158088703224912494e-3,
                                -.210264441724104883e-3,.217439618115212643e-3,-.164318106536763890e-3,
                                .844182239838527433e-4,-.261908384015814087e-4,.368991826595316234e-5};
    y = x = z;
    tmp = x + g + 0.5;
    tmp = (x + 0.5) * mw_log(tmp) - tmp;
    A_g = 0.999999999999997092; //this is c0
    
    for (int j = 0; j < 14; j++) 
    {
        A_g += coeff[j] / ++y;
    } //calculates the series approx sum
        
    //sqrt(2 * pi) = 2.5066282746310005
    tmp += mw_log(2.5066282746310005 * A_g / x);//returns the log of the gamma function
    
    return tmp;
}

/*
 * Returns the incomplete gamma function P(a, x), evaluated by its
 * series representation. Also returns ln(Gamma(a)) through gln.
 */
static real gser(const real a, const real x, real* gln)
{
    const int itmax = 100;
    const real eps = (real) 3.0e-7;
    real sum, del, ap;

    *gln = gammln(a);
    if (x <= 0.0)
    {
        return 0.0;
    }

    ap = a;
    del = sum = 1.0 / a;
    for (int n = 1; n <= itmax; ++n)
    {
        ap += 1.0;
        del *= x / ap;
        sum += del;
        if (mw_fabs(del) < mw_fabs(sum) * eps)
        {
            return sum * mw_exp(-x + a * mw_log(x) - (*gln));
        }
    }

    mw_printf("WARNING: gser did not converge (a=%f, x=%f)\n", a, x);
    return sum * mw_exp(-x + a * mw_log(x) - (*gln));
}

/*
 * Returns the incomplete gamma function Q(a, x), evaluated by its
 * continued-fraction representation (modified Lentz method).
 * Also returns ln(Gamma(a)) through gln.
 */
static real gcf(const real a, const real x, real* gln)
{
    const int itmax = 100;
    const real eps = (real) 3.0e-7;
    const real fpmin = (real) 1.0e-30;
    real an, b, c, d, del, h;

    *gln = gammln(a);
    b = x + 1.0 - a;
    if (mw_fabs(b) < fpmin)
    {
        b = fpmin;
    }
    c = 1.0 / fpmin;
    d = 1.0 / b;
    h = d;

    for (int i = 1; i <= itmax; ++i)
    {
        an = -(real) i * ((real) i - a);
        b += 2.0;
        d = an * d + b;
        if (mw_fabs(d) < fpmin)
        {
            d = fpmin;
        }
        c = b + an / c;
        if (mw_fabs(c) < fpmin)
        {
            c = fpmin;
        }
        d = 1.0 / d;
        del = d * c;
        h *= del;
        if (mw_fabs(del - 1.0) < eps)
        {
            return mw_exp(-x + a * mw_log(x) - (*gln)) * h;
        }
    }

    mw_printf("WARNING: gcf did not converge (a=%f, x=%f)\n", a, x);
    return mw_exp(-x + a * mw_log(x) - (*gln)) * h;
}

/*
 * Returns the regularized lower incomplete gamma function:
 * P(a, x) = gamma(a, x) / Gamma(a).
 * Uses series for x < a + 1, otherwise returns 1 - Q(a, x).
 */
real gammp(const real a, const real x)
{
    real gln;

    if (x < 0.0 || a <= 0.0)
    {
        mw_printf("WARNING: Invalid arguments in gammp (a=%f, x=%f)\n", a, x);
        return NAN;
    }

    if (x < (a + 1.0))
    {
        return gser(a, x, &gln);
    }
    return 1.0 - gcf(a, x, &gln);
}

/*
 * Returns the regularized upper incomplete gamma function:
 * Q(a, x) = Gamma(a, x) / Gamma(a) = 1 - P(a, x).
 * Uses series for x < a + 1, otherwise continued fraction.
 */
real gammq(const real a, const real x)
{
    real gln;

    if (x < 0.0 || a <= 0.0)
    {
        mw_printf("WARNING: Invalid arguments in gammq (a=%f, x=%f)\n", a, x);
        return NAN;
    }

    if (x < (a + 1.0))
    {
        return 1.0 - gser(a, x, &gln);
    }
    return gcf(a, x, &gln);
}

/* Returns the complete gamma function Gamma(z). */
real GammaFunc(const real z)
{
    return mw_exp(gammln(z));
}

/* Returns the upper incomplete gamma function Gamma(a, x). */
real UpperIncompleteGammaFunc(real a, real x)
{
    return GammaFunc(a) * gammq(a, x);
}

/* Returns the lower incomplete gamma function Gamma(a, x). */
real LowerIncompleteGammaFunc(real a, real x)
{
    return GammaFunc(a) * gammp(a, x);
}

/*
 * Numerical Recipes erff(x):
 * returns erf(x) using the regularized incomplete gamma P(1/2, x^2).
 */
real ErrorFunc(real x)
{
    const real xsq = x * x;
    return (x < 0.0) ? -gammp(0.5, xsq) : gammp(0.5, xsq);
}

/*
 * Numerical Recipes erffc(x):
 * returns erfc(x) using P(1/2, x^2) for x < 0 and Q(1/2, x^2) for x >= 0.
 */
real ComplementaryErrorFunc(real x)
{
    const real xsq = x * x;
    return (x < 0.0) ? 1.0 + gammp(0.5, xsq) : gammq(0.5, xsq);
}

/*
 * Numerical Recipes erfcc(x):
 * fast erfc approximation with fractional error < 1.2e-7.
 */
real ComplementaryErrorFuncApprox(real x)
{
    real z = mw_fabs(x);
    real t = 1.0 / (1.0 + 0.5 * z);
    real ans = t * mw_exp(-z * z - 1.26551223
                          + t * (1.00002368
                          + t * (0.37409196
                          + t * (0.09678418
                          + t * (-0.18628806
                          + t * (0.27886807
                          + t * (-1.13520398
                          + t * (1.48851587
                          + t * (-0.82215223
                          + t * 0.17087277)))))))));
    return (x >= 0.0) ? ans : 2.0 - ans;
}

// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // 

real calc_vLOS(const mwvector v, const mwvector p, real sunGCdist)
{
    real xsol = X(p) + sunGCdist;
    real mag = mw_sqrt( xsol * xsol + Y(p) * Y(p) + Z(p) * Z(p) );
    real vl = xsol * X(v) + Y(p) * Y(v) + Z(p) * Z(v);
    vl = vl / mag;
    
    return vl;
}

real calc_distance(const mwvector p, real sunGCdist)  /**Calculating the distance to each body **/
{
    real xsol = X(p) + sunGCdist;
    real distance = mw_sqrt(xsol * xsol + Y(p) * Y(p) + Z(p) * Z(p) );

    return distance;
}

/* Get the dispersion in each bin*/
void nbCalcDisp(NBodyHistogram* histogram, mwbool initial, real correction_factor)
{
    unsigned int i;
    unsigned int j;
    unsigned int Histindex;
    
    unsigned int lambdaBins = histogram->lambdaBins;
    unsigned int betaBins = histogram->betaBins;
    
    HistData* histData = histogram->data;

    real count;
    real n_ratio;
    real n_new;
    real sum, sq_sum, dispsq;
    
    
    for (i = 0; i < lambdaBins; ++i)
    {
        for(j = 0; j < betaBins; ++j)
        {
            Histindex = i * betaBins + j;
            count = (real) histData[Histindex].rawCount;
            count -= histData[Histindex].outliersRemoved;
            
            if(count > 10.0)//need enough counts so that bins with minimal bodies do not throw the vel disp off
            {
                n_new = count - 1.0; //because the mean is calculated from the same populations set
                n_ratio = count / (n_new); 
                
                sq_sum = histData[Histindex].sq_sum;
                sum = histData[Histindex].sum;
                
                 dispsq = (sq_sum / n_new) - n_ratio * sqr(sum / count);
                
                /* The following requires explanation. For the first calculation of dispersions, the bool initial 
                 * needs to be set to true. After that false.
                 * It will correct if there was no outliers removed because then the distribution does not have wings
                 * It will also correct if outliers were removed because then the wings were removed. 
                 * Does one correction everytime there was an outlier removed. Corrects once if no outliers were removed. 
                 */
                
                if(!initial)
                {
                    dispsq *= correction_factor;
                }//correcting for truncating the distribution when removing outliers.

                histData[Histindex].variable = mw_sqrt(dispsq);
                histData[Histindex].err =  mw_sqrt( (count + 1) /(count * n_new ) ) * histData[Histindex].variable ;
                
            }
        }
    }
    
}

void nbRemoveOutliers(const NBodyState* st, NBodyHistogram* histogram, real * use_body, real * var, real sigma_cutoff, real sunGCdist __attribute__((unused)), int histBins)
{
    unsigned int Histindex;
    Body* p;
    HistData* histData;
    const Body* endp = st->bodytab + st->nbody;

    unsigned int counter = 0;
    
    histData = histogram->data;

    real bin_sigma, new_count, this_var;

    /*Calculate old average and reset counters and sums for each histogram bin*/
    real bin_ave[histBins];
    real temp_sum[histBins];
    real temp_sqr[histBins];
    real temp_removed[histBins];

    for (int indx1 = 0; indx1 < histBins; ++indx1)
    {
        new_count = (real) (histData[indx1].rawCount - histData[indx1].outliersRemoved);
        bin_ave[indx1] = histData[indx1].sum / new_count;
        temp_sum[indx1] = 0.0;
        temp_sqr[indx1] = 0.0;
        temp_removed[indx1] = 0.0;
        //mw_printf("Cleared Bin %d\n",indx1);
    }
    /*------------------------------------------------------------------------*/
    
    for (p = st->bodytab; p < endp; ++p)
    {
        /* Only include bodies in models we aren't ignoring */
        if (!ignoreBody(p))
        {
            
            /* Check if the position is within the bounds of the histogram */
            if (use_body[counter] >= 0)//if it's not -1 then it was in the hist and set to the Histindex   
            {   
                Histindex = (int) use_body[counter];
                //mw_printf("Histogram Index = %d\n",Histindex);
                
                this_var = var[counter];
                
                /* Use old standard deviation calculated before */
                bin_sigma = histData[Histindex].variable;
                
                if(mw_fabs(bin_ave[Histindex] - this_var) < sigma_cutoff * bin_sigma)//if it is inside of the sigma limit
                {
                    temp_sum[Histindex] += this_var;
                    temp_sqr[Histindex] += this_var*this_var;
                }
                else
                {
                    temp_removed[Histindex]+=1.0;//keep track of how many are being removed
                }

                
            }
            counter++;
        }
    }
    for (int indx2 = 0; indx2 < histBins; ++indx2)
    {
        histData[indx2].sum = temp_sum[indx2];
        histData[indx2].sq_sum = temp_sqr[indx2];
        histData[indx2].outliersRemoved = temp_removed[indx2];
        //mw_printf("Outliers Removed @ Index %d = %.15f\n",indx2,histData[indx2].outliersRemoved);
    }
}

/*Removes outliers for momentum calculation. Requires a different function as it is not stored in a hist*/
void nbRemoveMomentumOutliers(const NBodyState* st, NBodyHistogram* histogram, int* in_hist, real sigma_cutoff, real IterMax, real correction_factor, real nbody, real counts)
{
    Body* p = NULL;

    mwbool initial = TRUE;
    mwbool corrected = FALSE;

    mwvector L_sum = ZERO_VECTOR; /*angular momentum vector sum*/
    mwvector L_sum2 = ZERO_VECTOR; /*angular momentum vector sum squared*/
    mwvector L = ZERO_VECTOR; /*angular momentum vector per particle*/
    mwvector L2 = ZERO_VECTOR; /*angular momentum vector per particle squared*/
    mwvector L_var = ZERO_VECTOR; /*variance in angular momentum vector of simulation*/
    mwvector r = ZERO_VECTOR; /*position vector*/
    mwvector v = ZERO_VECTOR; /*velocity vector*/
    real mass = 1.0; //histogram->massPerParticle; /*mass of each particle*/

    /*Read in old average*/
    mwvector L_avg = histogram->params.L;
    mwvector LErr = histogram->params.LErr;

    for(unsigned int i = 0; i < IterMax; i++)
    {
        for (unsigned int j = 0; j < nbody; ++j)
        {
            p = &st->bodytab[j];
            /* Only include bodies in models we aren't ignoring */
            if (!ignoreBody(p))
            {

                /* Check if the position is within the bounds of the histogram */
                if (in_hist[j] > 0)//if it's not 0 then it was in the hist and set to the Histindex   
                {   
                    /*Find momentum of particle*/
                    r = Pos(p);
                    v = Vel(p);
                    L = mw_crossv(r, v);
                    L = mw_mulvs(L, mass);

                    /*If outside of cutoff in any component, remove the particle*/
                    if(fabs(X(L_avg) - X(L)) < sigma_cutoff * X(LErr) &&
                       fabs(Y(L_avg) - Y(L)) < sigma_cutoff * Y(LErr) &&
                       fabs(Z(L_avg) - Z(L)) < sigma_cutoff * Z(LErr))//if it is inside of the sigma limit
                    {
                        L_sum = mw_addv(L_sum, L);
                        L2 = mw_mulv(L, L);
                        L_sum2 = mw_addv(L_sum2, L2);
                    }
                    else
                    {
                        in_hist[j] = 0; //mark as being removed
                        counts -= 1.0;
                        corrected = TRUE;
                    }
                }
            }
        }
        L_avg = mw_mulvs(L_sum, 1.0/(real)counts);
        if(corrected)
        {
            L_var = mw_divvs(mw_subv(L_sum2, mw_mulvs(mw_mulv(L_avg, L_avg), counts)), (real)(counts-1)); /*variance in angular momentum vector of simulation*/
            LErr.x = mw_sqrt(X(L_var));
            LErr.y = mw_sqrt(Y(L_var));
            LErr.z = mw_sqrt(Z(L_var));
        }
        else
        {
            LErr= histogram->params.LErr; /*use previous error if no outliers were removed*/
        }

        if(corrected || initial) //only correct if wings are removed or if this is the first calculation and there were no wings
        {
            LErr = mw_mulvs(LErr, correction_factor);
        }

        histogram->params.L.x = L_avg.x;
        histogram->params.L.y = L_avg.y;
        histogram->params.L.z = L_avg.z;
        histogram->params.LErr.x = LErr.x;
        histogram->params.LErr.y = LErr.y;
        histogram->params.LErr.z = LErr.z;

        L_avg = histogram->params.L;
        LErr = histogram->params.LErr;

        corrected = FALSE;
        initial = FALSE;
        L_sum = (mwvector)ZERO_VECTOR;
        L_sum2 = (mwvector)ZERO_VECTOR;
    }
    return;
}   

real nbCostComponent(const NBodyHistogram* data, const NBodyHistogram* histogram)
{
    unsigned int lambdaBins = data->lambdaBins;
    unsigned int betaBins = data->betaBins;
    unsigned int nbins = lambdaBins * betaBins;
    real n = (real) histogram->totalSimulated;
    real nSim_uncut = (real) histogram->totalNum;   /* Total simulated before dropping bins */
    real nData = (real) data->totalNum;
    real nDataVariance = 0.0;
    real histMass = histogram->massPerParticle;
    real dataMass = data->massPerParticle;
    real p; /* probability of observing an event */
    real rawCount;
    real nSim = nSim_uncut;
    
    if (data->lambdaBins != histogram->lambdaBins || data->betaBins != histogram->betaBins)
    {
        return NAN;
    }

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
    if (nSim == 0.0 || nData == 0.0)
    {
#pragma GCC diagnostic pop
        /* If the histogram is totally empty, it is worse than the worst case */
        return INFINITY;
    }

    if (histMass <= 0.0 || dataMass <= 0.0)
    {
        /*In order to calculate likelihood the masses are necessary*/
        return NAN;
    }
    
    
    /*
     * Correcting for bins in the comparison histogram that are not 
     * included in the comparison. Also calculating data errors.
     */
    for (unsigned int i = 0; i < nbins; ++i)
    {
        if(!data->data[i].useBin)
        {
            rawCount = mw_round(histogram->data[i].variable * nSim_uncut);
            nSim -= rawCount;
        }
        /*WARNING: These are NOT the errors in the normalized counts, but rather the errors in the
          counts divided by the total number of bodies within the histogram. There IS a difference!*/
        else
        {
            nDataVariance += sqr(data->data[i].err*nData);
        }
    }
    
    /* this is the newest version of the cost function
     * it uses a combination of the binomial error for sim 
     * and the poisson error for the data
     */
    p = ( nSim / n) ;

    /*Print statements for debugging likelihood*/
//    mw_printf("dataMass = %.15f\n",dataMass);
//    mw_printf("nData    = %.15f\n",nData);
//    mw_printf("histMass = %.15f\n",histMass);
//    mw_printf("nSim     = %.15f\n",nSim);
//    mw_printf("p        = %.15f\n",p);
//    mw_printf("Sim_Mass = %.15f\n",histMass*nSim);

    real num = - sqr(dataMass * nData - histMass * nSim);
    real denom = 2.0 * (sqr(dataMass) * nDataVariance + sqr(histMass) * nSim * p * (1.0 - p));
    real CostComponent = num / denom; //this is the log of the cost component

    /* the cost component is negative. Returning a postive value */
    return -CostComponent;
    
}

/* for use with velocity dispersion, beta dispersion, average vlos,
average beta, and average distance likelihood component calculations */
real nbLikelihood(const NBodyHistogram* data, const NBodyHistogram* histogram, int avgBins)
{
    unsigned int lambdaBins = data->lambdaBins;
    unsigned int betaBins = data->betaBins;
    unsigned int nbins = lambdaBins * betaBins;
    real Nsigma_sq = 0.0;
    real Data;
    real Hist;
    real err_data, err_hist;
    real probability;
    for (unsigned int i = 0; i < nbins; ++i)
    {
        if (data->data[i].useBin)
        {
            if (avgBins>1) /*calculate average of bins used for beta dispersion calculation, can only average over an odd number of bins*/
            {
                if (avgBins%2 == 0)
                {
                    printf("\tAveraging over an even number of bins is not currently supported, setting number to one less (%2d) \n", avgBins-1);
                    avgBins -= 1;
                }  
                unsigned int n = avgBins; /*number of bins used in average*/
                real varSum = 0.0;
                real valSum = 0.0; 
                for (int k = 0; k < avgBins; ++k)
                {
                    int index = i - (avgBins-1)/2 + k;
                    if (index<0 || (unsigned)index>=nbins) /*do not try to use bins that are off the range of the histogram*/
                    {
                        n -= 1;
                    }
                    else if (histogram->data[index].err<=0) /*do not try to use bins that have no data*/
                    {
                        n -= 1;
                    }
                    
                    else
                    {
                        varSum += sqr(histogram->data[index].err);
                        valSum += histogram->data[index].variable;
                    }
                }
                err_hist = mw_sqrt(varSum)/n; /*variance of the average is the sum of the variances/n^2 */ 
                err_data = data->data[i].err;
                Hist = valSum/n; /*average value across bins used */ 
                Data = data->data[i].variable;
            }
            else
            {
                err_data = data->data[i].err;
                err_hist = histogram->data[i].err;
            }
            if(err_data > 0)
            {
                if (avgBins<=1) 
                {
                    Data = data->data[i].variable;
                    Hist = histogram->data[i].variable;
                }
                if(err_hist > 0)
                {
                    Nsigma_sq += sqr( Data - Hist ) / ( sqr(err_data) + sqr(err_hist) );
                }
                else
                {
                    Nsigma_sq += 25;    /*Adding 5 sigma*/
                }
            }
        }

    }
    probability = (Nsigma_sq) / 2.0; //should be negative, but we return the negative of it anyway
    
    return probability;
}

/*Find momentum vector and momentum error of simulation. 
Calculations are done with light matter particles that fall within histogram range,
and outlier rejection is included as with all other calculations.
 */
void nbCalcMomentum(const NBodyState* st, const NBodyCtx* ctx, NBodyHistogram* data, NBodyHistogram* histogram)
{
    real nbody = st->nbody;
    HistogramParams* hp = &histogram->params;
    Body* p = NULL;
    int* in_hist = mwCalloc(nbody, sizeof(int)); /*array to mark which bodies are in the histogram*/
    real lambda = 0.0;
    real beta = 0.0;
    mwvector lambdaBetaR; /*array to hold lambda, beta, and distance*/
    NBHistTrig histTrig;
    nbGetHistTrig(&histTrig, hp);
    real counter = 0; /*number of bodies used*/
    mwvector L_sum = ZERO_VECTOR; /*angular momentum vector sum*/
    mwvector L_sum2 = ZERO_VECTOR; /*angular momentum vector sum squared*/
    mwvector L = ZERO_VECTOR; /*angular momentum vector per particle*/
    mwvector L2 = ZERO_VECTOR; /*angular momentum vector per particle squared*/
    mwvector r = ZERO_VECTOR; /*position vector*/
    mwvector v = ZERO_VECTOR; /*velocity vector*/
    /* Mass is not currently used in momentum calculation. Keeping it here in case someone needs it later*/
    real mass = 1.0; //histogram->massPerParticle; /*mass of each particle*/

    if(histogram->params.nRange >= 2) // Make sure any values given through lua are used
    {
        data->params.nRange = histogram->params.nRange;
        for(unsigned int j = 0; j < histogram->params.nRange; j++)
        {
            data->params.EMDRange[j] = histogram->params.EMDRange[j];
        }
    }
    if(data->params.nRange % 2 != 0)
    {
        data->params.nRange -= 1;
    }

    for (unsigned int i = 0; i < nbody; i++) /*sum over particles to find average momentum*/
    {
        p = &st->bodytab[i];
        if (!ignoreBody(p))
        {
            /* Get the position in lbr coorinates */
            lambdaBetaR = nbXYZToLambdaBeta(&histTrig, Pos(p), ctx->sunGCDist);
            lambda = L(lambdaBetaR);
            beta = B(lambdaBetaR);
            if (data->params.nRange == 0) /*If no EMD range is given, use entire hist*/
            {
                if ((lambda >= histogram->params.lambdaStart) && (lambda < histogram->params.lambdaEnd) &&
                    (beta >= histogram->params.betaStart) && (beta < histogram->params.betaEnd))
                {
                    in_hist[i] = 1.0; //mark that this body is in the histogram
                    r = Pos(p);
                    v = Vel(p);
                    L = mw_crossv(r, v); /*angular momentum vector of particle without conisderation of mass*/
                    L2 = mw_mulv(L, L); /*angular momentum vector squared*/
                    L_sum = mw_addv(L_sum, L);
                    L_sum2 = mw_addv(L_sum2, L2);
                    counter++;
                }
                else
                {
                    in_hist[i] = 0.0; //mark that this body is not in the histogram
                }
            }
            else /*If EMD range is given, only count particles in that range*/
            {
                mwbool counted = FALSE;
                for(unsigned int k = 0; k < data->params.nRange; k = k + 2)
                {
                    if ((lambda >= data->params.EMDRange[k]) && (lambda < data->params.EMDRange[k+1]) &&
                        (beta >= histogram->params.betaStart) && (beta < histogram->params.betaEnd))
                    {
                        counted = TRUE;
                    }
                }
                if(counted)
                {
                    in_hist[i] = 1.0; //mark that this body is in the histogram
                    r = Pos(p);
                    v = Vel(p);
                    L = mw_crossv(r, v); /*angular momentum vector of particle without conisderation of mass*/
                    L2 = mw_mulv(L, L); /*angular momentum vector squared*/
                    L_sum = mw_addv(L_sum, L);
                    L_sum2 = mw_addv(L_sum2, L2);
                    counter++;
                }
                else
                {
                    in_hist[i] = 0.0; //mark that this body is not in the histogram
                }
            }
        }
    }

    L_sum = mw_mulvs(L_sum, mass); /*adjust momentum for mass*/
    mwvector L_avg = mw_mulvs(L_sum, 1.0/(real)counter); /*average angular momentum vector of simulation*/
    L_sum2 = mw_mulvs(L_sum2, sqr(mass)); /*adjust momentum squared for mass*/
    mwvector L_var = mw_divvs(mw_subv(L_sum2, mw_mulvs(mw_mulv(L_avg, L_avg), counter)), (real)(counter-1)); /*variance in angular momentum vector of simulation*/

    mwvector LErr = ZERO_VECTOR;
    LErr.x = mw_sqrt(X(L_var));
    LErr.y = mw_sqrt(Y(L_var));
    LErr.z = mw_sqrt(Z(L_var));

    histogram->params.L.x = L_avg.x;
    histogram->params.L.y = L_avg.y;
    histogram->params.L.z = L_avg.z;
    histogram->params.LErr.x = LErr.x;
    histogram->params.LErr.y = LErr.y;
    histogram->params.LErr.z = LErr.z;

    nbRemoveMomentumOutliers(st, histogram, in_hist, ctx->MomentumSigma, ctx->IterMax, ctx->MomentumCorrect, nbody, counter); /*Remove outliers now that we have a standard deviation*/
    free(in_hist);
    return;
}

/*Actual likelihood calculation*/
real nbMomentumLikelihood(const NBodyHistogram* data, const NBodyHistogram* histogram)
{
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wfloat-equal"
    if(data->params.LErr.x == 0 || data->params.LErr.y == 0 || data->params.LErr.z == 0 || 
       histogram->params.LErr.x == 0 || histogram->params.LErr.y == 0 || histogram->params.LErr.z == 0)
    {
        mw_printf("WARNING: A momentum input is zero, it may not have been read in\n");
    }
    #pragma GCC diagnostic pop
    /* The likelihood only considers errors from the input data. This is so messy, unrealistic outputs with high errors will not return 
    reasonable scores, as the entire purpose of the momentum likelihood is to try to avoid these results*/
    real x_comp = (X(data->params.L) - X(histogram->params.L)) / X(data->params.LErr); //mw_sqrt(sqr(X(data->params.LErr)) + sqr(X(histogram->params.LErr)));
    real y_comp = (Y(data->params.L) - Y(histogram->params.L)) / Y(data->params.LErr); //mw_sqrt(sqr(Y(data->params.LErr)) + sqr(Y(histogram->params.LErr)));
    real z_comp = (Z(data->params.L) - Z(histogram->params.L)) / Z(data->params.LErr); //mw_sqrt(sqr(Z(data->params.LErr)) + sqr(Z(histogram->params.LErr)));

    real likelihood = 0.5 * (sqr(x_comp) + sqr(y_comp) + sqr(z_comp));

    return likelihood;
}
