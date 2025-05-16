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
real GammaFunc(const real z) 
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
    
    return mw_exp(tmp);
}

static real shift_factorial(real z, real n)
{
    int counter;
    real result = 0.0;
    for(counter = n; counter >= 1; counter--)
    {
        result += mw_log(z + (real) counter);
    }
    return mw_exp(result);
    
    
}


static real series_approx(real a, real x)
{

    real sum, del, ap;
    real pow_x = 1.0;
//     ap = a;
    ap = 0.0;
    del = sum = 1.0 / a;//starting: gammma(a) / gamma(a+1) = 1/a
    for (;;) 
    {
//         ++ap;
        ++ap;
//         del *= x / ap;

        pow_x *= x;
        del = pow_x / (a * shift_factorial(a, ap));
        
        sum += del;
        if (mw_fabs(del) < mw_fabs(sum) * 1.0e-15) 
        {
            return sum * exp(-x + a * log(x));
        }
    }
    
}
                            
real IncompleteGammaFunc(real a, real x)
{
    //the series approx returns gamma from 0 to X but we want from X to INF
    //Therefore, we subtract it from GammaFunc which is from 0 to INF
    //The continued frac approx is already from X to INF
    
//     static const real max_a = 100;
    real gamma = GammaFunc(a);

    if (x == 0.0) return gamma;
    // Use the series representation. 
    return gamma - series_approx(a,x);
    

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

void nbRemoveOutliers(const NBodyState* st, NBodyHistogram* histogram, real * use_body, real * var, real sigma_cutoff, real sunGCdist, int histBins)
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

    for (unsigned int indx1 = 0; indx1 < histBins; ++indx1)
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
    for (unsigned int indx2 = 0; indx2 < histBins; ++indx2)
    {
        histData[indx2].sum = temp_sum[indx2];
        histData[indx2].sq_sum = temp_sqr[indx2];
        histData[indx2].outliersRemoved = temp_removed[indx2];
        //mw_printf("Outliers Removed @ Index %d = %.15f\n",indx2,histData[indx2].outliersRemoved);
    }
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

    if (nSim == 0 || nData == 0)
    {
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
                for (unsigned int k = 0; k < avgBins; ++k)
                {
                    int index = i - (avgBins-1)/2 + k;
                    if (index<0 || index>=nbins) /*do not try to use bins that are off the range of the histogram*/
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
