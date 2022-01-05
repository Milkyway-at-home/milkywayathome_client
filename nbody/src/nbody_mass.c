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
static real_0 factorial(int n)
{
     int counter;
     real_0 result = 0.0;

     for (counter = n; counter >= 1; counter--)
       {
          result += mw_log_0((real_0) counter);
       }

     return result;
}


static real_0 choose(int n, int c)
{
    unsigned int i;
    real_0 result = 0.0;
    
    /* This for loop calulates log(n!/(n-c)!) */
    for (i = n - c + 1; i <= (unsigned int) n; ++i)
    {
        result += mw_log_0(i);
    }
    result -= factorial(c);
    return result;
}

real_0 probability_match(int n, real_0 ktmp, real_0 pobs)
{
    real_0 result = 0.0;

    /*
     * Previously, this function took in k as an int. Bad move.
     * This function was called twice, one of which sent a real valued k: (int) k1 and (real) k2
     * That real k2 was then converted to int. Could result in converted (int) k1 != (int) k2 when k1 = k2. 
     * Special result was poor likelihood for some histograms when check against themselves!
     * General results: unknown. But probably not good. (most likely caused different machines to report
     * different likelihood values).
     * 
     */
    int k = (int) mw_round_0(ktmp);    //patch. See above. 
    //The previous calculation does not return the right values.  Furthermore, we need a zeroed metric.                                                                                              
    result =  (real_0) choose(n, k);
    result += k * mw_log_0(pobs); 
    result += (n - k) * mw_log_0(1.0 - pobs);
    
    
    return mw_exp_0(result);
}
// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // 
// IMPLEMENTATION OF GAMMA FUNCTIONS. COMPLETE AND INCOMPLETE
real_0 GammaFunc(const real_0 z) 
{
    //Alogrithm for the calculation of the Lanczos Approx of the complete Gamma function 
    //as implemented in Numerical Recipes 3rd ed, 2007.
    real_0 g = 4.7421875; //g parameter for the gamma function
    real_0 x, tmp, y, A_g;
    
    //these are the cn's
    static const real_0 coeff[14] = {57.1562356658629235,-59.5979603554754912,
                                14.1360979747417471,-0.491913816097620199,.339946499848118887e-4,
                                .465236289270485756e-4,-.983744753048795646e-4,.158088703224912494e-3,
                                -.210264441724104883e-3,.217439618115212643e-3,-.164318106536763890e-3,
                                .844182239838527433e-4,-.261908384015814087e-4,.368991826595316234e-5};
    y = x = z;
    tmp = x + g + 0.5;
    tmp = (x + 0.5) * mw_log_0(tmp) - tmp;
    A_g = 0.999999999999997092; //this is c0
    
    for (int j = 0; j < 14; j++) 
    {
        A_g += coeff[j] / ++y;
    } //calculates the series approx sum
        
    //sqrt(2 * pi) = 2.5066282746310005
    tmp += mw_log_0(2.5066282746310005 * A_g / x);//returns the log of the gamma function
    
    return mw_exp_0(tmp);
}

static real_0 shift_factorial(real_0 z, real_0 n)
{
    int counter;
    real_0 result = 0.0;
    for(counter = n; counter >= 1; counter--)
    {
        result += mw_log_0(z + (real_0) counter);
    }
    return mw_exp_0(result);
    
    
}


static real_0 series_approx(real_0 a, real_0 x)
{

    real_0 sum, del, ap;
    real_0 pow_x = 1.0;
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
        if (mw_fabs_0(del) < mw_fabs_0(sum) * 1.0e-15) 
        {
            return sum * mw_exp_0(-x + a * mw_log_0(x));
        }
    }
    
}
                            
real_0 IncompleteGammaFunc(real_0 a, real_0 x)
{
    //the series approx returns gamma from 0 to X but we want from X to INF
    //Therefore, we subtract it from GammaFunc which is from 0 to INF
    //The continued frac approx is already from X to INF
    
//     static const real max_a = 100;
    real_0 gamma = GammaFunc(a);

    if (x == 0.0) return gamma;
    // Use the series representation. 
    return gamma - series_approx(a,x);
    

}    

// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // 

real calc_vLOS(const mwvector* v, const mwvector* p, real_0 sunGCdist)
{
    real tmp1, tmp2;
    real xsol = mw_add_s(&X(p), sunGCdist);

    tmp1 = mw_hypot(&xsol, &Y(p));
    real mag = mw_hypot(&tmp1, &Z(p));

    tmp1 = mw_mul(&xsol, &X(v));
    tmp2 = mw_mul(&Y(p), &Y(v));
    tmp1 = mw_add(&tmp1, &tmp2);
    tmp2 = mw_mul(&Z(p), &Z(v));
    real vl = mw_add(&tmp1, &tmp2);
    return mw_div(&vl, &mag);
}

real calc_distance(const mwvector* p, real_0 sunGCdist)  /**Calculating the distance to each body **/
{
    real xsol = mw_add_s(&X(p), sunGCdist);
    real tmp = mw_hypot(&xsol, &Y(p));
    return mw_hypot(&tmp, &Z(p));
}

/* Get the dispersion in each bin*/
void nbCalcDisp(NBodyHistogram* histogram, mwbool initial, real_0 correction_factor)
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
    real tmp1, tmp2;
    
    
    for (i = 0; i < lambdaBins; ++i)
    {
        for(j = 0; j < betaBins; ++j)
        {
            Histindex = i * betaBins + j;
            count = mw_sub(&histData[Histindex].rawCount, &histData[Histindex].outliersRemoved);
            
            if(showRealValue(&count) > 10.0)//need enough counts so that bins with minimal bodies do not throw the vel disp off
            {
                n_new = mw_add_s(&count, -1.0); //because the mean is calculated from the same populations set
                n_ratio = mw_div(&count, &n_new); 
                
                sq_sum = histData[Histindex].sq_sum;
                sum = histData[Histindex].sum;

                tmp1 = mw_div(&sq_sum, &n_new);
                tmp2 = mw_div(&sum, &count);
                tmp2 = sqr(&tmp2);
                tmp2 = mw_mul(&n_ratio, &tmp2);
                dispsq = mw_sub(&tmp1, &tmp2);
                
                /* The following requires explanation. For the first calculation of dispersions, the bool initial 
                 * needs to be set to true. After that false.
                 * It will correct if there was no outliers removed because then the distribution does not have wings
                 * It will also correct if outliers were removed because then the wings were removed. 
                 * Does one correction everytime there was an outlier removed. Corrects once if no outliers were removed. 
                 */
                
                if(!initial)
                {
                    dispsq = mw_mul_s(&dispsq, correction_factor);
                }//correcting for truncating the distribution when removing outliers.

                histData[Histindex].variable = mw_sqrt(&dispsq);
                tmp1 = mw_add_s(&count, 1.0);
                tmp2 = mw_mul(&count, &n_new);
                tmp1 = mw_div(&tmp1, &tmp2);
                tmp1 = mw_sqrt(&tmp1);
                histData[Histindex].err =  mw_mul(&tmp1, &histData[Histindex].variable) ;
                
            }
        }
    }
    
}

//FIXME:Need to determine how derivatives here propagate
void nbRemoveOutliers(const NBodyState* st, NBodyHistogram* histogram, real_0 * use_body, real * var, real_0 sigma_cutoff, int histBins)
{
    unsigned int Histindex;
    Body* p;
    HistData* histData;
    const Body* endp = st->bodytab + st->nbody;
    real tmp;

    unsigned int counter = 0;
    
    histData = histogram->data;

    real_0 bin_sigma;
    real new_count, this_var;

    /*Calculate old average and reset counters and sums for each histogram bin*/
    real bin_ave[histBins];
    real temp_sum[histBins];
    real temp_sqr[histBins];
    real temp_removed[histBins];

    for (unsigned int indx1 = 0; indx1 < histBins; ++indx1)
    {
        new_count = mw_sub(&histData[indx1].rawCount, &histData[indx1].outliersRemoved);
        bin_ave[indx1] = mw_div(&histData[indx1].sum, &new_count);
        temp_sum[indx1] = ZERO_REAL;
        temp_sqr[indx1] = ZERO_REAL;
        temp_removed[indx1] = ZERO_REAL;
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
                bin_sigma = showRealValue(&histData[Histindex].variable);
                
                if(mw_fabs_0(showRealValue(&bin_ave[Histindex]) - showRealValue(&this_var)) < sigma_cutoff * bin_sigma)//if it is inside of the sigma limit
                {
                    temp_sum[Histindex] = mw_add(&temp_sum[Histindex], &this_var);
                    tmp = sqr(&this_var);
                    temp_sqr[Histindex] = mw_add(&temp_sqr[Histindex], &tmp);
                }
                else
                {
                    temp_removed[Histindex] = mw_add_s(&temp_removed[Histindex], 1.0);//keep track of how many are being removed
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
    real_0 n = (real_0) histogram->totalSimulated;
    real nSim_uncut = (real) histogram->totalNum;   /* Total simulated before dropping bins */
    real nData = (real) data->totalNum;
    real histMass = histogram->massPerParticle;
    real dataMass = data->massPerParticle;
    real p; /* probability of observing an event */
    real rawCount;
    real nSim = nSim_uncut;
    real tmp1, tmp2, tmp3;
    
    if (data->lambdaBins != histogram->lambdaBins || data->betaBins != histogram->betaBins)
    {
        return mw_real_const(NAN);
    }

    if (showRealValue(&nSim) == 0 || showRealValue(&nData) == 0)
    {
        /* If the histogram is totally empty, it is worse than the worst case */
        return mw_real_const(INFINITY);
    }

    if (showRealValue(&histMass) <= 0.0 || showRealValue(&dataMass) <= 0.0)
    {
        /*In order to calculate likelihood the masses are necessary*/
        return mw_real_const(NAN);
    }
    
    
    /*
     * Correcting for bins in the comparison histogram that are not 
     * included in the comparison.
     */
    for (unsigned int i = 0; i < nbins; ++i)
    {
        if(!data->data[i].useBin)
        {
            tmp1 = mw_mul(&histogram->data[i].variable, &nSim_uncut);
            rawCount = mw_round(&tmp1);
            nSim = mw_sub(&nSim, &rawCount);
        }

    }
    
    /* this is the newest version of the cost function
     * it uses a combination of the binomial error for sim 
     * and the poisson error for the data
     */
    p = mw_mul_s(&nSim, inv_0(n));

    /*Print statements for debugging likelihood*/
//    mw_printf("dataMass = %.15f\n",dataMass);
//    mw_printf("nData    = %.15f\n",nData);
//    mw_printf("histMass = %.15f\n",histMass);
//    mw_printf("nSim     = %.15f\n",nSim);
//    mw_printf("p        = %.15f\n",p);
//    mw_printf("Sim_Mass = %.15f\n",histMass*nSim);

    tmp1 = mw_mul(&dataMass, &nData);
    tmp2 = mw_mul(&histMass, &nSim);
    tmp1 = mw_sub(&tmp1, &tmp2);
    real num = sqr(&tmp1);

    tmp1 = sqr(&dataMass);
    tmp1 = mw_mul(&tmp1, &nData);
    tmp2 = sqr(&histMass);
    tmp3 = sqr(&p);
    tmp3 = mw_sub(&p, &tmp3);
    tmp3 = mw_mul(&nSim, &tmp3);
    tmp2 = mw_mul(&tmp2, &tmp3);
    tmp1 = mw_add(&tmp1, &tmp2);
    real denom = mw_mul_s(&tmp1, 2.0);
    real CostComponent = mw_div(&num, &denom); //this is the log of the cost component

    /* the cost component is negative. Returning a postive value */
    return CostComponent;
    
}

/* for use with velocity dispersion, beta dispersion, average vlos,
average beta, and average distance likelihood component calculations */
real nbLikelihood(const NBodyHistogram* data, const NBodyHistogram* histogram)
{
    unsigned int lambdaBins = data->lambdaBins;
    unsigned int betaBins = data->betaBins;
    unsigned int nbins = lambdaBins * betaBins;
    real Nsigma_sq = ZERO_REAL;
    real Data;
    real Hist;
    real err_data, err_hist;
    real probability;
    real tmp1, tmp2, tmp3;
    for (unsigned int i = 0; i < nbins; ++i)
    {
        if (data->data[i].useBin)
        {
            err_data = data->data[i].err;
            err_hist = histogram->data[i].err;

            if(showRealValue(&err_data) > 0)
            {
                Data = data->data[i].variable;
                Hist = histogram->data[i].variable;

                if(showRealValue(&err_hist) > 0)
                {
                    tmp1 = mw_sub(&Data, &Hist);
                    tmp1 = sqr(&tmp1);
                    tmp2 = sqr(&err_data);
                    tmp3 = sqr(&err_hist);
                    tmp2 = mw_add(&tmp2, &tmp3);
                    tmp1 = mw_div(&tmp1, &tmp2);
                    Nsigma_sq = mw_add(&Nsigma_sq, &tmp1);
                }
                else
                {
                    Nsigma_sq = mw_add_s(&Nsigma_sq, 25.0);    /*Adding 5 sigma*/
                }
            }
        }

    }
        probability = mw_mul_s(&Nsigma_sq, 0.5); //should be negative, but we return the negative of it anyway
    
    return probability;
}
