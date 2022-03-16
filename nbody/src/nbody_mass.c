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
#include "milkyway_util.h"


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

mwvector LBtoY(real* lambda, real* beta)
{
    mwvector result;
    real tmp1, tmp2;

    tmp1 = d2r(beta);
    tmp1 = mw_cos(&tmp1);
    tmp2 = d2r(lambda);
    tmp2 = mw_cos(&tmp2);
    result.x = mw_mul(&tmp1, &tmp2);

    tmp1 = d2r(beta);
    tmp1 = mw_cos(&tmp1);
    tmp2 = d2r(lambda);
    tmp2 = mw_sin(&tmp2);
    result.y = mw_mul(&tmp1, &tmp2);

    tmp1 = d2r(beta);
    result.z = mw_sin(&tmp1);

    return result;
}

real ln_FB5_dist(mwvector* x, mwvector* gamma1, mwvector* gamma2, mwvector* gamma3, real* kappa, real* beta) //Returns natural log of FB5 since these quantities are really small
{
    if(2*mw_abs_0(showRealValue(beta)) > showRealValue(kappa))
    {
        mw_printf("WARNING!: INVALID KAPPA AND BETA DETECTED!\n");
        return ZERO_REAL;
    }
    real tmp1, tmp2, tmp3;

    tmp1 = sqr(kappa);
    tmp2 = sqr(beta);
    tmp2 = mw_mul_s(&tmp2, 4.0);
    tmp1 = mw_sub(&tmp1, &tmp2);
    tmp1 = mw_log(&tmp1);
    tmp1 = mw_mul_s(&tmp1, -0.5);
    tmp1 = mw_add_s(&tmp1, mw_log_0(2*M_PI));
    real ln_c = mw_add(kappa, &tmp1);

    tmp1 = mw_dotv(gamma1, x);
    tmp1 = mw_mul(kappa, &tmp1);
    tmp2 = mw_dotv(gamma2, x);
    tmp3 = mw_dotv(gamma3, x);
    tmp3 = sqr(&tmp3);
    tmp2 = sqr(&tmp2);
    tmp2 = mw_sub(&tmp2, &tmp3);
    tmp2 = mw_mul(beta, &tmp2);
    real expon = mw_add(&tmp1, &tmp2);

    real result = mw_sub(&expon, &ln_c);
    //mw_printf("FB5 = %.15f\n", showRealValue(&result));

   return result;
}

real add_logspace(real* r1, real* r2)
{
    real big, small, tmp, result;
    if (showRealValue(r1) < showRealValue(r2))
    {
        big = *r2;
        small = *r1;
    }
    else
    {
        big = *r1;
        small = *r2;
    }

    tmp = mw_sub(&small, &big);
    tmp = mw_exp(&tmp);
    tmp = mw_log1p(&tmp);
    result = mw_add(&big, &tmp);
    return result;
}

real getBodyBinFrac(const NBodyCtx* ctx, const HistogramParams* hp, const Body* p, int HistIndex)
{
    real tmp1, tmp2;
    mwvector tmp;
    NBHistTrig histTrig;
    nbGetHistTrig(&histTrig, hp, ctx->leftHanded);

    //mw_printf("Body POS = [%.15f, %.15f, %.15f]\n", showRealValue(&X(&Pos(p))), showRealValue(&Y(&Pos(p))), showRealValue(&Z(&Pos(p))));

    /* Calculate vlos and distance */
    real vlos = calc_vLOS(&Vel(p), &Pos(p), ctx->sunGCDist);
    real dist = calc_distance(&Pos(p), ctx->sunGCDist);

    /* Calculate gamma1 (direction from Sun to body) */
    mwvector bodyLBR = nbXYZToLambdaBeta(&histTrig, &Pos(p), ctx->sunGCDist);
    //mw_printf("Body L,B = %.15f, %.15f\n", showRealValue(&bodyLBR.x), showRealValue(&bodyLBR.y));
    mwvector gamma1 = LBtoY(&bodyLBR.x, &bodyLBR.y);

    mwvector y; //Cartesian coords of transformed LB
    y.x = mw_mul(&gamma1.x, &bodyLBR.z);
    y.y = mw_mul(&gamma1.y, &bodyLBR.z);
    y.z = mw_mul(&gamma1.z, &bodyLBR.z);
    //mw_printf("Transformed Body POS = [%.15f, %.15f, %.15f]\n", showRealValue(&y.x), showRealValue(&y.y), showRealValue(&y.z));

    /* Rotate body velocity into Lambda-Beta system by rotating two positions */
    mwvector shiftedPos = mw_addv(&Pos(p), &Vel(p));
    mwvector shiftedLBR = nbXYZToLambdaBeta(&histTrig, &shiftedPos, ctx->sunGCDist);

    tmp = LBtoY(&shiftedLBR.x, &shiftedLBR.y);
    tmp.x = mw_mul(&tmp.x, &shiftedLBR.z);
    tmp.y = mw_mul(&tmp.y, &shiftedLBR.z);
    tmp.z = mw_mul(&tmp.z, &shiftedLBR.z);

    mwvector rotVel = mw_subv(&tmp, &y);

    /* Calculate proper motion to get gamma2 */
    mwvector gamma2;
    mwvector vlos_vec;
    vlos_vec.x = mw_mul(&gamma1.x, &vlos);
    vlos_vec.y = mw_mul(&gamma1.y, &vlos);
    vlos_vec.z = mw_mul(&gamma1.z, &vlos);

    mwvector properMotion_vec = mw_subv(&rotVel, &vlos_vec);
    real properMotion = mw_length(&properMotion_vec);

    gamma2.x = mw_div(&properMotion_vec.x, &properMotion);
    gamma2.y = mw_div(&properMotion_vec.y, &properMotion);
    gamma2.z = mw_div(&properMotion_vec.z, &properMotion);

    properMotion = mw_div(&properMotion, &dist); //radians per Gyr
    //mw_printf("properMotion = %.15f\n", showRealValue(&properMotion));

    /* Calculate gamma3 by taking the cross product of gamma1 and gamma2 (parity doesn't matter) */
    mwvector gamma3 = mw_crossv(&gamma1, &gamma2);

    //mw_printf("gamma1 = [%.15f, %.15f, %.15f]\n", showRealValue(&gamma1.x), showRealValue(&gamma1.y), showRealValue(&gamma1.z));
    //mw_printf("gamma2 = [%.15f, %.15f, %.15f]\n", showRealValue(&gamma2.x), showRealValue(&gamma2.y), showRealValue(&gamma2.z));
    //mw_printf("gamma3 = [%.15f, %.15f, %.15f]\n", showRealValue(&gamma3.x), showRealValue(&gamma3.y), showRealValue(&gamma3.z));

    /* Calculate smei-major and semi-minor axes (a and b respectively) */
    real b = inv(&dist);
    b = mw_mul_s(&b, mw_sqrt_0(ctx->eps2));

    real a = mw_mul_s(&properMotion, ctx->timestep);
    a = mw_hypot(&a, &b);

    //mw_printf("a, b = %.15f, %.15f\n", showRealValue(&a), showRealValue(&b));

    /* Get kappa and beta estimators for FB5 from a and b */
    tmp1 = mw_tan(&a);
    tmp2 = mw_tan(&b);
    tmp1 = mw_div(&tmp1, &tmp2);
    tmp1 = mw_add_s(&tmp1, -1.0);
    tmp2 = mw_sin(&b);
    tmp2 = sqr(&tmp2);
    tmp1 = mw_mul(&tmp2, &tmp1);
    tmp1 = mw_mul_s(&tmp1, -0.5);
    tmp1 = mw_add_s(&tmp1, 1.0);
    tmp2 = mw_cos(&b);
    tmp1 = mw_mul(&tmp1, &tmp2);
    tmp1 = mw_neg(&tmp1);
    tmp1 = mw_add_s(&tmp1, 1.0);
    real solidAngle = mw_mul_s(&tmp1, 2.0*M_PI);

    tmp1 = mw_sin(&a);
    tmp2 = mw_sin(&b);
    tmp1 = mw_mul(&tmp1, &tmp2);
    tmp1 = mw_div(&tmp1, &solidAngle);
    real r1 = mw_mul_s(&tmp1, M_PI);

    tmp1 = mw_sin(&a);
    tmp1 = sqr(&tmp1);
    tmp2 = mw_sin(&b);
    tmp2 = sqr(&tmp2);
    tmp1 = mw_sub(&tmp1, &tmp2);
    tmp1 = mw_mul(&r1, &tmp1);
    real r2 = mw_mul_s(&tmp1, 0.25);

    tmp1 = mw_mul_s(&r1, -2.0);
    tmp1 = mw_add_s(&tmp1, 2.0);
    tmp1 = mw_sub(&tmp1, &r2);
    real part1 = inv(&tmp1);

    tmp1 = mw_mul_s(&r1, -2.0);
    tmp1 = mw_add_s(&tmp1, 2.0);
    tmp1 = mw_add(&tmp1, &r2);
    real part2 = inv(&tmp1);

    real kappa = mw_add(&part1, &part2);
    real beta = mw_sub(&part1, &part2);
    beta = mw_mul_s(&beta, 0.5);

    //mw_printf("kappa = %.15f\n", showRealValue(&kappa));
    //mw_printf("beta  = %.15f\n", showRealValue(&beta));

    /* Get the bounds of integration (Lambda & Beta bin edges) */
    unsigned int BetaIndex = HistIndex % (hp->betaBins);
    unsigned int LambdaIndex = (HistIndex - BetaIndex)/(hp->betaBins);
    real_0 lambdaSize = (hp->lambdaEnd - hp->lambdaStart) / (real_0) hp->lambdaBins;
    real_0 betaSize = (hp->betaEnd - hp->betaStart) / (real_0) hp->betaBins;

    real_0 LambdaStart = hp->lambdaStart + lambdaSize*(LambdaIndex);
    real_0 LambdaEnd   = hp->lambdaStart + lambdaSize*(LambdaIndex+1);
    real_0 BetaStart   = hp->betaStart + betaSize*(BetaIndex);
    real_0 BetaEnd     = hp->betaStart + betaSize*(BetaIndex+1);

    /* Integrate distribution over full bin */
    unsigned int intLambdaBins = mw_ceil_0(3.0 * d2r_0(lambdaSize) / (1.0*showRealValue(&b)));
    unsigned int intBetaBins = mw_ceil_0(3.0 * d2r_0(betaSize) / (1.0*showRealValue(&b)));
    //mw_printf("L_bins, B_bins = %u, %u\n",intLambdaBins, intBetaBins);
    real B_val, L_val, B_val_rad, val;
    mwvector x;
    real integral;
    /* Since the quantities of this integral are so small, we are required
       to perform this integral calculation in log space */
    for(int i = 0; i < intBetaBins; i++) //Integrate over Beta
    {
        B_val = mw_real_const(BetaStart + (i + 0.5)*betaSize/(1.0*intBetaBins));
        tmp1 = d2r(&B_val);
        tmp1 = mw_cos(&tmp1);
        tmp1 = mw_log(&tmp1);
        for(int j = 0; j < intLambdaBins; j++) //Integrate over Lambda
        {
            L_val = mw_real_const(LambdaStart + (j + 0.5)*lambdaSize/(1.0*intLambdaBins));
            x = LBtoY(&L_val, &B_val);
            val = ln_FB5_dist(&x, &gamma1, &gamma2, &gamma3, &kappa, &beta);
            val = mw_add(&val, &tmp1);
            val = mw_add_s(&val, mw_log_0(d2r_0(betaSize)*d2r_0(lambdaSize)/(1.0*intBetaBins*intLambdaBins)));
            //mw_printf("val = %.15f\n", showRealValue(&val));
            if(i==0)
            {
                integral = val;
            }
            else
            {
                integral = add_logspace(&integral, &val);
            }
        }
    }
    //mw_printf("integral = %.15f\n", showRealValue(&integral));
    return integral;
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

void nbRemoveOutliers(const NBodyState* st, NBodyHistogram* histogram, real_0 * use_body, real * var, real_0 sigma_cutoff, int histBins, mwbool contBins, real * bodyFraction)
{
    unsigned int Histindex;
    Body* p;
    HistData* histData;
    const Body* endp = st->bodytab + st->nbody;
    real tmp;

    unsigned int counter = 0;
    
    histData = histogram->data;

    real_0 bin_sigma;
    real new_count, this_var, weight;

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
            if(contBins)
            {
                this_var = var[counter];
                for(int i = 0; i <histBins; i++)
                {
                    weight = bodyFraction[counter*histBins + i];
                    /* Use old standard deviation calculated before */
                    bin_sigma = showRealValue(&histData[i].variable);
                
                    if(mw_fabs_0(showRealValue(&bin_ave[i]) - showRealValue(&this_var)) < sigma_cutoff * bin_sigma)//if it is inside of the sigma limit
                    {
                        tmp = mw_mul(&this_var, &weight);
                        temp_sum[i] = mw_add(&temp_sum[i], &tmp);
                        tmp = sqr(&this_var);
                        tmp = mw_mul(&tmp, &weight);
                        temp_sqr[i] = mw_add(&temp_sqr[i], &tmp);
                    }
                    else
                    {
                        temp_removed[i] = mw_add(&temp_removed[i], &weight);//keep track of how many are being removed
                    }
                }
            }
            /* Check if the position is within the bounds of the histogram */
            else if (use_body[counter] >= 0)//if it's not -1 then it was in the hist and set to the Histindex   
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
        mw_printf("DIFFERENT NUMBER OF BINS DETECTED!\n");
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
        mw_printf("NEGATIVE MASS DETECTED!\n");
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
    //mw_printf("dataMass = %.15f\n", showRealValue(&dataMass));
    //mw_printf("nData    = %.15f\n", showRealValue(&nData));
    //mw_printf("histMass = %.15f\n", showRealValue(&histMass));
    //mw_printf("nSim     = %.15f\n", showRealValue(&nSim));
    //mw_printf("p        = %.15f\n", showRealValue(&p));
    //mw_printf("Sim_Mass = %.15f\n", showRealValue(&histMass)*showRealValue(&nSim));

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
    //mw_printf("num = %.15f\n", showRealValue(&num));
    //mw_printf("denom = %.15f\n", showRealValue(&denom));
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
