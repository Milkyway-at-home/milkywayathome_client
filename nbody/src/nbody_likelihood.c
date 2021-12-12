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
#include "milkyway_util.h"
#include "nbody_show.h"
#include "nbody_lua.h"


/*
  Load information necessary to calculate the likelihood from input lua script.

  Load parameters where in the sky the histogram goes and likelihood type.

  Return TRUE on failure.
*/
int nbGetLikelihoodInfo(const NBodyFlags* nbf, HistogramParams* hp, NBodyLikelihoodMethod* method)
{
    lua_State* luaSt = NULL;

    luaSt = nbOpenLuaStateWithScript(nbf, NULL);
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

real nbMatchHistogramFiles(const char* datHist, const char* matchHist, mwbool use_veldisp, mwbool use_betadisp, mwbool use_betacomp, mwbool use_vlos, mwbool use_dist)
{
    MainStruct* dat;
    MainStruct* match;
    real emd = mw_real_const(NAN);
    real cost_component = mw_real_const(NAN);
    real vel_disp = mw_real_const(NAN);
    real beta_disp = mw_real_const(NAN);
    real beta_component = mw_real_const(NAN);
    real LOS_velocity_component = mw_real_const(NAN);
    real distance_component = mw_real_const(NAN);
    real likelihood = mw_real_const(NAN);
    dat = nbReadHistogram(datHist);
    match = nbReadHistogram(matchHist);

    if (dat && match)
    {
        emd = nbMatchEMD(dat, match);
        cost_component = nbCostComponent(dat->histograms[0], match->histograms[0]);
        likelihood = mw_add(emd, cost_component);
        
        if(use_betadisp)
        {
            beta_disp = nbLikelihood(dat->histograms[1], match->histograms[1]);
            likelihood = mw_add(likelihood, beta_disp);
        }
        if(use_veldisp)
        {
            vel_disp = nbLikelihood(dat->histograms[2], match->histograms[2]);
            likelihood = mw_add(likelihood, vel_disp);
        }
        if(use_betacomp)
        {
            if(!dat->usage[4] || !match->usage[4])
            {
                mw_printf("One of these files does not contain any info for average beta\n");
                return mw_real_const(NAN);
            }
            beta_component = nbLikelihood(dat->histograms[3], match->histograms[3]);
            likelihood = mw_add(likelihood, beta_component);
        }
        if(use_vlos)
        {
            if(!dat->usage[3] || !match->usage[3])
            {
                mw_printf("One of these files does not contain any info for average vlos\n");
                return mw_real_const(NAN);
            }
            LOS_velocity_component = nbLikelihood(dat->histograms[4], match->histograms[4]);
            likelihood = mw_add(likelihood, LOS_velocity_component);
        }
        if(use_dist)
        {
            if(!dat->usage[5] || !match->usage[5])
            {
                mw_printf("One of these files does not contain any info for average distance\n");
                return mw_real_const(NAN);
            }
            distance_component = nbLikelihood(dat->histograms[5], match->histograms[5]);
            likelihood = mw_add(likelihood, distance_component);
        }
        
    }

    free(dat);
    free(match);
    return likelihood;
}


/* Calculate the likelihood from the final state of the simulation */
real * nbSystemLikelihood(const NBodyState* st,
                     const MainStruct* data,
                     const MainStruct* histogram,
                     NBodyLikelihoodMethod method)
{
    
    real geometry_component;
    real cost_component;
    real velocity_dispersion_component = mw_real_const(NAN);
    real beta_dispersion_component = mw_real_const(NAN);
    real beta_component = mw_real_const(NAN);
    real LOS_velocity_component = mw_real_const(NAN);
    real distance_component = mw_real_const(NAN);
    real likelihood = mw_real_const(NAN);

    static real likelihoodArray[8];

    static real NANArray[8]; /*This array contains NANs for each element so that it may be properly passed from this function*/
    NANArray[0] = mw_real_const(NAN);
    NANArray[1] = mw_real_const(NAN);
    NANArray[2] = mw_real_const(NAN);
    NANArray[3] = mw_real_const(NAN);
    NANArray[4] = mw_real_const(NAN);
    NANArray[5] = mw_real_const(NAN);
    NANArray[6] = mw_real_const(NAN);
    NANArray[7] = mw_real_const(NAN);
    
    if (data->histograms[0]->lambdaBins != histogram->histograms[0]->lambdaBins)
    {
        mw_printf("Number of bins does not match those in histogram file. "
                  "Expected %u, got %u\n",
                  histogram->histograms[0]->lambdaBins,
                  data->histograms[0]->lambdaBins);
        return NANArray;
    }

    
    /* likelihood due to shape of the histograms */
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
        if (showRealValue(histogram->histograms[0]->totalNum) < 0.0001 * (real_0) st->nbody)
        {
            static real worstEMD_Array[8];
            real_0 worstEMD;

            mw_printf("Number of particles in bins is very small compared to total. "
                      "(%u << %u). Skipping distance calculation\n",
                      histogram->histograms[0]->totalNum,
                      st->nbody
                );
            worstEMD = nbWorstCaseEMD(histogram->histograms[0]);

            worstEMD_Array[0] = mw_real_const(worstEMD);
            worstEMD_Array[1] = mw_real_const(worstEMD);
            worstEMD_Array[2] = ZERO_REAL;
            worstEMD_Array[3] = ZERO_REAL;
            worstEMD_Array[4] = ZERO_REAL;
            worstEMD_Array[5] = ZERO_REAL;
            worstEMD_Array[6] = ZERO_REAL;
            worstEMD_Array[7] = ZERO_REAL;
            //return 2.0 * worstEMD;
            return worstEMD_Array; //Changed.  See above comment.
        }
        // this function has been changed to accept MainStruct
        geometry_component = nbMatchEMD(data, histogram);
        //mw_printf("EMD Calculated!\n");
    }
    else
    {
        geometry_component = nbCalcChisq(data->histograms[0], histogram->histograms[0], method);
    }
    
    /* likelihood due to the amount of mass in the histograms */
    
    cost_component = nbCostComponent(data->histograms[0], histogram->histograms[0]);

    likelihood = mw_add(geometry_component, cost_component);
    
    /* likelihood due to the vel dispersion per bin of the two hist */
    if(st->useBetaDisp)
    {
        beta_dispersion_component = nbLikelihood(data->histograms[1], histogram->histograms[1]);
        likelihood = mw_add(likelihood, beta_dispersion_component);
    }
    if(st->useVelDisp)
    {
        velocity_dispersion_component = nbLikelihood(data->histograms[2], histogram->histograms[2]);
        likelihood = mw_add(likelihood, velocity_dispersion_component);
    }
    if(st->useVlos)
    {
        if(!data->usage[3] || !histogram->usage[3])
        {
            mw_printf("One of these files does not contain any info for average velocity\n");
            return NANArray;
        }  
        LOS_velocity_component = nbLikelihood(data->histograms[3], histogram->histograms[3]);
        likelihood = mw_add(likelihood, LOS_velocity_component);
    }
    if(st->useBetaComp)
    {
        if(!data->usage[4] || !histogram->usage[4])
        {
            mw_printf("One of these files does not contain any info for average beta\n");
            return NANArray;
        }  
        beta_component = nbLikelihood(data->histograms[4], histogram->histograms[4]);
        likelihood = mw_add(likelihood, beta_component);
    }
    if(st->useDist)
    {
        if(!data->usage[5] || !histogram->usage[5])
        {
            mw_printf("One of these files does not contain any info for average distance\n");
            return NANArray;
        }  
        distance_component = nbLikelihood(data->histograms[5], histogram->histograms[5]);
        likelihood = mw_add(likelihood, distance_component);
    }

    likelihoodArray[0]=likelihood;
    likelihoodArray[1]=geometry_component;
    likelihoodArray[2]=cost_component;
    likelihoodArray[3]=beta_dispersion_component;
    likelihoodArray[4]=velocity_dispersion_component;
    likelihoodArray[5]=beta_component;
    likelihoodArray[6]=LOS_velocity_component;
    likelihoodArray[7]=distance_component;

    return likelihoodArray;   
} 
