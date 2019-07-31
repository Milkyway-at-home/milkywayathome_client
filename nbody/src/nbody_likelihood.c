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
    real emd = NAN;
    real cost_component = NAN;
    real vel_disp = NAN;
    real beta_disp = NAN;
    real beta_component = NAN;
    real LOS_velocity_component = NAN;
    real distance_component = NAN;
    real likelihood = NAN;
    dat = nbReadHistogram(datHist);
    match = nbReadHistogram(matchHist);

    if (dat && match)
    {
        emd = nbMatchEMD(dat, match);
        cost_component = nbCostComponent(dat->histograms[0], match->histograms[0]);
        likelihood = emd + cost_component;
        
        if(use_betadisp)
        {
            beta_disp = nbLikelihood(dat->histograms[1], match->histograms[1]);
            likelihood += beta_disp;
        }
        if(use_veldisp)
        {
            vel_disp = nbLikelihood(dat->histograms[2], match->histograms[2]);
            likelihood += vel_disp;
        }
         if(use_betacomp)
        {
            beta_component = nbLikelihood(dat->histograms[3], match->histograms[3]);
            likelihood += beta_component;
        }
        if(use_vlos)
        {
            LOS_velocity_component = nbLikelihood(dat->histograms[4], match->histograms[4]);
            likelihood += LOS_velocity_component;
        }
        if(use_dist)
        {
            distance_component = nbLikelihood(dat->histograms[5], match->histograms[5]);
            likelihood += distance_component;
        }
        
    }

    free(dat);
    free(match);
    return likelihood;
}


/* Calculate the likelihood from the final state of the simulation */
real nbSystemLikelihood(const NBodyState* st,
                     const MainStruct* data,
                     const MainStruct* histogram,
                     NBodyLikelihoodMethod method)
{
    
    real geometry_component;
    real cost_component;
    real velocity_dispersion_component = NAN;
    real beta_dispersion_component = NAN;
    real beta_component = NAN;
    real LOS_velocity_component = NAN;
    real distance_component = NAN;
    real likelihood = NAN;
    
    if (data->histograms[0]->lambdaBins != histogram->histograms[0]->lambdaBins)
    {
        mw_printf("Number of bins does not match those in histogram file. "
                  "Expected %u, got %u\n",
                  histogram->histograms[0]->lambdaBins,
                  data->histograms[0]->lambdaBins);
        return NAN;
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
        if (histogram->histograms[0]->totalNum < 0.0001 * (real) st->nbody)
        {
            real worstEMD;

            mw_printf("Number of particles in bins is very small compared to total. "
                      "(%u << %u). Skipping distance calculation\n",
                      histogram->histograms[0]->totalNum,
                      st->nbody
                );
            worstEMD = nbWorstCaseEMD(histogram->histograms[0]);
            //return 2.0 * worstEMD;
            return worstEMD; //Changed.  See above comment.
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

    likelihood = geometry_component + cost_component;
    
    /* likelihood due to the vel dispersion per bin of the two hist */
    if(st->useBetaDisp)
    {
        beta_dispersion_component = nbLikelihood(data->histograms[1], histogram->histograms[1]);
        likelihood += beta_dispersion_component;
    }
    if(st->useVelDisp)
    {
        velocity_dispersion_component = nbLikelihood(data->histograms[2], histogram->histograms[2]);
        likelihood += velocity_dispersion_component;
    }
    if(st->useBetaComp)
    {
        beta_component = nbLikelihood(data->histograms[3], histogram->histograms[3]);
        likelihood += beta_component;
    }
    if(st->useVlos)
    {
        LOS_velocity_component = nbLikelihood(data->histograms[4], histogram->histograms[4]);
        likelihood += LOS_velocity_component;
    }
    if(st->useDist)
    {
        distance_component = nbLikelihood(data->histograms[5], histogram->histograms[5]);
        likelihood += distance_component;
    }
    return likelihood;   
} 
