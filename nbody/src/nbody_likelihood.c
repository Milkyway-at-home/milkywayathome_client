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

real nbMatchHistogramFiles(const char* datHist, const char* matchHist, mwbool use_veldisp)
{
    NBodyHistogram* dat;
    NBodyHistogram* match;
    real emd = NAN;
    real cost_component = NAN;
    real vel_disp = NAN;
    real likelihood = NAN;
    dat = nbReadHistogram(datHist);
    match = nbReadHistogram(matchHist);

    if (dat && match)
    {
        emd = nbMatchEMD(dat, match);
        cost_component = nbCostComponent(dat, match);
        likelihood = emd + cost_component;
        
        if(use_veldisp)
        {
            vel_disp = nbVelocityDispersion(dat, match);
            likelihood += vel_disp;
        }
        
    }
    free(dat);
    free(match);
    return likelihood;
}


/* Calculate the likelihood from the final state of the simulation */
real nbSystemLikelihood(const NBodyState* st,
                     const NBodyHistogram* data,
                     const NBodyHistogram* histogram,
                     NBodyLikelihoodMethod method)
{
    
    real geometry_component;
    real cost_component;
    real velocity_dispersion_component;
    real likelihood = NAN;
    
    if (data->lambdaBins != histogram->lambdaBins)
    {
        mw_printf("Number of bins does not match those in histogram file. "
                  "Expected %u, got %u\n",
                  histogram->lambdaBins,
                  data->lambdaBins);
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
        if (histogram->totalNum < 0.0001 * (real) st->nbody)
        {
            real worstEMD;

            mw_printf("Number of particles in bins is very small compared to total. "
                      "(%u << %u). Skipping distance calculation\n",
                      histogram->totalNum,
                      st->nbody
                );
            worstEMD = nbWorstCaseEMD(histogram);
            //return 2.0 * worstEMD;
            return worstEMD; //Changed.  See above comment.
        }

        geometry_component = nbMatchEMD(data, histogram);
    }
    else
    {
        geometry_component = nbCalcChisq(data, histogram, method);
    }
    
    /* likelihood due to the amount of mass in the histograms */
    
    cost_component = nbCostComponent(data, histogram);

    likelihood = geometry_component + cost_component;
    
    /* likelihood due to the vel dispersion per bin of the two hist */
    if(st->useVelDisp)
    {
        velocity_dispersion_component = nbVelocityDispersion(data, histogram);
        likelihood += velocity_dispersion_component;
    }
    return likelihood;
    
}


