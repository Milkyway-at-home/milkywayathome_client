/*
 * Copyright (c) 2010-2011 Matthew Arsenault
 * Copyright (c) 2010-2011 Rensselaer Polytechnic Institute
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
#include "nbody.h"
#include "nbody_plain.h"
#include "nbody_shmem.h"
#include "nbody_curses.h"
#include "nbody_defaults.h"
#include "nbody_util.h"
#include "nbody_checkpoint.h"
#include "nbody_mass.h"
#include "nbody_grav.h"
#include "nbody_histogram.h"
#include "nbody_likelihood.h"
#include "nbody_devoptions.h"
#include "nbody_orbit_integrator.h"
#include "nbody_potential.h"
#include "nbody_friction.h"

#ifdef NBODY_BLENDER_OUTPUT
  #include "blender_visualizer.h"
#endif

static void nbReportProgress(const NBodyCtx* ctx, NBodyState* st)
{
    real frac = (real) st->step / (real) ctx->nStep;

    mw_fraction_done(frac);

    if (st->reportProgress)
    {
        mw_mvprintw(0, 0,
                    "Running: %f / %f (%f%%)\n",
                    frac * ctx->timeEvolve,
                    ctx->timeEvolve,
                    100.0 * frac
            );

        mw_refresh();
    }
}

static NBodyStatus nbCheckpoint(const NBodyCtx* ctx, NBodyState* st)
{
    if (nbTimeToCheckpoint(ctx, st))
    {
        if (nbWriteCheckpoint(ctx, st))
        {
            return NBODY_CHECKPOINT_ERROR;
        }

        mw_checkpoint_completed();
    }

    return NBODY_SUCCESS;
}

static inline int get_likelihood(const NBodyCtx* ctx, NBodyState* st, const NBodyFlags* nbf)
{
    MainStruct* data = NULL;
    MainStruct* histogram = NULL;
    real likelihood = NAN;
    real likelihood_EMD = NAN;
    real likelihood_Mass = NAN;
    real likelihood_Beta = NAN;
    real likelihood_Vel = NAN;
    real likelihood_BetaAvg = NAN;
    real likelihood_VelAvg = NAN;
    real likelihood_Dist = NAN;

    real *likelihoodArray;

    NBodyLikelihoodMethod method;
    HistogramParams hp;
    
    mwbool calculateLikelihood = (nbf->histogramFileName != NULL);
    
    if (calculateLikelihood || nbf->histoutFileName || nbf->printHistogram)
    {
        if (nbGetLikelihoodInfo(nbf, &hp, &method) || method == NBODY_INVALID_METHOD)
        {
            /* this would normally return a print statement 
             * but I do not want to overload the output since 
             * this would run every time step.
             */
            return 0;
        }
    }
    
    if (calculateLikelihood)
    {
        
        histogram = nbCreateHistogram(ctx, st, &hp);
 
        if (!histogram)
        {
            /* this would normally return a print statement 
             * but I do not want to overload the output since 
             * this would run every time step.
             */
            return 0;
        }
        
        data = nbReadHistogram(nbf->histogramFileName);
        
        if (!data)
        {
            free(histogram);
            /* if the input histogram does not exist, I do not want the 
             * simulation to terminate as you can still get the output file
             * from it. Therefore, this function will end here but with 0
             */
            return 0;
        }
        likelihoodArray = nbSystemLikelihood(st, data, histogram, method);
        likelihood         = likelihoodArray[0];
        likelihood_EMD     = likelihoodArray[1];
        likelihood_Mass    = likelihoodArray[2];
        likelihood_Beta    = likelihoodArray[3];
        likelihood_Vel     = likelihoodArray[4];
        likelihood_BetaAvg = likelihoodArray[5];
        likelihood_VelAvg  = likelihoodArray[6];
        likelihood_Dist    = likelihoodArray[7];

        /*
          Used to fix Windows platform issues.  Windows' infinity is expressed as:
          1.#INF00000, -1.#INF00000, or 0.#INF000000.  The server reads these as -1, 1, and 0
          respectively, accounting for the sign change.  Thus, I have changed overflow
          infinities (not errors) to be the worst case.  The worst case is now the actual
          worst thing that can happen.

        * It previous returned the worse case when the likelihood == 0. 
        * Changed it to be best case, 1e-9 which has been added in nbody_defaults.h
        */
        if (likelihood > DEFAULT_WORST_CASE || likelihood < (-1 * DEFAULT_WORST_CASE) || isnan(likelihood))
        {
            likelihood = DEFAULT_WORST_CASE;
        }
        else if(likelihood == 0.0)
        {
            likelihood = DEFAULT_BEST_CASE;
        }

        /* this checks to see if the likelihood is an improvement */
        if(mw_fabs(likelihood) < mw_fabs(st->bestLikelihood))
        {
            st->bestLikelihood = likelihood;

            st->bestLikelihood_EMD = likelihood_EMD;

            st->bestLikelihood_Mass = likelihood_Mass;

            if (st->useBetaDisp)
            {
                st->bestLikelihood_Beta = likelihood_Beta;
            }
            else st->bestLikelihood_Beta = 0.0;

            if (st->useVelDisp)
            {
                st->bestLikelihood_Vel = likelihood_Vel;
            }
            else st->bestLikelihood_Vel = 0.0;

            if (st->useBetaComp)
            {
                st->bestLikelihood_BetaAvg = likelihood_BetaAvg;
            }
            else st->bestLikelihood_BetaAvg = 0.0;

            if (st->useVlos)
            {
                st->bestLikelihood_VelAvg = likelihood_VelAvg;
            }
            else st->bestLikelihood_VelAvg = 0.0;

            if (st->useDist)
            {
                st->bestLikelihood_Dist = likelihood_Dist;
            }
            else st->bestLikelihood_Dist = 0.0;
            
            /* Calculating the time that the best likelihood occurred */
            st->bestLikelihood_time = ((real) st->step / (real) ctx->nStep) * ctx->timeEvolve;
            
            /* checking how many times the likelihood was improved */
            st->bestLikelihood_count++;
            
            /* if it is an improvement then write out this histogram */
            if (nbf->histoutFileName)
            {
                nbWriteHistogram(nbf->histoutFileName, ctx, st, histogram);
            }
            
            /* if we are creating an out file, store the body tab */
            if(nbf->outFileName)
            {
                memcpy(st->bestLikelihoodBodyTab, st->bodytab, st->nbody * sizeof(Body));
            }
        }
    }
    
    if(data != NULL)
    {
        for(int i = 0; i < 6; i++)
        {
            free(histogram->histograms[i]);
            free(data->histograms[i]);
        }
    }
    free(histogram);
    free(data);

    return NBODY_SUCCESS;
    
}

/* Advance velocity by half a timestep */
static inline void bodyAdvanceVel(Body* p, const mwvector a, const real dtHalf)
{
    mwvector dv;

    dv = mw_mulvs(a, dtHalf);   /* get velocity increment */
    mw_incaddv(Vel(p), dv);     /* advance v by 1/2 step */
}

/* Advance body position by 1 timestep */
static inline void bodyAdvancePos(Body* p, const real dt)
{
    mwvector dr;
    
    dr = mw_mulvs(Vel(p), dt);  /* get position increment */
    mw_incaddv(Pos(p), dr);     /* advance r by 1 step */
}

static inline void advancePosVel(NBodyState* st, const int nbody, const real dt, const mwvector acc_i)
{
    int i;
    real dtHalf = 0.5 * dt;
    Body* bodies = mw_assume_aligned(st->bodytab, 16);
    const mwvector* accs = mw_assume_aligned(st->acctab, 16);

  #ifdef _OPENMP
    #pragma omp parallel for private(i) shared(bodies, accs) schedule(dynamic, 4096 / sizeof(accs[0]))
  #endif
    for (i = 0; i < nbody; ++i)
    {
        bodyAdvanceVel(&bodies[i], mw_addv(accs[i], acc_i), dtHalf);
        bodyAdvancePos(&bodies[i], dt);
    }

}

static inline void advancePosVel_LMC(NBodyState* st, const real dt, const mwvector acc, const mwvector acc_i)
{
    real dtHalf = 0.5 * dt;
    mwvector dr;
    mwvector dv;

    mwvector acc_total = mw_addv(acc, acc_i);
    dv = mw_mulvs(acc_total, dtHalf);
    mw_incaddv(st->LMCvel,dv);

    dr = mw_mulvs(st->LMCvel,dt);
    mw_incaddv(st->LMCpos,dr);
    
}

static inline void advanceVelocities(NBodyState* st, const int nbody, const real dt, const mwvector acc_i1)
{
    int i;
    real dtHalf = 0.5 * dt;
    Body* bodies = mw_assume_aligned(st->bodytab, 16);
    const mwvector* accs = mw_assume_aligned(st->acctab, 16);

  #ifdef _OPENMP
    #pragma omp parallel for private(i) schedule(dynamic, 4096 / sizeof(accs[0]))
  #endif
    for (i = 0; i < nbody; ++i)      /* loop over all bodies */
    {
        bodyAdvanceVel(&bodies[i], mw_addv(accs[i], acc_i1), dtHalf);
    }
}

static inline void advanceVelocities_LMC(NBodyState* st, const real dt, const mwvector acc, const mwvector acc_i)
{
    real dtHalf = 0.5 * dt;
    mwvector dv;

    mwvector acc_total = mw_addv(acc, acc_i);
    dv = mw_mulvs(acc_total, dtHalf);
    mw_incaddv(st->LMCvel,dv);
}


/* stepSystem: advance N-body system one time-step. */
NBodyStatus nbStepSystemPlain(const NBodyCtx* ctx, NBodyState* st, const mwvector acc_i, const mwvector acc_i1)
{
    NBodyStatus rc;
    mwvector acc_LMC;
    
    const real dt = ctx->timestep;
    
    real barTime = st->step * dt - st->previousForwardTime;

    advancePosVel(st, st->nbody, dt, acc_i);   /* acc_i and acc_i1 are accelerations due to the shifting Milky Way */
    if(ctx->LMC){
	acc_LMC = mw_addv(nbExtAcceleration(&ctx->pot, st->LMCpos, barTime), dynamicalFriction_LMC(&ctx->pot, st->LMCpos, st->LMCvel, ctx->LMCmass, ctx->LMCscale, ctx->LMCDynaFric, barTime, ctx->coulomb_log));
        advancePosVel_LMC(st, dt, acc_LMC, acc_i);
    }
    //printf("LMC position: %f %f %f, LMC mass: %f, LMC scale: %f \n", X(st->LMCpos), Y(st->LMCpos), 
    //       Z(st->LMCpos), ctx->LMCmass, ctx->LMCscale);
    //mw_printf("LMC position: %f %f %f, LMC mass: %f, LMC scale: %f \n", X(st->LMCpos), Y(st->LMCpos), 
    //       Z(st->LMCpos), ctx->LMCmass, ctx->LMCscale);
    rc = nbGravMap(ctx, st);
    advanceVelocities(st, st->nbody, dt, acc_i1);
    if(ctx->LMC){
	acc_LMC = mw_addv(nbExtAcceleration(&ctx->pot, st->LMCpos, barTime), dynamicalFriction_LMC(&ctx->pot, st->LMCpos, st->LMCvel, ctx->LMCmass, ctx->LMCscale, ctx->LMCDynaFric, barTime, ctx->coulomb_log));
        advanceVelocities_LMC(st, dt, acc_LMC, acc_i1);
    }

//    mw_printf("(%.15f) LMC position: [%.15f,%.15f,%.15f] | ",(int)(st->step)*(ctx->timestep),X(st->LMCpos[0]),Y(st->LMCpos[0]),Z(st->LMCpos[0]));
//    mw_printf("LMC velocity: [%.15f,%.15f,%.15f]\n",X(st->LMCvel[0]),Y(st->LMCvel[0]),Z(st->LMCvel[0]));

    st->step++;
    #ifdef NBODY_BLENDER_OUTPUT
        blenderPrintBodies(st, ctx);
        printf("Frame: %d\n", (int)(st->step));
    #endif

    return rc;
}

NBodyStatus nbRunSystemPlain(const NBodyCtx* ctx, NBodyState* st, const NBodyFlags* nbf)
{
    NBodyLikelihoodMethod method;
    HistogramParams hp;
    nbGetLikelihoodInfo(nbf, &hp, &method);
   
    if (ctx->LMC){
        //These values are set in nbody_orbit_integrator.c. In the event of a checkpoint, these values are already stored, so running this code would reset them to NULL pointers.
        if (!st->shiftByLMC) {
            mwvector* shiftLMC;
            size_t sizeLMC;
            mwvector LMCx;
            mwvector LMCv;

            getLMCArray(&shiftLMC, &sizeLMC);
            setLMCShiftArray(st, shiftLMC, sizeLMC);
            getLMCPosVel(&LMCx, &LMCv);
            setLMCPosVel(st, LMCx, LMCv);
        }
    }

    NBodyStatus rc = NBODY_SUCCESS;
    rc |= nbGravMap(ctx, st); /* Calculate accelerations for 1st step this episode */
    if (nbStatusIsFatal(rc))
        return rc;

    #ifdef NBODY_BLENDER_OUTPUT
        if(mkdir("./frames", S_IRWXU | S_IRWXG) < 0)
        {
          return NBODY_ERROR;
        }
        deleteOldFiles(st);
        mwvector startCmPos;
        mwvector perpendicularCmPos;
        mwvector nextCmPos;
        nbFindCenterOfMass(&startCmPos, st);
        perpendicularCmPos=startCmPos;
    #endif
        
    real curStep = st->step;
    real Nstep = ctx->nStep;
    
    st->bestLikelihood = DEFAULT_WORST_CASE; //initializing it.

    while (st->step < ctx->nStep)
    {
        #ifdef NBODY_BLENDER_OUTPUT
            nbFindCenterOfMass(&nextCmPos, st);
            blenderPossiblyChangePerpendicularCmPos(&nextCmPos,&perpendicularCmPos,&startCmPos);
        #endif
            
        /* if one needs to add run time options, here is the place to do it
         * this will not run on client side. and provides a good environment
         * to add options without bogging down the client side application
         */    
        #ifdef NBODY_DEV_OPTIONS
            if(ctx->MultiOutput)
            {
                dev_write_outputs(ctx, st, nbf, ctx->OutputFreq);
            }
                
        #endif
        if(!ctx->LMC) {
            mwvector zero;
            SET_VECTOR(zero,0,0,0);
            rc |= nbStepSystemPlain(ctx, st, zero, zero); 
        } else {
            rc |= nbStepSystemPlain(ctx, st, st->shiftByLMC[st->step], st->shiftByLMC[st->step+1]);
        }

        curStep = st->step;
        
        if(curStep / Nstep >= ctx->BestLikeStart && ctx->useBestLike)
        {
            get_likelihood(ctx, st, nbf);
        }
    
        if (nbStatusIsFatal(rc))   /* advance N-body system */
            return rc;

        rc |= nbCheckpoint(ctx, st);
        if (nbStatusIsFatal(rc))
            return rc;
        /* We report the progress at step + 1. 0 is the original
           center of mass. */
        nbReportProgress(ctx, st);
        nbUpdateDisplayedBodies(ctx, st);
    }
    
    #ifdef NBODY_BLENDER_OUTPUT
        blenderPrintMisc(st, ctx, startCmPos, perpendicularCmPos);
    #endif


    return nbWriteFinalCheckpoint(ctx, st);
}

