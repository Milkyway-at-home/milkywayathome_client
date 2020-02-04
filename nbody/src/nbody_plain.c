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


static inline int get_likelihood(const NBodyCtx* ctx, NBodyState* st, const NBodyFlags* nbf)
{
    NBodyHistogram* data = NULL;
    NBodyHistogram* histogram = NULL;
    real likelihood = NAN;
    real likelihood_EMD = NAN;
    real likelihood_Mass = NAN;
    real likelihood_Beta = NAN;
    real likelihood_Vel = NAN;
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
        likelihood = nbSystemLikelihood(st, data, histogram, method);

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

            st->bestLikelihood_Mass = nbCostComponent(data, histogram);

            if (st->useBetaDisp)
            {
                st->bestLikelihood_Beta = nbBetaDispersion(data, histogram);
            }
            else st->bestLikelihood_Beta = 0.0;

            if (st->useVelDisp)
            {
                st->bestLikelihood_Vel = nbVelocityDispersion(data, histogram);
            }
            else st->bestLikelihood_Vel = 0.0;

            st->bestLikelihood_EMD = likelihood-(st->bestLikelihood_Mass)-(st->bestLikelihood_Beta)-(st->bestLikelihood_Vel);
            
            /* Calculating the time that the best likelihood occurred */
            st->bestLikelihood_time = ((real) st->step / (real) ctx->nStep) * ctx->timeEvolve;
            
            /* checking how many times the likelihood was improved */
            st->bestLikelihood_count++;
            
            /* if it is an improvement then write out this histogram */
            if (nbf->histoutFileName)
            {
                nbWriteHistogram(nbf->histoutFileName, ctx, st, histogram);
            }
        }
    }
    
    free(histogram);
    free(data);
    return NBODY_SUCCESS;
    
}


/* stepSystem: advance N-body system one time-step. */
NBodyStatus nbStepSystemPlain(const NBodyCtx* ctx, NBodyState* st, const mwvector acc_i, const mwvector acc_i1)
{
    NBodyStatus rc;
    
    const real dt = ctx->timestep;

    advancePosVel(st, st->nbody, dt, acc_i);

    rc = nbGravMap(ctx, st);
    advanceVelocities(st, st->nbody, dt, acc_i1);

    st->step++;
    #ifdef NBODY_BLENDER_OUTPUT
        blenderPrintBodies(st, ctx);
        printf("Frame: %d\n", (int)(st->step));
    #endif

    return rc;
}

NBodyStatus nbRunSystemPlain(const NBodyCtx* ctx, NBodyState* st, const NBodyFlags* nbf)
{   
    FILE *fp;
    mwvector array[ctx->nStep], shift;
    float ax,ay,az;
    int i = 0;
    if (ctx->LMC){
        for(int j=0; j < ctx->nStep; j++){
            real tempX = st->shiftByLMC[j][0];
            real tempY = st->shiftByLMC[j][1];
            real tempZ = st->shiftByLMC[j][2];
            SET_VECTOR(shift, tempX, tempY, tempZ);
            array[j] = shift;
        }

        /*shift.txt code
    	fp = fopen("shift.txt", "r");
    	for(int j=0; j<ctx->nStep; j++){
        	fscanf(fp, "%f %f %f", &ax, &ay, &az);
        	SET_VECTOR(shift,ax,ay,az);
        	array[j] = shift;
        }*/
    }
    else{
    	for( int j=0; j<ctx->nStep; j++){
        	SET_VECTOR(shift,0,0,0);
        	array[j] = shift;
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
        rc |= nbStepSystemPlain(ctx, st, array[i], array[i+1]);
        i ++;
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


