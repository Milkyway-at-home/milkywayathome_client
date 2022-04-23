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
#include "nbody_autodiff.h"

#include <time.h>

#ifdef NBODY_BLENDER_OUTPUT
  #include "blender_visualizer.h"
#endif

static void nbReportProgress(const NBodyCtx* ctx, NBodyState* st)
{
    real_0 boinc_frac = (real_0) st->step / (real_0) ctx->nStep;
    mw_fraction_done(boinc_frac);

    real_0 frac = (real_0) st->step / (real_0) ctx->nStepRev;

    if (st->reportProgress)
    {
        mw_mvprintw(0, 0,
                    "Running: %f / %f (%f%%)\n",
                    frac * ctx->timeBack,
                    ctx->timeEvolve,
                    100.0 * boinc_frac
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
    real likelihood = mw_real_const(NAN);
    real likelihood_EMD = mw_real_const(NAN);
    real likelihood_Mass = mw_real_const(NAN);
    real likelihood_Beta = mw_real_const(NAN);
    real likelihood_Vel = mw_real_const(NAN);
    real likelihood_BetaAvg = mw_real_const(NAN);
    real likelihood_VelAvg = mw_real_const(NAN);
    real likelihood_Dist = mw_real_const(NAN);

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
        //mw_printf("BEFORE CREATE HISTOGRAM\n");
        histogram = nbCreateHistogram(ctx, st, &hp);
        //mw_printf("AFTER CREATE HISTOGRAM\n");
 
        if (!histogram)
        {
            /* this would normally return a print statement 
             * but I do not want to overload the output since 
             * this would run every time step.
             */
            return 0;
        }
        //mw_printf("BEFORE READ HISTOGRAM\n");
        data = nbReadHistogram(nbf->histogramFileName);
        //mw_printf("AFTER READ HISTOGRAM\n");
        
        if (!data)
        {
            free(histogram);
            /* if the input histogram does not exist, I do not want the 
             * simulation to terminate as you can still get the output file
             * from it. Therefore, this function will end here but with 0
             */
            return 0;
        }
        //mw_printf("BEFORE SYSTEM LIKELIHOOD\n");
        likelihoodArray = nbSystemLikelihood(st, data, histogram, method);
        //mw_printf("AFTER SYSTEM LIKELIHOOD\n");
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
        if (showRealValue(&likelihood) > DEFAULT_WORST_CASE || showRealValue(&likelihood) < (-1 * DEFAULT_WORST_CASE) || isnan(showRealValue(&likelihood)))
        {
            likelihood = mw_real_const(DEFAULT_WORST_CASE);
        }
        else if(showRealValue(&likelihood) == 0.0)
        {
            likelihood = mw_real_const(DEFAULT_BEST_CASE);
        }

        /* this checks to see if the likelihood is an improvement */
        if(mw_fabs_0(showRealValue(&likelihood)) < mw_fabs_0(showRealValue(&st->bestLikelihood)))
        {
            st->bestLikelihood = likelihood;

            st->bestLikelihood_EMD = likelihood_EMD;

            st->bestLikelihood_Mass = likelihood_Mass;

            if (st->useBetaDisp)
            {
                st->bestLikelihood_Beta = likelihood_Beta;
            }
            else st->bestLikelihood_Beta = ZERO_REAL;

            if (st->useVelDisp)
            {
                st->bestLikelihood_Vel = likelihood_Vel;
            }
            else st->bestLikelihood_Vel = ZERO_REAL;

            if (st->useBetaComp)
            {
                st->bestLikelihood_BetaAvg = likelihood_BetaAvg;
            }
            else st->bestLikelihood_BetaAvg = ZERO_REAL;

            if (st->useVlos)
            {
                st->bestLikelihood_VelAvg = likelihood_VelAvg;
            }
            else st->bestLikelihood_VelAvg = ZERO_REAL;

            if (st->useDist)
            {
                st->bestLikelihood_Dist = likelihood_Dist;
            }
            else st->bestLikelihood_Dist = ZERO_REAL;
            
            /* Calculating the time that the best likelihood occurred */
            st->bestLikelihood_time = ((real_0) st->step / (real_0) ctx->nStepRev) * ctx->timeBack;
            
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
static inline void bodyAdvanceVel(Body* p, const mwvector* a_old, const mwvector* acc_i, const real_0 dtHalf, int i)
{
    mwvector a = mw_addv(a_old, acc_i);
    mwvector dv;
    //mw_printf("a = [%.15f,%.15f,%.15f]\n", showRealValue(&X(&a)), showRealValue(&Y(&a)), showRealValue(&Z(&a)));

    dv.x = mw_mul_s(&a.x, dtHalf);   /* get velocity increment */
    dv.y = mw_mul_s(&a.y, dtHalf);
    dv.z = mw_mul_s(&a.z, dtHalf);

    //mw_printf("DV = [%.15f,%.15f,%.15f]\n", showRealValue(&X(&dv)), showRealValue(&Y(&dv)), showRealValue(&Z(&dv)));
    //if (i==0) mw_printf("OLD VEL = [%.15f,%.15f,%.15f]\n", showRealValue(&X(&Vel(p))), showRealValue(&Y(&Vel(p))), showRealValue(&Z(&Vel(p))));
    Vel(p).x = mw_add(&Vel(p).x, &dv.x);     /* advance v by 1/2 step */
    Vel(p).y = mw_add(&Vel(p).y, &dv.y);
    Vel(p).z = mw_add(&Vel(p).z, &dv.z);
    //if (i==0) mw_printf("NEW VEL = [%.15f,%.15f,%.15f]\n", showRealValue(&X(&Vel(p))), showRealValue(&Y(&Vel(p))), showRealValue(&Z(&Vel(p))));
}

/* Advance body position by 1 timestep */
static inline void bodyAdvancePos(Body* p, const real_0 dt, int i)
{
    mwvector dr;
    
    dr.x = mw_mul_s(&X(&Vel(p)), dt);  /* get position increment */
    dr.y = mw_mul_s(&Y(&Vel(p)), dt);
    dr.z = mw_mul_s(&Z(&Vel(p)), dt);

    //mw_printf("DR = [%.15f,%.15f,%.15f]\n", showRealValue(&X(&dr)), showRealValue(&Y(&dr)), showRealValue(&Z(&dr)));
    //if (i==0) mw_printf("OLD POS = [%.15f,%.15f,%.15f]\n", showRealValue(&X(&Pos(p))), showRealValue(&Y(&Pos(p))), showRealValue(&Z(&Pos(p))));
    if (i==0) printVectorFull(&Pos(p), "POS");
    Pos(p).x = mw_add(&Pos(p).x, &dr.x);     /* advance r by 1 step */
    Pos(p).y = mw_add(&Pos(p).y, &dr.y);
    Pos(p).z = mw_add(&Pos(p).z, &dr.z);
    //if (i==0) mw_printf("NEW POS = [%.15f,%.15f,%.15f]\n", showRealValue(&X(&Pos(p))), showRealValue(&Y(&Pos(p))), showRealValue(&Z(&Pos(p))));
}

static inline void advancePosVel(NBodyState* st, const int nbody, const real_0 dt, const mwvector* acc_i)
{
    int i;
    real_0 dtHalf = 0.5 * dt;
    Body* bodies = mw_assume_aligned(st->bodytab, 16);
    const mwvector* accs = mw_assume_aligned(st->acctab, 16);

  #ifdef _OPENMP
    #pragma omp parallel for private(i) shared(bodies, accs) schedule(dynamic, (int) MAX(4096 / sizeof(accs[0]), 1))
  #endif
    for (i = 0; i < nbody; ++i)
    {
        bodyAdvanceVel(&bodies[i], &(accs[i]), acc_i, dtHalf, i);
        bodyAdvancePos(&bodies[i], dt, i);
    }

}

static inline void advancePosVel_LMC(NBodyState* st, const real_0 dt, const mwvector* acc, const mwvector* acc_i)
{
    real_0 dtHalf = 0.5 * dt;
    mwvector dr;
    mwvector dv;

    mwvector acc_total = mw_addv(acc, acc_i);

    dv.x = mw_mul_s(&acc_total.x, dtHalf);
    dv.y = mw_mul_s(&acc_total.y, dtHalf);
    dv.z = mw_mul_s(&acc_total.z, dtHalf);

    st->LMCvel.x = mw_add(&st->LMCvel.x, &dv.x);
    st->LMCvel.y = mw_add(&st->LMCvel.y, &dv.y);
    st->LMCvel.z = mw_add(&st->LMCvel.z, &dv.z);

    dr.x = mw_mul_s(&X(&st->LMCvel), dt);
    dr.y = mw_mul_s(&Y(&st->LMCvel), dt);
    dr.z = mw_mul_s(&Z(&st->LMCvel), dt);

    st->LMCpos.x = mw_add(&st->LMCpos.x, &dr.x);
    st->LMCpos.y = mw_add(&st->LMCpos.y, &dr.y);
    st->LMCpos.z = mw_add(&st->LMCpos.z, &dr.z);
}

static inline void advanceVelocities(NBodyState* st, const int nbody, const real_0 dt, const mwvector* acc_i1)
{
    int i;
    real_0 dtHalf = 0.5 * dt;
    Body* bodies = mw_assume_aligned(st->bodytab, 16);
    const mwvector* accs = mw_assume_aligned(st->acctab, 16);

  #ifdef _OPENMP
    #pragma omp parallel for private(i) schedule(dynamic, (int) MAX(4096 / sizeof(accs[0]), 1))
  #endif
    for (i = 0; i < nbody; ++i)      /* loop over all bodies */
    {
        bodyAdvanceVel(&bodies[i], &(accs[i]), acc_i1, dtHalf, i);
    }
}

static inline void advanceVelocities_LMC(NBodyState* st, const real_0 dt, const mwvector* acc, const mwvector* acc_i)
{
    real_0 dtHalf = 0.5 * dt;
    mwvector dv;

    mwvector acc_total = mw_addv(acc, acc_i);
    dv.x = mw_mul_s(&acc_total.x, dtHalf);
    dv.y = mw_mul_s(&acc_total.y, dtHalf);
    dv.z = mw_mul_s(&acc_total.z, dtHalf);

    st->LMCvel.x = mw_add(&st->LMCvel.x, &dv.x);
    st->LMCvel.y = mw_add(&st->LMCvel.y, &dv.y);
    st->LMCvel.z = mw_add(&st->LMCvel.z, &dv.z);
}


/* stepSystem: advance N-body system one time-step. */
NBodyStatus nbStepSystemPlain(const NBodyCtx* ctx, NBodyState* st, const mwvector* acc_i, const mwvector* acc_i1)
{
    NBodyStatus rc;
    mwvector acc_LMC, acc_DF;
    real massLMC, scaleLMC;
    //mw_printf("ACC_I  = [%.15f, %.15f, %.15f]\n", showRealValue(&X(acc_i)), showRealValue(&Y(acc_i)), showRealValue(&Z(acc_i)) );
    //mw_printf("ACC_I1 = [%.15f, %.15f, %.15f]\n", showRealValue(&X(acc_i1)), showRealValue(&Y(acc_i1)), showRealValue(&Z(acc_i1)) );
    
    const real_0 dt = ctx->timestep;
    
    real_0 barTime = st->step * dt - st->previousForwardTime;

    //mw_printf("   LMC_POS = [%.15f, %.15f, %.15f]\n", showRealValue(&X(&st->LMCpos)), showRealValue(&Y(&st->LMCpos)), showRealValue(&Z(&st->LMCpos)));
    //mw_printf("   LMC_VEL = [%.15f, %.15f, %.15f]\n", showRealValue(&X(&st->LMCvel)), showRealValue(&Y(&st->LMCvel)), showRealValue(&Z(&st->LMCvel)));

    advancePosVel(st, st->nbody, dt, acc_i);   /* acc_i and acc_i1 are accelerations due to the shifting Milky Way */
    if(ctx->LMC){
        massLMC = mw_real_var(ctx->LMCmass, LMC_MASS_POS);
        scaleLMC = mw_real_var(ctx->LMCscale, LMC_RADIUS_POS);
        acc_DF = dynamicalFriction_LMC(&ctx->pot, &st->LMCpos, &st->LMCvel, &massLMC, &scaleLMC, ctx->LMCDynaFric, barTime, ctx->coulomb_log);
        acc_LMC = nbExtAcceleration(&ctx->pot, &st->LMCpos, barTime);
        //mw_printf("  DF_ACC = [%.15f, %.15f, %.15f]\n", showRealValue(&X(&acc_DF)), showRealValue(&Y(&acc_DF)), showRealValue(&Z(&acc_DF)));
        //mw_printf(" LMC_ACC = [%.15f, %.15f, %.15f]\n", showRealValue(&X(&acc_LMC)), showRealValue(&Y(&acc_LMC)), showRealValue(&Z(&acc_LMC)));
	acc_LMC = mw_addv(&acc_LMC, &acc_DF);
        //mw_printf(" LMC_ACC = [%.15f, %.15f, %.15f]\n", showRealValue(&X(&acc_LMC)), showRealValue(&Y(&acc_LMC)), showRealValue(&Z(&acc_LMC)));
        advancePosVel_LMC(st, dt, &acc_LMC, acc_i);
    }
    rc = nbGravMap(ctx, st);
    advanceVelocities(st, st->nbody, dt, acc_i1);
    if(ctx->LMC){
        acc_DF = dynamicalFriction_LMC(&ctx->pot, &st->LMCpos, &st->LMCvel, &massLMC, &scaleLMC, ctx->LMCDynaFric, barTime, ctx->coulomb_log);
        acc_LMC = nbExtAcceleration(&ctx->pot, &st->LMCpos, barTime);
        //mw_printf("  DF_ACC = [%.15f, %.15f, %.15f]\n", showRealValue(&X(&acc_DF)), showRealValue(&Y(&acc_DF)), showRealValue(&Z(&acc_DF)));
        //mw_printf(" LMC_ACC = [%.15f, %.15f, %.15f]\n", showRealValue(&X(&acc_LMC)), showRealValue(&Y(&acc_LMC)), showRealValue(&Z(&acc_LMC)));
	acc_LMC = mw_addv(&acc_LMC, &acc_DF);
        //mw_printf(" LMC_ACC = [%.15f, %.15f, %.15f]\n", showRealValue(&X(&acc_LMC)), showRealValue(&Y(&acc_LMC)), showRealValue(&Z(&acc_LMC)));
        advanceVelocities_LMC(st, dt, &acc_LMC, acc_i1);
    }

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
    clock_t time_stamp;
    real_0 time_taken;

    //real_0 input0 = 2.0;
    //real input1 = mw_real_var(3.0, 0);
    //real input2 = mw_real_var(4.0, 1);
    //real input3 = mw_real_var(5.0, 2);
    //real output4 = mw_mad(&input1, &input2, &input3);
    //printReal(&output4);
   
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
    mw_printf("Calculating Initial Accelerations...\n");
    rc |= nbGravMap(ctx, st); /* Calculate accelerations for 1st step this episode */
    mw_printf("Found Initial Accelerations!\n");
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
        
    real_0 curStep = st->step;
    real_0 Nstep = ctx->nStep;
    
    st->bestLikelihood = mw_real_const(DEFAULT_WORST_CASE); //initializing it.

    mw_printf("Beginning N-Body Simulation...\n");
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

        //time_stamp = clock();
        //mw_printf("    BEFORE STEP SYSTEM\n");
        if(!ctx->LMC) {
            mwvector zero;
            SET_VECTOR(&zero,ZERO_REAL,ZERO_REAL,ZERO_REAL);
            //mw_printf("    ZERO VECTOR = [%.15f, %.15f, %.15f]\n", showRealValue(&X(&zero)), showRealValue(&Y(&zero)), showRealValue(&Z(&zero)) );
            rc |= nbStepSystemPlain(ctx, st, &zero, &zero); 
        } else {
            rc |= nbStepSystemPlain(ctx, st, &st->shiftByLMC[st->step], &st->shiftByLMC[st->step+1]);
        }
        //mw_printf("    AFTER STEP SYSTEM\n");
        //time_stamp = clock() - time_stamp;
        //time_taken = ((real_0) time_stamp) / CLOCKS_PER_SEC;
        //mw_printf("%.15f\n", time_taken);

        curStep = st->step;
        
        if(curStep / Nstep >= ctx->BestLikeStart && ctx->useBestLike)
        {
            //mw_printf(" TIME = %.15f\n", ctx->timestep*curStep);
            get_likelihood(ctx, st, nbf);
            //mw_printf("    AFTER LIKELIHOOD\n");
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
        blenderPrintMisc(st, ctx, &startCmPos, &perpendicularCmPos);
    #endif


    return nbWriteFinalCheckpoint(ctx, st);
}

