/*
 *  Copyright (c) 1993, 2001 Joshua E. Barnes, Honolulu, HI.
 *  Copyright (c) 2010-2011 Rensselaer Polytechnic Institute
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

#include "nbody_priv.h"
#include "nbody_util.h"
#include "nbody_grav.h"
#include "milkyway_util.h"

#ifdef _OPENMP
  #include <omp.h>
#endif /* _OPENMP */


/*
 * nbodyGravity: Walk the tree starting at the root to do force
 * calculations.
 *
 *
 * Random notes:
 *   - Not inlined without inline from multiple calls in
 *     mapForceBody(). Measurably better with the inline, but only
 *     slightly.
 */
static inline mwvector nbGravity(const NBodyCtx* ctx, NBodyState* st, const Body* p)
{
    mwbool showProblem = FALSE;
    NBodyNode* q;

    mwvector pos0 = Pos(p);
    mwvector acc0 = ZERO_VECTOR;
    real tmp1, tmp2;
    real drab, phii, mor3;

    /*Code here is for debugging tree-incest errors*/
    mwbool checkIncest = TRUE;
    mwbool skipSelf = FALSE;
    if (checkIncest)
    {
        q = (NBodyNode*) st->tree.root; 
        while (q != NULL)               /* while not at end of scan */
        {
            mwvector dr = mw_subv(&Pos(q), &pos0);   /* Then compute distance */
            real drSq = mw_sqrv(&dr);               /* and distance squared */

            if (isBody(q) || (showRealValue(&drSq) >= Rcrit2(q)))      /* If is a body */
            {
                if (mw_likely((Body*) q != p))   /* self-interaction? */
                {
                    if ((showProblem)&&(!isBody(q)))
                    {
                        //mw_printf("%u / %.15f >= %.15f / MASS = %.15f\n", isBody(q), showRealValue(&drSq), Rcrit2(q), showRealValue(&Mass(q)));
                        mw_printf("");
                    }
                }
                else
                {
                    skipSelf = TRUE;   /* Encountered self */
                }
                q = Next(q);  /* Follow next link */
            }
            else
            {
                 q = More(q); /* Follow to the next level if need to go deeper */
            }
        }

        if (!skipSelf)
        {
            /* If a body does not encounter itself in its traversal of the
             * tree, it is "tree incest" */
            //mw_printf("Incest found in first part of nbody_grav.c at %p = [%.15f, %.15f, %.15f]\n", p, showRealValue(&X(&Pos(p))), showRealValue(&Y(&Pos(p))), showRealValue(&Z(&Pos(p))) );
            mw_printf("");
        }
    }
    
    /*----------------------------------------------*/

    skipSelf = FALSE;
    showProblem = FALSE;
    q = (NBodyNode*) st->tree.root; /* Start at the root */

    while (q != NULL)               /* while not at end of scan */
    {
        mwvector dr = mw_subv(&Pos(q), &pos0);   /* Then compute distance */
        real drSq = mw_sqrv(&dr);               /* and distance squared */
//        if (showProblem)
//        {
//            mw_printf("    Q POINTER %p = [%.15f, %.15f, %.15f]\n", q, showRealValue(&X(&Pos(q))), showRealValue(&Y(&Pos(q))), showRealValue(&Z(&Pos(q))) );
//            mw_printf("    DR = [%.15f, %.15f, %.15f]\n", showRealValue(&X(&Pos(&dr))), showRealValue(&Y(&Pos(&dr))), showRealValue(&Z(&Pos(&dr))) );
//            mw_printf("    DRSQ = %.15f\n", showRealValue(&drSq) );
//        }

        if (isBody(q) || (showRealValue(&drSq) >= Rcrit2(q)))      /* If is a body or far enough away to approximate */
        {
            if (mw_likely((Body*) q != p))   /* self-interaction? */
            {
                /* Compute gravity */

                drSq = mw_add_s(&drSq, ctx->eps2);   /* use standard softening */
                drab = mw_sqrt(&drSq);
                phii = mw_div(&Mass(q), &drab);
                mor3 = mw_div(&phii, &drSq);

                acc0.x = mw_mad(&mor3, &dr.x, &acc0.x);
                acc0.y = mw_mad(&mor3, &dr.y, &acc0.y);
                acc0.z = mw_mad(&mor3, &dr.z, &acc0.z);

                if (ctx->useQuad && isCell(q))          /* if cell, add quad term */
                {
                    real dr5inv, drQdr, phiQ;
                    mwvector Qdr;

                    /* form Q * dr */
                    tmp1 = mw_mul(&Quad(q).xy, &dr.y);
                    tmp2 = mw_mul(&Quad(q).xz, &dr.z);
                    tmp2 = mw_add(&tmp1, &tmp2);
                    tmp1 = mw_mul(&Quad(q).xx, &dr.x);
                    Qdr.x = mw_add(&tmp1, &tmp2);

                    tmp1 = mw_mul(&Quad(q).yy, &dr.y);
                    tmp2 = mw_mul(&Quad(q).yz, &dr.z);
                    tmp2 = mw_add(&tmp1, &tmp2);
                    tmp1 = mw_mul(&Quad(q).xy, &dr.x);
                    Qdr.y = mw_add(&tmp1, &tmp2);

                    tmp1 = mw_mul(&Quad(q).yz, &dr.y);
                    tmp2 = mw_mul(&Quad(q).zz, &dr.z);
                    tmp2 = mw_add(&tmp1, &tmp2);
                    tmp1 = mw_mul(&Quad(q).xz, &dr.x);
                    Qdr.z = mw_add(&tmp1, &tmp2);


                    /* form dr * Q * dr */
                    tmp1 = mw_mul(&Qdr.y, &dr.y);
                    tmp2 = mw_mul(&Qdr.z, &dr.z);
                    tmp2 = mw_add(&tmp1, &tmp2);
                    tmp1 = mw_mul(&Qdr.x, &dr.x);
                    drQdr = mw_add(&tmp1, &tmp2);

                    tmp1 = sqr(&drSq);
                    tmp1 = mw_mul(&tmp1, &drab);
                    dr5inv = inv(&tmp1);  /* form dr^-5 */

                    /* get quad. part of phi */
                    tmp1 = mw_mul(&dr5inv, &drQdr);
                    tmp1 = mw_div(&tmp1, &drSq);
                    phiQ = mw_mul_s(&tmp1,2.5);

                    acc0.x = mw_mad(&phiQ, &dr.x, &acc0.x);
                    acc0.y = mw_mad(&phiQ, &dr.y, &acc0.y);
                    acc0.z = mw_mad(&phiQ, &dr.z, &acc0.z);

                    /* acceleration */
                    tmp1 = mw_mul(&dr5inv, &Qdr.x);
                    acc0.x = mw_sub(&acc0.x, &tmp1);

                    tmp1 = mw_mul(&dr5inv, &Qdr.y);
                    acc0.y = mw_sub(&acc0.y, &tmp1);

                    tmp1 = mw_mul(&dr5inv, &Qdr.z);
                    acc0.z = mw_sub(&acc0.z, &tmp1);
                }
            }
            else
            {
                skipSelf = TRUE;   /* Encountered self */
                //mw_printf("Encountered Self! ");
            }

            q = Next(q);  /* Follow next link */
        }
        else
        {
             q = More(q); /* Follow to the next level if need to go deeper */
        }
    }

    if (!skipSelf)
    {
        /* If a body does not encounter itself in its traversal of the
         * tree, it is "tree incest" */

        nbReportTreeIncest(ctx, st);
        //mw_printf("Incest found in second part of nbody_grav.c at %p = [%.15f, %.15f, %.15f]\n", p, showRealValue(&X(&Pos(p))), showRealValue(&Y(&Pos(p))), showRealValue(&Z(&Pos(p))) );
    }

    return acc0;
}

static inline void nbMapForceBody(const NBodyCtx* ctx, NBodyState* st)
{
    int i;
    const int nbody = st->nbody;  /* Prevent reload on each loop */
    mwvector LMCx;
    mwvector a, a_tmp;
    const Body* b;
    real lmcmass, lmcscale;

    const Body* bodies = mw_assume_aligned(st->bodytab, 16);
    mwvector* accels = mw_assume_aligned(st->acctab, 16);
    real_0 curTime = st->step * ctx->timestep;
    real_0 timeFromStart = (-1)*ctx->Ntsteps*ctx->timestep + curTime;

    //use previous calibration run to shift time and calibrate the bar
    real_0 barTime = st->step * ctx->timestep - st->previousForwardTime;

    if (ctx->LMC) {
        LMCx = st->LMCpos;
        lmcmass = mw_real_var(ctx->LMCmass, LMC_MASS_POS);      //Sets LMC mass and radius at positions 19 and 20 in the gradient
        lmcscale = mw_real_var(ctx->LMCscale, LMC_RADIUS_POS);
    }
    else {
        SET_VECTOR(&LMCx, ZERO_REAL, ZERO_REAL, ZERO_REAL);
        lmcmass = mw_real_var(0.0, LMC_MASS_POS);
        lmcscale = mw_real_var(1.0, LMC_RADIUS_POS);
    }

  //mw_printf("Size of Accels = %u\n", sizeof(accels[0]));
  #ifdef _OPENMP
    #pragma omp parallel for private(i, b, a, a_tmp) shared(bodies, accels) schedule(dynamic, (int) MAX(4096 / sizeof(accels[0]), 1))
  #endif
    for (i = 0; i < nbody; ++i)      /* get force on each body */
    {
        /* Repeat the base hackGrav part in each case or else GCC's
         * -funswitch-loops doesn't happen. Without that this constant
         * gets checked on every body on every step which is dumb.  */
        switch (ctx->potentialType)
        {
            case EXTERNAL_POTENTIAL_DEFAULT:
                //mw_printf("DEFAULT POTENTIAL - TREE\n");
                b = &bodies[i];
                a = nbGravity(ctx, st, b);
                //mw_printf("nbGravity = [%.15f, %.15f, %.15f]\n", showRealValue(&X(&a)), showRealValue(&Y(&a)), showRealValue(&Z(&a)) );
                a_tmp = nbExtAcceleration(&ctx->pot, &Pos(b), barTime);
                a = mw_addv(&a, &a_tmp); //Adding External Acceleration
                //mw_printf("nbExtAcc = [%.15f, %.15f, %.15f]\n", showRealValue(&X(&a_tmp)), showRealValue(&Y(&a_tmp)), showRealValue(&Z(&a_tmp)) );
                a_tmp = plummerAccel(&Pos(b), &LMCx, &lmcmass, &lmcscale);
                a = mw_addv(&a, &a_tmp); //Adding LMC Acceleration
                //mw_printf("plummerAcc = [%.15f, %.15f, %.15f]\n", showRealValue(&X(&a_tmp)), showRealValue(&Y(&a_tmp)), showRealValue(&Z(&a_tmp)) );

                //mw_printf("ACC = [%.15f, %.15f, %.15f]\n", showRealValue(&X(&a)), showRealValue(&Y(&a)), showRealValue(&Z(&a)) );
                accels[i] = a;
                break;

            case EXTERNAL_POTENTIAL_NONE:
                //mw_printf("NULL POTENTIAL - TREE\n");
                accels[i] = nbGravity(ctx, st, &bodies[i]);
                break;

            case EXTERNAL_POTENTIAL_CUSTOM_LUA:
                //mw_printf("CUSTOM POTENTIAL - TREE\n");
                a = nbGravity(ctx, st, &bodies[i]);
                nbEvalPotentialClosure(st, &Pos(&bodies[i]), &a_tmp);
                a = mw_addv(&a, &a_tmp); //Adding External Acceleration
                a_tmp = plummerAccel(&Pos(&bodies[i]), &LMCx, &lmcmass, &lmcscale);
                a = mw_addv(&a, &a_tmp); //Adding LMC Acceleration

                accels[i] = a;
                break;

            default:
                mw_fail("Bad external potential type: %d\n", ctx->potentialType);
        }
    }
}

static mwvector nbGravity_Exact(const NBodyCtx* ctx, NBodyState* st, const Body* p)
{
    int i;
    const int nbody = st->nbody;
    mwvector a = ZERO_VECTOR;
    const real_0 eps2 = ctx->eps2;
    real tmp;
    mwvector tmp_vec;

    //mw_printf("P POS = [%.15f, %.15f, %.15f]\n", showRealValue(&X(&Pos(p))), showRealValue(&Y(&Pos(p))), showRealValue(&Z(&Pos(p))) );
    for (i = 0; i < nbody; ++i)
    {
        const Body* b = &st->bodytab[i];
        //mw_printf(" B POS = [%.15f, %.15f, %.15f]\n", showRealValue(&X(&Pos(b))), showRealValue(&Y(&Pos(b))), showRealValue(&Z(&Pos(b))) );

        mwvector dr = mw_subv(&Pos(b), &Pos(p));
        //mw_printf(" DR = [%.15f, %.15f, %.15f]\n", showRealValue(&X(&dr)), showRealValue(&Y(&dr)), showRealValue(&Z(&dr)) );
        tmp = mw_sqrv(&dr);
        //mw_printf(" tmp = %.15f\n", showRealValue(&tmp));
        real drSq = mw_add_s(&tmp, eps2);
        //mw_printf(" drSq = %.15f\n", showRealValue(&drSq));

        real drab = mw_sqrt(&drSq);
        real phii = mw_div(&Mass(b), &drab);
        real mor3 = mw_div(&phii, &drSq);

        tmp_vec = mw_mulvs(&dr, &mor3);
        //mw_printf("tmp_vec = [%.15f, %.15f, %.15f]\n", showRealValue(&X(&tmp_vec)), showRealValue(&Y(&tmp_vec)), showRealValue(&Z(&tmp_vec)) );
        a = mw_addv(&a, &tmp_vec);
        //mw_printf("ACC = [%.15f, %.15f, %.15f]\n", showRealValue(&X(&a)), showRealValue(&Y(&a)), showRealValue(&Z(&a)) );
    }
    //mw_printf("ACC = [%.15f, %.15f, %.15f]\n", showRealValue(&X(&a)), showRealValue(&Y(&a)), showRealValue(&Z(&a)) );
    return a;
}

static inline void nbMapForceBody_Exact(const NBodyCtx* ctx, NBodyState* st)
{
    int i;
    const int nbody = st->nbody;  /* Prevent reload on each loop */
    mwvector LMCx;
    mwvector a, a_tmp;
    const Body* b;
    real lmcmass, lmcscale;

    Body* bodies = mw_assume_aligned(st->bodytab, 16);
    mwvector* accels = mw_assume_aligned(st->acctab, 16);
    real_0 curTime = st->step * ctx->timestep;
    real_0 timeFromStart = -ctx->Ntsteps*ctx->timestep + curTime;
    real_0 barTime = st->step * ctx->timestep - st->previousForwardTime;

    if (ctx->LMC) {
        LMCx = st->LMCpos;
        lmcmass = mw_real_var(ctx->LMCmass, LMC_MASS_POS);      //Sets LMC mass and radius at positions 19 and 20 in the gradient
        lmcscale = mw_real_var(ctx->LMCscale, LMC_RADIUS_POS);
    }
    else {
        SET_VECTOR(&LMCx, ZERO_REAL, ZERO_REAL, ZERO_REAL);
        lmcmass = mw_real_var(0.0, LMC_MASS_POS);
        lmcscale = mw_real_var(1.0, LMC_RADIUS_POS);
    }

  #ifdef _OPENMP
    #pragma omp parallel for private(i, b, a, a_tmp) shared(bodies, accels) schedule(dynamic, (int) MAX(4096 / sizeof(accels[0]), 1))
  #endif

    for (i = 0; i < nbody; ++i)      /* get force on each body */
    {
        switch (ctx->potentialType)
        {
            case EXTERNAL_POTENTIAL_DEFAULT:
                //mw_printf("DEFAULT POTENTIAL - EXACT\n");
                b = &bodies[i];
                a = nbGravity_Exact(ctx, st, b);
                a_tmp = nbExtAcceleration(&ctx->pot, &Pos(b), barTime);
                a = mw_addv(&a, &a_tmp); //Adding External Acceleration
                a_tmp = plummerAccel(&Pos(b), &LMCx, &lmcmass, &lmcscale);
                a = mw_addv(&a, &a_tmp); //Adding LMC Acceleration
                
                accels[i] = a;
                break;

            case EXTERNAL_POTENTIAL_NONE:
                //mw_printf("NULL POTENTIAL - EXACT\n");
                accels[i] = nbGravity_Exact(ctx, st, &bodies[i]);
                break;

            case EXTERNAL_POTENTIAL_CUSTOM_LUA:
                //mw_printf("CUSTOM POTENTIAL - EXACT\n");
                a = nbGravity_Exact(ctx, st, &bodies[i]);
                nbEvalPotentialClosure(st, &Pos(&bodies[i]), &a_tmp);
                a = mw_addv(&a, &a_tmp); //Adding External Acceleration
                a_tmp = plummerAccel(&Pos(&bodies[i]), &LMCx, &lmcmass, &lmcscale);
                a = mw_addv(&a, &a_tmp); //Adding LMC Acceleration

                accels[i] = a;
                break;

            default:
                mw_fail("Bad external potential type: %d\n", ctx->potentialType);
        }
    }
}

static inline NBodyStatus nbIncestStatusCheck(const NBodyCtx* ctx, const NBodyState* st)
{
    if (st->treeIncest)
    {
        return ctx->allowIncest ? NBODY_TREE_INCEST_NONFATAL : NBODY_TREE_INCEST_FATAL;
    }

    return NBODY_SUCCESS;
}

NBodyStatus nbGravMap(const NBodyCtx* ctx, NBodyState* st)
{
    NBodyStatus rc;

    if (mw_likely(ctx->criterion != Exact))
    {
        rc = nbMakeTree(ctx, st);
        if (nbStatusIsFatal(rc))
            return rc;

        nbMapForceBody(ctx, st);
    }
    else
    {
        nbMapForceBody_Exact(ctx, st);
    }

    if (st->potentialEvalError)
    {
        return NBODY_LUA_POTENTIAL_ERROR;
    }

    return nbIncestStatusCheck(ctx, st); /* Check if incest occured during step */
}

