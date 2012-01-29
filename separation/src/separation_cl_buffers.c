/*
 *  Copyright (c) 2010-2011 Matthew Arsenault
 *  Copyright (c) 2010-2011 Rensselaer Polytechnic Institute
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

#include "milkyway_util.h"
#include "separation_cl_buffers.h"
#include "r_points.h"
#include "calculated_constants.h"

static cl_int createSummarizationBuffers(CLInfo* ci,
                                         SeparationCLMem* cm,
                                         const SeparationSizes* sizes)
{
    cm->summarizationBufs[0] = mwCreateZeroReadWriteBuffer(ci, sizes->summarizationBufs[0]);
    if (!cm->summarizationBufs[0])
    {
        mw_printf("Error creating summarizaton buffer of size "ZU"\n", sizes->summarizationBufs[0]);
        return MW_CL_ERROR;
    }

    cm->summarizationBufs[1] = mwCreateZeroReadWriteBuffer(ci, sizes->summarizationBufs[1]);
    if (!cm->summarizationBufs[1])
    {
        mw_printf("Error creating summarizaton buffer of size "ZU"\n", sizes->summarizationBufs[1]);
        return MW_CL_ERROR;
    }

    return CL_SUCCESS;
}

static cl_int createOutBgBuffer(CLInfo* ci,
                                SeparationCLMem* cm,
                                const SeparationSizes* sizes)
{
    cm->outBg = mwCreateZeroReadWriteBuffer(ci, sizes->outBg);
    if (!cm->outBg)
    {
        mw_printf("Error creating out bg buffer of size "ZU"\n", sizes->outBg);
        return MW_CL_ERROR;
    }

    return CL_SUCCESS;
}

static cl_int createOutStreamsBuffer(CLInfo* ci, SeparationCLMem* cm, const SeparationSizes* sizes)
{
    cm->outStreams = mwCreateZeroReadWriteBuffer(ci, sizes->outStreams);
    if (!cm->outStreams)
    {
        mw_printf("Error creating out probs buffer of size "ZU"\n", sizes->outStreams);
        return MW_CL_ERROR;
    }

    return CL_SUCCESS;
}

static cl_int createSCBuffer(CLInfo* ci,
                             SeparationCLMem* cm,
                             const StreamConstants* sc,
                             const SeparationSizes* sizes,
                             const cl_mem_flags constBufFlags)
{
    cl_int err;
    real* buf;
    cl_int i;

    buf = mwCallocA(sizes->nStream * 8, sizeof(real));

    /* Pack into format used by kernel */
    for (i = 0; i < sizes->nStream; ++i)
    {
        buf[8 * i + 0] = X(sc[i].a);
        buf[8 * i + 1] = X(sc[i].c);

        buf[8 * i + 2] = Y(sc[i].a);
        buf[8 * i + 3] = Y(sc[i].c);

        buf[8 * i + 4] = Z(sc[i].a);
        buf[8 * i + 5] = Z(sc[i].c);

        buf[8 * i + 6] = sc[i].sigma_sq2_inv;
        buf[8 * i + 7] = 0.0;
    }

    cm->sc = clCreateBuffer(ci->clctx, constBufFlags, sizes->sc, (void*) buf, &err);
    mwFreeA(buf);

    if (err != CL_SUCCESS)
    {
        mwPerrorCL(err, "Error creating stream constants buffer of size "ZU, sizes->sc);
        return err;
    }

    return CL_SUCCESS;
}

static cl_int createRBuffers(CLInfo* ci,
                             SeparationCLMem* cm,
                             const AstronomyParameters* ap,
                             const IntegralArea* ia,
                             const StreamGauss sg,
                             const SeparationSizes* sizes,
                             cl_mem_flags constBufFlags)
{
    cl_int err;
    RPoints* r_pts;
    RConsts* rc;

    r_pts = precalculateRPts(ap, ia, sg, &rc, FALSE);

    cm->rPts = clCreateBuffer(ci->clctx, constBufFlags, sizes->rPts, r_pts, &err);

    if (err != CL_SUCCESS)
    {
        mwPerrorCL(err, "Error creating stream r points buffer of size "ZU, sizes->rPts);
        return err;
    }

    cm->rc = clCreateBuffer(ci->clctx, constBufFlags, sizes->rc, rc, &err);
    if (err != CL_SUCCESS)
    {
        mwPerrorCL(err, "Error creating stream r consts buffer of size "ZU, sizes->rc);
        return err;
    }

    cm->sg_dx = clCreateBuffer(ci->clctx, constBufFlags, sizes->sg_dx, sg.dx, &err);
    if (err != CL_SUCCESS)
    {
        mwPerrorCL(err, "Error creating stream sg_dx buffer of size "ZU, sizes->sg_dx);
        return err;
    }

    mwFreeA(r_pts);
    mwFreeA(rc);

    return CL_SUCCESS;
}

static void getSplitLBTrig(const AstronomyParameters* ap,
                           const IntegralArea* ia,
                           LTrigPair** lTrigOut,
                           real** bSinOut)
{
    cl_uint i, j;
    LTrigPair* lTrig;
    real* bSin;
    LBTrig* lbts;
    size_t idx;

    lbts = precalculateLBTrig(ap, ia, FALSE);

    lTrig = (LTrigPair*) mwMallocA(ia->mu_steps * ia->nu_steps * sizeof(LTrigPair));
    bSin = (real*) mwMallocA(ia->mu_steps * ia->nu_steps * sizeof(real));

    for (i = 0; i < ia->nu_steps; ++i)
    {
        for (j = 0; j < ia->mu_steps; ++j)
        {
            idx = i * ia->mu_steps + j;

            lTrig[idx].lCosBCos = lbts[idx].lCosBCos;
            lTrig[idx].lSinBCos = lbts[idx].lSinBCos;

            bSin[idx] = lbts[idx].bSin;
        }
    }

    mwFreeA(lbts);

    *lTrigOut = lTrig;
    *bSinOut = bSin;
}

static cl_int createLBTrigBuffer(CLInfo* ci,
                                 SeparationCLMem* cm,
                                 const AstronomyParameters* ap,
                                 const IntegralArea* ia,
                                 const SeparationSizes* sizes,
                                 const cl_mem_flags constBufFlags)
{
    cl_int err = CL_SUCCESS;
    LTrigPair* lTrig = NULL;
    real* bSin = NULL;

    getSplitLBTrig(ap, ia, &lTrig, &bSin);

    cm->lTrig = clCreateBuffer(ci->clctx, constBufFlags, sizes->lTrig, lTrig, &err);
    if (err != CL_SUCCESS)
    {
        mwPerrorCL(err, "Error creating lTrig buffer of size "ZU, sizes->lTrig);
        return err;
    }

    cm->bSin = clCreateBuffer(ci->clctx, constBufFlags, sizes->bSin, bSin, &err);
    if (err != CL_SUCCESS)
    {
        mwPerrorCL(err, "Error creating bSin buffer of size "ZU, sizes->bSin);
        return err;
    }

    mwFreeA(lTrig);
    mwFreeA(bSin);

    return CL_SUCCESS;
}

static cl_int createAPBuffer(CLInfo* ci,
                             SeparationCLMem* cm,
                             const AstronomyParameters* ap,
                             const SeparationSizes* sizes,
                             const cl_mem_flags constBufFlags)
{
    cl_int err = CL_SUCCESS;
    double buf[16];
    union
    {
        double d;
        cl_uint i[2];
    } item;

    memset(buf, 0, sizeof(buf));

    buf[0] = 0.0;
    buf[1] = 0.0;

    item.i[0] = ap->convolve;
    item.i[1] = ap->number_streams;
    buf[2] = item.d;
    buf[3] = 0.0;

    buf[4] = ap->m_sun_r0;
    buf[5] = ap->r0;

    buf[6] = ap->q_inv_sqr;
    buf[7] = 0.0;

    buf[8] = ap->bg_a;
    buf[9] = ap->bg_b;

    buf[10] = ap->bg_c;
    buf[11] = 0.0;

    cm->ap = clCreateBuffer(ci->clctx, constBufFlags, sizeof(buf), (void*) buf, &err);
    if (err != CL_SUCCESS)
    {
        mwPerrorCL(err, "Error creating astronomy parameters buffer of size "ZU, sizes->ap);
        return err;
    }

    return CL_SUCCESS;
}

void calculateSizes(SeparationSizes* sizes, const AstronomyParameters* ap, const IntegralArea* ia)
{
    assert(_summarizationWorkgroupSize != 0);

    sizes->nStream = ap->number_streams;

    /* globals */
    sizes->summarizationBufs[0] = 2 * sizeof(real) * mwDivRoundup(ia->mu_steps * ia->r_steps, _summarizationWorkgroupSize);
    sizes->summarizationBufs[1] = sizes->summarizationBufs[0]; /* Could make it smaller I guess */
    sizes->outBg = 2 * sizeof(real) * ia->mu_steps * ia->r_steps;
    sizes->outStreams = 2 * sizeof(real) * ia->mu_steps * ia->r_steps * ap->number_streams;

    sizes->rPts = sizeof(RPoints) * ap->convolve * ia->r_steps;
    sizes->lTrig = sizeof(LTrigPair) * ia->mu_steps * ia->nu_steps;
    sizes->bSin = sizeof(real) * ia->mu_steps * ia->nu_steps;

    /* Constant buffer things */
    sizes->ap = 16 * sizeof(double);
    sizes->ia = sizeof(IntegralArea);
    sizes->sc = 8 * sizeof(real) * ap->number_streams;
    sizes->rc = sizeof(RConsts) * ia->r_steps;
    sizes->sg_dx = sizeof(real) * ap->convolve;
}

cl_int createSeparationBuffers(CLInfo* ci,
                               SeparationCLMem* cm,
                               const AstronomyParameters* ap,
                               const IntegralArea* ia,
                               const StreamConstants* sc,
                               const StreamGauss sg,
                               const SeparationSizes* sizes)
{
    cl_int err = CL_SUCCESS;
    cl_mem_flags constBufFlags = CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR;

    err |= createSummarizationBuffers(ci, cm, sizes);
    err |= createOutBgBuffer(ci, cm, sizes);
    err |= createOutStreamsBuffer(ci, cm, sizes);

    err |= createRBuffers(ci, cm, ap, ia, sg, sizes, constBufFlags);
    err |= createLBTrigBuffer(ci, cm, ap, ia, sizes, constBufFlags);

    err |= createSCBuffer(ci, cm, sc, sizes, constBufFlags);
    err |= createAPBuffer(ci, cm, ap, sizes, constBufFlags);

    return err;
}

void releaseSeparationBuffers(SeparationCLMem* cm)
{
    clReleaseMemObject(cm->summarizationBufs[0]);
    clReleaseMemObject(cm->summarizationBufs[1]);
    clReleaseMemObject(cm->outStreams);
    clReleaseMemObject(cm->outBg);

    clReleaseMemObject(cm->rc);
    clReleaseMemObject(cm->rPts);
    clReleaseMemObject(cm->lTrig);
    clReleaseMemObject(cm->bSin);

    clReleaseMemObject(cm->ap);
    clReleaseMemObject(cm->sc);
    clReleaseMemObject(cm->sg_dx);
}

