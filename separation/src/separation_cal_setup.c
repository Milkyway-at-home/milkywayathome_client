/*
Copyright (C) 2010, 2011  Matthew Arsenault

This file is part of Milkway@Home.

Milkyway@Home is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Milkyway@Home is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Milkyway@Home.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "milkyway_util.h"
#include "separation_types.h"
#include "calculated_constants.h"
#include "r_points.h"
#include "show_cal_types.h"
#include "separation_cal_setup.h"
#include "separation_cal_types.h"
#include "separation_cal_kernelgen.h"

#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wswitch-enum"
#endif



PFNCALCTXWAITFOREVENTS mw_calCtxWaitForEvents = NULL;

static FILE* _isaLogFunctionFile = NULL;

static void isaLogFunction(const char* msg)
{
    if (_isaLogFunctionFile)
        fputs(msg, _isaLogFunctionFile);
}

static CALresult printISA(FILE* f, CALimage image)
{
    if (!image)
        return CAL_RESULT_ERROR;

    _isaLogFunctionFile = f;
    calclDisassembleImage(image, isaLogFunction);
    _isaLogFunctionFile = NULL;

    return CAL_RESULT_OK;
}

static CALresult releaseMWMemRes(CALcontext ctx, MWMemRes* mr)
{
    CALresult err = CAL_RESULT_OK;

    if (mr->mem)
    {
        err = calCtxReleaseMem(ctx, mr->mem);
        if (err != CAL_RESULT_OK)
            cal_warn("Failed to release CALmem", err);
        mr->mem = 0;
    }

    if (mr->res)
    {
        err = calResFree(mr->res);
        if (err != CAL_RESULT_OK)
            cal_warn("Failed to release CAL resource", err);
        mr->res = 0;
    }

    return err;
}

CALresult mwUnloadKernel(MWCALInfo* ci)
{
    CALresult err = CAL_RESULT_OK;
    CALresult erri;

    if (ci->module)
    {
        erri = calModuleUnload(ci->calctx, ci->module);
        if (erri != CAL_RESULT_OK)
            cal_warn("Failed to unload module", erri);
        ci->module = 0;
        err |= erri;
    }

    if (ci->image)
    {
        erri = calImageFree(ci->image);
        if (erri != CAL_RESULT_OK)
            cal_warn("Failed free image", erri);
        ci->image = 0;
    }

    return err;
}

static CALresult mwDestroyCALInfo(MWCALInfo* ci)
{
    CALresult err = CAL_RESULT_OK;
    CALresult erri;

    err |= mwUnloadKernel(ci);
    if (ci->calctx)
    {
        erri = calCtxDestroy(ci->calctx);
        if (erri != CAL_RESULT_OK)
            cal_warn("Failed to destroy CAL context", erri);
        ci->calctx = 0;
        err |= erri;
    }

    if (ci->dev)
    {
        erri = calDeviceClose(ci->dev);
        if (erri != CAL_RESULT_OK)
            cal_warn("Failed to close device", erri);
        ci->dev = 0;
        err |= erri;
    }

    if (err != CAL_RESULT_OK)
        cal_warn("Failed to cleanup CAL info", err);

    return err;
}

static CALresult mwGetDevice(MWCALInfo* ci, CALuint devID)
{
    CALresult err;

    err = calDeviceGetCount(&ci->numDevices);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to get CAL device count", err);
        return err;
    }

    if (ci->numDevices == 0)
    {
        warn("Didn't find any CAL devices\n");
        return CAL_RESULT_ERROR;
    }

    if (devID + 1 > ci->numDevices)
    {
        warn("Requested device ID %u > found number of devices (%u)\n",
             devID, ci->numDevices);
        return CAL_RESULT_ERROR;
    }

    ci->devID = devID;
    err = calDeviceOpen(&ci->dev, ci->devID);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to open CAL device", err);
        return err;
    }

    err = calDeviceGetInfo(&ci->devInfo, ci->devID);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to get CAL device information", err);
        return err;
    }

    ci->devAttribs.struct_size = sizeof(struct CALdeviceattribsRec);
    err = calDeviceGetAttribs(&ci->devAttribs, ci->devID);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to get CAL device attributes", err);
        return err;
    }

    return CAL_RESULT_OK;
}

/* Find devices and create context */
static CALresult mwGetCALInfo(MWCALInfo* ci, CALuint devID)
{
    CALresult err;

    err = calGetVersion(&ci->version.major, &ci->version.minor, &ci->version.patchLevel);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to get CAL version", err);
        return err;
    }

    err = mwGetDevice(ci, devID);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Error getting device information", err);
        return err;
    }

    err = calCtxCreate(&ci->calctx, ci->dev);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to create CAL context", err);
        return err;
    }

    return CAL_RESULT_OK;
}

/* Try to get memory handle and cleanup resource if that fails */
static CALresult getMemoryHandle(MWMemRes* mr, MWCALInfo* ci)
{
    CALresult err = CAL_RESULT_OK;

    err = calCtxGetMem(&mr->mem, ci->calctx, mr->res);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to get memory handle", err);
        if (calResFree(mr->res) != CAL_RESULT_OK)
            warn("Failed to release CAL resource\n");
        else
            mr->res = 0;
    }

    return err;
}

/* Try to map the resource and free it on failure */
CALresult mapMWMemRes(MWMemRes* mr, CALvoid** pPtr, CALuint* pitch)
{
    CALresult err;

    err = calResMap(pPtr, pitch, mr->res, 0);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to map resource", err);
        if (calResFree(mr->res) != CAL_RESULT_OK)
            warn("Failed to release CAL resource\n");
        else
            mr->res = 0;
    }

    return err;
}

/* Try to unmap resource and free it on failure */
CALresult unmapMWMemRes(MWMemRes* mr)
{
    CALresult err;

    err = calResUnmap(mr->res);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to unmap resource", err);
        if (calResFree(mr->res) != CAL_RESULT_OK)
            warn("Failed to release CAL resource\n");
        else
            mr->res = 0;
    }

    return err;
}

static CALresult printBufferDouble(MWMemRes* mr,
                                   CALuint numberElements,
                                   CALuint width,
                                   CALuint height)
{
    CALuint i, j, k, pitch;
    CALdouble* bufPtr;
    CALdouble* tmp;
    CALresult err;

    err = mapMWMemRes(mr, (CALvoid**) &bufPtr, &pitch);
    if (err != CAL_RESULT_OK)
        return err;

    for (i = 0; i < height; ++i)
    {
        tmp = &bufPtr[i * numberElements * pitch];
        for (j = 0; j < width; ++j)
        {
            for (k = 0; k < numberElements; ++k)
                warn("%22.16lf", tmp[numberElements * j + k]);
            warn(" \n");
        }
    }

    return unmapMWMemRes(mr);
}

/* both arguments set on nu step share 1 single element buffer */
static CALresult createNuCB(MWMemRes* mr, MWCALInfo* ci)
{
    CALresult err;

    err = calResAllocRemote1D(&mr->res, &ci->dev, 1, 1, constantFormatReal2, CAL_RESALLOC_CACHEABLE);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to allocate nu buffer", err);
        releaseMWMemRes(ci->calctx, mr);
        return err;
    }

    return getMemoryHandle(mr, ci);
}


static CALresult setConstants(SeparationCALMem* cm,
                              const AstronomyParameters* ap,
                              const IntegralArea* ia,
                              const StreamConstants* sc)
{
    CALresult err;
    CALdouble* cPtr;
    CALuint pitch = 0;

    union
    {
        CALdouble d;
        CALuint size[2];
    } area;

    const CALuint extra = mwNextMultiple(ia->r_steps * ia->mu_steps, 64) - ia->r_steps * ia->mu_steps;

    warn("Extra = %u\n", extra);


    err = mapMWMemRes(&cm->consts, (CALvoid**) &cPtr, &pitch);
    if (err != CAL_RESULT_OK)
        return err;

    area.size[0] = ia->mu_steps;
    area.size[1] = ia->r_steps;
    cPtr[0] = area.d;        // cb1[0].xy

    area.size[0] = ia->nu_steps;
    area.size[1] = extra;
    cPtr[1] = area.d;        // cb1[0].zw


    area.size[0] = ap->convolve;
    area.size[1] = ap->number_streams;
    cPtr[2] = area.d;        // cb1[1].xy
    cPtr[3] = 0.0;           // cb1[1].zw


    cPtr[4] = ap->m_sun_r0;  // cb1[2].xy
    cPtr[5] = ap->r0;        // cb1[2].zw

    cPtr[6] = ap->q_inv_sqr; // cb1[3].xy
    cPtr[7] = 0.0;           // cb1[3].zw

    cPtr[8] = ap->bg_a;      // cb1[4].xy
    cPtr[9] = ap->bg_b;      // cb1[4].zw

    cPtr[10] = ap->bg_c;     // cb1[5].xy
    cPtr[11] = 0.0;          // cb1[5].zw

    cPtr[12] = 0.0;
    cPtr[13] = 0.0;

    cPtr[14] = 0.0;
    cPtr[15] = 0.0;

    err = unmapMWMemRes(&cm->consts);
    if (err != CAL_RESULT_OK)
        return err;

    return CAL_RESULT_OK;
}

static CALresult createConstsCB(MWMemRes* mr, MWCALInfo* ci)
{
    CALresult err;

    err = calResAllocLocal1D(&mr->res, ci->dev, 8, constantFormatReal2, CAL_RESALLOC_CACHEABLE);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to allocate constant buffer", err);
        releaseMWMemRes(ci->calctx, mr);
        return err;
    }

    return getMemoryHandle(mr, ci);
}

static CALuint formatToNumElements(CALformat x)
{
    switch (x)
    {
        case formatReal1:
        case constantFormatReal1:
            return 1;
        case formatReal2:
        case constantFormatReal2:
            return 2;
        default:
            warn("Unhandled format to number elements: %d\n", x);
            return 0;
    }
}

static size_t formatToSize(CALformat x)
{
    switch (x)
    {
        case formatReal1:
        case constantFormatReal1:
            return sizeof(real);
        case formatReal2:
        case constantFormatReal2:
            return 2 * sizeof(real);
        default:
            warn("Unhandled format to size: %d\n", x);
            return 0;
    }
}

CALresult createConstantBuffer1D(MWMemRes* mr,
                                 MWCALInfo* ci,
                                 const CALvoid* src,
                                 CALformat format,
                                 CALuint width,
                                 CALuint flags)
{
    CALresult err;
    CALvoid* bufPtr;
    CALuint pitch;

    /* Create buffer */
    err = calResAllocLocal1D(&mr->res, ci->dev,
                             formatToSize(format) * formatToNumElements(format) * width,
                             CAL_FORMAT_UINT_1, flags);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to create constant 1D resource", err);
        releaseMWMemRes(ci->calctx, mr);
        return err;
    }

    /* Map and write to the buffer */
    err = mapMWMemRes(mr, &bufPtr, &pitch);
    if (err != CAL_RESULT_OK)
    {
        releaseMWMemRes(ci->calctx, mr);
        return err;
    }

    memcpy(bufPtr, src, formatToSize(format) * width);

    err = unmapMWMemRes(mr);
    if (err != CAL_RESULT_OK)
    {
        releaseMWMemRes(ci->calctx, mr);
        return err;
    }

    err = getMemoryHandle(mr, ci);
    if (err != CAL_RESULT_OK)
    {
        releaseMWMemRes(ci->calctx, mr);
        return err;
    }

    return CAL_RESULT_OK;
}

static CALresult zeroBuffer(MWMemRes* mr, CALuint numberElements, CALuint width, CALuint height)
{
    CALresult err;
    CALdouble* ptr;
    CALuint pitch;

    err = mapMWMemRes(mr, (CALvoid**) &ptr, &pitch);
    if (err != CAL_RESULT_OK)
        return err;

    memset(ptr, 0, width * height * numberElements * sizeof(real));

    err = unmapMWMemRes(mr);
    if (err != CAL_RESULT_OK)
        return err;

    return CAL_RESULT_OK;
}

/* Output appropriate for width * height real1 elements */
static CALresult createOutputBuffer2D(MWMemRes* mr, MWCALInfo* ci, CALuint width, CALuint height)
{
    CALresult err;
    const CALuint nElem = USE_KAHAN ? 2 : 1;

#if USE_KAHAN
    #warning SETUP KAHAN
#else
    #warning SETUP NO KAHAN
#endif

    /* 2*sizeof(uint1) = sizeof(double) */
    err = calResAllocLocal1D(&mr->res, ci->dev, 2 * nElem * width * height, CAL_FORMAT_UINT_1, CAL_RESALLOC_GLOBAL_BUFFER);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to create output resource", err);
        releaseMWMemRes(ci->calctx, mr);
        return err;
    }

    err = zeroBuffer(mr, nElem, width, height);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to zero output buffer", err);
        return err;
    }

    /* Get the handle for the context */
    err = getMemoryHandle(mr, ci);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to create handle for output buffer", err);
        return err;
    }

    return CAL_RESULT_OK;
}

static CALresult createOutBgBuffer(MWCALInfo* ci,
                                   SeparationCALMem* cm,
                                   const CALSeparationSizes* sizes)
{
    CALresult err;

    err = createOutputBuffer2D(&cm->outBg, ci, sizes->rSteps, sizes->muSteps);
    if (err != CAL_RESULT_OK)
        cal_warn("Failed to create output buffer", err);

    return err;
}

/* Create a separate output buffer for each stream */
static CALresult createOutStreamBuffers(MWCALInfo* ci,
                                        SeparationCALMem* cm,
                                        const CALSeparationSizes* sizes)
{
    CALresult err = CAL_RESULT_OK;

    err = createOutputBuffer2D(&cm->outStreams, ci, sizes->rSteps, cm->numberStreams * sizes->muSteps);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to create out streams buffer", err);
    }

    return err;
}

static CALresult createRBuffers(MWCALInfo* ci,
                                SeparationCALMem* cm,
                                const AstronomyParameters* ap,
                                const IntegralArea* ia,
                                const StreamGauss sg)
{
    RPoints* r_pts;
    RConsts* rc;
    CALresult err = CAL_RESULT_OK;

    r_pts = precalculateRPts(ap, ia, sg, &rc, FALSE);

    err = createConstantBuffer1D(&cm->rPts, ci, (CALdouble*) r_pts,
                                 formatReal2, ap->convolve * ia->r_steps, CAL_RESALLOC_GLOBAL_BUFFER);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to create r_pts buffer", err);
        goto fail;
    }

    err = createConstantBuffer1D(&cm->rc, ci, (CALdouble*) rc,
                                 formatReal2, ia->r_steps, CAL_RESALLOC_GLOBAL_BUFFER);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to create rc buffer", err);
        goto fail;
    }

fail:
    mwFreeA(r_pts);
    mwFreeA(rc);

    return err;
}

/* Might be more convenient to split l and b stuff for CAL */
static void getSplitLBTrig(const AstronomyParameters* ap,
                           const IntegralArea* ia,
                           LTrigPair** lTrigBCosOut,
                           real** bTrigOut)
{
    CALuint i, j;
    LTrigPair* lTrigBCos;
    real* bTrig;
    LBTrig* lbts;
    size_t idx;
    CALboolean transpose = CAL_FALSE;

    lTrigBCos = (LTrigPair*) mwMallocA(ia->mu_steps * ia->nu_steps * sizeof(LTrigPair));
    bTrig = (real*) mwMallocA(ia->mu_steps * ia->nu_steps * sizeof(real));

    lbts = precalculateLBTrig(ap, ia, transpose);

    for (i = 0; i < ia->nu_steps; ++i)
    {
        for (j = 0; j < ia->mu_steps; ++j)
        {
            idx = transpose ? j * ia->nu_steps + i : i * ia->mu_steps + j;

            lTrigBCos[idx].lCosBCos = lbts[idx].lCosBCos;
            lTrigBCos[idx].lSinBCos = lbts[idx].lSinBCos;

            bTrig[idx] = lbts[idx].bSin;
        }
    }

    mwFreeA(lbts);

    *lTrigBCosOut = lTrigBCos;
    *bTrigOut = bTrig;
}

static CALresult createLBTrigBuffers(MWCALInfo* ci,
                                     SeparationCALMem* cm,
                                     const AstronomyParameters* ap,
                                     const IntegralArea* ia)
{
    CALresult err = CAL_RESULT_OK;
    LTrigPair* lTrig;
    real* bTrig;

    getSplitLBTrig(ap, ia, &lTrig, &bTrig);

    err = createConstantBuffer1D(&cm->lTrig, ci, (CALdouble*) lTrig,
                                 formatReal2, ia->nu_steps * ia->mu_steps, CAL_RESALLOC_GLOBAL_BUFFER);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to create l trig buffer", err);
        goto fail;
    }

    err = createConstantBuffer1D(&cm->bTrig, ci, (CALdouble*) bTrig,
                                 formatReal1, ia->nu_steps * ia->mu_steps, CAL_RESALLOC_GLOBAL_BUFFER);
    if (err != CAL_RESULT_OK)
        cal_warn("Failed to create b trig buffer", err);

fail:
    mwFreeA(lTrig);
    mwFreeA(bTrig);

    return err;
}

CALresult releaseSeparationBuffers(MWCALInfo* ci, SeparationCALMem* cm)
{
    CALresult err = CAL_RESULT_OK;

    err |= releaseMWMemRes(ci->calctx, &cm->outBg);
    err |= releaseMWMemRes(ci->calctx, &cm->outStreams);
    err |= releaseMWMemRes(ci->calctx, &cm->rPts);
    err |= releaseMWMemRes(ci->calctx, &cm->rc);
    err |= releaseMWMemRes(ci->calctx, &cm->lTrig);
    err |= releaseMWMemRes(ci->calctx, &cm->bTrig);
    err |= releaseMWMemRes(ci->calctx, &cm->nuBuf);
    err |= releaseMWMemRes(ci->calctx, &cm->starsXY);
    err |= releaseMWMemRes(ci->calctx, &cm->starsZ);
    err |= releaseMWMemRes(ci->calctx, &cm->consts);
    err |= releaseMWMemRes(ci->calctx, &cm->streamConsts);

    if (err != CAL_RESULT_OK)
        cal_warn("Failed to release buffers", err);

    return err;
}

CALresult createStreamConstantBuffers(MWCALInfo* ci, SeparationCALMem* cm, const StreamConstants* sc)
{
    CALdouble* buf;
    CALresult err = CAL_RESULT_OK;
    CALuint i;

    buf = mwCallocA(5 * 8, sizeof(real));

    /* Need 7 constants per stream + 1 extra to round to 8 */
    for (i = 0; i < cm->numberStreams; ++i)
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

    err = createConstantBuffer1D(&cm->streamConsts, ci, buf, formatReal2, 4 * 5, CAL_RESALLOC_GLOBAL_BUFFER);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to create rc buffer", err);
    }

    mwFreeA(buf);

    return err;
}

CALresult createSeparationBuffers(MWCALInfo* ci,
                                  SeparationCALMem* cm,
                                  const AstronomyParameters* ap,
                                  const IntegralArea* ia,
                                  const StreamConstants* sc,
                                  const StreamGauss sg,
                                  const CALSeparationSizes* sizes)
{
    CALresult err = CAL_RESULT_OK;

    cm->numberStreams = ap->number_streams;

    err |= createOutBgBuffer(ci, cm, sizes);
    err |= createOutStreamBuffers(ci, cm, sizes);

    err |= createRBuffers(ci, cm, ap, ia, sg);
    err |= createLBTrigBuffers(ci, cm, ap, ia);

    err |= createNuCB(&cm->nuBuf, ci);

    err |= createStreamConstantBuffers(ci, cm, sc);
    err |= createConstsCB(&cm->consts, ci);
    err |= setConstants(cm, ap, ia, sc);

    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to create buffers", err);
        releaseSeparationBuffers(ci, cm);
    }

    return err;
}


static CALboolean checkDeviceCapabilities(const struct CALdeviceattribsRec* attrs)
{

  #if DOUBLEPREC
    if (!attrs->doublePrecision)
    {
        warn("Device does not support double precision\n");
        return CAL_FALSE;
    }
  #endif

    /* TODO: Memory */

    return CAL_TRUE;
}

static void printCALInfo(const MWCALInfo* ci)
{
    warn("Found %u CAL devices\n"
         "Chose device %u\n"
         "\n"
         "Device target:         %s\n"
         "Revision:              %u\n"
         "CAL Version:           %u.%u.%u\n"
         "Engine clock:          %u Mhz\n"
         "Memory clock:          %u Mhz\n"
         "GPU RAM:               %u\n"
         "Wavefront size:        %u\n"
         "Double precision:      %s\n"
         "Compute shader:        %s\n"
         "Number SIMD:           %u\n"
         "Number shader engines: %u\n"
         "Pitch alignment:       %u\n"
         "Surface alignment:     %u\n"
         "Max size 2D:           { %u, %u }\n"
         "\n"
         ,
         ci->numDevices,
         ci->devID,
         showCALtargetEnum(ci->devInfo.target),
         ci->devAttribs.targetRevision,
         ci->version.major, ci->version.minor, ci->version.patchLevel,
         ci->devAttribs.engineClock,
         ci->devAttribs.memoryClock,
         ci->devAttribs.localRAM,
         ci->devAttribs.wavefrontSize,
         showCALboolean(ci->devAttribs.doublePrecision),
         showCALboolean(ci->devAttribs.computeShader),
         ci->devAttribs.numberOfSIMD,
         ci->devAttribs.numberOfShaderEngines,
         ci->devAttribs.pitch_alignment,
         ci->devAttribs.surface_alignment,
         ci->devInfo.maxResource2DWidth, ci->devInfo.maxResource2DHeight
        );
}

static CALobject createCALBinary(const char* srcIL, CALtarget target)
{
    CALresult err;
    CALobject obj = NULL;

    if (!srcIL)
        return NULL;

    err = calclCompile(&obj, CAL_LANGUAGE_IL, srcIL, target);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Error compiling kernel for target '%s'", err, showCALtargetEnum(target));
        return NULL;
    }

    return obj;
}

static CALimage createCALImage(const char* src, CALtarget target)
{
    CALobject obj;
    CALresult rc;
    CALimage img;

    if (!src)
        return NULL;

    obj = createCALBinary(src, target);
    if (!obj)
        return NULL;

    rc = calclLink(&img, &obj, 1);
    calclFreeObject(obj);
    if (rc != CAL_RESULT_OK)
    {
        warn("Error linking image (%d) : %s\n", rc, calclGetErrorString());
        return NULL;
    }

    return img;
}

static CALimage createCALImageFromFile(const char* filename, CALtarget target)
{
    char* src;
    CALimage img;

    src = mwReadFile(filename);
    if (!src)
    {
        perror("IL source file");
        return NULL;
    }

    img = createCALImage(src, target);
    free(src);

    return img;
}

static CALresult printISAToFile(const char* filename, CALimage img)
{
    CALresult err;
    FILE* srcLog;

    srcLog = fopen(filename, "w");
    if (!srcLog)
    {
        perror("ISA output file");
        return CAL_RESULT_ERROR;
    }

    err = printISA(srcLog, img);
    fclose(srcLog);

    return err;
}

#define IL_INTEGRAL_OUT_FILE "calpp_integral_kernel.il"
#define ISA_INTEGRAL_OUT_FILE "calpp_integral_kernel.isa"

#ifndef NDEBUG
static void writeDebugKernels(const char* src, CALimage img, const char* ilOut, const char* isaOut)
{
    char buf[1024];

    mwWriteFile(ilOut, src);
    if (printISAToFile(isaOut, img) == CAL_RESULT_OK)
    {
        sprintf(buf, "grep GPR \"%s\"", isaOut);
        system(buf);
    }
}
#endif /* NDEBUG */


static CALimage createCALImageFromGeneratedKernel(const MWCALInfo* ci,
                                                  const AstronomyParameters* ap,
                                                  const StreamConstants* sc)

{
    CALimage img;
    char* src;

    //src = separationIntegralKernelSrc(ci->devID, IL_MAX_STREAMS);
    src = separationIntegralKernelSrc(ci->devID, ap->number_streams);
    img = createCALImage(src, ci->devInfo.target);

  #ifndef NDEBUG
    writeDebugKernels(src, img, IL_INTEGRAL_OUT_FILE, ISA_INTEGRAL_OUT_FILE);
  #endif /* NDEBUG */

    free(src);

    return img;
}


static inline CALresult getNameMWCALInfo(MWCALInfo* ci, CALname* name, const CALchar* varName)
{
    return calModuleGetName(name, ci->calctx, ci->module, varName);
}

CALresult getModuleNames(MWCALInfo* ci, SeparationCALNames* cn, CALuint numberStreams)
{
    CALresult err = CAL_RESULT_OK;
    CALuint i;
    char buf[20] = "";

    err |= getNameMWCALInfo(ci, &cn->nuBuf, "cb0");
    err |= getNameMWCALInfo(ci, &cn->consts, "cb1");
    err |= getNameMWCALInfo(ci, &cn->streamConsts, "cb2");

    err |= getNameMWCALInfo(ci, &cn->outBg, "uav0");
    err |= getNameMWCALInfo(ci, &cn->outStreams, "uav1");
    err |= getNameMWCALInfo(ci, &cn->rPts,  "uav2");
    err |= getNameMWCALInfo(ci, &cn->rc,    "uav3");
    err |= getNameMWCALInfo(ci, &cn->lTrig, "uav4");
    err |= getNameMWCALInfo(ci, &cn->bTrig, "uav5");

    if (err != CAL_RESULT_OK)
        cal_warn("Failed to get module names", err);

    return err;
}

CALresult setKernelArguments(MWCALInfo* ci, SeparationCALMem* cm, SeparationCALNames* cn)
{
    CALresult err = CAL_RESULT_OK;
    CALuint i;


    err |= calCtxSetMem(ci->calctx, cn->nuBuf, cm->nuBuf.mem);
    if (err != CAL_RESULT_OK)
        cal_warn("nubuf", err);

    err |= calCtxSetMem(ci->calctx, cn->consts, cm->consts.mem);
    if (err != CAL_RESULT_OK)
        cal_warn("consts", err);

    err |= calCtxSetMem(ci->calctx, cn->streamConsts, cm->streamConsts.mem);
    if (err != CAL_RESULT_OK)
        cal_warn("STREAM consts", err);

    if (err != CAL_RESULT_OK)
        cal_warn("OK?", err);

    err |= calCtxSetMem(ci->calctx, cn->outBg, cm->outBg.mem);
    if (err != CAL_RESULT_OK)
        cal_warn("SHOOP", err);

    err |= calCtxSetMem(ci->calctx, cn->outStreams, cm->outStreams.mem);
    if (err != CAL_RESULT_OK)
        cal_warn("DA", err);


    err |= calCtxSetMem(ci->calctx, cn->rPts, cm->rPts.mem);

    if (err != CAL_RESULT_OK)
        cal_warn("WOOP", err);


    err |= calCtxSetMem(ci->calctx, cn->rc, cm->rc.mem);
    if (err != CAL_RESULT_OK)
        cal_warn("OK 1", err);

    err |= calCtxSetMem(ci->calctx, cn->lTrig, cm->lTrig.mem);
    err |= calCtxSetMem(ci->calctx, cn->bTrig, cm->bTrig.mem);
    if (err != CAL_RESULT_OK)
        cal_warn("OK 2", err);


    if (err != CAL_RESULT_OK)
        cal_warn("Failed to set kernel arguments", err);

    return err;
}

CALresult separationLoadKernel(MWCALInfo* ci,
                               const AstronomyParameters* ap,
                               const StreamConstants* sc)
{
    CALresult err;

    ci->image = createCALImageFromGeneratedKernel(ci, ap, sc);
    if (!ci->image)
    {
        warn("Failed to load image\n");
        return CAL_RESULT_ERROR;
    }

    err = calModuleLoad(&ci->module, ci->calctx, ci->image);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to load module", err);
        return err;
    }

    err = calModuleGetEntry(&ci->func, ci->calctx, ci->module, "main");
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to find main in module", err);
        return err;
    }

    return CAL_RESULT_OK;
}


CALresult mwCALShutdown(MWCALInfo* ci)
{
    CALresult err = CAL_RESULT_OK;

    err |= mwDestroyCALInfo(ci);
    err |= calShutdown();
    if (err != CAL_RESULT_OK)
        cal_warn("Failed to shutdown CAL", err);

    return err;
}


/* Init CAL, check device capabilities, and prepare new kernel from workunit parameters */
CALresult separationCALInit(MWCALInfo* ci, const CLRequest* clr)
{
    CALresult err;

    err = calInit();
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to init CAL", err);
        return err;
    }

    err = calExtSupported(0x8009);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("calCtxWaitForEvents not supported\n", err);
    }
    else
    {
        err = calExtGetProc((CALextproc*) &mw_calCtxWaitForEvents, 0x8009, "calCtxWaitForEvents");
        if (err != CAL_RESULT_OK)
        {
            cal_warn("Error getting calCtxWaitForEvents", err);
        }
    }

    err = mwGetCALInfo(ci, clr->devNum);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to get CAL info", err);
        calShutdown();
        return err;
    }

    printCALInfo(ci);

    if (!checkDeviceCapabilities(&ci->devAttribs))
    {
        warn("Device failed capability check\n");
        mwCALShutdown(ci);
        return CAL_RESULT_ERROR;
    }

    return CAL_RESULT_OK;
}

