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

static CALresult writeConstantBufferDouble(MWMemRes* mr,
                                           const CALdouble* dataPtr,
                                           CALuint numberElements,
                                           CALuint width,
                                           CALuint height)
{
    CALuint i, pitch;
    CALdouble* bufPtr;
    size_t rowWidth;
    CALresult err = CAL_RESULT_OK;

    err = mapMWMemRes(mr, (CALvoid**) &bufPtr, &pitch);
    if (err != CAL_RESULT_OK)
        return err;

    rowWidth = sizeof(real) * numberElements * width;
    for (i = 0; i < height; ++i)
    {
        memcpy(&bufPtr[i * numberElements * pitch],
               &dataPtr[i * numberElements * width],
               rowWidth);
    }

    return unmapMWMemRes(mr);
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

static CALresult createConstantBuffer2D(MWMemRes* mr,
                                        MWCALInfo* ci,
                                        const CALdouble* dataPtr,
                                        CALformat format,
                                        CALboolean linearBuffer,
                                        CALuint width,
                                        CALuint height)
{
    CALresult err;
    CALuint flags = linearBuffer ? CAL_RESALLOC_GLOBAL_BUFFER : 0;

    err = calResAllocLocal2D(&mr->res, ci->dev, width, height, format, flags);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to create 2D constant resource", err);
        releaseMWMemRes(ci->calctx, mr);
        return err;
    }

    err = writeConstantBufferDouble(mr, dataPtr, formatToNumElements(format), width, height);
    if (err != CAL_RESULT_OK)
        return err;

    err = getMemoryHandle(mr, ci);
    if (err != CAL_RESULT_OK)
        return err;

    return CAL_RESULT_OK;
}


CALresult createConstantBuffer1D(MWMemRes* mr,
                                 MWCALInfo* ci,
                                 const CALvoid* src,
                                 CALformat format,
                                 CALuint width)
{
    CALresult err;
    CALvoid* bufPtr;
    CALuint pitch;

    /* Create buffer */
    err = calResAllocLocal1D(&mr->res, ci->dev, width, format, 0);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to create constant 1D resource", err);
        releaseMWMemRes(ci->calctx, mr);
        return err;
    }

    /* Map and write to the buffer */
    err = mapMWMemRes(mr, &bufPtr, &pitch);
    if (err != CAL_RESULT_OK)
        return err;

    memcpy(bufPtr, src, formatToSize(format) * width);

    err = unmapMWMemRes(mr);
    if (err != CAL_RESULT_OK)
        return err;

    err = getMemoryHandle(mr, ci);
    if (err != CAL_RESULT_OK)
        return err;

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

    memset(ptr, 0, height * pitch * numberElements * sizeof(real));

    err = unmapMWMemRes(mr);
    if (err != CAL_RESULT_OK)
        return err;

    return CAL_RESULT_OK;
}

/* Output appropriate for width * height real1 elements */
static CALresult createOutputBuffer2D(MWMemRes* mr, MWCALInfo* ci, CALuint width, CALuint height)
{
    CALresult err;
    CALuint flags = 0;
    //CALuint flags = CAL_RESALLOC_GLOBAL_BUFFER;

    err = calResAllocLocal2D(&mr->res, ci->dev, width, height, formatReal2, flags);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to create output resource", err);
        releaseMWMemRes(ci->calctx, mr);
        return err;
    }

    err = zeroBuffer(mr, 2, width, height);
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
    CALuint i, numberAllocated;
    CALresult err = CAL_RESULT_OK;

    cm->outStreams = (MWMemRes*) mwCalloc(cm->numberStreams, sizeof(MWMemRes));
    for (i = 0; i < cm->numberStreams; ++i)
    {
        err = createOutputBuffer2D(&cm->outStreams[i], ci, sizes->rSteps, sizes->muSteps);
        if (err != CAL_RESULT_OK)
        {
            cal_warn("Failed to create out streams buffer", err);
            break;
        }
    }

    /* Clean up if any failed before */
    numberAllocated = i;
    if (numberAllocated < cm->numberStreams)
    {
        for (i = 0; i < numberAllocated; ++i)
            releaseMWMemRes(ci->calctx, &cm->outStreams[i]);
        free(cm->outStreams);
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

    r_pts = precalculateRPts(ap, ia, sg, &rc, TRUE);

    err = createConstantBuffer2D(&cm->rPts, ci, (CALdouble*) r_pts, formatReal2, CAL_FALSE, ia->r_steps, ap->convolve);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to create r_pts buffer", err);
        goto fail;
    }

    err = createConstantBuffer1D(&cm->rc, ci, (CALdouble*) rc, formatReal2, ia->r_steps);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to create rc buffer", err);
        goto fail;
    }

    err = createConstantBuffer1D(&cm->sg_dx, ci, sg.dx, constantFormatReal1, ap->convolve);
    if (err != CAL_RESULT_OK)
        cal_warn("Failed to create sg_dx buffer", err);

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
    CALboolean transpose = CAL_TRUE;

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

    err = createConstantBuffer2D(&cm->lTrig, ci, (CALdouble*) lTrig, formatReal2, CAL_FALSE, ia->nu_steps, ia->mu_steps);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to create l trig buffer", err);
        goto fail;
    }

    err = createConstantBuffer2D(&cm->bTrig, ci, (CALdouble*) bTrig, formatReal1, CAL_FALSE, ia->nu_steps, ia->mu_steps);
    if (err != CAL_RESULT_OK)
        cal_warn("Failed to create b trig buffer", err);

fail:
    mwFreeA(lTrig);
    mwFreeA(bTrig);

    return err;
}

static CALresult releaseStreamBuffers(MWCALInfo* ci, SeparationCALMem* cm)
{
    CALuint i;
    CALresult err = CAL_RESULT_OK;

    if (!cm->outStreams) /* Nothing to free */
        return err;

    for (i = 0; i < cm->numberStreams; ++i)
        err |= releaseMWMemRes(ci->calctx, &cm->outStreams[i]);

    free(cm->outStreams);
    cm->outStreams = NULL;

    return err;
}

CALresult releaseSeparationBuffers(MWCALInfo* ci, SeparationCALMem* cm)
{
    CALresult err = CAL_RESULT_OK;

    err |= releaseMWMemRes(ci->calctx, &cm->outBg);
    err |= releaseStreamBuffers(ci, cm);
    err |= releaseMWMemRes(ci->calctx, &cm->rPts);
    err |= releaseMWMemRes(ci->calctx, &cm->rc);
    err |= releaseMWMemRes(ci->calctx, &cm->sg_dx);
    err |= releaseMWMemRes(ci->calctx, &cm->lTrig);
    err |= releaseMWMemRes(ci->calctx, &cm->bTrig);
    err |= releaseMWMemRes(ci->calctx, &cm->nuBuf);

    err |= releaseMWMemRes(ci->calctx, &cm->sg_qgauss_W);
    err |= releaseMWMemRes(ci->calctx, &cm->starsXY);
    err |= releaseMWMemRes(ci->calctx, &cm->starsZ);

    if (err != CAL_RESULT_OK)
        cal_warn("Failed to release buffers", err);

    return err;
}

CALresult createSeparationBuffers(MWCALInfo* ci,
                                  SeparationCALMem* cm,
                                  const AstronomyParameters* ap,
                                  const IntegralArea* ia,
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

    src = separationIntegralKernelSrc(ap, sc, ci->devID);
    img = createCALImage(src, ci->devInfo.target);
    free(src);

  #ifndef NDEBUG
    writeDebugKernels(src, img, IL_INTEGRAL_OUT_FILE, ISA_INTEGRAL_OUT_FILE);
  #endif /* NDEBUG */

    return img;
}


static inline CALresult getNameMWCALInfo(MWCALInfo* ci, CALname* name, const CALchar* varName)
{
    return calModuleGetName(name, ci->calctx, ci->module, varName);
}

void destroyModuleNames(SeparationCALNames* cn)
{
    free(cn->outStreams);
    cn->outStreams = NULL;

    free(cn->inStreams);
    cn->inStreams = NULL;
}

CALresult getModuleNames(MWCALInfo* ci, SeparationCALNames* cn, CALuint numberStreams)
{
    CALresult err = CAL_RESULT_OK;
    CALuint i;
    char buf[20] = "";

    cn->outStreams = (CALname*) mwCalloc(numberStreams, sizeof(CALname));
    cn->inStreams = (CALname*) mwCalloc(numberStreams, sizeof(CALname));

    err |= getNameMWCALInfo(ci, &cn->nuBuf, "cb0");
    err |= getNameMWCALInfo(ci, &cn->sg_dx, "cb1");

    err |= getNameMWCALInfo(ci, &cn->rPts,  "i0");
    err |= getNameMWCALInfo(ci, &cn->rc,    "i1");
    err |= getNameMWCALInfo(ci, &cn->lTrig, "i2");
    err |= getNameMWCALInfo(ci, &cn->bTrig, "i3");

    err |= getNameMWCALInfo(ci, &cn->outBg, "o0");
    err |= getNameMWCALInfo(ci, &cn->inMu, "i4");
    for (i = 0; i < numberStreams; ++i)
    {
        sprintf(buf, "i%u", i + 5);
        err |= getNameMWCALInfo(ci, &cn->inStreams[i], buf);

        sprintf(buf, "o%u", i + 1);
        err |= getNameMWCALInfo(ci, &cn->outStreams[i], buf);
    }

    if (err != CAL_RESULT_OK)
        cal_warn("Failed to get module names", err);

    return err;
}

CALresult setKernelArguments(MWCALInfo* ci, SeparationCALMem* cm, SeparationCALNames* cn)
{
    CALresult err = CAL_RESULT_OK;
    CALuint i;

    err |= calCtxSetMem(ci->calctx, cn->outBg, cm->outBg.mem);
    err |= calCtxSetMem(ci->calctx, cn->inMu, cm->outBg.mem);
    for (i = 0; i < cm->numberStreams; ++i)
    {
        err |= calCtxSetMem(ci->calctx, cn->outStreams[i], cm->outStreams[i].mem);
        err |= calCtxSetMem(ci->calctx, cn->inStreams[i],  cm->outStreams[i].mem);
    }

    err |= calCtxSetMem(ci->calctx, cn->rPts, cm->rPts.mem);
    err |= calCtxSetMem(ci->calctx, cn->rc, cm->rc.mem);
    err |= calCtxSetMem(ci->calctx, cn->lTrig, cm->lTrig.mem);
    err |= calCtxSetMem(ci->calctx, cn->bTrig, cm->bTrig.mem);
    err |= calCtxSetMem(ci->calctx, cn->sg_dx, cm->sg_dx.mem);
    err |= calCtxSetMem(ci->calctx, cn->nuBuf, cm->nuBuf.mem);

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
        return err;
    }

    err = calExtGetProc((CALextproc*) &mw_calCtxWaitForEvents, 0x8009, "calCtxWaitForEvents");
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Error getting calCtxWaitForEvents", err);
        return err;
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

