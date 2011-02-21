/*
Copyright (C) 2010  Matthew Arsenault

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
#include "cal_binary.h"
#include "separation_cal_types.h"
#include "separation_cal_kernelgen.h"


#define cal_warn(msg, err, ...) fprintf(stderr, msg ": %s (%s)\n", ##__VA_ARGS__, calGetErrorString(), showCALresult(err))

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

static CALresult mwDestroyCALInfo(MWCALInfo* ci)
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
        return -1;
    }

    if (devID > ci->numDevices)
    {
        warn("Requested device ID %u > found number of devices (%u)\n",
             devID,
             ci->numDevices);
        return -1;
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

#if DOUBLEPREC
  /* VERY IMPORTANT:
     http://developer.amd.com/support/KnowledgeBase/Lists/KnowledgeBase/DispForm.aspx?ID=92
   */
  /* For some reason it doesn't work if you try to use the uint32 ones
   * with a cb[]. */
  #define constantFormatReal1 CAL_FORMAT_FLOAT_2
  #define constantFormatReal2 CAL_FORMAT_FLOAT_4

  #define formatReal1 CAL_FORMAT_UNSIGNED_INT32_2
  #define formatReal2 CAL_FORMAT_UNSIGNED_INT32_4
#else
  #define constantFormatReal1 CAL_FORMAT_FLOAT_2
  #define constantFormatReal2 CAL_FORMAT_FLOAT_4

  #define formatReal1 CAL_FORMAT_UNSIGNED_INT32_1
  #define formatReal2 CAL_FORMAT_UNSIGNED_INT32_2
#endif /* DOUBLEPREC */

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
static CALresult mapMWMemRes(MWMemRes* mr, CALvoid** pPtr, CALuint* pitch)
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
static CALresult unmapMWMemRes(MWMemRes* mr)
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
               &dataPtr[i * width * numberElements],
               rowWidth);
    }

    return unmapMWMemRes(mr);
}

static CALresult printBufferDouble(const char* name,
                                   MWMemRes* mr,
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

    err = calResAllocRemote1D(&mr->res, &ci->dev, 1, 1, constantFormatReal2, 0);
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
                                        CALuint width,
                                        CALuint height)
{
    CALresult err;

    err = calResAllocLocal2D(&mr->res, ci->dev, width, height, format, 0);
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


static CALresult createConstantBuffer1D(MWMemRes* mr,
                                        MWCALInfo* ci,
                                        SeparationCALMem* cm,
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
    size_t rowWidth;
    CALuint i, pitch;

    err = mapMWMemRes(mr, (CALvoid**) &ptr, &pitch);
    if (err != CAL_RESULT_OK)
        return err;

    rowWidth = sizeof(real) * width * numberElements;
    for (i = 0; i < height; ++i)
        memset(&ptr[i * numberElements * pitch], 0, rowWidth);

    err = unmapMWMemRes(mr);
    if (err != CAL_RESULT_OK)
        return err;

    return CAL_RESULT_OK;
}

/* Output appropriate for width * height real1 elements */
static CALresult createOutputBuffer2D(MWMemRes* mr, MWCALInfo* ci, CALuint width, CALuint height)
{
    CALresult err;

    err = calResAllocLocal2D(&mr->res, ci->dev, width, height, formatReal2, 0);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to create output resource", err);
        releaseMWMemRes(ci->calctx, mr);
        return err;
    }

    /* Get the handle for the context */
    err = getMemoryHandle(mr, ci);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to create handle for output buffer", err);
        return err;
    }

    err = zeroBuffer(mr, 2, width, height);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to zero output buffer", err);
        return err;
    }

    return CAL_RESULT_OK;
}

static CALresult createOutMuBuffer(MWCALInfo* ci,
                                   SeparationCALMem* cm,
                                   const CALSeparationSizes* sizes)
{
    CALresult err;

    err = createOutputBuffer2D(&cm->outBg, ci, sizes->muSteps, sizes->rSteps);
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
        err = createOutputBuffer2D(&cm->outStreams[i], ci, sizes->muSteps, sizes->rSteps);
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
                                const StreamGauss sg,
                                const CALSeparationSizes* sizes)
{
    RPoints* r_pts;
    RConsts* rc;
    CALresult err = CAL_RESULT_OK;

    r_pts = precalculateRPts(ap, ia, sg, &rc, FALSE);

    err = createConstantBuffer2D(&cm->rPts, ci, (CALdouble*) r_pts, formatReal2, ap->convolve, ia->r_steps);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to create r_pts buffer", err);
        goto fail;
    }

    err = createConstantBuffer1D(&cm->rc, ci, cm, (CALdouble*) rc, formatReal2, ia->r_steps);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to create rc buffer", err);
        goto fail;
    }

    err = createConstantBuffer1D(&cm->sg_dx, ci, cm, sg.dx, constantFormatReal1, ap->convolve);
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
                                     const IntegralArea* ia,
                                     const CALSeparationSizes* sizes)
{
    CALresult err = CAL_RESULT_OK;
    LTrigPair* lTrig;
    real* bTrig;

    getSplitLBTrig(ap, ia, &lTrig, &bTrig);

    err = createConstantBuffer2D(&cm->lTrig, ci, (CALdouble*) lTrig, formatReal2, ia->nu_steps, ia->mu_steps);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to create l trig buffer", err);
        goto fail;
    }

    err = createConstantBuffer2D(&cm->bTrig, ci, (CALdouble*) bTrig, formatReal1, ia->nu_steps, ia->mu_steps);
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

    if (err != CAL_RESULT_OK)
        cal_warn("Failed to release buffers", err);

    return err;
}

static CALresult createSeparationBuffers(MWCALInfo* ci,
                                         SeparationCALMem* cm,
                                         const AstronomyParameters* ap,
                                         const IntegralArea* ia,
                                         const StreamConstants* sc,
                                         const StreamGauss sg,
                                         const CALSeparationSizes* sizes)
{
    CALresult err = CAL_RESULT_OK;

    cm->numberStreams = ap->number_streams;

    err |= createOutMuBuffer(ci, cm, sizes);
    err |= createOutStreamBuffers(ci, cm, sizes);

    err |= createRBuffers(ci, cm, ap, ia, sg, sizes);
    err |= createLBTrigBuffers(ci, cm, ap, ia, sizes);

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
         "Compute shader:        %s\n"
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
         showCALboolean(ci->devAttribs.computeShader),

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

static CALobject createCALBinary(const char* srcIL)
{
    CALresult err;
    CALobject obj = NULL;

    if (!srcIL)
        return NULL;

    err = calclCompile(&obj, CAL_LANGUAGE_IL, srcIL, CAL_TARGET_CYPRESS);
    if (err != CAL_RESULT_OK)
    {
        warn("Error compiling kernel (%d) : %s\n", err, calclGetErrorString());
        return NULL;
    }

    return obj;
}

static CALimage createCALImage(const char* src)
{
    CALobject obj;
    CALresult rc;
    CALimage img;

    if (!src)
        return NULL;

    obj = createCALBinary(src);
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

static CALimage createCALImageFromFile(const char* filename)
{
    char* src;
    CALimage img;

    src = mwReadFile(filename);
    if (!src)
    {
        perror("IL source file");
        return NULL;
    }

    img = createCALImage(src);
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

static CALimage createCALImageFromGeneratedKernel(const AstronomyParameters* ap,
                                                  const IntegralArea* ia,
                                                  const StreamConstants* sc)

{
    CALimage img;
    char* src;
    char buf[512];
    const char* ilFile = "calpp_kernel.il";
    const char* isaFile = "calpp_kernel_Cypress.isa";

    src = separationKernelSrc(ap, ia, sc);
    mwWriteFile(ilFile, src);

    img = createCALImage(src);
    free(src);

    if (printISAToFile(isaFile, img) == CAL_RESULT_OK)
    {
        sprintf(buf, "grep GPR %s", isaFile);
        system(buf);
    }

    return img;
}


static inline CALresult getNameMWCALInfo(MWCALInfo* ci, CALname* name, const CALchar* varName)
{
    return calModuleGetName(name, ci->calctx, ci->module, varName);
}

static void destroyModuleNames(SeparationCALNames* cn)
{
    free(cn->outStreams);
    cn->outStreams = NULL;

    free(cn->inStreams);
    cn->inStreams = NULL;
}

static CALresult getModuleNames(MWCALInfo* ci, SeparationCALNames* cn, CALuint numberStreams)
{
    CALresult err = CAL_RESULT_OK;
    CALuint i;
    char buf[20] = "";

    cn->outStreams = mwCalloc(numberStreams, sizeof(CALname));
    cn->inStreams = mwCalloc(numberStreams, sizeof(CALname));

    err |= getNameMWCALInfo(ci, &cn->sg_dx, "cb0");
    err |= getNameMWCALInfo(ci, &cn->nuBuf, "cb1");

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

static CALresult setKernelArguments(MWCALInfo* ci, SeparationCALMem* cm, SeparationCALNames* cn)
{
    CALresult err = CAL_RESULT_OK;
    CALuint i;

    /* CHECKME: Bind the same output buffer to the input OK? */
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

static CALresult separationSetupCAL(MWCALInfo* ci,
                                    const AstronomyParameters* ap,
                                    const IntegralArea* ia,
                                    const StreamConstants* sc)
{
    CALresult err;

    ci->image = createCALImageFromGeneratedKernel(ap, ia, sc);
    if (!ci->image)
    {
        warn("Failed to load image\n");
        return -1;
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

static CALresult runKernel(MWCALInfo* ci, SeparationCALMem* cm, const IntegralArea* ia)
{
    CALresult err;
    CALevent ev = 0;
    CALdomain domain = { 0, 0, ia->mu_steps, ia->r_steps };

    err = calCtxRunProgram(&ev, ci->calctx, ci->func, &domain);

#if 0
    CALdomain3D global = { ia->mu_steps, ia->r_steps, 1 };
    //CALdomain3D local = { 400, ia->r_steps, 1 };
    //CALdomain3D local = { ia->mu_steps, 2, 1 };
    //CALdomain3D local = { 400, ia->r_steps, 1 };
    CALdomain3D local = { 64, 28, 1 };
    CALprogramGrid grid;

    grid.func = ci->func;
    grid.flags = 0;

    grid.gridBlock = local;

    grid.gridSize.width  = (global.width + local.width - 1) / local.width;
    grid.gridSize.height = (global.height + local.height - 1) / local.height;
    grid.gridSize.depth  = (global.depth + local.depth - 1) / local.depth;

    warn("arst %u %u %u -> { %u %u }\n", grid.gridSize.width, grid.gridSize.height, grid.gridSize.depth,

         grid.gridSize.width * local.width,
         grid.gridSize.height * local.height
        );

    err = calCtxRunProgramGrid(&ev, ci->calctx, &grid);
#endif

    if (err != CAL_RESULT_OK)
    {
        cal_warn("Error running kernel", err);
        return err;
    }

    while (calCtxIsEventDone(ci->calctx, ev) == CAL_RESULT_PENDING);

    return CAL_RESULT_OK;
}
static CALresult setNuKernelArgs(MWCALInfo* ci,
                                 SeparationCALMem* cm,
                                 SeparationCALNames* cn,
                                 const IntegralArea* ia,
                                 CALuint nuStep)
{
    CALresult err;
    CALdouble* nuBufPtr;
    CALfloat* nuStepPtr;
    CALuint pitch = 0;
    NuId nuid;

    err = mapMWMemRes(&cm->nuBuf, (CALvoid**) &nuBufPtr, &pitch);
    if (err != CAL_RESULT_OK)
        return err;

    nuid = calcNuStep(ia, nuStep);

    nuStepPtr = (CALfloat*) nuBufPtr;
    *nuStepPtr = (CALfloat) nuStep;
    nuBufPtr[1] = nuid.id;

    err = unmapMWMemRes(&cm->nuBuf);
    if (err != CAL_RESULT_OK)
        return err;

    return CAL_RESULT_OK;
}

static real sumResults(MWMemRes* mr, const IntegralArea* ia)
{
    CALuint i, j, pitch;
    Kahan* bufPtr;
    Kahan* tmp;
    Kahan ksum = ZERO_KAHAN;
    CALresult err = CAL_RESULT_OK;

    err = mapMWMemRes(mr, (CALvoid**) &bufPtr, &pitch);
    if (err != CAL_RESULT_OK)
        return NAN;

    for (i = 0; i < ia->r_steps; ++i)
    {
        tmp = &bufPtr[i * pitch];
        for (j = 0; j < ia->mu_steps; ++j)
        {
            KAHAN_ADD(ksum, tmp[j].sum);
        }
    }

    err = unmapMWMemRes(mr);
    if (err != CAL_RESULT_OK)
        return NAN;

    return ksum.sum + ksum.correction;
}

static real readResults(MWCALInfo* ci,
                        SeparationCALMem* cm,
                        const IntegralArea* ia,
                        real* probs_results,
                        CALuint numberStreams)
{
    CALuint i;
    real result;

    result = sumResults(&cm->outBg, ia);

    #if 1
    for (i = 0; i < numberStreams; ++i)
        probs_results[i] = sumResults(&cm->outStreams[i], ia);
    #endif

    return result;
}


static real runIntegral(MWCALInfo* ci,
                        SeparationCALMem* cm,
                        const IntegralArea* ia,
                        real* probs_results)
{
    CALresult err;
    unsigned int i;
    SeparationCALNames cn;
    double t1, t2;

    memset(&cn, 0, sizeof(SeparationCALNames));
    err = getModuleNames(ci, &cn, cm->numberStreams);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to get module names", err);
        return NAN;
    }

    err = setKernelArguments(ci, cm, &cn);
    if (err != CAL_RESULT_OK)
    {
        destroyModuleNames(&cn);
        return NAN;
    }

    for (i = 0; i < ia->nu_steps; ++i)
    //for (i = 0; i < 10; ++i)
    {
        warn("Trying to run step: %u\n", i);
        t1 = mwGetTime();
        err = setNuKernelArgs(ci, cm, &cn, ia, i);
        if (err != CAL_RESULT_OK)
            break;

        err = runKernel(ci, cm, ia);
        if (err != CAL_RESULT_OK)
            break;

        t2 = mwGetTime();
        warn("Time = %fms\n", 1000.0 * (t2 - t1));
    }

    destroyModuleNames(&cn);

    return (err != CAL_RESULT_OK) ? NAN : readResults(ci, cm, ia, probs_results, cm->numberStreams);
}

static void calculateCALSeparationSizes(CALSeparationSizes* sizes,
                                        const AstronomyParameters* ap,
                                        const IntegralArea* ia)
{
    sizes->outBg = sizeof(Kahan) * ia->mu_steps * ia->r_steps;
    sizes->outStreams = sizeof(Kahan) * ia->mu_steps * ia->r_steps * ap->number_streams;
    sizes->rPts = sizeof(RPoints) * ap->convolve * ia->r_steps;
    sizes->rc = sizeof(RConsts) * ia->r_steps;
    sizes->sg_dx = sizeof(real) * ap->convolve;
    sizes->lTrig = sizeof(LTrigPair) * ia->mu_steps * ia->nu_steps;
    sizes->bTrig = sizeof(real) * ia->mu_steps * ia->nu_steps;

    sizes->nuSteps = ia->nu_steps;
    sizes->muSteps = ia->mu_steps;
    sizes->rSteps = ia->r_steps;
}

/* Init CAL, and prepare separation kernel from workunit parameters */
CALresult separationCALInit(MWCALInfo* ci,
                            const AstronomyParameters* ap,
                            const IntegralArea* ia,
                            const StreamConstants* sc)
{
    CALresult err;

    err = calInit();
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to init CAL", err);
        return err;
    }

    err = mwGetCALInfo(ci, 0);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to get CAL info", err);
        return err;
    }

    printCALInfo(ci);

    err = separationSetupCAL(ci, ap, ia, sc);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to setup CAL", err);
        mwDestroyCALInfo(ci);
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

real integrateCAL(const AstronomyParameters* ap,
                  const IntegralArea* ia,
                  const StreamConstants* sc,
                  const StreamGauss sg,
                  real* st_probs,
                  EvaluationState* es,
                  const CLRequest* clr,
                  MWCALInfo* ci)
{
    CALresult err;
    SeparationCALMem cm;
    CALSeparationSizes sizes;
    real result = NAN;

    calculateCALSeparationSizes(&sizes, ap, ia);

    memset(&cm, 0, sizeof(SeparationCALMem));
    err = createSeparationBuffers(ci, &cm, ap, ia, sc, sg, &sizes);
    if (err != CAL_RESULT_OK)
        return NAN;

    result = runIntegral(ci, &cm, ia, st_probs);

    err = releaseSeparationBuffers(ci, &cm);
    if (err != CAL_RESULT_OK)
        result = NAN;

    return result;
}

