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

#include <cal.h>
#include <calcl.h>
#include "milkyway_util.h"
#include "separation_types.h"
#include "calculated_constants.h"
#include "show_cal_types.h"
#include "cal_binary.h"

#include "r_points.h"

/* FIXME: Also defined in CL version */
typedef struct
{
    size_t outMu;
    size_t outProbs;

    size_t ap;        /* Constants */
    size_t sc;
    size_t ia;
    size_t ic;
    size_t rc;
    size_t rPts;
    size_t sg_dx;
    size_t lbts;
} SeparationSizes;



typedef struct
{
    CALuint major, minor, patchLevel;
} MWCALVersion;

typedef struct
{
    MWCALVersion version;
    CALuint numDevices;
    CALdevice devID;    /* Index of device chosen */
    CALdevice dev;
    CALdeviceinfo devInfo;
    CALdeviceattribs devAttribs;
    CALcontext calctx;
    CALmodule module;
    CALimage image;
    CALfunc func;
} MWCALInfo;

#define EMPTY_CAL_INFO { 0 }

/* Pair of resource and associated CALmem */
typedef struct
{
    CALresource res;
    CALmem mem;
} MWMemRes;

typedef struct
{
    MWMemRes outMu;
    MWMemRes outProbs;

    /* constant, read only buffers */
    MWMemRes ap;
    MWMemRes ia;
    MWMemRes sc;
    MWMemRes rc;        /* r constants */
    MWMemRes rPts;
    MWMemRes sg_dx;
    MWMemRes lbts;      /* sin, cos of l, b */
} SeparationCALMem;

#define EMPTY_SEPARATION_CAL_MEM { 0 }

#define cal_warn(str, err) fprintf(stderr, str ": %s\n", showCALresult(err))


static CALresult releaseMWMemRes(CALcontext ctx, MWMemRes* mr)
{
    CALresult err;

    err = calCtxReleaseMem(ctx, mr->mem);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to release CALmem", err);
        return err;
    }

    err = calResFree(mr->res);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to release CAL resource", err);
        return err;
    }

    return CAL_RESULT_OK;
}

static CALresult mwDestroyCALInfo(MWCALInfo* ci)
{
    CALresult err;

    err = calCtxDestroy(ci->calctx);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to destroy CAL context", err);
        return err;
    }

    err = calDeviceClose(ci->dev);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to close device", err);
        return err;
    }

    err = calShutdown();
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to shutdown CAL", err);
        return err;
    }

    return CAL_RESULT_OK;
}

static CALresult mwInitCAL(MWCALInfo* ci, CALuint devID)
{
    CALresult err;

    err = calInit();
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to init CAL", err);
        return err;
    }

    err = calGetVersion(&ci->version.major, &ci->version.minor, &ci->version.patchLevel);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to get CAL version", err);
        return err;
    }

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

    err = calDeviceOpen(&ci->dev, ci->devID);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to open CAL device", err);
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
static const CALformat formatReal1 = CAL_FORMAT_FLOAT_1;
static const CALformat formatReal2 = CAL_FORMAT_FLOAT_2;
#else
static const CALformat formatReal1 = CAL_FORMAT_DOUBLE_1;
static const CALformat formatReal2 = CAL_FORMAT_DOUBLE_2;
#endif

/* Try to get memory handle and cleanup resource if that fails */
static CALresult getMemoryHandle(MWMemRes* mr, MWCALInfo* ci)
{
    CALresult err = CAL_RESULT_OK;

    err = calCtxGetMem(&mr->mem, ci->calctx, mr->res);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to get memory handle handle", err);
        if (calResFree(mr->res) != CAL_RESULT_OK)
            warn("Failed to release CAL resource\n");
    }

    return err;
}

/* Try to map the resource and free it on failure */
static CALresult mapMWMemRes(MWMemRes* mr, MWCALInfo* ci, CALvoid** pPtr, CALuint* pitch)
{
    CALresult err;

    err = calResMap(pPtr, pitch, mr->res, 0);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to map ap resource", err);
        if (calResFree(mr->res) != CAL_RESULT_OK)
            warn("Failed to release CAL resource\n");
    }
    return err;
}

/* Try to unmap resource and free it on failure */
static CALresult unmapMWMemRes(MWMemRes* mr, MWCALInfo* ci)
{
    CALresult err;

    err = calResUnmap(mr->res);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to unmap resource", err);
        if (calResFree(mr->res) != CAL_RESULT_OK)
            warn("Failed to release CAL resource\n");
    }

    return err;
}

static CALresult createConstantBuffer(MWMemRes* mr,
                                      MWCALInfo* ci,
                                      SeparationCALMem* cm,
                                      const CALvoid* src,
                                      size_t size)
{
    CALresult err;
    CALvoid* cbPtr;
    CALuint cbPitch;
    CALuint width = size / 4; /* width = totalsize / formatsize */

    if (size % 4 != 0)
        warn("Constant size %zu is not divisible by sizeof(uint32)\n", size);

    /* Create buffer */
    err = calResAllocLocal1D(&mr->res, ci->dev,
                             size, CAL_FORMAT_UNSIGNED_INT32_1, 0);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to create ap resource", err);
        return err;
    }

    /* Map and write to the buffer */
    err = mapMWMemRes(&cm->ap, ci, &cbPtr, &cbPitch);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to map ap resource", err);
        return err;
    }

    memcpy(cbPtr, src, size);

    err = unmapMWMemRes(mr, ci);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to unmap resource", err);
        return err;
    }

    err = getMemoryHandle(mr, ci);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to get memory handle", err);
        return err;
    }

    return CAL_RESULT_OK;
}

static CALresult zeroBuffer1D(MWMemRes* mr, MWCALInfo* ci, size_t size)
{
    CALresult err;
    CALuint pitch;
    CALvoid* ptr;

    err = mapMWMemRes(mr, ci, &ptr, &pitch);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to map ap resource", err);
        return err;
    }

    memset(ptr, 0, size);

    err = unmapMWMemRes(mr, ci);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to unmap resource", err);
        return err;
    }

    return CAL_RESULT_OK;
}

static CALresult createOutputBuffer1D(MWMemRes* mr, MWCALInfo* ci, size_t size)
{
    CALresult err;
    CALuint width = size / 4; /* width = totalsize / formatsize */

    if (size % 4 != 0)
        warn("Constant size %zu is not divisible by sizeof(uint32)\n", size);

    err = calResAllocLocal1D(&mr->res, ci->dev, width, CAL_FORMAT_UNSIGNED_INT32_1, 0);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to create output resource", err);
        return err;
    }

    err = zeroBuffer1D(mr, ci, size);
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

static CALresult createOutMuBuffer(MWCALInfo* ci,
                                   SeparationCALMem* cm,
                                   const SeparationSizes* sizes)
{
    CALresult err;

    err = createOutputBuffer1D(&cm->outMu, ci, sizes->outMu);
    if (err != CAL_RESULT_OK)
        cal_warn("Failed to create output buffer", err);

    return err;
}

static CALresult createOutProbsBuffer(MWCALInfo* ci,
                                      SeparationCALMem* cm,
                                      const SeparationSizes* sizes)
{
    CALresult err;

    err = createOutputBuffer1D(&cm->outProbs, ci, sizes->outProbs);
    if (err != CAL_RESULT_OK)
        cal_warn("Failed to create out probs buffer", err);

    return err;
}

static CALresult createSCBuffer(MWCALInfo* ci,
                                SeparationCALMem* cm,
                                const StreamConstants* sc,
                                const SeparationSizes* sizes)
{
    CALresult err = CAL_RESULT_OK;

    err = createConstantBuffer(&cm->sc, ci, cm, sc, sizes->sc);
    if (err != CAL_RESULT_OK)
        cal_warn("Failed to create sc buffer", err);

    return err;
}

static CALresult createAPBuffer(MWCALInfo* ci,
                                SeparationCALMem* cm,
                                const AstronomyParameters* ap,
                                const SeparationSizes* sizes)
{
    CALresult err = CAL_RESULT_OK;

    err = createConstantBuffer(&cm->ap, ci, cm,
                               ap, sizeof(AstronomyParameters));
    if (err != CAL_RESULT_OK)
        cal_warn("Failed to create ap buffer", err);

    return err;
}

static CALresult createIABuffer(MWCALInfo* ci,
                                SeparationCALMem* cm,
                                const IntegralArea* ia,
                                const SeparationSizes* sizes)
{
    CALresult err = CAL_RESULT_OK;

    err = createConstantBuffer(&cm->ia, ci, cm,
                               ia, sizeof(IntegralArea));
    if (err != CAL_RESULT_OK)
        cal_warn("Failed to create ia buffer", err);

    return err;
}

static CALresult createRBuffers(MWCALInfo* ci,
                                SeparationCALMem* cm,
                                const AstronomyParameters* ap,
                                const IntegralArea* ia,
                                const StreamGauss sg,
                                const SeparationSizes* sizes)
{
    RPoints* r_pts;
    RConsts* rc;
    CALresult err = CAL_RESULT_OK;

    r_pts = precalculateRPts(ap, ia, sg, &rc, TRUE);

    err = createConstantBuffer(&cm->rPts, ci, cm, r_pts, sizes->rPts);
    if (err != CAL_RESULT_OK)
        cal_warn("Failed to create r_pts buffer", err);

    err = createConstantBuffer(&cm->rc, ci, cm, rc, sizes->rc);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to create rc buffer", err);
        releaseMWMemRes(ci->calctx, &cm->rPts);
    }

    err = createConstantBuffer(&cm->sg_dx, ci, cm, sg.dx, sizes->sg_dx);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to create sg_dx buffer", err);
        releaseMWMemRes(ci->calctx, &cm->rPts);
        releaseMWMemRes(ci->calctx, &cm->rc);
    }

    mwFreeA(r_pts);
    mwFreeA(rc);

    return err;
}


static CALresult createLBTrigBuffer(MWCALInfo* ci,
                                    SeparationCALMem* cm,
                                    const AstronomyParameters* ap,
                                    const IntegralArea* ia,
                                    const SeparationSizes* sizes)
{
    CALresult err = CAL_RESULT_OK;
    LBTrig* lbts;

    lbts = precalculateLBTrig(ap, ia, TRUE);
    err = createConstantBuffer(&cm->lbts, ci, cm, lbts, sizes->lbts);
    if (err != CAL_RESULT_OK)
        cal_warn("Failed to create lb trig buffer", err);

    mwFreeA(lbts);

    return err;
}

CALresult createSeparationBuffers(MWCALInfo* ci,
                                  SeparationCALMem* cm,
                                  const AstronomyParameters* ap,
                                  const IntegralArea* ia,
                                  const StreamConstants* sc,
                                  const StreamGauss sg,
                                  const SeparationSizes* sizes)
{
    CALresult err = CAL_RESULT_OK;

    err |= createOutMuBuffer(ci, cm, sizes);
    err |= createOutProbsBuffer(ci, cm, sizes);
    err |= createAPBuffer(ci, cm, ap, sizes);
    err |= createIABuffer(ci, cm, ia, sizes);
    err |= createSCBuffer(ci, cm, sc, sizes);
    err |= createRBuffers(ci, cm, ap, ia, sg, sizes);
    err |= createLBTrigBuffer(ci, cm, ap, ia, sizes);

    return err;
}

CALresult releaseSeparationBuffers(MWCALInfo* ci, SeparationCALMem* cm)
{
    CALresult err = CAL_RESULT_OK;

    err |= releaseMWMemRes(ci->calctx, &cm->outProbs);
    err |= releaseMWMemRes(ci->calctx, &cm->outMu);

    err |= releaseMWMemRes(ci->calctx, &cm->ap);
    err |= releaseMWMemRes(ci->calctx, &cm->ia);
    err |= releaseMWMemRes(ci->calctx, &cm->sc);
    err |= releaseMWMemRes(ci->calctx, &cm->rPts);
    err |= releaseMWMemRes(ci->calctx, &cm->rc);
    err |= releaseMWMemRes(ci->calctx, &cm->sg_dx);
    err |= releaseMWMemRes(ci->calctx, &cm->lbts);

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
         "Wavefront size:        %u Mhz\n"
         "Double precision:      %s\n"
         "Number SIMD:           %u\n"
         "Number shader engines: %u\n"
         "GPU RAM:               %u\n",
         ci->numDevices,
         ci->devID,
         showCALtargetEnum(ci->devInfo.target),
         ci->devAttribs.targetRevision,
         showCALboolean(ci->devAttribs.computeShader),

         ci->devAttribs.engineClock,
         ci->devAttribs.memoryClock,
         ci->devAttribs.wavefrontSize,
         showCALboolean(ci->devAttribs.doublePrecision),
         ci->devAttribs.numberOfSIMD,
         ci->devAttribs.numberOfShaderEngines,
         ci->devAttribs.localRAM);
}

static CALobject createCALBinary(const char* srcIL)
{
    CALobject obj;
    CALresult err;

    err = calclCompile(&obj, CAL_LANGUAGE_IL, srcIL, CAL_TARGET_CYPRESS);
    if (err != CAL_RESULT_OK)
    {
        warn("Error compiling kernel (%d) : %s\n", err, calclGetErrorString());
        return NULL;
    }

    return obj;
}

CALimage readCALImageFromFile(const char* filename)
{
    char* src;
    CALobject obj;
    CALuint size;
    CALresult rc;
    CALimage img;

    src = mwReadFile(filename);
    if (!src)
    {
        perror("IL source file");
        return NULL;
    }

    obj = createCALBinary(src);
    free(src);

    rc = calclLink(&img, &obj, 1);
    calclFreeObject(obj);

    if (rc != CAL_RESULT_OK)
    {
        warn("Error linking image (%d) : %s\n", rc, calclGetErrorString());
        return NULL;
    }

    return img;
}

static void isaLogFunction(const char* msg)
{
    fputs(msg, stdout);
}

CALresult getISA(CALimage image)
{
    calclDisassembleImage(image, isaLogFunction);
    return CAL_RESULT_OK;
}

typedef struct
{
    CALname outMu;
    CALname outProbs;

    CALname ap;
    CALname ia;
    CALname sc;
    CALname rc;
    CALname rPts;
    CALname sg_dx;
    CALname lbts;
} SeparationCALNames;

static inline CALresult getNameMWCALInfo(MWCALInfo* ci, CALname* name, const CALchar* varName)
{
    return calModuleGetName(name, ci->calctx, ci->module, varName);
}

static CALresult getModuleNames(SeparationCALNames* cn, MWCALInfo* ci)
{
    CALresult err = CAL_RESULT_OK;

    err |= getNameMWCALInfo(ci, &cn->outMu, "o0");
    err |= getNameMWCALInfo(ci, &cn->outProbs, "o1");
    err |= getNameMWCALInfo(ci, &cn->ap, "cb0");
    err |= getNameMWCALInfo(ci, &cn->ia, "cb1");
    err |= getNameMWCALInfo(ci, &cn->sc, "cb2");
    err |= getNameMWCALInfo(ci, &cn->rc, "cb3");
    err |= getNameMWCALInfo(ci, &cn->rPts, "cb4");
    err |= getNameMWCALInfo(ci, &cn->sg_dx, "cb5");
    err |= getNameMWCALInfo(ci, &cn->lbts, "cb6");
    err |= getNameMWCALInfo(ci, &cn->lbts, "cb7");

    return err;
}

static CALresult setKernelArguments(MWCALInfo* ci, SeparationCALMem* cm, SeparationCALNames* cn)
{
    CALresult err = CAL_RESULT_OK;

     #if 0

    err |= calCtxSetMem(ci->calctx, cn->outMu, cm->outMu.mem);
    err |= calCtxSetMem(ci->calctx, cn->outProbs, cm->outProbs.mem);

    err |= calCtxSetMem(ci->calctx, cn->ap, cm->ap.mem);
    err |= calCtxSetMem(ci->calctx, cn->ia, cm->ia.mem);
    err |= calCtxSetMem(ci->calctx, cn->sc, cm->sc.mem);
    err |= calCtxSetMem(ci->calctx, cn->rc, cm->rc.mem);
    err |= calCtxSetMem(ci->calctx, cn->rPts, cm->rPts.mem);
    err |= calCtxSetMem(ci->calctx, cn->sg_dx, cm->sg_dx.mem);
    err |= calCtxSetMem(ci->calctx, cn->lbts, cm->lbts.mem);
    #endif

    return err;
}

static CALresult setupCAL(MWCALInfo* ci, SeparationCALMem* cm)
{
    CALresult err;
    SeparationCALNames cn;

    //ci->image = readCALImageFromFile("/home/matt/RawATIAppKernel.il");
    ci->image = readCALImageFromFile("mu_sum_kernel_Cypress.il");
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

    err = getModuleNames(&cn, ci);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to get module names", err);
        return err;
    }

    err = calModuleGetEntry(&ci->func, ci->calctx, ci->module, "main");
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to get module entry", err);
        return err;
    }

    err = setKernelArguments(ci, cm, &cn);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to set kernel arguments", err);
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
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Error running kernel", err);
        return err;
    }

    while (calCtxIsEventDone(ci->calctx, ev) == CAL_RESULT_PENDING);

    return CAL_RESULT_OK;
}

static void sumMuResults(Kahan* bg_prob,
                         const real* mu_results,
                         const CALuint mu_steps,
                         const CALuint r_steps)
{
    CALuint i, j;

    for (i = 0; i < mu_steps; ++i)
    {
        for (j = 0; j < r_steps; ++j)
        {
            KAHAN_ADD(*bg_prob, mu_results[r_steps * i + j]);
        }
    }
}

static void sumStreamResults(Kahan* probs_results,
                             const real* probs_V_reff_xr_rp3,
                             const CALuint number_streams)
{
    CALuint i;

    for (i = 0; i < number_streams; ++i)
        KAHAN_ADD(probs_results[i], probs_V_reff_xr_rp3[i]);
}

static void sumProbsResults(Kahan* probs_results,
                            const real* st_probs_V_reff_xr_rp3_mu_r,
                            const CALuint mu_steps,
                            const CALuint r_steps,
                            const CALuint number_streams)
{
    CALuint i, j, idx;

    for (i = 0; i < mu_steps; ++i)
    {
        for (j = 0; j < r_steps; ++j)
        {
            idx = (i * r_steps * number_streams) + (j * number_streams);
            sumStreamResults(probs_results, &st_probs_V_reff_xr_rp3_mu_r[idx], number_streams);
        }
    }
}

static CALresult readMuResults(MWCALInfo* ci,
                               SeparationCALMem* cm,
                               const IntegralArea* ia,
                               Kahan* bg_prob)
{
    CALresult err;
    CALuint pitch;
    real* mu_results;

    err = mapMWMemRes(&cm->outMu, ci, (CALvoid**) &mu_results, &pitch);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to map result buffer", err);
        return err;
    }

    sumMuResults(bg_prob, mu_results, ia->mu_steps, ia->r_steps);

    err = unmapMWMemRes(&cm->outMu, ci);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to unmap result buffer", err);
        return err;
    }

    return CAL_RESULT_OK;
}

static CALresult runIntegral(MWCALInfo* ci,
                             SeparationCALMem* cm,
                             const IntegralArea* ia,
                             Kahan* probs_results)
{
    CALresult err;
    Kahan bg_sum = ZERO_KAHAN;
    unsigned int i;

    for (i = 0; i < ia->nu_steps; ++i)
    {
        err = runKernel(ci, cm, ia);
        if (err != CAL_RESULT_OK)
        {
            cal_warn("Failed to run kernel", err);
            return err;
        }
    }

    err = readMuResults(ci, cm, ia, &bg_sum);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to read output buffer", err);
        return err;
    }

    //readProbsResults(ci, cm, ia, probs_results);

    return CAL_RESULT_OK;
}
real integrateCAL(const AstronomyParameters* ap,
                  const IntegralArea* ia,
                  const StreamConstants* sc,
                  const StreamGauss sg,
                  real* st_probs,
                  EvaluationState* es,
                  const CLRequest* clr)
{
    real result = NAN;
    MWCALInfo ci = EMPTY_CAL_INFO;
    SeparationCALMem cm = EMPTY_SEPARATION_CAL_MEM;
    CALresult err;

    err = mwInitCAL(&ci, 0);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to init CAL", err);
        return NAN;
    }

    err = setupCAL(&ci, &cm);
    if (err != CAL_RESULT_OK)
        cal_warn("Failed to setup CAL", err);
    else
        result = runIntegral(&ci, &cm, ia, st_probs);

    err = releaseSeparationBuffers(&ci, &cm);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to release buffers", err);
        result = NAN;
    }

    err = mwDestroyCALInfo(&ci);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to cleanup CAL", err);
        result = NAN;
    }

    return result;
}

