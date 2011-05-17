/*
Copyright (C) 2011  Matthew Arsenault

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

#ifndef _SEPARATION_CAL_TYPES_H_
#define _SEPARATION_CAL_TYPES_H_

#include <CAL/cal.h>
#include <CAL/calcl.h>


#include "milkyway_util.h"
#include "separation_types.h"
#include "evaluation_state.h"

#ifdef __cplusplus
extern "C" {
#endif

/* FIXME: Also defined in CL version */
typedef struct
{
    size_t outBg;
    size_t outStreams;

    size_t nuSteps, muSteps, rSteps;

    size_t rPts;
    size_t rc;
    size_t sg_dx;
    size_t lTrig;
    size_t bTrig;
} CALSeparationSizes;


typedef struct
{
    CALname outBg;
    CALname* outStreams;

    CALname inMu;
    CALname* inStreams;

    CALname nuBuf;
    CALname rPts;
    CALname rc;
    CALname sg_dx;
    CALname lTrig;
    CALname bTrig;
} SeparationCALNames;


typedef struct
{
    CALuint major, minor, patchLevel;
} MWCALVersion;

#define EMPTY_CAL_VERSION { 0, 0, 0 }

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


/* Pair of resource and associated CALmem */
typedef struct
{
    CALresource res;
    CALmem mem;
} MWMemRes;

#define EMPTY_MEM_RES { 0, 0 }

typedef struct
{
    MWMemRes outBg;
    MWMemRes* outStreams;

    /* constant, read only buffers */
    MWMemRes rc;        /* r constants */
    MWMemRes rPts;

    MWMemRes sg_dx;
    MWMemRes sg_qgauss_W;

    MWMemRes starsXY;
    MWMemRes starsZ;

    MWMemRes lTrig;      /* sin, cos of l */
    MWMemRes bTrig;      /* sin, cos of b */
    MWMemRes nuBuf;
    CALuint numberStreams;
} SeparationCALMem;

typedef struct
{
    CALuint nChunkMu;
    //CALuint nChunkR;

    CALuint chunkWaitTime;   /* Estimated time (ms) per chunk for waiting */

    CALuint* chunkMuBorders;
    //CALuint* chunkRBorders;
} SeparationCALChunks;

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


#ifdef __cplusplus
}
#endif

#endif /* _SEPARATION_CAL_TYPES_H_ */

