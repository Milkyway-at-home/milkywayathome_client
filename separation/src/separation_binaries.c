/* Copyright 2010 Matthew Arsenault, Travis Desell, Boleslaw
Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail and
Rensselaer Polytechnic Institute.

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
#include "milkyway_cl.h"
#include "mw_cl.h"
#include "setup_cl.h"
#include "separation_cl_buffers.h"
#include "separation_cl_defs.h"
#include "separation_binaries.h"

#define SEPARATION_BINARY_HEADER "milkyway_separation_CL_kernel"
#define SEPARATION_BINARY_TAIL   "end_CL_kernel"


static void setSeparationBinaryHeader(SeparationBinaryHeader* hdr,
                                      const AstronomyParameters* ap,
                                      const DevInfo* di,
                                      size_t binSize)
{
    hdr->versionMajor   = SEPARATION_VERSION_MAJOR;
    hdr->versionMinor   = SEPARATION_VERSION_MINOR;
    hdr->doublePrec     = DOUBLEPREC;

    hdr->number_streams = ap->number_streams;
    hdr->fast_h_prob    = ap->fast_h_prob;
    hdr->aux_bg_profile = ap->aux_bg_profile;

    hdr->binSize        = binSize;

    hdr->devType        = di->devType;
    hdr->vendorID       = di->vendorID;

    strncpy(hdr->devName, di->devName, sizeof(hdr->devName));
    strncpy(hdr->deviceVersion, di->version, sizeof(hdr->deviceVersion));
    strncpy(hdr->driverVersion, di->driver, sizeof(hdr->driverVersion));
    memset(hdr->_reserved, 0, sizeof(hdr->_reserved));
}

static cl_bool checkBinaryHeader(const AstronomyParameters* ap,
                                 const DevInfo* di,
                                 const SeparationBinaryHeader* hdr)
{
    if (   hdr->versionMajor != SEPARATION_VERSION_MAJOR
        || hdr->versionMinor != SEPARATION_VERSION_MINOR)
    {
        warn("Version of precompiled kernel doesn't match\n");
        return CL_FALSE;
    }

    if (hdr->doublePrec != DOUBLEPREC)
    {
        warn("Precompiled kernel precision does not match\n");
        return CL_FALSE;
    }

    if (   hdr->number_streams != ap->number_streams
        || hdr->fast_h_prob    != ap->fast_h_prob
        || hdr->aux_bg_profile != ap->aux_bg_profile)
    {
        warn("Kernel compiled parameters do not match\n");
        return CL_FALSE;
    }

    if (hdr->vendorID != di->vendorID)
    {
        warn("Vendor ID of binary does not match\n");
        return CL_FALSE;
    }

    if (strncmp(hdr->devName, di->devName, sizeof(di->devName)))
    {
        warn("Device does not match\n");
        return CL_FALSE;
    }

    if (strncmp(hdr->deviceVersion, di->version, sizeof(di->version)))
    {
        warn("Device version does not match\n");
        return CL_FALSE;
    }

    if (strncmp(hdr->driverVersion, di->driver, sizeof(di->driver)))
    {
        warn("Device driver version does not match\n");
        return CL_FALSE;
    }

    return CL_TRUE;
}

static unsigned char* readCoreBinary(FILE* f, SeparationBinaryHeader* hdr)
{
    unsigned char* bin;
    size_t readn;

    bin = (unsigned char*) mwMalloc(hdr->binSize);

    readn = fread(bin, sizeof(unsigned char), hdr->binSize, f);
    if (readn != hdr->binSize)
    {
        warn("Error reading program binary header: read "ZU", expected "ZU"\n",
             readn, hdr->binSize);
        hdr->binSize = 0;
        free(bin);
        bin = NULL;
    }

    return bin;
}

static mwbool freadCheckedStr(char* buf, const char* str, size_t len, FILE* f)
{
    size_t readn;

    readn = fread(buf, sizeof(char), len, f);
    if (readn != len)
    {
        warn("Failed to read '%s' from file: read "ZU", expected "ZU"\n", str, readn, len);
        return TRUE;
    }

    if (strncmp(buf, str, len))
    {
        warn("Read string doesn't match '%s'\n", str);
        return TRUE;
    }

    return FALSE;
}

static unsigned char* separationLoadBinaryFile(FILE* f,
                                               SeparationBinaryHeader* hdr,
                                               size_t* binSizeOut)
{
    char buf[64] = "";
    size_t readn;
    unsigned char* bin;

    if (freadCheckedStr(buf, SEPARATION_BINARY_HEADER, sizeof(SEPARATION_BINARY_HEADER), f))
    {
        warn("Failed to find binary prefix\n");
        return NULL;
    }

    readn = fread(hdr, sizeof(SeparationBinaryHeader), 1, f);
    if (readn != 1)
    {
        warn("Error reading program binary header: read "ZU", expected %d\n", readn, 1);
        return NULL;
    }

    bin = readCoreBinary(f, hdr);

    if (freadCheckedStr(buf, SEPARATION_BINARY_TAIL, sizeof(SEPARATION_BINARY_TAIL), f))
    {
        warn("Failed to find end marker of program binary\n");
        free(bin);
        bin = NULL;
        hdr->binSize = 0;
    }

    if (binSizeOut)
        *binSizeOut = hdr->binSize;
    return bin;
}

unsigned char* separationLoadBinary(const AstronomyParameters* ap,
                                    const DevInfo* di,
                                    const char* filename,
                                    size_t* binSizeOut)
{
    unsigned char* bin;
    FILE* f;
    SeparationBinaryHeader hdr;

    f = mwOpenResolved(filename, "rb");
    if (!f)
    {
        perror("Failed to open file to read program binary");
        return NULL;
    }

    bin = separationLoadBinaryFile(f, &hdr, binSizeOut);
    if (fclose(f))
        perror("Failed to close program binary");

    if (!checkBinaryHeader(ap, di, &hdr))
    {
        warn("Binary header invalid for this device\n");
        free(bin);
        return NULL;
    }

    return bin;
}

cl_bool separationSaveBinary(const AstronomyParameters* ap,
                             const DevInfo* di,
                             const unsigned char* bin,
                             const size_t binSize,
                             const char* filename)
{
    SeparationBinaryHeader hdr;
    FILE* f;

    f = mwOpenResolved(filename, "wb");
    if (!f)
    {
        perror("Failed to open file to save program binary");
        return CL_TRUE;
    }

    setSeparationBinaryHeader(&hdr, ap, di, binSize);

    fwrite(SEPARATION_BINARY_HEADER, sizeof(SEPARATION_BINARY_HEADER), 1, f);
    fwrite(&hdr, sizeof(SeparationBinaryHeader), 1, f);
    fwrite(bin, binSize, 1, f);
    fwrite(SEPARATION_BINARY_TAIL, sizeof(SEPARATION_BINARY_TAIL), 1, f);

    if (fclose(f))
    {
        perror("Failed to close program binary");
        return CL_TRUE;
    }

    return CL_FALSE;
}

