/*
 * Copyright (c) 2011 Matthew Arsenault
 * Copyright (c) 2011 Rensselaer Polytechnic Institute
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

#include <sys/stat.h>
#include <string.h>
#include <stdarg.h>
#include <libelf.h>

#include "replace_amd_il.h"
#include "milkyway_util.h"
#include "milkyway_cl.h"
#include "il_kernels.h"


static char* replaceUAVIds(const char* ilSrc, size_t* lenOut, ...)
{
    char* buf;
    va_list argPtr;
    size_t len;
    int rc;

    len = strlen(ilSrc);
    buf = (char*) mwMalloc(len + 1);
    buf[len] = '\0';

    va_start(argPtr, lenOut);
    rc = vsprintf(buf, ilSrc, argPtr);
    va_end(argPtr);

    /* Should be == len when uavid = 2 digits, slighly less when uavid = 1 digit */
    if ((size_t) rc > len)
    {
        free(buf);
        return NULL;
    }

    if (lenOut)
    {
        *lenOut = len;
    }

    return buf;
}

/* Different UAV IDs are used for different GPUs, and it must match */
static char* getILSrc(int nStream, MWCALtargetEnum target, cl_int uavGuess, size_t* len)
{
    char* ilSrc = NULL;
    cl_int u = (uavGuess >= 0) ? uavGuess : mwUAVIdFromMWCALtargetEnum(target);

    /* Should we be checking which UAV is used from the binary? */
    if (u > 99)
    {
        /* We rely on the size of a 2 digit UAV id being the same as the '%d' format,
           but there are only a few UAVs anyway
         */
        mw_printf("UAV id %u is absurd\n", u);
        return NULL;
    }

    /* This is pretty special */
    switch (nStream)
    {
        case 1:   /* 9 */
            ilSrc = replaceUAVIds(ilKernelSrc1, len, u, u, u, u, u, u, u, u, u);
            break;

        case 2:  /* 11 */
            ilSrc = replaceUAVIds(ilKernelSrc2, len, u, u, u, u, u, u, u, u, u, u, u);
            break;

        case 3:  /* 13 */
            ilSrc = replaceUAVIds(ilKernelSrc3, len, u, u, u, u, u, u, u, u, u, u, u, u, u);
            break;

        case 4:  /* 15 */
            ilSrc = replaceUAVIds(ilKernelSrc4, len, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u);
            break;

        default:
            mw_unreachable();
    }

    if (!ilSrc)
    {
        mw_printf("Error getting processed IL kernel source\n");
    }
    else
    {
        if (target < MW_CAL_TARGET_CYPRESS)
        {
            char* dclPos;

            /* Stupid hack.

               The OpenCL compiler generates a dcl_arena_uav_id(8) on Evergreen and later,
               but it's then unused. This isn't supported on R700, and the build
               fails if it's there. On Evergreen, if we remove the declaration the build is OK,
               but then when it silently fails / runs instantly.
               There is probably a better way of dealing with this.

               We could probably avoid this and a lot of this othero
               junk if we reverse engineered the .rodata section, but
               I'm too lazy to do that.

               Remove this declaration by finding the line and overwriting with spaces.
            */
            dclPos = strstr(ilSrc, "dcl_arena_uav_id(8)\n");
            if (dclPos)
            {
                while (*dclPos != '\n')
                {
                    *dclPos++ = ' ';
                }
            }
        }
    }

    return ilSrc;
}

static char* ilSrc = NULL;
static size_t ilSrcLen = 0;

static void freeILSrc()
{
    free(ilSrc);
    ilSrc = NULL;
    ilSrcLen = 0;
}

static int replaceAMDILSection(Elf* e, int nStream, MWCALtargetEnum target)
{
    Elf_Scn* scn = NULL;
    size_t shstrndx = 0;
    static const int verbose = 1;
    static const int verboseDebug = 0;

    /* Get section index of section containing the string table of section names */
    if (elf_getshdrstrndx(e, &shstrndx) != 0)
    {
        mw_printf("elf_getshdrstrndx failed: %s\n", elf_errmsg(-1));
        return 1;
    }

    /* Iterate through all the sections */
    while ((scn = elf_nextscn(e, scn)))
    {
        Elf32_Shdr* shdr;
        const char* name;

        /* Get the header for this section */
        shdr = elf32_getshdr(scn);
        if (!shdr)
        {
            mw_printf("elf32_getshdr() failed: %s\n", elf_errmsg(-1));
            return 1;
        }

        /* Look up the name of the section in the string table */
        name = elf_strptr(e, shstrndx, shdr->sh_name);
        if (!name)
        {
            mw_printf("elf_strptr() failed: %s\n", elf_errmsg(-1));
            return 1;
        }

        /*
        if (strstr(name, ".rodata") != NULL)
        {
            Elf_Data* data = elf_getdata(scn, NULL);


            FILE* f = fopen("rodata_section.bin", "wb");
            if (f)
            {
                fwrite(data->d_buf, 1, data->d_size, f);
                fclose(f);
            }
            else
            {
                perror("Failed to open file");
            }

            size_t roSize;
            char* r770RO = mwReadFileWithSize("rodata_section_RV770.bin", &roSize);
            assert(r770RO);

            data->d_buf = r770RO;
            data->d_size = roSize;


            if (!elf_flagdata(data, ELF_C_SET, ELF_F_DIRTY))
            {
                mw_printf("elf_flagdata() failed: %s\n", elf_errmsg(-1));
                return 1;
            }

            if (!elf_flagscn(scn, ELF_C_SET, ELF_F_DIRTY))
            {
                mw_printf("elf_flagscn() failed: %s\n", elf_errmsg(-1));
                return 1;
            }

            if (elf_update(e, ELF_C_NULL) < 0)
            {
                mw_printf("elf_update(NULL) failed: %s\n", elf_errmsg(-1));
                return 1;
            }
        }
        */

        if (strstr(name, ".amdil") != NULL)
        {
            int uavId;
            const char* uavComment;
            Elf_Data* data;

            data = elf_getdata(scn, NULL);
            if (!data)
            {
                mw_printf("Failed to get data for .amdil section: %s\n", elf_errmsg(-1));
                return 1;
            }

            if (verbose)
            {
                mw_printf("Replacing section data of type %d, off %d align "ZU"\n",
                          data->d_type,
                          (int) data->d_off,
                          data->d_align);
            }


            /* Before we overwrite it, there is information we would like to extract */
            uavComment = strstr((const char*) data->d_buf, ";uavid:");
            if (!uavComment || (sscanf(uavComment, ";uavid:%d\n", &uavId) != 1))
            {
                mw_printf("Error reading uavid from IL comment\n");
                uavId = -1;
            }

            ilSrc = getILSrc(nStream, target, uavId, &ilSrcLen);
            if (!ilSrc || (ilSrcLen == 0))
            {
                mw_printf("Failed to get IL source\n");
                return 1;
            }

            data->d_buf = (void*) ilSrc;
            data->d_size = ilSrcLen;

            if (!elf_flagdata(data, ELF_C_SET, ELF_F_DIRTY))
            {
                mw_printf("elf_flagdata() failed: %s\n", elf_errmsg(-1));
                return 1;
            }

            if (!elf_flagscn(scn, ELF_C_SET, ELF_F_DIRTY))
            {
                mw_printf("elf_flagscn() failed: %s\n", elf_errmsg(-1));
                return 1;
            }

            /* Don't let libelf rearrange the sections when writing.

               clBuildProgram() crashes on Windows if you don't do
               this with some(?) Catalyst versions
             */
            if (!elf_flagelf(e, ELF_C_SET, ELF_F_LAYOUT))
            {
                mw_printf("elf_flagelf() failed: %s\n", elf_errmsg(-1));
                return 1;
            }

            if (elf_update(e, ELF_C_NULL) < 0)
            {
                mw_printf("elf_update(NULL) failed: %s\n", elf_errmsg(-1));
                return 1;
            }
        }

        if (verboseDebug)
        {
            printf("Section %u %s\n", (unsigned int) elf_ndxscn(scn), name);
        }
    }

    if (elf_update(e, ELF_C_WRITE) < 0)
    {
        mw_printf("elf_update(ELF_C_WRITE) failed: %s\n", elf_errmsg(-1));
        return 1;
    }

    return 0;
}

static char* readFD(int fd, size_t* lenOut)
{
    char* strBuf = NULL;
    struct stat props;
    size_t len = 0;

    if (fd < 0)
    {
        return NULL;
    }

    if (fstat(fd, &props) < 0)
    {
        mwPerror("fstat on temporary AMD binary file");
    }

    if (props.st_size <= 0)
    {
        mw_printf("Modified binary file is empty\n");
        return NULL;
    }

    len = props.st_size + 1;

    strBuf = (char*) malloc(len);
    if (!strBuf)
    {
        mwPerror("Failed to allocate "ZU" for AMD binary file", len);
        return NULL;
    }
    strBuf[props.st_size] = '\0';

    if (read(fd, strBuf, props.st_size) < 0)
    {
        mwPerror("Error reading from AMD Binary file");
        free(strBuf);
        strBuf = NULL;
    }

    if (lenOut)
    {
        *lenOut = props.st_size;
    }

    return strBuf;
}

static const char* showElfKind(Elf_Kind ek)
{
    switch (ek)
    {
        case ELF_K_NONE:
            return "ELF_K_NONE";
        case ELF_K_AR:
            return "ELF_K_AR";
        case ELF_K_COFF:
            return "ELF_K_COFF";
        case ELF_K_ELF:
            return "ELF_K_ELF";
        case ELF_K_NUM:
            return "ELF_K_NUM";
        default:
            return "Unknown Elf_Kind";
    }
}

static int processElf(int fd, int nStream, MWCALtargetEnum target)
{
    int rc;
    Elf* e;
    Elf_Kind ek = ELF_K_NONE;

    if (elf_version(EV_CURRENT) == EV_NONE)
    {
        mw_printf("ELF library initialization failed: %s", elf_errmsg(-1));
        return 1;
    }

    e = elf_begin(fd, ELF_C_RDWR, NULL);
    if (!e)
    {
        mw_printf("elf_begin() failed: %s\n", elf_errmsg(-1));
        return 1;
    }

    ek = elf_kind(e);
    if (ek != ELF_K_ELF)
    {
        mw_printf("ELF object is wrong kind: %s\n", showElfKind(ek));
        return 1;
    }

    rc = replaceAMDILSection(e, nStream, target);
    freeILSrc();
    if (rc < 0)
    {
        mw_printf("Failed to replace .amdil section\n");
        return 1;
    }

    elf_end(e);

    return 0;
}


static int getTmpBinaryName(char* buf, size_t size)
{
    return snprintf(buf, size, "tmp_cl_binary_%d.bin", (int) getpid());
}

#ifdef _WIN32
static const int openPermMode = _S_IREAD | _S_IWRITE;
static const int openMode = _O_RDWR | _O_TRUNC | _O_CREAT | _O_BINARY | _O_SHORT_LIVED;
#else
static const int openPermMode = S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH;
static const int openMode = O_RDWR | O_TRUNC | O_CREAT;
#endif /* _WIN32 */

/* Take a program binary from clGetPropramInfo() and a replacement IL source as a string.
   Replace the .amdil section in the ELF image and return a new copy of the binary.
 */
unsigned char* getModifiedAMDBinary(unsigned char* bin, size_t binSize, int nStream, MWCALtargetEnum target, size_t* newBinSizeOut)
{
    int fd = -1;
    int rc = 0;
    char tmpBinFile[128];
    unsigned char* newBin = NULL;

    if (!bin)
        return NULL;

    getTmpBinaryName(tmpBinFile, sizeof(tmpBinFile));

    /* Write binary to a temporary file since we need a file descriptor for libelf */
    fd = open(tmpBinFile, openMode, openPermMode);
    if (fd < 0)
    {
        mwPerror("Failed to open AMD binary file '%s", tmpBinFile);
        return NULL;
    }

    if (write(fd, bin, binSize) <= 0)
    {
        mwPerror("Failed to write temporary binary file '%s'", tmpBinFile);
        return NULL;
    }

    rc = processElf(fd, nStream, target);
    if (rc == 0)
    {
        if (lseek(fd, 0, SEEK_SET) != 0)
        {
            mwPerror("Failed to seek temporary binary file '%s'", tmpBinFile);
            return NULL;
        }

        newBin = (unsigned char*) readFD(fd, newBinSizeOut);
    }

    if (close(fd) < 0)
    {
        mwPerror("Failed to close binary file '%s'", tmpBinFile);
        free(newBin);
        return NULL;
    }

    mw_remove(tmpBinFile);

    return newBin;
}

