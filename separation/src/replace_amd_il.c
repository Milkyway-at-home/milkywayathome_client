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


#include <err.h>
#include <errno.h>
#include <sysexits.h>
#include <sys/stat.h>
#include <inttypes.h>
#include <string.h>

#include <libelf.h>

#include "replace_amd_il.h"
#include "milkyway_util.h"



static int replaceAMDILSection(Elf* e, const char* ilBuf, size_t ilLen)
{
    Elf_Scn* scn = NULL;
    size_t shstrndx = 0;
    static const int verbose = 1;
    static const int verboseDebug = 0;


    /* Get section index of section containing the string table of section names */
    if (elf_getshdrstrndx(e, &shstrndx) != 0)
    {
        warnx("elf_getshdrstrndx failed: %s.", elf_errmsg(-1));
        return EX_SOFTWARE;
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
            warnx("elf32_getshdr() failed: %s.", elf_errmsg(-1));
            return EX_SOFTWARE;
        }

        /* Look up the name of the section in the string table */
        name = elf_strptr(e, shstrndx, shdr->sh_name);
        if (!name)
        {
            warnx("elf_strptr() failed: %s.", elf_errmsg(-1));
            return EX_SOFTWARE;
        }

        if (strstr(name, ".amdil") != NULL)
        {
            Elf_Data* data = elf_getdata(scn, NULL);

            if (!data)
            {
                warnx("Failed to get data for section: %s.", elf_errmsg(-1));
                return EX_SOFTWARE;
            }

            if (verbose)
            {
                printf("Replacing section data of type %d, off %d align %zu\n", data->d_type, (int) data->d_off, data->d_align);
            }

            //memset(data->d_buf, 0, data->d_size);
            data->d_buf = (void*) ilBuf;
            data->d_size = ilLen;

            if (!elf_flagdata(data, ELF_C_SET, ELF_F_DIRTY))
            {
                warnx("elf_flagdata() failed: %s.", elf_errmsg(-1));
                return EX_SOFTWARE;
            }

            if (!elf_flagscn(scn, ELF_C_SET, ELF_F_DIRTY))
            {
                warnx("elf_flagscn() failed: %s.", elf_errmsg(-1));
                return EX_SOFTWARE;
            }

            /* Don't let libelf rearrange the sections when writing. Not sure if necessary. */
            if (!elf_flagelf(e, ELF_C_SET, ELF_F_LAYOUT))
            {
                warnx("elf_flagelf() failed: %s.", elf_errmsg(-1));
                return EX_SOFTWARE;
            }

            if (elf_update(e, ELF_C_NULL) < 0)
            {
                warnx("elf_update(NULL) failed: %s.", elf_errmsg(-1));
                return EX_SOFTWARE;
            }
        }

        if (verboseDebug)
        {
            printf("Section %-4.4jd %s\n", (uintmax_t) elf_ndxscn(scn), name);
        }
    }

    if (elf_update(e, ELF_C_WRITE) < 0)
    {
        warnx("elf_update(ELF_C_WRITE) failed: %s.", elf_errmsg(-1));
        return EX_SOFTWARE;
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
        warn("fstat on temporary AMD binary file");
    }

    if (props.st_size <= 0)
    {
        warn("Modified binary file is empty\n");
        return NULL;
    }

    len = props.st_size + 1;

    warn("SIZE IS %zu\n", len);
    strBuf = malloc(len);
    if (!strBuf)
    {
        err(errno, "Failed to allocate space for AMD binary file\n");
    }
    strBuf[props.st_size] = '\0';

    if (read(fd, strBuf, props.st_size) < 0)
    {
        warn("Error reading from AMD Binary file\n");
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

static int processElf(int fd, const char* ilBuf, size_t ilLen)
{
    Elf* e;
    Elf_Kind ek = ELF_K_NONE;

    if (elf_version(EV_CURRENT) == EV_NONE)
    {
        warnx("ELF library initialization failed: %s", elf_errmsg(-1));
        return EX_SOFTWARE;
    }

    e = elf_begin(fd, ELF_C_RDWR, NULL);
    if (!e)
    {
        warnx("elf_begin() failed: %s.", elf_errmsg(-1));
        return EX_SOFTWARE;
    }

    ek = elf_kind(e);
    if (ek != ELF_K_ELF)
    {
        warnx("ELF object is wrong kind: %s", showElfKind(ek));
        return EX_SOFTWARE;
    }

    if (replaceAMDILSection(e, ilBuf, ilLen) < 0)
    {
        warnx("Failed to replace .amdil section");
        return EX_SOFTWARE;
    }

    elf_end(e);

    return 0;
}


static int getTmpBinaryName(char* buf, size_t size)
{
    return snprintf(buf, size, "tmp_cl_binary_%d.bin", (int) getpid());
}

/* Take a program binary from clGetPropramInfo() and a replacement IL source as a string.
   Replace the .amdil section in the ELF image and return a new copy of the binary.
 */
unsigned char* getModifiedAMDBinary(unsigned char* bin, size_t binSize, const char* ilSrc, size_t ilLen, size_t* newBinSizeOut)
{
    int fd = -1;
    int rc = 0;
    char tmpBinFile[128];
    unsigned char* newBin = NULL;

    if (!bin || !ilSrc)
        return NULL;

    getTmpBinaryName(tmpBinFile, sizeof(tmpBinFile));

    /* Write binary to a temporary file since we need a file descriptor for libelf */
    fd = open(tmpBinFile, O_RDWR | O_TRUNC | O_CREAT, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
    if (fd < 0)
    {
        warn("Failed to open AMD binary file '%s'", tmpBinFile);
        return NULL;
    }

    if (write(fd, bin, binSize) <= 0)
    {
        warn("Failed to write temporary binary file '%s'", tmpBinFile);
        return NULL;
    }

    rc = processElf(fd, ilSrc, ilLen);

    if (rc == 0)
    {
        if (lseek(fd, 0, SEEK_SET) != 0)
        {

            warn("Failed to seek temporary binary file");
            return NULL;
        }

        newBin = (unsigned char*) readFD(fd, newBinSizeOut);
    }

    if (close(fd) < 0)
    {
        warn("Failed to close binary file '%s'", tmpBinFile);
        free(newBin);
        return NULL;
    }

    return newBin;
}

