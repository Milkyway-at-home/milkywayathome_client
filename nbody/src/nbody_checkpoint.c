/* Copyright 2010 Matthew Arsenault, Travis Desell, Dave Przybylo,
Nathan Cole, Boleslaw Szymanski, Heidi Newberg, Carlos Varela, Malik
Magdon-Ismail and Rensselaer Polytechnic Institute.

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

#ifndef _WIN32
  #include <unistd.h>
  #include <sys/mman.h>
  #include <sys/types.h>
  #include <sys/stat.h>
  #include <fcntl.h>
#else
  #include <windows.h>
#endif /* _WIN32 */

#include <string.h>
#include "nbody_types.h"
#include "nbody_checkpoint.h"
#include "milkyway_util.h"
#include "nbody_defaults.h"



/* Checkpoint file: Very simple binary "format"
   Name     Type    Values     Notes
-------------------------------------------------------
   NBodyCheckpointHeader
   bodytab       Body*    anything   Array of bodies
   ending        string   "end"      Kind of dumb and pointless
 */

static const char hdr[] = "mwnbody";
static const char tail[] = "end";

typedef struct __attribute__((packed, aligned(512)))
{
    char header[sizeof(hdr)];        /* "mwnbody" */
    size_t realSize;                 /* Does the checkpoint use float or double */
    size_t ptrSize;
    int majorVersion, minorVersion;  /* Version check */
    unsigned int nbody;              /* Saved copies of state */
    real tnow;
    real rsize;
    int treeIncest;
    NBodyCtx ctx;
} NBodyCheckpointHeader;

static const size_t hdrSize = sizeof(NBodyCheckpointHeader) + sizeof(tail);



static void prepareWriteCheckpointHeader(NBodyCheckpointHeader* cp, const NBodyCtx* ctx, const NBodyState* st)
{
    strcpy(cp->header, hdr);
    cp->realSize = sizeof(real);
    cp->ptrSize = sizeof(void*);

    cp->majorVersion = MILKYWAY_NBODY_VERSION_MAJOR;
    cp->minorVersion= MILKYWAY_NBODY_VERSION_MINOR;

    memcpy(&cp->ctx, ctx, sizeof(cp->ctx));
    cp->nbody = st->nbody;
    cp->tnow = st->tnow;
    cp->rsize = st->tree.rsize;
    cp->treeIncest = st->treeIncest;
}

static void readCheckpointHeader(NBodyCheckpointHeader* cp, NBodyCtx* ctx, NBodyState* st)
{
    memcpy(ctx, &cp->ctx, sizeof(*ctx));
    st->nbody = cp->nbody;
    st->tnow = cp->tnow;
    st->tree.rsize = cp->rsize;
    st->treeIncest = cp->treeIncest;
}

static int verifyCheckpointHeader(const NBodyCheckpointHeader* cpHdr,
                                  const CheckpointHandle* cp,
                                  const NBodyState* st,
                                  size_t supposedCheckpointSize)
{
    if (strncmp(cpHdr->header, hdr, sizeof(cpHdr->header)))
    {
        return warn1("Didn't find header for checkpoint file.\n");
    }

    /* Make sure the file isn't lying about how many bodies there are */
    if (supposedCheckpointSize != cp->cpFileSize)
    {
        return warn1("Expected checkpoint file size ("ZU") is incorrect for expected number of bodies "
                     "(%u bodies, real size "ZU")\n",
                     supposedCheckpointSize,
                     st->nbody,
                     (size_t) cp->cpFileSize);
    }

    if (cpHdr->realSize != sizeof(real))
    {
        return warn1("Got checkpoint file for wrong type. "
                     "Expected sizeof(real) = "ZU", got "ZU"\n",
                     sizeof(real), cpHdr->realSize);
    }

    if (cpHdr->ptrSize != sizeof(void*))
    {
        return warn1("Got checkpoint file for wrong architecture. "
                     "Expected sizeof(void*) = "ZU", got "ZU"\n", sizeof(void*), cpHdr->ptrSize);
    }

    if (   cpHdr->majorVersion != MILKYWAY_NBODY_VERSION_MAJOR
        || cpHdr->minorVersion != MILKYWAY_NBODY_VERSION_MINOR)
    {
        return warn1("Version mismatch in checkpoint file. File is for %u.%u, But version is %u.%u\n",
                     cpHdr->majorVersion, cpHdr->minorVersion,
                     MILKYWAY_NBODY_VERSION_MAJOR, MILKYWAY_NBODY_VERSION_MINOR);
    }

    return 0;
}


#ifndef _WIN32

static int openCheckpointHandle(const NBodyState* st,
                                CheckpointHandle* cp,
                                const char* filename, int writing)
{
    struct stat sb;

    cp->fd = open(filename, O_RDWR | O_CREAT, S_IWUSR | S_IRUSR);
    if (cp->fd == -1)
    {
        perror("open checkpoint tmp");
        return TRUE;
    }

    if (fstat(cp->fd, &sb) == -1)
    {
        perror("checkpoint fstat");
        return TRUE;
    }

    if (!S_ISREG(sb.st_mode))
        return warn1("checkpoint file is not a file\n");

    if (writing)
    {
        cp->cpFileSize = hdrSize + st->nbody * sizeof(Body);
        /* Make the file the right size in case it's a new file */
        if (ftruncate(cp->fd, cp->cpFileSize) < 0)
        {
            perror("ftruncate checkpoint");
            return TRUE;
        }
    }
    else
    {
        cp->cpFileSize = sb.st_size;
        if (cp->cpFileSize == 0)
            return warn1("checkpoint file is empty\n");
    }

    cp->mptr = mmap(0, cp->cpFileSize, PROT_READ | PROT_WRITE, MAP_SHARED, cp->fd, 0);
    if (cp->mptr == MAP_FAILED)
    {
        perror("mmap: Failed to open checkpoint file for writing");
        return TRUE;
    }

    return FALSE;
}

static int closeCheckpointHandle(CheckpointHandle* cp)
{
    struct stat sb;

    /* Clean up the checkpointing */
    if (cp->fd != -1)
    {
        if (fstat(cp->fd, &sb) == -1)
        {
            perror("fstat on closing checkpoint");
            return TRUE;
        }

        if (close(cp->fd) == -1)
        {
            perror("closing checkpoint file");
            return TRUE;
        }

        if (cp->mptr && cp->mptr != MAP_FAILED && munmap(cp->mptr, sb.st_size) == -1)
        {
            perror("munmap checkpoint");
            return TRUE;
        }
    }

    return FALSE;
}

#else  /* Windows version */

/* Relevant: http://msdn.microsoft.com/en-us/library/aa366556(VS.85).aspx
             http://msdn.microsoft.com/en-us/library/aa366548(VS.85).aspx
             http://msdn.microsoft.com/en-us/library/ms810613.aspx

             Flushing:
             http://msdn.microsoft.com/en-us/library/aa366563(v=VS.85).aspx
 */
static int openCheckpointHandle(const NBodyState* st, CheckpointHandle* cp, const char* filename, int writing)
{
    SYSTEM_INFO si;
    DWORD sysGran;
    DWORD mapViewSize, fileMapStart, fileMapSize;

    /* Try to create a new file */
    cp->file = CreateFile(filename,
                          GENERIC_READ | GENERIC_WRITE,
                          0,     /* Other processes can't touch this */
                          NULL,
                          CREATE_NEW,
                          FILE_FLAG_SEQUENTIAL_SCAN,
                          NULL);

    /* If the checkpoint already exists, open it */
    if (GetLastError() == ERROR_FILE_EXISTS)
    {
        cp->file = CreateFile(filename,
                              GENERIC_READ | GENERIC_WRITE,
                              0,     /* Other processes can't touch this */
                              NULL,
                              OPEN_EXISTING,
                              FILE_FLAG_SEQUENTIAL_SCAN,
                              NULL);
    }

    /* TODO: More filetype checking and stuff */

    if (cp->file == INVALID_HANDLE_VALUE)
        return warn1("Failed to open checkpoint file '%s': %ld\n", filename, GetLastError());

    if (writing)
    {
        cp->cpFileSize = hdrSize + st->nbody * sizeof(Body);
    }
    else
    {
        /* We don't know how much to expect when reading the file */
        cp->cpFileSize = GetFileSize(cp->file, NULL);
        if (cp->cpFileSize == INVALID_FILE_SIZE || cp->cpFileSize == 0)
        {
            warn("Invalid checkpoint file size (%ld) or empty checkpoint file '%s': %ld\n",
                 cp->cpFileSize, filename, GetLastError());
            CloseHandle(cp->file);
            return TRUE;
        }
    }

    GetSystemInfo(&si);
    sysGran = si.dwAllocationGranularity;
    fileMapStart = 0;
    mapViewSize = cp->cpFileSize;
    fileMapSize = cp->cpFileSize;

    cp->mapFile = CreateFileMapping(cp->file,
                                    NULL,
                                    PAGE_READWRITE,
                                    0,
                                    fileMapSize,
                                    NULL);
    if (cp->mapFile == NULL)
    {
        return warn1("Failed to create mapping for checkpoint file '%s': %ld\n",
                     filename, GetLastError());
    }

    cp->mptr = (char*) MapViewOfFile(cp->mapFile,
                                     FILE_MAP_ALL_ACCESS,
                                     0,
                                     fileMapStart,
                                     mapViewSize);
    if (cp->mptr == NULL)
    {
        return warn1("Failed to open checkpoint file view for file '%s': %ld\n",
                     filename, GetLastError());
    }

    return FALSE;
}

static int closeCheckpointHandle(CheckpointHandle* cp)
{
    if (cp->file != INVALID_HANDLE_VALUE)
    {
        if (cp->mptr && !UnmapViewOfFile((LPVOID) cp->mptr))
            return warn1("Error %ld occurred unmapping the checkpoint view object!\n", GetLastError());

        if (cp->mapFile && !CloseHandle(cp->mapFile))
            return warn1("Error %ld occurred closing the checkpoint mapping!\n", GetLastError());

        if (cp->file && !CloseHandle(cp->file))
            return warn1("Error %ld occurred closing checkpoint file\n", GetLastError());
    }

    return FALSE;
}

#endif /* _WIN32 */

/* Should be given the same context as the dump. Returns nonzero if the state failed to be thawed */
static int thawState(NBodyCtx* ctx, NBodyState* st, CheckpointHandle* cp)
{
    size_t bodySize, supposedCheckpointSize;
    NBodyCheckpointHeader cpHdr;
    char* p = cp->mptr;

    memset(&cpHdr, 0, sizeof(cpHdr));
    memcpy(&cpHdr, p, sizeof(cpHdr));
    p += sizeof(cpHdr);

    readCheckpointHeader(&cpHdr, ctx, st);

    assert(cp->cpFileSize != 0);
    bodySize = st->nbody * sizeof(Body);
    supposedCheckpointSize = hdrSize + bodySize;

    verifyCheckpointHeader(&cpHdr, cp, st, supposedCheckpointSize);

    /* Read the bodies */
    st->bodytab = (Body*) mwMallocA(bodySize);
    memcpy(st->bodytab, p, bodySize);
    p += bodySize;

    if (strncmp(p, tail, sizeof(tail)))
    {
        free(st->bodytab);
        st->bodytab = NULL;
        return warn1("Failed to find end marker in checkpoint file.\n");
    }

    return FALSE;
}

static void freezeState(const NBodyCtx* ctx, const NBodyState* st, CheckpointHandle* cp)
{
    const size_t bodySize = sizeof(Body) * st->nbody;
    char* p = cp->mptr;
    NBodyCheckpointHeader cpHdr;

    memset(&cpHdr, 0, sizeof(cpHdr));
    prepareWriteCheckpointHeader(&cpHdr, ctx, st);

    memcpy(p, &cpHdr, sizeof(cpHdr));
    p += sizeof(cpHdr);

    /* The main piece of state*/
    memcpy(p, st->bodytab, bodySize);
    p += bodySize;

    strcpy(p, tail);
}

/* Open the temporary checkpoint file for writing */
int resolveCheckpoint(NBodyState* st, const char* checkpointFileName)
{
    int rc = 0;

    st->checkpointResolved = (char*) mwCalloc(2048, sizeof(char));

    rc = mw_resolve_filename(checkpointFileName, st->checkpointResolved, 2048 * sizeof(char));
    if (rc)
    {
        warn("Failed to resolve checkpoint file '%s': %d\n", checkpointFileName, rc);
        free(st->checkpointResolved);
        st->checkpointResolved = NULL;
    }

    return rc;
}

int resolvedCheckpointExists(const NBodyState* st)
{
    if (!st->checkpointResolved)
        mw_panic("Checking if checkpoint exists, but haven't resolved yet\n");

    return mw_file_exists(st->checkpointResolved);
}


/* Read the actual checkpoint file to resume */
int readCheckpoint(NBodyCtx* ctx, NBodyState* st)
{
    CheckpointHandle cp = EMPTY_CHECKPOINT_HANDLE;

    if (openCheckpointHandle(st, &cp, st->checkpointResolved, FALSE))
    {
        warn("Opening checkpoint for resuming failed\n");
        closeCheckpointHandle(&cp);
        return TRUE;
    }

    if (thawState(ctx, st, &cp))
    {
        warn("Thawing state failed\n");
        return TRUE;
    }

    if (closeCheckpointHandle(&cp))
        warn("Failed to close checkpoint properly\n");

    /* Make sure state is ready to use */
    st->acctab = (mwvector*) mwCallocA(st->nbody, sizeof(mwvector));

    return FALSE;
}

/* Use specified temporary file to avoid bad things happening if
 * multiple tests running at a time */
int writeCheckpointWithTmpFile(const NBodyCtx* ctx, const NBodyState* st, const char* tmpFile)
{
    int failed = FALSE;
    CheckpointHandle cp = EMPTY_CHECKPOINT_HANDLE;

    assert(st->checkpointResolved);
    if (openCheckpointHandle(st, &cp, tmpFile, TRUE))
    {
        closeCheckpointHandle(&cp);
        return warn1("Failed to open temporary checkpoint file\n"
                     "Failed to write checkpoint\n");
    }

    freezeState(ctx, st, &cp);

    if (closeCheckpointHandle(&cp))
    {
        warn("Failed to properly close temporary checkpoint file\n");
        failed = TRUE;
    }

    /* Swap the real checkpoint with the temporary atomically. This
     * should avoid corruption in the event the file write is
     * interrupted. */
    /* Don't update if the file was not closed properly; it can't be trusted. */
    if (!failed && mw_rename(tmpFile, st->checkpointResolved))
    {
        mw_win_perror("Failed to update checkpoint with temporary");
        failed = TRUE;
    }

    if (failed)
        warn("Failed to write checkpoint\n");

    return failed;
}

int writeCheckpoint(const NBodyCtx* ctx, const NBodyState* st)
{
    return writeCheckpointWithTmpFile(ctx, st, CHECKPOINT_TMP_FILE);
}

