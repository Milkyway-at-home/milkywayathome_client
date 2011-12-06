/*
 *  Copyright (c) 2010-2011 Rensselaer Polytechnic Institute
 *  Copyright (c) 2010-2011 Matthew Arsenault
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

#include "nbody_types.h"
#include "nbody_checkpoint.h"
#include "milkyway_util.h"
#include "nbody_defaults.h"

#if HAVE_FCNTL_H
  #include <fcntl.h>
#endif

#if HAVE_WINDOWS_H
  #include <windows.h>
#endif

#if HAVE_SYS_MMAN_H
  #include <sys/mman.h>
#endif

#if HAVE_SYS_TYPES_H
  #include <sys/types.h>
#endif

#if HAVE_SYS_STAT_H
  #include <sys/stat.h>
#endif

#ifndef _WIN32

typedef struct
{
    int fd;            /* File descriptor for checkpoint file */
    char* mptr;        /* mmap'd pointer for checkpoint file */
    size_t cpFileSize; /* For checking how big the file should be for expected bodies */
} CheckpointHandle;

#define EMPTY_CHECKPOINT_HANDLE { 1, NULL, 0 }

#else

typedef struct
{
    HANDLE file;
    HANDLE mapFile;
    char* mptr;
    DWORD cpFileSize;
} CheckpointHandle;

#define EMPTY_CHECKPOINT_HANDLE { INVALID_HANDLE_VALUE, INVALID_HANDLE_VALUE, NULL, 0 }

#endif /* _WIN32 */


/* Checkpoint file: Very simple binary "format"
   Name        Type         Values     Notes
-------------------------------------------------------
   NBodyCheckpointHeader
   bodytab       Body[]     anything   Array of bodies
   orbitTrace    mwvector[] anything   Array of center of mass history
   ending        string     "end"      Kind of dumb and pointless
 */

static const char hdr[] = "mwnbody";
static const char tail[] = "end";

typedef struct
{
    char header[128];                     /* "mwnbody" */
    uint32_t majorVersion, minorVersion;  /* Version check */
    uint32_t nbody;
    uint32_t step;
    uint32_t realSize;                   /* Does the checkpoint use float or double */
    uint32_t ptrSize;
    real rsize;
    uint32_t treeIncest;
    uint32_t nOrbitTrace;
    NBodyCtx ctx;
} NBodyCheckpointHeader;

static const size_t hdrSize = sizeof(NBodyCheckpointHeader) + sizeof(tail);



static void prepareWriteCheckpointHeader(NBodyCheckpointHeader* cp, const NBodyCtx* ctx, const NBodyState* st)
{
    strcpy(cp->header, hdr);
    cp->realSize = sizeof(real);
    cp->ptrSize = sizeof(void*);

    cp->majorVersion = NBODY_VERSION_MAJOR;
    cp->minorVersion= NBODY_VERSION_MINOR;

    memcpy(&cp->ctx, ctx, sizeof(cp->ctx));
    cp->nbody = st->nbody;
    cp->step = st->step;
    cp->rsize = st->tree.rsize;
    cp->treeIncest = st->treeIncest;
    cp->nOrbitTrace = N_ORBIT_TRACE_POINTS;
}

static void readCheckpointHeader(NBodyCheckpointHeader* cp, NBodyCtx* ctx, NBodyState* st)
{
    memcpy(ctx, &cp->ctx, sizeof(*ctx));
    st->nbody = cp->nbody;
    st->step = cp->step;
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
        mw_printf("Didn't find header for checkpoint file.\n");
        return 1;
    }

    /* Make sure the file isn't lying about how many bodies there are */
    if (supposedCheckpointSize != cp->cpFileSize)
    {
        mw_printf("Expected checkpoint file size ("ZU") is incorrect for expected number of bodies "
                  "(%u bodies, real size "ZU")\n",
                  supposedCheckpointSize,
                  st->nbody,
                  (size_t) cp->cpFileSize);
        return 1;
    }

    if (cpHdr->realSize != sizeof(real))
    {
        mw_printf("Got checkpoint file for wrong type. "
                  "Expected sizeof(real) = "ZU", got "ZU"\n",
                  sizeof(real), (size_t) cpHdr->realSize);
        return 1;
    }

    if (cpHdr->ptrSize != sizeof(void*))
    {
        mw_printf("Got checkpoint file for wrong architecture. "
                  "Expected sizeof(void*) = "ZU", got "ZU"\n", sizeof(void*), (size_t) cpHdr->ptrSize);
        return 1;
    }

    if (   cpHdr->majorVersion != NBODY_VERSION_MAJOR
        || cpHdr->minorVersion != NBODY_VERSION_MINOR)
    {
        mw_printf("Version mismatch in checkpoint file. File is for %u.%u, But version is %u.%u\n",
                  cpHdr->majorVersion, cpHdr->minorVersion,
                  NBODY_VERSION_MAJOR, NBODY_VERSION_MINOR);
        return 1;
    }

    return 0;
}


#ifndef _WIN32

static int nbOpenCheckpointHandle(const NBodyState* st,
                                  CheckpointHandle* cp,
                                  const char* filename,
                                  int writing)
{
    struct stat sb;

    cp->fd = open(filename, O_RDWR | O_CREAT, S_IWUSR | S_IRUSR);
    if (cp->fd == -1)
    {
        mwPerror("Error opening checkpoint '%s'", filename);
        return TRUE;
    }

    if (fstat(cp->fd, &sb) == -1)
    {
        mwPerror("Error on fstat() of checkpoint '%s'", filename);
        return TRUE;
    }

    if (!S_ISREG(sb.st_mode))
    {
        mw_printf("Checkpoint '%s' is not a file\n", filename);
        return TRUE;
    }

    if (writing)
    {
        cp->cpFileSize = hdrSize + st->nbody * sizeof(Body) + N_ORBIT_TRACE_POINTS * sizeof(mwvector);
        /* Make the file the right size in case it's a new file */
        if (ftruncate(cp->fd, cp->cpFileSize) < 0)
        {
            mwPerror("Error ftruncate() on checkpoint '%s'", filename);
            return TRUE;
        }
    }
    else
    {
        cp->cpFileSize = sb.st_size;
        if (cp->cpFileSize == 0)
        {
            mw_printf("Checkpoint '%s' is empty\n", filename);
            return 1;
        }
    }

    cp->mptr = mmap(NULL, cp->cpFileSize, PROT_READ | PROT_WRITE, MAP_SHARED, cp->fd, 0);
    if (cp->mptr == MAP_FAILED)
    {
        mwPerror("Error mmap()ing checkpoint '%s'", filename);
        return TRUE;
    }

    return FALSE;
}

static int nbCloseCheckpointHandle(CheckpointHandle* cp)
{
    struct stat sb;

    /* Clean up the checkpointing */
    if (cp->fd != -1)
    {
        if (fstat(cp->fd, &sb) == -1)
        {
            mwPerror("Error on fstat() closing checkpoint");
            return TRUE;
        }

        if (close(cp->fd) == -1)
        {
            mwPerror("closing checkpoint file");
            return TRUE;
        }

        if (cp->mptr && (cp->mptr != MAP_FAILED) && (munmap(cp->mptr, sb.st_size) == -1))
        {
            mwPerror("munmap() checkpoint");
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
static int nbOpenCheckpointHandle(const NBodyState* st,
                                  CheckpointHandle* cp,
                                  const char* filename,
                                  int writing)
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
    {
        mwPerrorW32("Failed to open checkpoint file '%s'\n", filename);
        return TRUE;
    }

    if (writing)
    {
        cp->cpFileSize = hdrSize + st->nbody * sizeof(Body) + N_ORBIT_TRACE_POINTS * sizeof(mwvector);
    }
    else
    {
        /* We don't know how much to expect when reading the file */
        cp->cpFileSize = GetFileSize(cp->file, NULL);
        if (cp->cpFileSize == INVALID_FILE_SIZE || cp->cpFileSize == 0)
        {
            mwPerrorW32("Invalid checkpoint file size (%ld) or empty checkpoint file '%s'",
                        cp->cpFileSize, filename);
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
        mwPerrorW32("Failed to create mapping for checkpoint file '%s'", filename);
        return TRUE;
    }

    cp->mptr = (char*) MapViewOfFile(cp->mapFile,
                                     FILE_MAP_ALL_ACCESS,
                                     0,
                                     fileMapStart,
                                     mapViewSize);
    if (cp->mptr == NULL)
    {
        mwPerrorW32("Failed to open checkpoint file view for file '%s'", filename);
        return TRUE;
    }

    return FALSE;
}

static int nbCloseCheckpointHandle(CheckpointHandle* cp)
{
    if (cp->file != INVALID_HANDLE_VALUE)
    {
        if (cp->mptr && !UnmapViewOfFile((LPVOID) cp->mptr))
        {
            mwPerrorW32("Error unmapping the checkpoint view object");
            return TRUE;
        }

        if (cp->mapFile && !CloseHandle(cp->mapFile))
        {
            mwPerrorW32("Error closing the checkpoint mapping");
            return TRUE;
        }

        if (cp->file && !CloseHandle(cp->file))
        {
            mwPerrorW32("Error closing checkpoint file");
            return TRUE;
        }
    }

    return FALSE;
}

#endif /* _WIN32 */

/* Should be given the same context as the dump. Returns nonzero if the state failed to be thawed */
static int nbThawState(NBodyCtx* ctx, NBodyState* st, CheckpointHandle* cp)
{
    size_t bodySize, traceSize, supposedCheckpointSize;
    NBodyCheckpointHeader cpHdr;
    char* p = cp->mptr;

    memset(&cpHdr, 0, sizeof(cpHdr));
    memcpy(&cpHdr, p, sizeof(cpHdr));
    p += sizeof(cpHdr);

    readCheckpointHeader(&cpHdr, ctx, st);

    assert(cp->cpFileSize != 0);
    bodySize = st->nbody * sizeof(Body);
    traceSize = cpHdr.nOrbitTrace * sizeof(mwvector);
    supposedCheckpointSize = hdrSize + bodySize + traceSize;

    if (verifyCheckpointHeader(&cpHdr, cp, st, supposedCheckpointSize))
    {
        return TRUE;
    }


    /* Read the bodies */
    st->bodytab = (Body*) mwMallocA(bodySize);
    memcpy(st->bodytab, p, bodySize);
    p += bodySize;

    st->orbitTrace = (mwvector*) mwMallocA(traceSize);
    memcpy(st->orbitTrace, p, traceSize);
    p += traceSize;

    if (strncmp(p, tail, sizeof(tail)))
    {
        free(st->bodytab);
        st->bodytab = NULL;

        free(st->orbitTrace);
        st->orbitTrace = NULL;

        mw_printf("Failed to find end marker in checkpoint file.\n");
        return TRUE;
    }

    return FALSE;
}

static void nbFreezeState(const NBodyCtx* ctx, const NBodyState* st, CheckpointHandle* cp)
{
    const size_t bodySize =  st->nbody * sizeof(Body);
    const size_t traceSize = N_ORBIT_TRACE_POINTS * sizeof(mwvector);
    char* p = cp->mptr;
    NBodyCheckpointHeader cpHdr;

    memset(&cpHdr, 0, sizeof(cpHdr));
    prepareWriteCheckpointHeader(&cpHdr, ctx, st);

    memcpy(p, &cpHdr, sizeof(cpHdr));
    p += sizeof(cpHdr);

    /* The main piece of state*/
    memcpy(p, st->bodytab, bodySize);
    p += bodySize;

    memcpy(p, st->orbitTrace, traceSize);
    p += traceSize;

    strcpy(p, tail);
}

/* Open the temporary checkpoint file for writing */
int nbResolveCheckpoint(NBodyState* st, const char* checkpointFileName)
{
    int rc = 0;

    st->checkpointResolved = (char*) mwCalloc(4096, sizeof(char));

    rc = mw_resolve_filename(checkpointFileName, st->checkpointResolved, 4096 * sizeof(char));
    if (rc)
    {
        mw_printf("Failed to resolve checkpoint file '%s': %d\n", checkpointFileName, rc);
        free(st->checkpointResolved);
        st->checkpointResolved = NULL;
    }

    return rc;
}

int nbResolvedCheckpointExists(const NBodyState* st)
{
    if (!st->checkpointResolved)
        mw_panic("Checking if checkpoint exists, but haven't resolved yet\n");

    return mw_file_exists(st->checkpointResolved);
}

/* Try to open a checkpoint with a few tries if the open fails.
   This is in case of weird/rare failures like interrupted system calls.
 */
static int nbOpenCheckpointHandleWithAttempts(const NBodyState* st,
                                              CheckpointHandle* cp,
                                              const char* filename,
                                              int writing)
{
    unsigned int tries = 0;
    const unsigned int maxTries = 5;

    do
    {
        if (!nbOpenCheckpointHandle(st, cp, filename, writing))
            break;

        if (nbCloseCheckpointHandle(cp))
            return TRUE;

        ++tries;
        mwMilliSleep(10);
    }
    while (tries < maxTries);

    if (tries >= maxTries)
    {
        mw_printf("Failed to open checkpoint '%s' after %d tries\n", filename, tries);
        return TRUE;
    }

    return FALSE;
}

/* Read the actual checkpoint file to resume */
int nbReadCheckpoint(NBodyCtx* ctx, NBodyState* st)
{
    CheckpointHandle cp = EMPTY_CHECKPOINT_HANDLE;

    if (nbOpenCheckpointHandleWithAttempts(st, &cp, st->checkpointResolved, FALSE))
    {
        mw_printf("Opening checkpoint '%s' for resuming failed\n", st->checkpointResolved);
        nbCloseCheckpointHandle(&cp);
        return TRUE;
    }

    if (nbThawState(ctx, st, &cp))
    {
        nbCloseCheckpointHandle(&cp);
        return TRUE;
    }

    if (nbCloseCheckpointHandle(&cp))
    {
        return TRUE;
    }

    /* Make sure state is ready to use */
    st->acctab = (mwvector*) mwCallocA(st->nbody, sizeof(mwvector));

    return FALSE;
}

/* Use specified temporary file to avoid bad things happening if
 * multiple tests running at a time */
int nbWriteCheckpointWithTmpFile(const NBodyCtx* ctx, const NBodyState* st, const char* tmpFile)
{
    int failed = FALSE;
    CheckpointHandle cp = EMPTY_CHECKPOINT_HANDLE;

    assert(st->checkpointResolved);

    if (nbOpenCheckpointHandleWithAttempts(st, &cp, tmpFile, TRUE))
    {
        return TRUE;
    }

    nbFreezeState(ctx, st, &cp);

    if (nbCloseCheckpointHandle(&cp))
    {
        mw_printf("Failed to properly close temporary checkpoint file\n");
        failed = TRUE;
    }

    /* Swap the real checkpoint with the temporary atomically. This
     * should avoid corruption in the event the file write is
     * interrupted. */
    /* Don't update if the file was not closed properly; it can't be trusted. */
    if (!failed && mw_rename(tmpFile, st->checkpointResolved))
    {
        mwPerror("Failed to update checkpoint '%s' with temporary", st->checkpointResolved);
        failed = TRUE;
    }

    if (failed)
    {
        mw_printf("Failed to write checkpoint\n");
    }

    return failed;
}

int nbWriteCheckpoint(const NBodyCtx* ctx, const NBodyState* st)
{
    char path[256];

    snprintf(path, sizeof(path), "nbody_checkpoint_tmp_%d", (int) getpid());

    return nbWriteCheckpointWithTmpFile(ctx, st, path);
}

int nbTimeToCheckpoint(const NBodyCtx* ctx, NBodyState* st)
{
    time_t now;

    if (BOINC_APPLICATION)
    {
        return mw_time_to_checkpoint();
    }

    if (ctx->checkpointT < 0)
    {
        return FALSE;
    }

    now = time(NULL);
    if ((now - st->lastCheckpoint) > ctx->checkpointT)
    {
        st->lastCheckpoint = now;
        return TRUE;
    }

    return FALSE;
}

