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

static const char hdr[] = "mwnbody";
static const char tail[] = "end";

/* Everything except the size of all the bodies */
static const size_t hdrSize = sizeof(size_t)                              /* size of real */
                            + sizeof(void*)                               /* Check pointer size */
                            + sizeof(NBodyCtx)
                            + 2 * sizeof(int)                    /* Major, minor version number */
                            + sizeof(char) * (sizeof(tail) + sizeof(hdr) - 2) /* error checking tags */
                            + 1 * sizeof(unsigned int)                        /* nbody count */
                            + 2 * sizeof(real)                                /* tnow, rsize */
                            + sizeof(int);                                    /* tree incest */

/* Macros to read/write the buffer and advance the pointer the correct size */
#define DUMP_CTX(p, x) { *((NBodyCtx*) (p)) = *(x); (p) += sizeof(NBodyCtx); }
#define DUMP_REAL(p, x) { *((real*) (p)) = (x); (p) += sizeof(real); }
#define DUMP_INT(p, x) { *((int*) (p)) = (x); (p) += sizeof(int); }
#define DUMP_SIZE_T(p, x) { *((size_t*) (p)) = (x); (p) += sizeof(size_t); }
#define DUMP_STR(p, x, size) { memcpy((p), (x), (size)); (p) += (size); }

#define READ_CTX(x, p) { *(x) = *((NBodyCtx*) (p)); (p) += sizeof(NBodyCtx); }
#define READ_REAL(x, p) { (x) = *((real*) (p)); (p) += sizeof(real); }
#define READ_INT(x, p) { (x) = *((int*) (p)); (p) += sizeof(int); }
#define READ_SIZE_T(x, p) { (x) = *((size_t*) (p)); (p) += sizeof(size_t); }
#define READ_STR(x, p, size) { memcpy((x), (p), (size)); (p) += (size); }


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


static int openCheckpointHandle(const NBodyCtx* ctx, CheckpointHandle* cp, const char* filename, int writing)
{
    SYSTEM_INFO si;
    DWORD sysGran;
    DWORD fSize;
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
        cp->cpFileSize = hdrSize + nbody * sizeof(Body);
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
static inline int thawState(NBodyCtx* ctx, NBodyState* st, CheckpointHandle* cp)
{
    size_t realSize, ptrSize;
    unsigned int majorVersion, minorVersion;
    char buf[sizeof(hdr)];
    char tailBuf[sizeof(tail)];
    size_t bodySize, supposedCheckpointSize;
    char* p = cp->mptr;

    warn("Thawing state\n");

    READ_STR(buf, p, sizeof(hdr) - 1);

    READ_SIZE_T(realSize, p);
    READ_SIZE_T(ptrSize, p);
    READ_INT(majorVersion, p);
    READ_INT(minorVersion, p);

    READ_CTX(ctx, p);

    READ_INT(st->nbody, p);
    bodySize = st->nbody * sizeof(Body);

    READ_REAL(st->tnow, p);
    READ_REAL(st->tree.rsize, p);
    READ_INT(st->treeIncest, p);

    if (strncmp(hdr, buf, sizeof(hdr) - 1))
    {
        return warn1("Didn't find header for checkpoint file.\n");
    }

    assert(cp->cpFileSize != 0);
    /* Make sure the file isn't lying about how many bodies there are */
    supposedCheckpointSize = hdrSize + st->nbody * sizeof(Body);
    if (supposedCheckpointSize != cp->cpFileSize)
    {
        return warn1("Expected checkpoint file size ("ZU") is incorrect for expected number of bodies "
                     "(%u bodies, real size "ZU")\n",
                     supposedCheckpointSize,
                     st->nbody,
                     cp->cpFileSize);
    }

    if (realSize != sizeof(real))
    {
        return warn1("Got checkpoint file for wrong type. "
                     "Expected sizeof(real) = "ZU", got "ZU"\n",
                     sizeof(real), realSize);
    }

    if (ptrSize != sizeof(void*))
    {
        return warn1("Got checkpoint file for wrong architecture. "
                     "Expected sizeof(void*) = "ZU", got "ZU"\n", sizeof(void*), ptrSize);
    }

    if (majorVersion != MILKYWAY_NBODY_VERSION_MAJOR || minorVersion != MILKYWAY_NBODY_VERSION_MINOR)
    {
        return warn1("Version mismatch in checkpoint file. File is for %u.%u, But version is %u.%u\n",
                     majorVersion, minorVersion,
                     MILKYWAY_NBODY_VERSION_MAJOR, MILKYWAY_NBODY_VERSION_MINOR);
    }

    /* Read the bodies */
    st->bodytab = (Body*) mwMallocA(bodySize);
    memcpy(st->bodytab, p, bodySize);
    p += bodySize;

    READ_STR(tailBuf, p, sizeof(tailBuf) - 1);
    if (strncmp(tail, tailBuf, sizeof(tailBuf) - 1))
    {
        free(st->bodytab);
        st->bodytab = NULL;
        return warn1("Failed to find end marker in checkpoint file.\n");

    }

    return FALSE;
}

/* Checkpoint file: Very simple binary "format"
   Name     Type    Values     Notes
-------------------------------------------------------
   header        string   "mwnbody"  No null terminator
   sizeof(real)  size_t   4, 8       Does the checkpoint use float or double
   sizeof(void*) size_t   4, 8
   version major int      4
   version minor int      4
                                        Saved parts of the program state
   ctx           NBodyCtx sizeof(NBodyCtx)
   nbody         uint     anything   Num. of bodies expected. Error if doesn't match nbody in context.
   tnow          real     anything
   rsize         real     anything
   treeIncest    int      0, 1       If incest has occured
   bodytab       Body*    anything   Array of bodies
   ending        string   "end"      No null terminator
 */

/* Use a very simple flag to mark when writing the checkpoint file
 * begins and ends. I think this should always be good enough, unless
 * something really weird happens. If the read is interrupted, the
 * checkpoint file is garbage and we lose everything. Uses the boinc
 * critical sections, so it hopefully won't be interrupted.
 */
static inline void freezeState(const NBodyCtx* ctx, const NBodyState* st, CheckpointHandle* cp)
{
    const size_t bodySize = sizeof(Body) * st->nbody;
    char* p = cp->mptr;

    /* -1 so we don't bother with the null terminator. It's slightly
        annoying since the strcmps use it, but memcpy doesn't. We
        don't need it anyway  */
    DUMP_STR(p, hdr, sizeof(hdr) - 1);  /* Simple marker for a checkpoint file */
    DUMP_SIZE_T(p, sizeof(real));   /* Make sure we don't confuse double and float checkpoints */
    DUMP_SIZE_T(p, sizeof(void*));
    DUMP_INT(p, MILKYWAY_NBODY_VERSION_MAJOR);
    DUMP_INT(p, MILKYWAY_NBODY_VERSION_MINOR);

    /* Now that we have some basic check stuff written, dump the state */
    DUMP_CTX(p, ctx);

    /* Little state pieces */
    DUMP_INT(p, st->nbody);         /* Make sure we get the right number of bodies */
    DUMP_REAL(p, st->tnow);
    DUMP_REAL(p, st->tree.rsize);
    DUMP_INT(p, st->treeIncest);

    /* The main piece of state*/
    memcpy(p, st->bodytab, bodySize);
    p += bodySize;

    DUMP_STR(p, tail, sizeof(tail) - 1);
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

