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
#include "nbody_priv.h"
#include "milkyway_util.h"
#include "io.h"

#define CHECKPOINT_TMP_FILE "nbody_checkpoint_tmp"


static const char hdr[] = "mwnbody";
static const char tail[] = "end";

/* Everything except the size of all the bodies */
static const size_t hdrSize = sizeof(size_t)                                  /* size of real */
                            + sizeof(char) * (sizeof(tail) + sizeof(hdr) - 2) /* error checking tags */
                            + 1 * sizeof(unsigned int)                        /* nbody count */
                            + 2 * sizeof(real);                               /* tnow, rsize */

/* Macros to read/write the buffer and advance the pointer the correct size */
#define DUMP_REAL(p, x) { *((real*) (p)) = (x); (p) += sizeof(real); }
#define DUMP_INT(p, x) { *((int*) (p)) = (x); (p) += sizeof(int); }
#define DUMP_SIZE_T(p, x) { *((size_t*) (p)) = (x); (p) += sizeof(size_t); }
#define DUMP_STR(p, x, size) { memcpy((p), (x), (size)); (p) += (size); }

#define READ_REAL(x, p) { (x) = *((real*) (p)); (p) += sizeof(real); }
#define READ_INT(x, p) { (x) = *((int*) (p)); (p) += sizeof(int); }
#define READ_SIZE_T(x, p) { (x) = *((size_t*) (p)); (p) += sizeof(size_t); }
#define READ_STR(x, p, size) { memcpy((x), (p), (size)); (p) += (size); }


#ifndef _WIN32

static int openCheckpointHandle(const NBodyCtx* ctx, CheckpointHandle* cp, const char* filename)
{
    struct stat sb;
    const size_t checkpointFileSize = hdrSize + ctx->nbody * sizeof(body);

    cp->fd = open(filename, O_RDWR | O_CREAT, S_IWUSR | S_IRUSR);
    if (cp->fd == -1)
    {
        perror("open checkpoint tmp");
        return TRUE;
    }

    /* Make the file the right size in case it's a new file */
    if (ftruncate(cp->fd, checkpointFileSize) < 0)
    {
        perror("ftruncate checkpoint");
        return TRUE;
    }

    if (fstat(cp->fd, &sb) == -1)
    {
        perror("fstat");
        return TRUE;
    }

    if (!S_ISREG(sb.st_mode))
        return warn("checkpoint file is not a file\n");

    cp->mptr = mmap(0, sb.st_size, PROT_READ | PROT_WRITE, MAP_SHARED, cp->fd, 0);
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

        if (munmap(cp->mptr, sb.st_size) == -1)
        {
            perror("munmap");
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


static int openCheckpointHandle(const NBodyCtx* ctx, CheckpointHandle* cp, const char* filename)
{
    SYSTEM_INFO si;
    DWORD sysGran;
    DWORD mapViewSize;
    DWORD fileMapStart;
    DWORD fileMapSize;
    const DWORD checkpointFileSize = hdrSize + ctx->nbody * sizeof(body);

    /* Try to create a new file */
    cp->file = CreateFile(filename,
                          GENERIC_READ | GENERIC_WRITE,
                          0,     /* Other processes can't touch this */
                          NULL,
                          CREATE_NEW,
                          FILE_FLAG_SEQUENTIAL_SCAN,
                          NULL);

    /* If the checkpoint already exists, open it */
    if ( GetLastError() == ERROR_FILE_EXISTS )
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

    if ( cp->file == INVALID_HANDLE_VALUE )
        return warn("Failed to open checkpoint file '%s': %ld\n", filename, GetLastError());

    GetSystemInfo(&si);
    sysGran = si.dwAllocationGranularity;
    fileMapStart = 0;
    mapViewSize = checkpointFileSize;
    fileMapSize = checkpointFileSize;

    cp->mapFile = CreateFileMapping(cp->file,
                                    NULL,
                                    PAGE_READWRITE,
                                    0,
                                    fileMapSize,
                                    NULL);

    if (cp->mapFile == NULL)
    {
        return warn("Failed to create mapping for checkpoint file '%s': %ld\n",
                    filename, GetLastError());
    }

    cp->mptr = (char*) MapViewOfFile(cp->mapFile,
                                     FILE_MAP_ALL_ACCESS,
                                     0,
                                     fileMapStart,
                                     mapViewSize);
    if (cp->mptr == NULL)
    {
        return warn("Failed to open checkpoint file view for file '%s': %ld\n",
                    filename, GetLastError());
    }

    return FALSE;
}

static int closeCheckpointHandle(CheckpointHandle* cp)
{
    if ( cp->file != INVALID_HANDLE_VALUE )
    {
        if (!UnmapViewOfFile((LPVOID) cp->mptr))
        {
            warn("Error %ld occurred unmapping the checkpoint view object!\n", GetLastError());
            return TRUE;
        }

        if (!CloseHandle(cp->mapFile))
        {
            warn("Error %ld occurred closing the checkpoint mapping!\n", GetLastError());
            return TRUE;
        }

        if (!CloseHandle(cp->file))
        {
            warn("Error %ld occurred closing checkpoint file\n", GetLastError());
            return TRUE;
        }
    }

    return FALSE;
}

#endif /* _WIN32 */

/* Should be given the same context as the dump. Returns nonzero if the state failed to be thawed */
static inline int thawState(const NBodyCtx* ctx, NBodyState* st, CheckpointHandle* cp)
{
    unsigned int nbody;
    size_t realSize;
    char buf[sizeof(hdr)];
    char tailBuf[sizeof(tail)];
    const size_t bodySize = ctx->nbody * sizeof(body);
    char* p = cp->mptr;
    int failed = FALSE;

    warn("Thawing state\n");

    READ_STR(buf, p, sizeof(hdr) - 1);

    READ_INT(nbody, p);
    READ_SIZE_T(realSize, p);

    READ_REAL(st->tnow, p);
    READ_REAL(st->tree.rsize, p);

    if (strncmp(hdr, buf, sizeof(hdr) - 1))
    {
        warn("Didn't find header for checkpoint file.\n");
        failed = TRUE;
    }

    if (ctx->nbody != nbody)
    {
        warn("Number of bodies in checkpoint file (%u) "
             "does not match number expected by context (%u).\n",
             nbody,
             ctx->nbody);
        failed = TRUE;
    }

    if (realSize != sizeof(real))
    {
        warn("Got checkpoint file for wrong type. "
             "Expected sizeof(real) = %lu, got %lu\n",
             (long unsigned int) sizeof(real),
             (long unsigned int) realSize);
        failed = TRUE;
    }

    /* Read the bodies */
    st->bodytab = (bodyptr) mwMalloc(bodySize);
    memcpy(st->bodytab, p, bodySize);
    p += bodySize;

    READ_STR(tailBuf, p, sizeof(tailBuf) - 1);
    if (strncmp(tail, tailBuf, sizeof(tailBuf) - 1))
    {
        warn("Failed to find end marker in checkpoint file.\n");
        failed = TRUE;
    }

    return failed;
}

/* Checkpoint file: Very simple binary "format"
   Name     Type    Values     Notes
-------------------------------------------------------
   header       string   "mwnbody"  No null terminator
   nbody        uint     anything   Num. of bodies expected. Error if doesn't match nbody in context.
   sizeof(real) size_t   4, 8       Does the checkpoint use float or double
                                    Saved parts of the program state
   tnow         real     anything
   rsize        real     anything
   bodytab      bodyptr  anything   Array of bodies
   ending       string   "end"      No null terminator
 */

/* Use a very simple flag to mark when writing the checkpoint file
 * begins and ends. I think this should always be good enough, unless
 * something really weird happens. If the read is interrupted, the
 * checkpoint file is garbage and we lose everything. Uses the boinc
 * critical sections, so it hopefully won't be interrupted.
 */
static inline void freezeState(const NBodyCtx* ctx, const NBodyState* st, CheckpointHandle* cp)
{
    const size_t bodySize = sizeof(body) * ctx->nbody;
    char* p = cp->mptr;

    /* -1 so we don't bother with the null terminator. It's slightly
        annoying since the strcmps use it, but memcpy doesn't. We
        don't need it anyway  */
    DUMP_STR(p, hdr, sizeof(hdr) - 1);  /* Simple marker for a checkpoint file */
    DUMP_INT(p, ctx->nbody);        /* Make sure we get the right number of bodies */
    DUMP_SIZE_T(p, sizeof(real));   /* Make sure we don't confuse double and float checkpoints */

    /* Now that we have some basic check stuff written, dump the state */

    /* Little state pieces */
    DUMP_REAL(p, st->tnow);
    DUMP_REAL(p, st->tree.rsize);

    /* The main piece of state*/
    memcpy(p, st->bodytab, bodySize);
    p += bodySize;

    DUMP_STR(p, tail, sizeof(tail) - 1);
}

/* Open the temporary checkpoint file for writing */
int resolveCheckpoint(NBodyCtx* ctx)
{
    int rc;
    rc = boinc_resolve_filename(ctx->cp_filename,
                                ctx->cp_resolved,
                                sizeof(ctx->cp_resolved));
    if (rc)
    {
        warn("Failed to resolve checkpoint file '%s': %d\n", ctx->cp_filename, rc);
        return rc;
    }

    return FALSE;
}


/* Read the actual checkpoint file to resume */
int readCheckpoint(const NBodyCtx* ctx, NBodyState* st)
{
    CheckpointHandle cp;

    if (openCheckpointHandle(ctx, &cp, ctx->cp_resolved))
    {
        warn("Opening checkpoint for resuming failed\n");
        return TRUE;
    }

    if (thawState(ctx, st, &cp))
    {
        warn("Thawing state failed\n");
        return TRUE;
    }

    if (closeCheckpointHandle(&cp))
        warn("Failed to close checkpoint properly\n");

    return FALSE;
}

int writeCheckpoint(const NBodyCtx* ctx, const NBodyState* st)
{
    int failed = FALSE;
    CheckpointHandle cp;

    if (openCheckpointHandle(ctx, &cp, CHECKPOINT_TMP_FILE))
    {
        warn("Failed to open temporary checkpoint file\n"
             "Failed to write checkpoint\n");
        return TRUE;
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
    if (!failed && mw_rename(CHECKPOINT_TMP_FILE, ctx->cp_resolved))
    {
        mw_win_perror("Failed to update checkpoint with temporary");
        failed = TRUE;
    }

    if (failed)
        warn("Failed to write checkpoint\n");

    return failed;
}

