/* ************************************************************************** */
/* IO.C: I/O routines for export version of hierarchical N-body code. */
/* Public routines: inputdata(), initoutput(), stopoutput(), output(). */
/* */
/* Copyright (c) 1993 by Joshua E. Barnes, Honolulu, HI. */
/* It's free because it's yours. */
/* ************************************************************************** */

#include <string.h>
#include "nbody_priv.h"
#include "nbody_util.h"
#include "io.h"
#include "nbody_boinc.h"

#ifndef _WIN32
  #include <unistd.h>
  #include <sys/mman.h>
  #include <sys/types.h>
  #include <sys/stat.h>
  #include <fcntl.h>
#else
  #include <windows.h>
#endif /* _WIN32 */


static const char hdr[] = "mwnbody";
static const char tail[] = "end";

/* Everything except the size of all the bodies */
const size_t hdrSize =   sizeof(size_t)                                  /* size of real */
                       + sizeof(char) * (sizeof(tail) + sizeof(hdr) - 2) /* error checking tags */
                       + 2 * sizeof(int)                                 /* nbody count + valid flag */
                       + 2 * sizeof(real);                               /* tout and tnow */

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

/* The ctx is ignored here. It's only there for windows. */
#define SET_LOCK(lock, x, ctx) { *((int*) (lock)) = (x); msync((lock), sizeof(int), MS_SYNC);}

#define SYNC_WRITE(ctx, size)  msync(ctx->cp.mptr, (size), MS_SYNC)

void openCheckpoint(NBodyCtx* ctx)
{
    int rc;
    char resolvedPath[1024];
    struct stat sb;
    const size_t checkpointFileSize = hdrSize + ctx->model.nbody * sizeof(body);

    rc = boinc_resolve_filename(ctx->cp.filename, resolvedPath, sizeof(resolvedPath));
    if (rc)
        fail("Error resolving checkpoint file '%s': %d\n", ctx->cp.filename, rc);

    ctx->cp.fd = open(resolvedPath, O_RDWR | O_CREAT, S_IWUSR | S_IRUSR);
    if (ctx->cp.fd == -1)
    {
        perror("open checkpoint");
        nbody_finish(EXIT_FAILURE);
    }

    /* Make the file the right size in case it's a new file */
    ftruncate(ctx->cp.fd, checkpointFileSize);

    if (fstat(ctx->cp.fd, &sb) == -1)
    {
        perror("fstat");
        nbody_finish(EXIT_FAILURE);
    }

    if (!S_ISREG(sb.st_mode))
        fail("checkpoint file is not a file\n");

    ctx->cp.mptr = mmap(0, sb.st_size, PROT_READ | PROT_WRITE, MAP_SHARED, ctx->cp.fd, 0);
    if (ctx->cp.mptr == MAP_FAILED)
    {
        perror("mmap: Failed to open checkpoint file for writing");
        nbody_finish(EXIT_FAILURE);
    }

}

void closeCheckpoint(NBodyCtx* ctx)
{
    struct stat sb;

    /* Clean up the checkpointing */
    if (ctx->cp.fd != -1)
    {
        if (fstat(ctx->cp.fd, &sb) == -1)
        {
            perror("fstat on closing checkpoint");
            nbody_finish(EXIT_FAILURE);
        }

        if (close(ctx->cp.fd) == -1)
        {
            perror("closing checkpoint file");
            nbody_finish(EXIT_FAILURE);
        }

        if (munmap(ctx->cp.mptr, sb.st_size) == -1)
        {
            perror("munmap");
            nbody_finish(EXIT_FAILURE);
        }
    }

    warn("Removing checkpoint file '%s'\n", ctx->cp.filename);
    nbody_remove(ctx->cp.filename);
}

#else

/* Relevant: http://msdn.microsoft.com/en-us/library/aa366556(VS.85).aspx
             http://msdn.microsoft.com/en-us/library/aa366548(VS.85).aspx
             http://msdn.microsoft.com/en-us/library/ms810613.aspx

             Flushing:
             http://msdn.microsoft.com/en-us/library/aa366563(v=VS.85).aspx
 */


/* CHECKME: Do we need to do the file handle, the mapFile handle, or both? */
/* TODO: Check that these actually succeed */
#define SET_LOCK(lock, x, ctx)                          \
    {                                                   \
        *((int*) (lock)) = (x);                         \
        FlushViewOfFile((lock), sizeof(int));           \
        FlushFileBuffers(((NBodyCtx*) ctx)->cp.file);   \
    }

#define SYNC_WRITE(ctx, size)                                   \
    {                                                           \
        FlushViewOfFile(((NBodyCtx*) ctx)->cp.mptr, (size));    \
        FlushFileBuffers(((NBodyCtx*) ctx)->cp.file);           \
    }

void openCheckpoint(NBodyCtx* ctx)
{
    const DWORD checkpointFileSize = hdrSize + ctx->model.nbody * sizeof(body);

    SYSTEM_INFO si;
    DWORD sysGran;
    DWORD mapViewSize;
    DWORD fileMapStart;
    DWORD fileMapSize;

    /* Try to create a new file */
    ctx->cp.file = CreateFile(ctx->cp.filename,
                              GENERIC_READ | GENERIC_WRITE,
                              0,     /* Other processes can't touch this */
                              NULL,
                              CREATE_NEW,
                              FILE_ATTRIBUTE_TEMPORARY | FILE_FLAG_SEQUENTIAL_SCAN,
                              NULL);

    /* If the checkpoint already exists, open it */
    if ( GetLastError() == ERROR_FILE_EXISTS )
    {
        ctx->cp.file = CreateFile(ctx->cp.filename,
                                  GENERIC_READ | GENERIC_WRITE,
                                  0,     /* Other processes can't touch this */
                                  NULL,
                                  OPEN_EXISTING,
                                  FILE_ATTRIBUTE_TEMPORARY | FILE_FLAG_SEQUENTIAL_SCAN,
                                  NULL);
    }

    /* TODO: More filetype checking and stuff */

    if ( ctx->cp.file == INVALID_HANDLE_VALUE )
        fail( "Failed to open checkpoint file '%s': %ld\n", ctx->cp.filename, GetLastError());

    GetSystemInfo(&si);
    sysGran = si.dwAllocationGranularity;
    fileMapStart = 0;
    mapViewSize = checkpointFileSize;
    fileMapSize = checkpointFileSize;

    ctx->cp.mapFile = CreateFileMapping(ctx->cp.file,
                                        NULL,
                                        PAGE_READWRITE,
                                        0,
                                        fileMapSize,
                                        "nbody checkpoint file");

    if ( ctx->cp.mapFile == NULL )
    {
        fail("Failed to creating mapping for checkpoint file '%s': %ld\n",
             ctx->cp.filename, GetLastError());
    }

    ctx->cp.mptr = (char*) MapViewOfFile(ctx->cp.mapFile,
                                         FILE_MAP_ALL_ACCESS,
                                         0,
                                         fileMapStart,
                                         mapViewSize);
    if ( ctx->cp.mptr == NULL )
    {
        fail("Failed to open checkpoint file view for file '%s': %ld\n",
             ctx->cp.filename, GetLastError());
    }

}

void closeCheckpoint(NBodyCtx* ctx)
{
    if ( ctx->cp.file != INVALID_HANDLE_VALUE )
    {
        if (!UnmapViewOfFile((LPVOID) ctx->cp.mptr))
        {
            fail("Error %ld occurred unmapping the checkpoint view object!\n",
                 GetLastError());
        }

        if (!CloseHandle(ctx->cp.mapFile))
        {
            fail("Error %ld occurred closing the checkpoint mapping!\n",
                 GetLastError());
        }

        if (!CloseHandle(ctx->cp.file))
        {
            fail("Error %ld occurred closing the checkpoint file '%s'\n",
                 GetLastError(), ctx->cp.filename);
        }
    }

    warn("Removing checkpoint file '%s'\n", ctx->cp.filename);
    nbody_remove(ctx->cp.filename);
}

#endif /* _WIN32 */

void initOutput(NBodyCtx* ctx)
{
    if (ctx->outfilename)                       /* output file specified? */
    {
        ctx->outfile = nbody_fopen(ctx->outfilename, "w");           /* setup output FILE* */
        if (ctx->outfile == NULL)
            fail("initOutput: cannot open file %s\n", ctx->outfilename);
    }
    else
    {
        /* Using stderr annoys me, so only do this if we have to BOINC */
      #if BOINC_APPLICATION && !BOINC_DEBUG
        ctx->outfile = stderr;
      #else
        ctx->outfile = stdout;
      #endif /* BOINC_APPLICATION */

    }
}

/* Low-level input and output operations. */

static void out_2vectors(FILE* str, vector vec1, vector vec2)
{
    fprintf(str, " %21.14E %21.14E %21.14E %21.14E %21.14E %21.14E\n", vec1[0], vec1[1], vec1[2], vec2[0], vec2[1], vec2[2]);
}

#if BOINC_APPLICATION
/* Should be given the same context as the dump. Returns nonzero if the state failed to be thawed */
int thawState(const NBodyCtx* ctx, NBodyState* st)
{
    const size_t bodySize = ctx->model.nbody * sizeof(body);

    int failed = FALSE;

    int nbody;
    size_t realSize;
    char buf[sizeof(hdr)];
    char tailBuf[sizeof(tail)];
    char* p = ctx->cp.mptr;
    int valid;

    warn("Thawing state\n");

    READ_STR(buf, p, sizeof(hdr) - 1);
    READ_INT(valid, p);

    READ_INT(nbody, p);
    READ_SIZE_T(realSize, p);

    READ_REAL(st->tout, p);
    READ_REAL(st->tnow, p);
    READ_REAL(st->tree.rsize, p);

    /* TODO: Better checking of things */
    if (strncmp(hdr, buf, sizeof(hdr) - 1))
    {
        warn("Didn't find header for checkpoint file.\n");
        failed = TRUE;
    }

    if (ctx->model.nbody != nbody)
    {
        warn("Number of bodies in checkpoint file does not match number expected by context.\n");
        failed = TRUE;
    }

    if (realSize != sizeof(real))
    {
        warn("Got checkpoint file for wrong type. "
                "Expected sizeof(real) = %lu, got %lu\n",
                sizeof(real),
                realSize);
        failed = TRUE;
    }

    if (!valid)
    {
        warn("Trying to read interrupted checkpoint file\n");
        failed = TRUE;
    }

    /* Read the bodies */
    st->bodytab = mallocSafe(bodySize);

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
   header  string   "mwnbody"  No null terminator
   lock    int      0 or 1     If 0, the checkpoint file is in the middle of a write and cannot be used.
   nbody   int      anything   Number of bodies expected in the file. Error if doesn't match nbody in reading context.
   tout    real     anything   Saved parts of the program state
   tnow    real     anything
   rsize   real     anything
   bodytab bodyptr  anything   Array of bodies
   ending  string   "end"      No null terminator
 */

/* Use a very simple flag to mark when writing the checkpoint file
 * begins and ends. I think this should always be good enough, unless
 * something really weird happens. If the read is interrupted, the
 * checkpoint file is garbage and we lose everything. Uses the boinc
 * critical sections, so it hopefully won't be interrupted. I'm not
 * entirely sure this is the proper use of the boinc critical
 * sections. */
inline static void freezeState(const NBodyCtx* ctx, const NBodyState* st)
{
    const size_t bodySize = sizeof(body) * ctx->model.nbody;
    char* p = ctx->cp.mptr;
    char* lock;

    warn("System freezing. tnow = %g\n", st->tnow);

    /* TODO: Better error checking */

    /* -1 so we don't bother with the null terminator. It's slightly
        annoying since the strcmps use it, but memcpy doesn't. We
        don't need it anyway  */
    DUMP_STR(p, hdr, sizeof(hdr) - 1);  /* Simple marker for a checkpoint file */

    lock = p;        /* We keep the lock here */
    p += sizeof(int);

    boinc_begin_critical_section();

    SET_LOCK(lock, 0, ctx);         /* Mark the file as in the middle of writing */

    DUMP_INT(p, ctx->model.nbody);  /* Make sure we get the right number of bodies */
    DUMP_SIZE_T(p, sizeof(real));   /* Make sure we don't confuse double and float checkpoints */

    /* Now that we have some basic check stuff written, dump the state */

    /* Little state pieces */
    DUMP_REAL(p, st->tout);
    DUMP_REAL(p, st->tnow);
    DUMP_REAL(p, st->tree.rsize);

    /* The main piece of state*/
    memcpy(p, st->bodytab, bodySize);
    p += bodySize;

    DUMP_STR(p, tail, sizeof(tail) - 1);

    SYNC_WRITE(ctx, hdrSize + bodySize);

    SET_LOCK(lock, 1, ctx);   /* Done writing, flag file as valid  */

    boinc_end_critical_section();
}

void nbodyCheckpoint(const NBodyCtx* ctx, const NBodyState* st)
{
    if (boinc_time_to_checkpoint())
    {
        freezeState(ctx, st);
        boinc_checkpoint_completed();
    }

    boinc_fraction_done(st->tnow / ctx->model.time_dwarf);
}
#endif

/* Output with the silly xml stuff that BOINC uses */
void boincOutput(const NBodyCtx* ctx, const NBodyState* st, const real chisq)
{
//    not printing out the bodies because it will food the server.
//
//    fprintf(ctx->outfile, "<bodies>\n");
//    output(ctx, st);
//    fprintf(ctx->outfile, "</bodies>\n");

    fprintf(ctx->outfile, "<search_likelihood>%.20g</search_likelihood>\n", chisq);
    fprintf(ctx->outfile, "<search_application>" BOINC_NBODY_APP_VERSION "</search_application>\n");
}

/* output: compute diagnostics and output data. */
void output(const NBodyCtx* ctx, const NBodyState* st)
{
    bodyptr p;
    vector lbR;
    const bodyptr endp = st->bodytab + ctx->model.nbody;

    for (p = st->bodytab; p < endp; p++)
    {
        if (ctx->outputCartesian)     /* Probably useful for making movies and such */
            out_2vectors(ctx->outfile, Pos(p), Vel(p));
        else
        {
            cartesianToLbr(ctx, lbR, Pos(p));
            out_2vectors(ctx->outfile, lbR, Vel(p));
        }
    }

    fflush(ctx->outfile);             /* drain output buffer */
}

void nbodyCtxDestroy(NBodyCtx* ctx)
{
    free(ctx->headline);
    if (ctx->outfile && ctx->outfile != stdout)
        fclose(ctx->outfile);
}

static void freeTree(Tree* t)
{
    nodeptr p, tmp;

    p = (nodeptr) t->root;
    while (p != NULL)
    {
        if (Type(p) == CELL)
        {
            tmp = More(p);
            free(p);
            p = tmp;
        }
        else                        /* skip over bodies */
            p = Next(p);
    }

    t->root = NULL;
    t->cellused = 0;
}

void nbodyStateDestroy(NBodyState* st)
{
    freeTree(&st->tree);
    free(st->bodytab);
    free(st->acctab);

  #if NBODY_OPENCL
    cleanupNBodyCL(st);
  #endif /* NBODY_OPENCL */
}


