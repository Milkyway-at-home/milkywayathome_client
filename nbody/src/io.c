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

#if BOINC_APPLICATION
  #include <boinc_api.h>
  #include <filesys.h>
#endif

#include <unistd.h>
#include <sys/mman.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

static const char hdr[] = "mwnbody";
static const char tail[] = "end";

/* Everything except the size of all the bodies */
const size_t hdrSize =   sizeof(size_t)                                  /* size of real */
                       + sizeof(char) * (sizeof(tail) + sizeof(hdr) - 2) /* error checking tags */
                       + 2 * sizeof(int)                                 /* nbody count + valid flag */
                       + 2 * sizeof(real);                               /* tout and tnow */

/* Macros to rea d/ write the buffer and advance the pointer the correct size */
#define DUMP_REAL(p, x) { *((real*) (p)) = (x); (p) += sizeof(real); }
#define DUMP_INT(p, x) { *((int*) (p)) = (x); (p) += sizeof(int); }
#define DUMP_SIZE_T(p, x) { *((size_t*) (p)) = (x); (p) += sizeof(size_t); }
#define DUMP_STR(p, x, size) { memcpy((p), (x), (size)); (p) += (size); }
#define SET_LOCK(lock, x) { *((int*) (lock)) = (x); msync((lock), sizeof(int), MS_SYNC);}

#define READ_REAL(x, p) { (x) = *((real*) (p)); (p) += sizeof(real); }
#define READ_INT(x, p) { (x) = *((int*) (p)); (p) += sizeof(int); }
#define READ_SIZE_T(x, p) { (x) = *((size_t*) (p)); (p) += sizeof(size_t); }
#define READ_STR(x, p, size) { memcpy((x), (p), (size)); (p) += (size); }

#ifndef _WIN32

void openCheckpoint(NBodyCtx* ctx)
{
    struct stat sb;
    const size_t checkpointFileSize = hdrSize + ctx->model.nbody * sizeof(body);

    /* TODO: Wuh wuh windows:
       http://msdn.microsoft.com/en-us/library/aa366556(VS.85).aspx
    */

    ctx->cp.fd = open(ctx->cp.file, O_RDWR | O_CREAT, S_IWUSR | S_IRUSR);
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
    {
        fprintf(stderr, "checkpoint file is not a file\n");
        nbody_finish(EXIT_FAILURE);
    }

    ctx->cp.mptr = mmap(0, sb.st_size, PROT_READ | PROT_WRITE, MAP_SHARED, ctx->cp.fd, 0);
    if (ctx->cp.mptr == MAP_FAILED)
    {
        perror("mmap: Failed to open checkpoint file for writing");
        nbody_finish(EXIT_FAILURE);
    }

}

#else

void openCheckpoint(NBodyCtx* ctx)
{
    const size_t checkpointFileSize = hdrSize + ctx->model.nbody * sizeof(body);

    HANDLE h;
    HANDLE map;

    /* Try to create a new file */
    ctx->cp.fd = CreateFile(ctx->cp.file,
                            GENERIC_READ | GENERIC_WRITE,
                            0,     /* Other processes can't touch this */
                            NULL,
                            CREATE_NEW,
                            FILE_ATTRIBUTE_TEMPORARY | FILE_FLAG_SEQUENTIAL_SCAN,
                            NULL);

    /* If the checkpoint already exists, open it */
    if ( GetLastError() == ERROR_FILE_EXISTS )
    {
        ctx->cp.fd = CreateFile(ctx->cp.file,
                                GENERIC_READ | GENERIC_WRITE,
                                0,     /* Other processes can't touch this */
                                NULL,
                                OPEN_EXISTING,
                                FILE_ATTRIBUTE_TEMPORARY | FILE_FLAG_SEQUENTIAL_SCAN,
                                NULL);
    }

    /* TODO: More filetype checking and stuff */

    if ( h == INVALID_HANDLE_VALUE )
    {
        fprintf(stderr, "Failed to open checkpoint file: %d\n", GetLastError());
        nbody_finish(EXIT_FAILURE);
    }

    ctx->cp.mptr = CreateFileMapping(h,
                                     NULL,
                                     PAGE_READWRITE,
                                     checkpointFileSize,
                                     checkpointFileSize,
                                     "checkpoint file");

    if ( ctx->cp.mptr == NULL )
    {
        fprintf(stderr, "Failed to creating mapping for checkpoint file: %d\n", GetLastError());
        nbody_finish(EXIT_FAILURE);
    }
}

#endif /* _WIN32 */

void initoutput(NBodyCtx* ctx)
{
    if (ctx->outfilename)                       /* output file specified? */
    {
        ctx->outfile = nbody_fopen(ctx->outfilename, "w");           /* setup output FILE* */
        if (ctx->outfile == NULL)
            fail("initoutput: cannot open file %s\n", ctx->outfilename);
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
                "Expected sizeof(real) = %zd, got %zd\n",
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
    st->bodytab = allocate(bodySize);

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

    printf("System freezing. tnow = %g\n", st->tnow);

    /* TODO: Better error checking */

    /* -1 so we don't bother with the null terminator. It's slightly
        annoying since the strcmps use it, but memcpy doesn't. We
        don't need it anyway  */
    DUMP_STR(p, hdr, sizeof(hdr) - 1);  /* Simple marker for a checkpoint file */

    lock = p;        /* We keep the lock here */
    p += sizeof(int);

    boinc_begin_critical_section();

    SET_LOCK(lock, 0);    /* Mark the file as in the middle of writing */

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

    msync(ctx->cp.mptr, hdrSize + bodySize, MS_SYNC);

    SET_LOCK(lock, 1);   /* Done writing, flag file as valid  */

    boinc_end_critical_section();
}

void nbody_boinc_output(const NBodyCtx* ctx, const NBodyState* st)
{
    if (boinc_time_to_checkpoint())
    {
        freezeState(ctx, st);
        boinc_checkpoint_completed();
    }

    boinc_fraction_done(st->tnow / ctx->model.time_dwarf);
}
#endif

void cartesianToLbr_rad(const NBodyCtx* ctx, vectorptr restrict lbR, const vectorptr restrict r)
{
    const real r0p = r[0] + ctx->sunGCDist;
    lbR[0] = ratan2(r[1], r0p);
    lbR[1] = ratan2(r[2], rsqrt(sqr(r0p) + sqr(r[1])));
    lbR[2] = rsqrt(sqr(r0p) + sqr(r[1]) + sqr(r[2]));

    if (lbR[0] < 0)
        lbR[0] += 2 * M_PI;
}

void cartesianToLbr(const NBodyCtx* ctx, vectorptr restrict lbR, const vectorptr restrict r)
{
    cartesianToLbr_rad(ctx, lbR, r);
    L(lbR) = r2d(L(lbR));
    B(lbR) = r2d(B(lbR));
}

#if DOUBLEPREC
  #define PRECSTRING "double"
#else
  #define PRECSTRING "float"
#endif

/* Output with the silly xml stuff that BOINC uses */
void boincOutput(const NBodyCtx* ctx, const NBodyState* st)
{
    fprintf(ctx->outfile, "<something>\n");
    output(ctx, st);
    fprintf(ctx->outfile, "</something>\n");
    fprintf(ctx->outfile, "<nbody_version>%s %s</nbody_version>\n", BOINC_NBODY_APP_VERSION, PRECSTRING);
}

/* OUTPUT: compute diagnostics and output data. */
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


/* A bunch of boilerplate for debug printing */

const char* showBool(const bool x)
{
    switch (x)
    {
        case FALSE:
            return "false";
        case TRUE:
            return "true";
        default:
            return "invalid boolean (but true)";
    }
}

const char* showCriterionT(const criterion_t x)
{
    static const char* table[] =
        {
            [EXACT]        = "exact",
            [NEWCRITERION] = "new criterion",
            [BH86]         = "bh86",
            [SW93]         = "sw93"
        };

    if (x > SW93)
        return "invalid criterion_t";
    else
        return table[x];

}

const char* showSphericalT(const spherical_t x)
{
    static const char* table[] =
        {
            [SphericalPotential]        = "SphericalPotential",
        };

    if (x > SphericalPotential)
        return "invalid spherical_t";
    else
        return table[x];
}

const char* showDiskT(const disk_t x)
{
    static const char* table[] =
        {
            [MiyamotoNagaiDisk] = "MiyamotoNagaiDisk",
            [ExponentialDisk]   = "ExponentialDisk"
        };

    if (x > ExponentialDisk)
        return "invalid disk_t";
    else
        return table[x];
}

const char* showHaloT(const halo_t x)
{
    static const char* table[] =
        {
            [LogarithmicHalo] = "LogarithmicHalo",
            [NFWHalo]         = "NFWHalo",
            [TriaxialHalo]    = "TriaxialHalo"
        };

    if (x > TriaxialHalo)
        return "invalid halo_t";
    else
        return table[x];
}

const char* showDwarfModelT(const dwarf_model_t x)
{
    static const char* table[] =
        {
            [DwarfModelPlummer] = "DwarfModelPlummer",
            [DwarfModelKing]    = "DwarfModelKing",
            [DwarfModelDehnen]  = "DwarfModelDehnen"
        };

    if (x > DwarfModelDehnen)
        return "invalid dwarf_model_t";
    else
        return table[x];
}

char* showSpherical(const Spherical* s)
{
    char* buf;

    if (0 > asprintf(&buf,
                     "{\n"
                     "      type  = %s\n"
                     "      mass  = %g\n"
                     "      scale = %g\n"
                     "    };\n",
                     showSphericalT(s->type),
                     s->mass,
                     s->scale))
    {
        fail("asprintf() failed\n");
    }

    return buf;
}

char* showHalo(const Halo* h)
{
    char* buf;

    if (0 > asprintf(&buf,
                     "{ \n"
                     "      type         = %s\n"
                     "      vhalo        = %g\n"
                     "      scale_length = %g\n"
                     "      flattenX     = %g\n"
                     "      flattenY     = %g\n"
                     "      flattenZ     = %g\n"
                     "      triaxAngle   = %g\n"
                     "    };\n",
                     showHaloT(h->type),
                     h->vhalo,
                     h->scale_length,
                     h->flattenX,
                     h->flattenY,
                     h->flattenZ,
                     h->triaxAngle))
    {
        fail("asprintf() failed\n");
    }

    return buf;
}

char* showDisk(const Disk* d)
{
    char* buf;

    if (0 > asprintf(&buf,
                     "{ \n"
                     "      type         = %s\n"
                     "      mass         = %g\n"
                     "      scale_length = %g\n"
                     "      scale_height = %g\n"
                     "    };\n",
                     showDiskT(d->type),
                     d->mass,
                     d->scale_length,
                     d->scale_height))
    {
        fail("asprintf() failed\n");
    }

    return buf;
}

/* For debugging. Need to make this go away for release since it uses
 * GNU extensions */
char* showPotential(const Potential* p)
{
    int rc;
    char* buf;
    char* sphBuf;
    char* diskBuf;
    char* haloBuf;

    sphBuf  = showSpherical(&p->sphere[0]);
    diskBuf = showDisk(&p->disk);
    haloBuf = showHalo(&p->halo);

    rc = asprintf(&buf,
                  "{\n"
                  "    sphere = %s\n"
                  "    disk = %s\n"
                  "    halo = %s\n"
                  "    rings  = { unused pointer %p }\n"
                  "  };\n",
                  sphBuf,
                  diskBuf,
                  haloBuf,
                  p->rings);

    if (rc < 0)
        fail("asprintf() failed\n");

    free(sphBuf);
    free(diskBuf);
    free(haloBuf);

    return buf;
}

char* showDwarfModel(const DwarfModel* d)
{
    char* buf;

    if (0 > asprintf(&buf,
                     "{ \n"
                     "      type           = %s\n"
                     "      nbody          = %d\n"
                     "      mass           = %g\n"
                     "      scale_radius   = %g\n"
                     "      timestep       = %g\n"
                     "      orbit_timestep = %g\n"
                     "      time_dwarf     = %g\n"
                     "      time_orbit     = %g\n"
                     "      eps            = %g\n"
                     "    };\n",
                     showDwarfModelT(d->type),
                     d->nbody,
                     d->mass,
                     d->scale_radius,
                     d->timestep,
                     d->orbit_timestep,
                     d->time_orbit,
                     d->time_dwarf,
                     d->eps))
    {
        fail("asprintf() failed\n");
    }

    return buf;
}

char* showInitialConditions(const InitialConditions* ic)
{
    char* buf;
    if (0 > asprintf(&buf,
                     "initial-conditions = { \n"
                     "  useGalC    = %s\n"
                     "  useRadians = %s\n"
                     "  position   = { %g, %g, %g }\n"
                     "  velocity   = { %g, %g, %g }\n"
                     "};\n",
                     showBool(ic->useGalC),
                     showBool(ic->useRadians),
                     ic->position[0],
                     ic->position[1],
                     ic->position[2],
                     ic->velocity[0],
                     ic->velocity[1],
                     ic->velocity[2]))
    {
        fail("asprintf() failed\n");
    }

    return buf;
}

char* showContext(const NBodyCtx* ctx)
{
    char* buf;
    char* potBuf;
    char* modelBuf;

    potBuf   = showPotential(&ctx->pot);
    modelBuf = showDwarfModel(&ctx->model);

    if (0 > asprintf(&buf,
                     "ctx = { \n"
                     "  pot             = %s\n"
                     "  model           = %s\n"
                     "  headline        = %s\n"
                     "  outfilename     = %s\n"
                     "  outfile         = %p\n"
                     "  sunGCDist       = %g\n"
                     "  criterion       = %s\n"
                     "  usequad         = %s\n"
                     "  allowIncest     = %s\n"
                     "  outputCartesian = %s\n"
                     "  seed            = %ld\n"
                     "  tree_rsize      = %g\n"
                     "  theta           = %g\n"
                     "  freq            = %g\n"
                     "  freqout         = %g\n"
                     "};\n",
                     potBuf,
                     modelBuf,
                     ctx->headline,
                     ctx->outfilename,
                     ctx->outfile,
                     ctx->sunGCDist,
                     showCriterionT(ctx->criterion),
                     showBool(ctx->usequad),
                     showBool(ctx->allowIncest),
                     showBool(ctx->outputCartesian),
                     ctx->seed,
                     ctx->tree_rsize,
                     ctx->theta,
                     ctx->freq,
                     ctx->freqout))
    {
        fail("asprintf() failed\n");
    }

    free(potBuf);
    free(modelBuf);

    return buf;
}

void printContext(const NBodyCtx* ctx)
{
    char* buf = showContext(ctx);
    puts(buf);
    free(buf);
}

void printInitialConditions(const InitialConditions* ic)
{
    char* buf = showInitialConditions(ic);
    puts(buf);
    free(buf);
}

char* showVector(const vector v)
{
    char* buf;

    if (asprintf(&buf, "{ %g, %g, %g }", v[0], v[1], v[2]) < 0)
        fail("asprintf() failed\n");

    return buf;

}

