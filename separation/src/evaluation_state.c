/*
Copyright 2008-2010 Travis Desell, Matthew Arsenault, Dave Przybylo,
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

#include "milkyway_util.h"
#include "separation.h"
#include "evaluation_state.h"


static char resolvedCheckpointPath[2048];


void initializeCut(Cut* integral, unsigned int number_streams)
{
    integral->bgIntegral = 0.0;
    integral->streamIntegrals = (real*) mwCallocA(number_streams, sizeof(real));
}

static void initializeState(const AstronomyParameters* ap, EvaluationState* es)
{
    unsigned int i;

    es->currentCut = 0;
    es->cut = &es->cuts[0];
    es->numberStreams = ap->number_streams;

    es->numberCuts = ap->number_integrals;
    es->cuts = (Cut*) mwCallocA(ap->number_integrals, sizeof(Cut));
    es->streamSums = (Kahan*) mwCallocA(ap->number_streams, sizeof(Kahan));
    es->streamTmps = (real*) mwCallocA(ap->number_streams, sizeof(real));

    for (i = 0; i < ap->number_integrals; i++)
        initializeCut(&es->cuts[i], ap->number_streams);
}

EvaluationState* newEvaluationState(const AstronomyParameters* ap)
{
    EvaluationState* es;

    es = mwCallocA(1, sizeof(EvaluationState));
    initializeState(ap, es);

    return es;
}

void copyEvaluationState(EvaluationState* esDest, const EvaluationState* esSrc)
{
    unsigned int i;

    assert(esDest && esSrc);

    *esDest = *esSrc;
    for (i = 0; i < esSrc->numberCuts; ++i)
        esDest->cuts[i] = esSrc->cuts[i];
}

static void freeCut(Cut* i)
{
    mwFreeA(i->streamIntegrals);
}

void freeEvaluationState(EvaluationState* es)
{
    unsigned int i;

    for (i = 0; i < es->numberCuts; ++i)
        freeCut(&es->cuts[i]);
    mwFreeA(es->cuts);
    mwFreeA(es->streamSums);
    mwFreeA(es->streamTmps);
    mwFreeA(es);
}

void clearEvaluationStateTmpSums(EvaluationState* es)
{
    unsigned int i;

    CLEAR_KAHAN(es->bgSum);
    for (i = 0; i < es->numberStreams; ++i)
        CLEAR_KAHAN(es->streamSums[i]);
}

void printEvaluationState(const EvaluationState* es)
{
    Cut* c;
    unsigned int j;

    printf("evaluation-state {\n"
           "  nu_step          = %u\n"
           "  mu_step          = %u\n"
           "  currentCut       = %u\n",
           es->nu_step,
           es->mu_step,
           es->currentCut);

    printf("  bgSum = { %25.15f, %25.15f }\n"
           "  bgTmp = %25.15f\n",
           es->bgSum.sum,
           es->bgSum.correction,
           es->bgTmp);

    for (j = 0; j < es->numberStreams; ++j)
    {
        printf("  streamSums[%u] = { %25.15f, %25.15f }\n"
               "  streamTmps[%u] = %25.15f\n",
               j, es->streamSums[j].sum, es->streamSums[j].correction,
               j, es->streamTmps[j]);
    }

    for (c = es->cuts; c < es->cuts + es->numberCuts; ++c)
    {
        printf("integral: bgIntegral = %20.15f\n", c->bgIntegral);
        printf("Stream integrals = ");
        for (j = 0; j < es->numberCuts; ++j)
            printf("  %20.15f, ", c->streamIntegrals[j]);
    }
    printf("\n");
}

static const char checkpoint_header[] = "separation_checkpoint";
static const char checkpoint_tail[] = "end_checkpoint";

typedef struct
{
    int major, minor, cl, cal;
} SeparationVersionHeader;

static const SeparationVersionHeader versionHeader =
{
    SEPARATION_VERSION_MAJOR,
    SEPARATION_VERSION_MINOR,
    SEPARATION_OPENCL,
    SEPARATION_CAL
};


static int versionMismatch(const SeparationVersionHeader* v)
{
    if (   v->major != SEPARATION_VERSION_MAJOR
        || v->minor != SEPARATION_VERSION_MINOR
        || v->cl    != SEPARATION_OPENCL
        || v->cal   != SEPARATION_CAL)
    {
        warn("Checkpoint version does not match:\n"
             "  Expected %d.%d, OpenCL = %d, CAL++ = %d,\n"
             "  Got %d.%d, OpenCL = %d, CAL++ = %d\n",
             SEPARATION_VERSION_MAJOR, SEPARATION_VERSION_MINOR, SEPARATION_OPENCL, SEPARATION_CAL,
             v->major, v->minor, v->cl, v->cal);
        return 1;
    }

    return 0;
}

static int readState(FILE* f, EvaluationState* es)
{
    SeparationVersionHeader version;
    Cut* c;
    char str_buf[sizeof(checkpoint_header) + 1];

    fread(str_buf, sizeof(checkpoint_header), 1, f);
    if (strncmp(str_buf, checkpoint_header, sizeof(str_buf)))
    {
        warn("Failed to find header in checkpoint file\n");
        return 1;
    }

    fread(&version, sizeof(version), 1, f);
    if (versionMismatch(&version))
        return 1;

    fread(&es->currentCut, sizeof(es->currentCut), 1, f);
    fread(&es->nu_step, sizeof(es->nu_step), 1, f);
    fread(&es->mu_step, sizeof(es->mu_step), 1, f);
    fread(&es->lastCheckpointNuStep, sizeof(es->lastCheckpointNuStep), 1, f);

    fread(&es->bgSum, sizeof(es->bgSum), 1, f);
    fread(es->streamSums, sizeof(es->streamSums[0]), es->numberStreams, f);

    fread(&es->bgTmp, sizeof(es->bgTmp), 1, f);
    fread(es->streamTmps, sizeof(es->streamTmps[0]), es->numberStreams, f);

    for (c = es->cuts; c < es->cuts + es->numberCuts; ++c)
    {
        fread(&c->bgIntegral, sizeof(c->bgIntegral), 1, f);
        fread(c->streamIntegrals, sizeof(c->streamIntegrals[0]), es->numberStreams, f);
    }

    fread(str_buf, sizeof(checkpoint_tail), 1, f);
    if (strncmp(str_buf, checkpoint_tail, sizeof(str_buf)))
    {
        warn("Failed to find tail in checkpoint file\n");
        return 1;
    }

    return 0;
}


/* The GPU checkpointing saves the sum of all resumes in the real integral.
   The temporary area holds the checkpointed sum of the last episode.
   Update the real integral from the results of the last episode.
*/
void addTmpSums(EvaluationState* es)
{
    unsigned int i;

    KAHAN_ADD(es->bgSum, es->bgTmp);
    for (i = 0; i < es->numberStreams; ++i)
        KAHAN_ADD(es->streamSums[i], es->streamTmps[i]);

    es->bgTmp = 0.0;
    memset(es->streamTmps, 0, sizeof(es->streamTmps[0]) * es->numberStreams);
}

int readCheckpoint(EvaluationState* es)
{
    int rc;
    FILE* f;

    f = mwOpenResolved(CHECKPOINT_FILE, "rb");
    if (!f)
    {
        perror("Opening checkpoint");
        return 1;
    }

    rc = readState(f, es);
    if (rc)
        warn("Failed to read state\n");

    fclose(f);

    /* The GPU checkpoint saved the last checkpoint results in temporaries.
       Add to the total from previous episodes. */
    addTmpSums(es);

    return rc;
}

static inline void writeState(FILE* f, const EvaluationState* es)
{
    Cut* c;
    const Cut* endc = es->cuts + es->numberCuts;

    fwrite(checkpoint_header, sizeof(checkpoint_header), 1, f);
    fwrite(&versionHeader, sizeof(versionHeader), 1, f);

    fwrite(&es->currentCut, sizeof(es->currentCut), 1, f);
    fwrite(&es->nu_step, sizeof(es->nu_step), 1, f);
    fwrite(&es->mu_step, sizeof(es->mu_step), 1, f);
    fwrite(&es->lastCheckpointNuStep, sizeof(es->lastCheckpointNuStep), 1, f);

    fwrite(&es->bgSum, sizeof(es->bgSum), 1, f);
    fwrite(es->streamSums, sizeof(es->streamSums[0]), es->numberStreams, f);

    /* Really only needs to be saved for GPU checkpointing */
    fwrite(&es->bgTmp, sizeof(es->bgTmp), 1, f);
    fwrite(es->streamTmps, sizeof(es->streamTmps[0]), es->numberStreams, f);

    for (c = es->cuts; c < endc; ++c)
    {
        fwrite(&c->bgIntegral, sizeof(c->bgIntegral), 1, f);
        fwrite(c->streamIntegrals, sizeof(c->streamIntegrals[0]), es->numberStreams, f);
    }

    fwrite(checkpoint_tail, sizeof(checkpoint_tail), 1, f);
}


#if BOINC_APPLICATION

/* Each checkpoint we introduce more errors from summing the entire
 * buffer. If we didn't we would have checkpoints hundreds of
 * megabytes in size. Make sure at least 10% has progressed (limiting
 * to ~10 max checkpoints over the summation). This keeps the
 * introduced error below acceptable levels. Additionally this
 * checkpointing isn't exactly cheap (a checkpoint is taking nearly as
 * long as an entire step on the 5870 for me), so we really want to
 * avoid doing it too often.*/
int timeToCheckpointGPU(const EvaluationState* es, const IntegralArea* ia)
{
    return (es->nu_step - es->lastCheckpointNuStep >= ia->nu_steps / 10) && boinc_time_to_checkpoint();
}

#else

int timeToCheckpointGPU(const EvaluationState* es, const IntegralArea* ia)
{
    return FALSE;
}

#endif /* BOINC_APPLICATION */


int resolveCheckpoint()
{
    int rc;

    rc = mw_resolve_filename(CHECKPOINT_FILE, resolvedCheckpointPath, sizeof(resolvedCheckpointPath));
    if (rc)
        warn("Error resolving checkpoint file '%s': %d\n", CHECKPOINT_FILE, rc);
    return rc;
}

int writeCheckpoint(EvaluationState* es)
{
    FILE* f;

    /* Avoid corrupting the checkpoint file by writing to a temporary file, and moving that */
    f = mw_fopen(CHECKPOINT_FILE_TMP, "wb");
    if (!f)
    {
        perror("Opening checkpoint temp");
        return 1;
    }

    es->lastCheckpointNuStep = es->nu_step;
    writeState(f, es);
    fclose(f);

    if (mw_rename(CHECKPOINT_FILE_TMP, resolvedCheckpointPath))
    {
        perror("Failed to update checkpoint file");
        return 1;
    }

    return 0;
}

int maybeResume(EvaluationState* es)
{
    if (mw_file_exists(resolvedCheckpointPath))
    {
        mw_report("Checkpoint exists. Attempting to resume from it\n");

        if (readCheckpoint(es))
        {
            mw_report("Reading checkpoint failed\n");
            mw_remove(CHECKPOINT_FILE);
            return 1;
        }
        else
            mw_report("Successfully resumed checkpoint\n");
    }

    return 0;
}

