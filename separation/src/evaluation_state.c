/*
 *  Copyright (c) 2008-2010 Travis Desell, Nathan Cole, Dave Przybylo
 *  Copyright (c) 2008-2010 Boleslaw Szymanski, Heidi Newberg
 *  Copyright (c) 2008-2010 Carlos Varela, Malik Magdon-Ismail
 *  Copyright (c) 2008-2011 Rensselaer Polytechnic Institute
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

#include "milkyway_util.h"
#include "separation.h"
#include "evaluation_state.h"


static char resolvedCheckpointPath[4096];

int integralsAreDone(const EvaluationState* es)
{
    return (es->currentCut >= es->numberCuts);
}


void initializeCut(Cut* integral, unsigned int number_streams)
{
    integral->bgIntegral = 0.0;
    integral->streamIntegrals = (real*) mwCallocA(number_streams, sizeof(real));
}

static void initializeState(const AstronomyParameters* ap, EvaluationState* es)
{
    int i;
    
    es->currentWU = ap->currentWU;
    es->WUPrinted = 0;
    es->currentCut = 0;
    es->cut = &es->cuts[0];
    es->numberStreams = ap->number_streams;

    es->numberCuts = ap->number_integrals;
    es->cuts = (Cut*) mwCallocA(ap->number_integrals, sizeof(Cut));
    es->streamSums = (Kahan*) mwCallocA(ap->number_streams, sizeof(Kahan));
    es->streamTmps = (real*) mwCallocA(ap->number_streams, sizeof(real));

    es->streamSumsCheckpoint = (Kahan*) mwCallocA(ap->number_streams, sizeof(Kahan));

    for (i = 0; i < ap->number_integrals; i++)
    {
        initializeCut(&es->cuts[i], ap->number_streams);
    }
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
    int i;

    assert(esDest && esSrc);

    *esDest = *esSrc;
    for (i = 0; i < esSrc->numberCuts; ++i)
    {
        esDest->cuts[i] = esSrc->cuts[i];
    }
}

static void freeCut(Cut* i)
{
    mwFreeA(i->streamIntegrals);
}

void freeEvaluationState(EvaluationState* es)
{
    int i;

    for (i = 0; i < es->numberCuts; ++i)
    {
        freeCut(&es->cuts[i]);
    }

    mwFreeA(es->cuts);
    mwFreeA(es->streamSums);
    mwFreeA(es->streamTmps);
    mwFreeA(es->streamSumsCheckpoint);
    mwFreeA(es);
}

void clearEvaluationStateTmpSums(EvaluationState* es)
{
    int i;

    CLEAR_KAHAN(es->bgSum);
    for (i = 0; i < es->numberStreams; ++i)
    {
        CLEAR_KAHAN(es->streamSums[i]);
    }
}

void printEvaluationState(const EvaluationState* es)
{
    Cut* c;
    int j;

    printf("evaluation-state {\n"
           "  currentWU        = %u\n"
           "  WUPrinted        = %u\n"
           "  nu_step          = %u\n"
           "  mu_step          = %u\n"
           "  currentCut       = %u\n",
           es->currentWU,
           es->WUPrinted,
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

void addTmpCheckpointSums(EvaluationState* es)
{
    int i;

    /* Add the GPU temporary results from previous episodes */
    KAHAN_REDUCTION(es->bgSum, es->bgSumCheckpoint);

    for (i = 0; i < es->numberStreams; ++i)
    {
        KAHAN_REDUCTION(es->streamSums[i], es->streamSumsCheckpoint[i]);
    }

    memset(&es->bgSumCheckpoint, 0, sizeof(Kahan));
    memset(es->streamSumsCheckpoint, 0, es->numberStreams * sizeof(Kahan));
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
    FALSE
};


static int versionMismatch(const SeparationVersionHeader* v)
{
    if (v->major != SEPARATION_VERSION_MAJOR || v->minor != SEPARATION_VERSION_MINOR)
    {
        mw_printf("Checkpoint version does not match:\n"
                  "  Expected %d.%d, got %d.%d, \n",
                  SEPARATION_VERSION_MAJOR, SEPARATION_VERSION_MINOR,
                  v->major, v->minor);
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
        mw_printf("Failed to find header in checkpoint file\n");
        return 1;
    }

    fread(&version, sizeof(version), 1, f);
    if (versionMismatch(&version))
        return 1;

    fread(&es->currentWU, sizeof(es->currentWU), 1, f);
    fread(&es->WUPrinted, sizeof(es->WUPrinted), 1, f);
    fread(&es->currentCut, sizeof(es->currentCut), 1, f);
    fread(&es->nu_step, sizeof(es->nu_step), 1, f);
    fread(&es->mu_step, sizeof(es->mu_step), 1, f);

    fread(&es->bgSum, sizeof(es->bgSum), 1, f);
    fread(es->streamSums, sizeof(es->streamSums[0]), es->numberStreams, f);

    fread(&es->lastCheckpointNuStep, sizeof(es->lastCheckpointNuStep), 1, f);
    fread(&es->bgSumCheckpoint, sizeof(es->bgSumCheckpoint), 1, f);
    fread(es->streamSumsCheckpoint, sizeof(es->streamSumsCheckpoint[0]), es->numberStreams, f);

    for (c = es->cuts; c < es->cuts + es->numberCuts; ++c)
    {
        fread(&c->bgIntegral, sizeof(c->bgIntegral), 1, f);
        fread(c->streamIntegrals, sizeof(c->streamIntegrals[0]), es->numberStreams, f);
    }

    fread(str_buf, sizeof(checkpoint_tail), 1, f);
    if (strncmp(str_buf, checkpoint_tail, sizeof(str_buf)))
    {
        mw_printf("Failed to find tail in checkpoint file\n");
        return 1;
    }

    return 0;
}

int readCheckpoint(EvaluationState* es)
{
    int rc;
    FILE* f;

    f = mwOpenResolved(CHECKPOINT_FILE, "rb");
    if (!f)
    {
        mwPerror("Opening checkpoint '%s'", CHECKPOINT_FILE);
        return 1;
    }

    rc = readState(f, es);
    if (rc)
        mw_printf("Failed to read state\n");

    fclose(f);

    addTmpCheckpointSums(es);

    return rc;
}

static inline void writeState(FILE* f, const EvaluationState* es)
{
    Cut* c;
    const Cut* endc = es->cuts + es->numberCuts;

    fwrite(checkpoint_header, sizeof(checkpoint_header), 1, f);
    fwrite(&versionHeader, sizeof(versionHeader), 1, f);

    fwrite(&es->currentWU, sizeof(es->currentWU), 1, f);
    fwrite(&es->WUPrinted, sizeof(es->WUPrinted), 1, f);
    fwrite(&es->currentCut, sizeof(es->currentCut), 1, f);
    fwrite(&es->nu_step, sizeof(es->nu_step), 1, f);
    fwrite(&es->mu_step, sizeof(es->mu_step), 1, f);

    fwrite(&es->bgSum, sizeof(es->bgSum), 1, f);
    fwrite(es->streamSums, sizeof(es->streamSums[0]), es->numberStreams, f);

    fwrite(&es->lastCheckpointNuStep, sizeof(es->lastCheckpointNuStep), 1, f);
    fwrite(&es->bgSumCheckpoint, sizeof(es->bgSumCheckpoint), 1, f);
    fwrite(es->streamSumsCheckpoint, sizeof(es->streamSumsCheckpoint[0]), es->numberStreams, f);

    for (c = es->cuts; c < endc; ++c)
    {
        fwrite(&c->bgIntegral, sizeof(c->bgIntegral), 1, f);
        fwrite(c->streamIntegrals, sizeof(c->streamIntegrals[0]), es->numberStreams, f);
    }

    fwrite(checkpoint_tail, sizeof(checkpoint_tail), 1, f);
}

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
    if (BOINC_APPLICATION)
    {
        return (es->nu_step - es->lastCheckpointNuStep >= ia->nu_steps / 10) && mw_time_to_checkpoint();
    }
    else
    {
        return FALSE;
    }
}

int resolveCheckpoint(void)
{
    int rc;

    rc = mw_resolve_filename(CHECKPOINT_FILE, resolvedCheckpointPath, sizeof(resolvedCheckpointPath));
    if (rc)
        mw_printf("Error resolving checkpoint file '%s': %d\n", CHECKPOINT_FILE, rc);
    return rc;
}

int writeCheckpoint(EvaluationState* es)
{
    FILE* f;

    /* Avoid corrupting the checkpoint file by writing to a temporary file, and moving that */
    f = mw_fopen(CHECKPOINT_FILE_TMP, "wb");
    if (!f)
    {
        mwPerror("Opening checkpoint '%s'", CHECKPOINT_FILE_TMP);
        return 1;
    }

    es->lastCheckpointNuStep = es->nu_step;
    writeState(f, es);
    fclose(f);

    if (mw_rename(CHECKPOINT_FILE_TMP, resolvedCheckpointPath))
    {
        mwPerror("Failed to update checkpoint file ('%s' to '%s')",
                 CHECKPOINT_FILE_TMP,
                 resolvedCheckpointPath
            );
        return 1;
    }

    return 0;
}

int deleteCheckpoint(void)
{
    return mw_remove(resolvedCheckpointPath);
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
        {
            mw_report("Successfully resumed checkpoint\n");
        }
    }

    return 0;
}

