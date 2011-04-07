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

    for (c = es->cuts; c < es->cuts + es->numberCuts; ++c)
    {
        printf("integral: bgIntegral = %g\n", c->bgIntegral);
        printf("Stream integrals = ");
        for (j = 0; j < es->numberCuts; ++j)
            printf("  %g, ", c->streamIntegrals[j]);
    }
    printf("\n");
}

static const char checkpoint_header[] = "separation_checkpoint";
static const char checkpoint_tail[] = "end_checkpoint";


static int readState(FILE* f, EvaluationState* es)
{
    Cut* c;
    char str_buf[sizeof(checkpoint_header) + 1];

    fread(str_buf, sizeof(checkpoint_header), 1, f);
    if (strncmp(str_buf, checkpoint_header, sizeof(str_buf)))
    {
        warn("Failed to find header in checkpoint file\n");
        return 1;
    }

    fread(&es->currentCut, sizeof(es->currentCut), 1, f);
    fread(&es->nu_step, sizeof(es->nu_step), 1, f);
    fread(&es->mu_step, sizeof(es->mu_step), 1, f);
    fread(&es->bgSum, sizeof(es->bgSum), 1, f);
    fread(es->streamSums, sizeof(es->streamSums[0]), es->numberStreams, f);

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

    return rc;
}

static inline void writeState(FILE* f, const EvaluationState* es)
{
    Cut* c;
    const Cut* endc = es->cuts + es->numberCuts;

    fwrite(checkpoint_header, sizeof(checkpoint_header), 1, f);

    fwrite(&es->currentCut, sizeof(es->currentCut), 1, f);
    fwrite(&es->nu_step, sizeof(es->nu_step), 1, f);
    fwrite(&es->mu_step, sizeof(es->mu_step), 1, f);
    fwrite(&es->bgSum, sizeof(es->bgSum), 1, f);
    fwrite(es->streamSums, sizeof(es->streamSums[0]), es->numberStreams, f);

    for (c = es->cuts; c < endc; ++c)
    {
        fwrite(&c->bgIntegral, sizeof(c->bgIntegral), 1, f);
        fwrite(c->streamIntegrals, sizeof(c->streamIntegrals[0]), es->numberStreams, f);
    }

    fwrite(checkpoint_tail, sizeof(checkpoint_tail), 1, f);
}



#if !SEPARATION_OPENCL

int resolveCheckpoint()
{
    int rc;

    rc = mw_resolve_filename(CHECKPOINT_FILE, resolvedCheckpointPath, sizeof(resolvedCheckpointPath));
    if (rc)
        warn("Error resolving checkpoint file '%s': %d\n", CHECKPOINT_FILE, rc);
    return rc;
}

int writeCheckpoint(const EvaluationState* es)
{
    FILE* f;

    /* Avoid corrupting the checkpoint file by writing to a temporary file, and moving that */
    f = mw_fopen(CHECKPOINT_FILE_TMP, "wb");
    if (!f)
    {
        perror("Opening checkpoint temp");
        return 1;
    }

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

#else /* SEPARATION_OPENCL */

int resolveCheckpoint()
{
    return 0;
}

int writeCheckpoint(const EvaluationState* es)
{
    return 0;
}

int maybeResume(EvaluationState* es)
{
    return 0;
}

#endif /* !SEPARATION_OPENCL */

