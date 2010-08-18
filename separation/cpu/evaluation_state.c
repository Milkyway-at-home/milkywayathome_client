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

#include "separation_types.h"
#include "evaluation_state.h"
#include "milkyway_util.h"

void initialize_integral(INTEGRAL* integral, unsigned int number_streams)
{
    integral->background_integral = 0.0;
    integral->stream_integrals = callocSafe(number_streams, sizeof(double));
    integral->probs = callocSafe(number_streams, sizeof(ST_PROBS));
}

void initialize_state(const ASTRONOMY_PARAMETERS* ap, EVALUATION_STATE* es)
{
    unsigned int i;

    es->current_integral = 0;
    es->number_streams = ap->number_streams;

    es->number_integrals = ap->number_integrals;
    es->integrals = mallocSafe(sizeof(INTEGRAL) * ap->number_integrals);

    for (i = 0; i < ap->number_integrals; i++)
        initialize_integral(&es->integrals[i], ap->number_streams);
}

void free_integral(INTEGRAL* i)
{
    free(i->stream_integrals);
    free(i->probs);
}

void free_evaluation_state(EVALUATION_STATE* es)
{
    unsigned int i;

    for (i = 0; i < es->number_integrals; ++i)
        free_integral(&es->integrals[i]);
    free(es->integrals);
}

#if BOINC_APPLICATION

void print_evaluation_state(const EVALUATION_STATE* es)
{
    INTEGRAL* i;
    unsigned int j;

    printf("evaluation-state {\n"
           "  r_step           = %u\n"
           "  nu_step          = %u\n"
           "  current_integral = %u\n"
           "  r_acc            = { %.20g, %.20g }\n"
           "  nu_acc           = { %.20g, %.20g }\n",
           es->r_step,
           es->nu_step,
           es->current_integral,
           es->r_acc.bg_int,
           es->r_acc.correction,
           es->nu_acc.bg_int,
           es->nu_acc.correction);

    for (i = es->integrals; i < es->integrals + es->number_integrals; ++i)
    {
        printf("integral: background_integral = %g\n", i->background_integral);
        printf("Stream integrals = ");
        for (j = 0; j < es->number_streams; ++j)
            printf("  %g, ", i->stream_integrals[j]);
        printf("Probs = ");
        for (j = 0; j < es->number_streams; ++j)
            printf(" { %g, %g },", i->probs[j].st_prob_int, i->probs[j].st_prob_int_c);
        printf("\n");
    }
    printf("\n");
}

static const char checkpoint_header[] = "separation_checkpoint";
static const char checkpoint_tail[] = "end_checkpoint";

int read_checkpoint(EVALUATION_STATE* es)
{
    INTEGRAL* i;
    char input_path[512];
    char str_buf[sizeof(checkpoint_tail)];

    if (boinc_resolve_filename(CHECKPOINT_FILE, input_path, sizeof(input_path)))
        return 1;

    FILE* f = boinc_fopen(input_path, "r");
    if (!f)
        return 1;

    fread(str_buf, sizeof(checkpoint_header), 1, f);
    if (strcmp(str_buf, checkpoint_header))
    {
        fprintf(stderr, "Failed to find header in checkpoint file\n");
        return 1;
    }

    fread(&es->current_integral, sizeof(es->current_integral), 1, f);
    fread(&es->r_step, sizeof(es->r_step), 1, f);
    fread(&es->nu_step, sizeof(es->nu_step), 1, f);
    fread(&es->r_acc, sizeof(es->r_acc), 1, f);
    fread(&es->nu_acc, sizeof(es->nu_acc), 1, f);

    for (i = es->integrals; i < es->integrals + es->number_integrals; ++i)
    {
        fread(&i->background_integral, sizeof(i->background_integral), 1, f);
        fread(i->stream_integrals, sizeof(double), es->number_streams, f);
        fread(i->probs, sizeof(ST_PROBS), es->number_streams, f);
    }

    fread(str_buf, sizeof(checkpoint_tail), 1, f);
    if (strcmp(str_buf, checkpoint_tail))
    {
        fprintf(stderr, "Failed to find tail in checkpoint file\n");
        return 1;
    }

    fclose(f);

    return 0;
}

int write_checkpoint(const EVALUATION_STATE* es)
{
    INTEGRAL* i;
    char output_path[512];
    FILE* f;
    const INTEGRAL* endi = es->integrals + es->number_integrals;

    /* Avoid corrupting the checkpoint file by writing to a temporary file, and moving that */
    boinc_resolve_filename(CHECKPOINT_FILE_TMP, output_path, sizeof(output_path));

    f = boinc_fopen(output_path, "wb");
    if (!f)
    {
        fprintf(stderr, "APP: error writing checkpoint (opening checkpoint file)\n");
        return 1;
    }

    fwrite(checkpoint_header, sizeof(checkpoint_header), 1, f);

    fwrite(&es->current_integral, sizeof(es->current_integral), 1, f);
    fwrite(&es->r_step, sizeof(es->r_step), 1, f);
    fwrite(&es->nu_step, sizeof(es->nu_step), 1, f);
    fwrite(&es->r_acc, sizeof(es->r_acc), 1, f);
    fwrite(&es->nu_acc, sizeof(es->nu_acc), 1, f);

    for (i = es->integrals; i < endi; ++i)
    {
        fwrite(&i->background_integral, sizeof(i->background_integral), 1, f);
        fwrite(i->stream_integrals, sizeof(double), es->number_streams, f);
        fwrite(i->probs, sizeof(ST_PROBS), es->number_streams, f);
    }

    fwrite(checkpoint_tail, sizeof(checkpoint_tail), 1, f);

    fclose(f);

    if (boinc_rename(CHECKPOINT_FILE_TMP, CHECKPOINT_FILE))
    {
        fprintf(stderr, "Failed to update checkpoint file\n");
        mw_finish(EXIT_FAILURE);
    }


    return 0;
}

#endif /* BOINC_APPLICATION */


