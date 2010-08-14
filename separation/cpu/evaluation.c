/*
Copyright 2008, 2009 Travis Desell, Dave Przybylo, Nathan Cole,
Boleslaw Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail
and Rensselaer Polytechnic Institute.

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

#include "separation.h"
#include "evaluation.h"
#include "evaluation_state.h"
#include "integral_constants.h"
#include "integrals_likelihood.c"

static void print_stream_integrals(EVALUATION_STATE* es, const unsigned int number_streams)
{
    unsigned int i;
    fprintf(stderr, "<background_integral> %.20lf </background_integral>\n", es->background_integral);
    fprintf(stderr, "<stream_integrals>");
    for (i = 0; i < number_streams; i++)
        fprintf(stderr, " %.20lf", es->stream_integrals[i]);
    fprintf(stderr, " </stream_integrals>\n");
}

static void final_stream_integrals(EVALUATION_STATE* es,
                                   const unsigned int number_streams,
                                   const unsigned int number_integrals)
{
    unsigned int i, j;

    es->background_integral = es->integrals[0].background_integral;
    for (i = 0; i < number_streams; ++i)
        es->stream_integrals[i] = es->integrals[0].stream_integrals[i];

    for (i = 1; i < number_integrals; ++i)
    {
        es->background_integral -= es->integrals[i].background_integral;
        for (j = 0; j < number_streams; j++)
            es->stream_integrals[j] -= es->integrals[i].stream_integrals[j];
    }

}

double evaluate(const ASTRONOMY_PARAMETERS* ap,
                const STAR_POINTS* sp,
                const STREAMS* streams,
                const STREAM_CONSTANTS* sc)
{
    double likelihood_val;
    EVALUATION_STATE es = EMPTY_EVALUATION_STATE;
    STREAM_GAUSS sg;

    initialize_state(ap, &es);
    get_stream_gauss(ap->convolve, &sg);

#if BOINC_APPLICATION
    if (boinc_file_exists(CHECKPOINT_FILE))
    {
        fprintf(stderr, "Checkpoint exists. Attempting to resume from it\n");

        if (read_checkpoint(&es))
        {
            fprintf(stderr, "Reading checkpoint failed\n");
            boinc_delete_file(CHECKPOINT_FILE);
            mw_finish(EXIT_FAILURE);
        }
    }
#endif

    vector* xyz = malloc(sizeof(vector) * ap->convolve);

    calculate_integrals(ap, sc, &sg, &es, xyz);

    /* Final checkpoint. */
    write_checkpoint(&es);

    final_stream_integrals(&es, ap->number_streams, ap->number_integrals);
    print_stream_integrals(&es, ap->number_streams);

    likelihood_val = likelihood(ap, sc, streams, &es, &sg, xyz, sp);

    free(xyz);
    free_evaluation_state(&es);
    free_stream_gauss(&sg);

  #if BOINC_APPLICATION
    boinc_delete_file(CHECKPOINT_FILE);
  #endif

    return likelihood_val;
}


