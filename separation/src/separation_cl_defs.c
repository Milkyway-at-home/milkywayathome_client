/* Copyright 2010 Matthew Arsenault, Travis Desell, Boleslaw
Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail and
Rensselaer Polytechnic Institute.

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

#include <stdio.h>
#include "separation_cl_defs.h"

/* We can define these as compile time constants in the CL kernel,
 * giving us loop unrolling, avoiding branches etc. */

/* TODO: Stream constants? */
char* separationCLDefs(const ASTRONOMY_PARAMETERS* ap, const char* extra)
{
    char* buf;

    asprintf(&buf,
             " %s "
             "-D AP_ALPHA=%.20e "
             "-D AP_Q=%.20e "
             "-D AP_sn=%.20e "
             "-D AP_r0=%.20e "
             "-D AP_DELTA=%.20e "
             "-D AP_COEFF=%.20e "
             "-D AP_ALPHA_DELTA3=%.20e "
             "-D AP_BG_A=%.20e "
             "-D AP_BG_B=%.20e "
             "-D AP_BG_C=%.20e "
             "-D AP_WEDGE=%d "
             "-D AP_AUX_BG_PROFILE=%d "
             "-D AP_CONVOLVE=%u "
             "-D AP_NUMBER_STREAMS=%u "
             "-D AP_SGR_COORDINATES=%d "
             "-D AP_TOTAL_CALC_PROBS=%.20e "
             "-D AP_BACKGROUND_WEIGHT=%.20e ",
             extra,

             ap->alpha,
             ap->q,
             ap->sn,
             ap->r0,
             ap->delta,
             ap->coeff,
             ap->alpha_delta3,
             ap->bg_a,
             ap->bg_b,
             ap->bg_c,
             ap->wedge,
             ap->aux_bg_profile,
             ap->convolve,
             ap->number_streams,
             ap->sgr_coordinates,
             ap->total_calc_probs,
             ap->background_weight
        );

    return buf;
}

