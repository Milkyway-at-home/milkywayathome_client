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

#include "separation_types.h"
#include "evaluation_state.h"
#include "separation_graphics.h"
#include "milkyway_util.h"

#include <assert.h>

/* There's more broken in graphics2 header which makes this be C++ */
#include <boinc_api.h>
#include <graphics2.h>

/* Work around lack of user_data in callback function used by BOINC */

/* Reference to actual evaluation state */
static const EvaluationState* ges = NULL;

static EvaluationState* sharedEs = NULL;


/* Update shared copy of evaluation state */
static void separationGraphicsCallback()
{
    *sharedEs = *ges;
    sharedEs->integrals = NULL;
    /* copyEvaluationState(sharedEs, ges); */
}

static EvaluationState* createSharedEvaluationState()
{
    void* es;

    es = boinc_graphics_make_shmem(SEPARATION_SHMEM_NAME, sizeof(EvaluationState));

    /* FIXME: need 2nd shared memory segment for integral state, but
     * do we really need it? */

    return (EvaluationState*) es;
}

int separationInitSharedEvaluationState(EvaluationState* es)
{
    sharedEs = createSharedEvaluationState();
    if (!sharedEs)
        return 1;

    ges = es;    /* should be user data to callback */
    separationGraphicsCallback();
    boinc_register_timer_callback((FUNC_PTR) separationGraphicsCallback);

    return 0;
}

void separationFreeShared()
{
    /* FIXME: ????? Seems to be missing from BOINC API */
    ges = NULL; /* Remove reference */
}




