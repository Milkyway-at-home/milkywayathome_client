/*
 * Copyright (c) 2012 Rensselaer Polytechnic Institute
 * Copyright (c) 2016-2018 Siddhartha Shelton
 * 
 * This file is part of Milkway@Home.
 *
 * Milkyway@Home is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Milkyway@Home is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Milkyway@Home.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "nbody_devoptions.h"
#include "milkyway_math.h"
#include "nbody_types.h"
#include "nbody.h"
#include "nbody_io.h"

int dev_write_outputs(const NBodyCtx* ctx, const NBodyState* st, const NBodyFlags* nbf, real freq)
{
    int rc = 0;
    freq = (freq == 0.0)*1.0 + freq;
    if((st->step + 1) % (int) freq == 0)
    {
        
        FILE* f;
        char output_file_name[1024];
        sprintf(output_file_name, "%d", st->step);
        
        f = mwOpenResolved(output_file_name, "w+");
        mw_boinc_print(f, "<bodies>\n");
        rc = nbOutputBodies(f, ctx, st, nbf);
        mw_boinc_print(f, "</bodies>\n");
        fclose(f);
        
    }
    
    return rc;
}

int dev_write_initial_output(const NBodyCtx* ctx, const NBodyState* st, const NBodyFlags* nbf)
{

    int rc = 0;
    FILE* f;
    f = mwOpenResolved("initial.out", "w+");
    if (!f)
    {
        mw_printf("Error opening 'initial.out' for writing\n");
        return 1;
    }

    mw_boinc_print(f, "<bodies>\n");
    rc = nbOutputBodies(f, ctx, st, nbf);
    mw_boinc_print(f, "</bodies>\n");

    fclose(f);
    
    mw_printf("Initial body data written to 'initial.out'\n");
    return rc;
}



