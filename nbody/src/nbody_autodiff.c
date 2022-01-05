/*
 *  Copyright (c) 2010-2018 Rensselaer Polytechnic Institute
 *  Copyright (c) 2018-2022 Eric Mendelsohn
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

#include "nbody_config.h"

#include "milkyway_util.h"

#if AUTODIFF

static void nbPrintLikelihood(FILE* f, real* likelihood)
{
    int i,j,k;

    mw_boinc_print(f, "<gradient>\n");
    fprintf(f, "[");
    for (i = 0; i < NumberOfModelParameters; i++)
    {
        fprintf(f, " %.15f",likelihood->gradient[i]);
    }
    fprintf(f, " ]\n");
    mw_boinc_print(f, "</gradient>\n");

    mw_boinc_print(f, "<hessian>\n");
    fprintf(f, "[");
    for (i = 0; i < NumberOfModelParameters; i++)
    {
        if (i != 0)
        {
            fprintf(f, " ");
        }
        fprintf(f, " [");
        for (j = 0; j < NumberOfModelParameters; j++)
        {
            if(i<j)
            {
                k=j;
                j=i;
                i=k;
            }
            k = (int) (i*(i+1)/2 + j);
            fprintf(f, " %.15f",likelihood->hessian[k]);
        }
        fprintf(f, " ]");
        if (i != NumberOfModelParameters)
        {
            fprintf(f, "\n");
        }
    }
    fprintf(f, " ]\n");
    mw_boinc_print(f, "</hessian>\n");
}

void printReal(real* a) //For debugging purposes
{
    int i;
    mw_printf("Value = %.15f\n",a->value);

    mw_printf("Gradient = [");
    for(i=0;i<NumberOfModelParameters;i++)
    {
        mw_printf(" %.15f", a->gradient[i]);
    }
    mw_printf(" ]\n");

    mw_printf("Hessian = [");
    for(i=0;i<HessianLength;i++)
    {
        mw_printf(" %.15f", a->hessian[i]);
    }
    mw_printf(" ]\n");
}

#else

static void nbPrintLikelihood(FILE* f, real* likelihood)
{
    mw_boinc_print(f, "<gradient>\n");
    fprintf(f, "[ NAN ]\n");
    mw_boinc_print(f, "</gradient>\n");

    mw_boinc_print(f, "<hessian>\n");
    fprintf(f, "[ NAN ]\n");
    mw_boinc_print(f, "</hessian>\n");
}

void printReal(real* a) //For debugging purposes
{
    mw_printf("Value = %.15f\n",*a);
}

#endif /*AUTODIFF*/


/* Write AUTODIFF file to given file name */
void nbWriteAutoDiff(const char* autoDiffFileName, real* likelihood)
{
    FILE* f = DEFAULT_OUTPUT_FILE;

    if (autoDiffFileName && strcmp(autoDiffFileName, ""))  /* If file specified, try to open it */
    {
        f = mwOpenResolved(autoDiffFileName, "w+");
        if (f == NULL) 
        {
            mwPerror("Error opening Derivative file '%s'. Using default output instead.", autoDiffFileName);
            f = DEFAULT_OUTPUT_FILE;
        }
    }

    nbPrintLikelihood(f, likelihood);

    if (f != DEFAULT_OUTPUT_FILE)
        fclose(f);
}
