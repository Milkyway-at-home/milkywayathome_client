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

static void nbPrintLikelihood(FILE* f, real* likelihood) //NOTE: THIS IS THE NEGATIVE LIKELIHOOD! MUST FLIP SIGNS HERE!
{
    int i,j,eff_i, eff_j,k;

    mw_boinc_print(f, "<gradient>\n");
    fprintf(f, "[");
    for (i = 0; i < NumberOfModelParameters; i++)
    {
        fprintf(f, " %.15e",-likelihood->gradient[i]);
        if(i!=NumberOfModelParameters-1) fprintf(f, ",");
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
                eff_i = j;
                eff_j = i;
            }
            else
            {
                eff_i = i;
                eff_j = j;
            }
            k = (int) (eff_i*(eff_i+1)/2 + eff_j);
            fprintf(f, " %.15e",-likelihood->hessian[k]);
            if(j!=NumberOfModelParameters-1) fprintf(f, ",");
        }
        fprintf(f, " ]");
        if (i != NumberOfModelParameters-1)
        {
            fprintf(f, "\n");
        }
    }
    fprintf(f, " ]\n");
    mw_boinc_print(f, "</hessian>\n");
}

void printRealFull(real* a, char var_name[]) //For debugging purposes
{
    int i;
    mw_printf("%s = %.15f\n", var_name, a->value);

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

void printRealGradient(real* a, char var_name[]) //For debugging purposes
{
    int i;
    mw_printf("NABLA %s = [", var_name);
    for(i=0;i<NumberOfModelParameters;i++)
    {
        mw_printf(" %.15f", a->gradient[i]);
    }
    mw_printf(" ]\n");
}

#else

static void nbPrintLikelihood(FILE* f, real* likelihood)
{
    mw_boinc_print(f, "<gradient>\n");
    fprintf(f, "[ No Derivative Information Calculated ]\n");
    mw_boinc_print(f, "</gradient>\n");

    mw_boinc_print(f, "<hessian>\n");
    fprintf(f, "[ No Derivative Information Calculated ]\n");
    mw_boinc_print(f, "</hessian>\n");
}

void printRealFull(real* a, char var_name[]) //For debugging purposes
{
    mw_printf("%s = %.15f\n", var_name, *a);
}

void printRealGradient(real* a, char var_name[]) //For debugging purposes
{
    mw_printf("%s = %.15f\n", var_name, *a);
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
