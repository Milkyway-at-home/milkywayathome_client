// A BOINC compatible fitness calculation routine for the Orphan FILE* density profile project
// B. Willett Feb. 25, 2010
// Adapted for use with B&H treecode May 11, 2010


// Takes a treecode position, converts it to (l,b), then to (lambda, beta), and then constructs a histogram of the density in lambda.

// Then calculates the (model density - data density)/sigma and adds them all up over all bins

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "nbody_priv.h"

#define r0 8.0
#define phi d2r(128.79)
#define theta d2r(54.39)
#define psi d2r(90.70)
#define beginning -50
#define end 50
#define binsize 3 

/* basic find max in array, O(n) */
static size_t findmax(real* arr, size_t n)
{
    size_t i, maxIdx = 0;

    for (i = 1; i < n; ++i)
    {
        if (arr[i] > arr[maxIdx])
            maxIdx = i;
    }

    return maxIdx;
}


real chisq(const NBodyCtx* ctx, NBodyState* st)
{

    printf("Initializing likelihood calculator...");
    real chisqval = 0.0;
    int i, j;
    int count = 0;
    int counttest = 1;
    int stillzero = 1;  /* For diagnosing zero chisq */

    const int nbody = ctx->model.nbody;
    const bodyptr endp = st->bodytab + nbody;

    FILE* f;

    real x[nbody+1];
    real y[nbody+1];
    real z[nbody+1];

    real l[nbody+1];
    real b[nbody+1];

    real lambda[nbody+1];
    real beta[nbody+1];

    bodyptr p;

    // Histogram prep
    int index1, index2;
    int largestbin1 = 0, largestbin2 = 0;
    int largestbin=0;
    int maxindex1 = (int)end / binsize;
    int maxindex2 = (int)abs(beginning) / binsize;

    real histodata1[maxindex1 + 1], histodata2[maxindex2 + 1];

    // Zero all of the histogram arrays
    memset(histodata1, 0, sizeof(real) * maxindex1);
    memset(histodata2, 0, sizeof(real) * maxindex2);

    printf("done\n");

    printf("Importing simulation results...");
    for (p = st->bodytab; p < endp; p++)
    {
        x[counttest-1] = Pos(p)[0];
        y[counttest-1] = Pos(p)[1];
        z[counttest-1] = Pos(p)[2];

        ++counttest;
    }
    printf("done\n");

    printf("Transforming simulation results...");
    for (i = 0; i < counttest - 1; i++)
    {
        count++;

        // Convert to (l,b) (involves convert x to Sun-centered)
        real r;

        x[count-1] += r0;

        r = rsqrt(x[count-1] * x[count-1] + y[count-1] * y[count-1] + z[count-1] * z[count-1]);

        // Leave in radians to make rotation easier
        b[count-1] = ratan2(z[count-1], rsqrt(x[count-1] * x[count-1] + y[count-1] * y[count-1]));
        l[count-1] = ratan2(y[count-1], x[count-1]);

        // Convert to (lambda, beta) (involves a rotation using the Newberg et al (2009) rotation matrices)

        beta[count-1] = r2d(rasin( rsin(theta) * rsin(phi) * rcos(b[count-1]) * rcos(l[count-1]) - rsin(theta) * rcos(phi) * rcos(b[count-1]) * rsin(l[count-1]) + rcos(theta) * rsin(b[count-1]) ));

        lambda[count-1] = r2d(atan2( (-rsin(psi) * rcos(phi) - rcos(theta) * rsin(phi) * rcos(psi)) * rcos(b[count-1]) * rcos(l[count-1]) + (-rsin(psi) * rsin(phi) + rcos(theta) * rcos(phi) * rcos(psi)) * rcos(b[count-1]) * rsin(l[count-1]) + rcos(psi) * rsin(theta) * rsin(b[count-1]), (rcos(psi) * rcos(phi) - rcos(theta) * rsin(phi) * rsin(psi)) * rcos(b[count-1]) * rcos(l[count-1]) + (rcos(psi) * rsin(phi) + rcos(theta) * rcos(phi) * rsin(psi)) * rcos(b[count-1]) * rsin(l[count-1]) + rsin(psi) * rsin(theta) * rsin(b[count-1]) ));

        // Create the histogram
        if (lambda[count-1] >= 0 && lambda[count-1] < end)
        {
            index1 = (int)(lambda[count-1] / binsize);
            histodata1[index1]++;

            if((int)histodata1[index1] > largestbin1) { largestbin1 = (int)histodata1[index1]; }
        }
        else if (lambda[count-1] > beginning && lambda[count-1] < 0)
        {
            index2 = abs((int)(lambda[count-1] / binsize));
            histodata2[abs(index2)]++;

            if((int)histodata2[abs(index2)] > largestbin2) { largestbin2 = (int)histodata2[abs(index2)]; }
        }
    }
    printf("done\n");

    // Find the largest entry
    // Get the single largest bin so we can normalize over it
    //CHECKME: Why the +1's in the sizes of histodatas?

    //largestbin1 = findmax(histodata1, maxindex1);
    //largestbin2 = findmax(histodata2, maxindex2);

    printf("Largestbin1: %i Largestbin2: %i\n", largestbin1, largestbin2);

    //largestbin = MAX(findmax(histodata1, maxindex1), findmax(histodata2, maxindex2));
    largestbin = MAX(largestbin1, largestbin2);

    printf("Largest bin: %i\n", largestbin);

    printf("Outputting to disk...");

    f = fopen("histout", "w");

    if (f == NULL)
    {
        printf("There was an error writing to the file\n");
        exit(EXIT_FAILURE);
    }

    printf("...file open...");

    // Print out the histogram
    real foo;
    for (i = 0, foo = -binsize; foo > beginning; foo -= binsize, ++i)
    {
        fprintf(f, "%f %f\n", foo + (binsize / 2.0) , histodata2[i] / ((real)largestbin));
    }

    for (i = 0, foo = 0; foo < end; foo += binsize, ++i)
    {
        fprintf(f, "%f %f\n", foo + (binsize / 2.0) , histodata1[i] / ((real)largestbin));
    }
    fclose(f);

    printf("done\n");

    // Calculate the chisq value by reading a file called "histogram" (the real data histogram should already be normalized)

    f = fopen("histogram", "r");

    if (f == NULL)
    {
        printf("histogram file not found...exiting\n");
        exit(EXIT_FAILURE);
    }

    int pos, fsize;
    real* fileLambda = NULL, *fileCount = NULL, *fileCountErr = NULL;
    int filecount = 0;

    pos = fseek(f, 0L, SEEK_END);
    fsize = ceil((double) (ftell(f) + 1) / 3);  /* Make sure it's big enough, avoid fun with integer division */
    fseek(f, 0L, SEEK_SET);

    fileLambda   = (real*) calloc(sizeof(real), fsize);
    fileCount    = (real*) calloc(sizeof(real), fsize);
    fileCountErr = (real*) calloc(sizeof(real), fsize);

    while (fscanf(f,
                  "%g %g %g\n",
                  &fileLambda[filecount],
                  &fileCount[filecount],
                  &fileCountErr[filecount]) != EOF)
    {
        ++filecount;
    }

    fclose(f);

    // Calculate the chisq
    for (i = 0, foo = -binsize; foo >= beginning; foo -= binsize, ++i)
    {
        for (j = 0; j < filecount; ++j)
        {
            // Don't include bins with zero, or zero error, this means there were no counts in that bin to begin with
            if (fileLambda[j] == foo + (binsize / 2.0) && histodata2[i] != 0 && fileCountErr[j] != 0)
            {
                chisqval += ((fileCount[j] - (histodata2[i] / largestbin)) / fileCountErr[j])
                              * ((fileCount[j] - (histodata2[i] / largestbin)) / fileCountErr[j]);

            }

            /* In principle, a zero chisq is the absolute best. But there's 2 ways that can happen: 1) all the bins perfectly lining up, or 2) the chisq having no contributions to it at all */
	    /* We need to check if the chisq is still zero after this addition AND the histogram data is also zero */
	    if(chisqval == 0.0 && histodata2[i] == 0.0) { stillzero = 1; } else { stillzero = 0; }
        }
    }

    for (i = 0, foo = 0.0; foo <= end; foo += binsize, ++i)
    {
        for (j = 0; j < filecount; ++j)
        {
            // Don't include bins with zero, or zero error, this means there were no counts in that bin to begin with
            if (fileLambda[j] == foo + (binsize / 2.0) && histodata1[i] != 0 && fileCountErr[j] != 0)
            {
                chisqval += ((fileCount[j] - (histodata1[i] / largestbin)) / fileCountErr[j])
                             * ((fileCount[j] - (histodata1[i] / largestbin)) / fileCountErr[j]);
            }

            /* In principle, a zero chisq is the absolute best. But there's 2 ways that can happen: 1) all the bins perfectly lining up, or 2) the chisq having no contributions to it at all */
	    /* We need to check if the chisq is still zero after this addition AND the histogram data is also zero */
	    if(chisqval == 0.0 && histodata1[i] == 0.0) { stillzero = 1; } else { stillzero = 0; }
        }
    }

    free(fileLambda);
    free(fileCount);
    free(fileCountErr);

    /* If stillzero = 1, then no contributions were ever added to the chisq, so there's actually no data in the range */
    /* Set it to a bad chisq */
    if(stillzero == 1) { chisqval = 998.00;}

    /* If the largest bin is zero, then this whole process breaks down and chisq blows up */
    /* Make a very bad chisq */ 
    if(largestbin == 0) { chisqval = 999.00;}

    printf("CHISQ = %f\n", chisqval);


    // MAXIMUM likelihood, multiply by -1
    real likelihood = -chisqval;
    printf("likelihood = %f\n", likelihood);
    return likelihood;
}

int compare (const void* a, const void* b)
{
    return ( *(int*)a - *(int*)b );
}

