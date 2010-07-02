// A BOINC compatible fitness calculation routine for the Orphan FILE* density profile project
// B. Willett Feb. 25, 2010
// Adapted for use with B&H treecode May 11, 2010


// Takes a treecode position, converts it to (l,b), then to (lambda, beta), and then constructs a histogram of the density in lambda.

// Then calculates the (model density - data density)/sigma and adds them all up over all bins

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "nbody_priv.h"

/* FIXME: Magic numbers */
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
    real chisqval = 0.0;
    int i, j;
    int stillzero = TRUE;  /* For diagnosing zero chisq */

    const int nbody = ctx->model.nbody;
    const bodyptr endp = st->bodytab + nbody;
    vector lbr;
    real lambda;
    bodyptr p;
    FILE* f;

    // Histogram prep
    int index1, index2;
    int largestbin1 = 0, largestbin2 = 0;
    int largestbin = 0;
    int maxindex1 = (int)end / binsize;
    int maxindex2 = (int)abs(beginning) / binsize;

    real* histodata1 = calloc(sizeof(real), maxindex1 + 1);
    real* histodata2 = calloc(sizeof(real), maxindex2 + 1);

    printf("Transforming simulation results...");
    for (p = st->bodytab; p < endp; ++p)
    {
        // Convert to (l,b) (involves convert x to Sun-centered)
        // Leave in radians to make rotation easier
        cartesianToLbr_rad(ctx, lbr, Pos(p));

        // Convert to (lambda, beta) (involves a rotation using the Newberg et al (2009) rotation matrices)

        //CHECKME: beta never actually used?
        //beta = r2d(rasin( rsin(theta) * rsin(phi) * rcos(B(lbr)) * rcos(L(lbr)) - rsin(theta) * rcos(phi) * rcos(B(lbr)) * rsin(L(lbr)) + rcos(theta) * rsin(B(lbr)) ));

        lambda = r2d(atan2( (-rsin(psi) * rcos(phi) - rcos(theta) * rsin(phi) * rcos(psi)) * rcos(B(lbr)) * rcos(L(lbr)) + (-rsin(psi) * rsin(phi) + rcos(theta) * rcos(phi) * rcos(psi)) * rcos(B(lbr)) * rsin(L(lbr)) + rcos(psi) * rsin(theta) * rsin(B(lbr)), (rcos(psi) * rcos(phi) - rcos(theta) * rsin(phi) * rsin(psi)) * rcos(B(lbr)) * rcos(L(lbr)) + (rcos(psi) * rsin(phi) + rcos(theta) * rcos(phi) * rsin(psi)) * rcos(B(lbr)) * rsin(L(lbr)) + rsin(psi) * rsin(theta) * rsin(B(lbr)) ));

        // Create the histogram
        if (lambda >= 0 && lambda < end)
        {
            index1 = (int)(lambda / binsize);
            histodata1[index1]++;

            if ((int)histodata1[index1] > largestbin1)
                largestbin1 = (int)histodata1[index1];
        }
        else if (lambda > beginning && lambda < 0)
        {
            index2 = abs((int)(lambda / binsize));
            histodata2[abs(index2)]++;

            if ((int)histodata2[abs(index2)] > largestbin2)
                largestbin2 = (int)histodata2[abs(index2)];
        }
        /* CHECKME: else? */
    }
    printf("done\n");

    // Find the largest entry
    // Get the single largest bin so we can normalize over it
    //CHECKME: Why the +1's in the sizes of histodatas?

    printf("Largestbin1: %i Largestbin2: %i\n", largestbin1, largestbin2);

    largestbin = MAX(largestbin1, largestbin2);

    printf("Largest bin: %i\n", largestbin);

    f = fopen("histout", "w");
    if (f == NULL)
    {
        perror("Writing histout");
        nbody_finish(EXIT_FAILURE);
    }

    printf("...file open...");

    // Print out the histogram
    real foo;
    for (i = 0, foo = -binsize; foo >= beginning; foo -= binsize, ++i)
    {
        fprintf(f, "%f %f\n", foo + (binsize / 2.0) , histodata2[i] / ((real)largestbin));
    }

    for (i = 0, foo = 0; foo <= end; foo += binsize, ++i)
    {
        fprintf(f, "%f %f\n", foo + (binsize / 2.0) , histodata1[i] / ((real)largestbin));
    }
    fclose(f);

    printf("done\n");

    // Calculate the chisq value by reading a file called "histogram" (the real data histogram should already be normalized)

    f = fopen("histogram", "r");
    if (f == NULL)
    {
        perror("histogram");
        nbody_finish(EXIT_FAILURE);
    }

    int fpos, fsize;
    real* fileLambda = NULL, *fileCount = NULL, *fileCountErr = NULL;
    int filecount = 0;

    fpos = fseek(f, 0L, SEEK_END);
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

                /* In principle, a zero chisq is the absolute best. But there's 2 ways that can happen: 1) all the bins perfectly lining up, or 2) the chisq having no contributions to it at all */
		/* We need to check if the chisq is still zero after this addition */
                if (chisqval != 0.0)
                    stillzero = FALSE;

            }
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
                /* In principle, a zero chisq is the absolute best. But there's 2 ways that can happen: 1) all the bins perfectly lining up, or 2) the chisq having no contributions to it at all */
		/* We need to check if the chisq is still zero after this addition */
		if(chisqval != 0.0) { stillzero = 0; };
            }
        }
    }

    free(fileLambda);
    free(fileCount);
    free(fileCountErr);

    free(histodata1);
    free(histodata2);

    /* If stillzero, then no contributions were ever added to the
     * chisq, so there's actually no data in the range */
    /* Set it to a bad chisq */
    if (stillzero)
        chisqval = 9999.99;

    printf("CHISQ = %f\n", chisqval);

    // MAXIMUM likelihood, multiply by -1
    real likelihood = -chisqval;
    printf("likelihood = %f\n", likelihood);
    return likelihood;
}

