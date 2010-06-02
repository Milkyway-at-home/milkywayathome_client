// A BOINC compatible fitness calculation routine for the Orphan FILE* density profile project
// B. Willett Feb. 25, 2010
// Adapted for use with B&H treecode May 11, 2010

// Takes a treecode position, converts it to (l,b), then to (lambda, beta), and then constructs a histogram of the density in lambda.

// Then calculates the (model density - data density)/sigma and adds them all up over all bins

#include <stdio.h>
#include <stdlib.h>
#include "code.h"

#define r0 8.0
#define phi 128.79*(3.141592654/180.0)
#define theta 54.39*(3.141592654/180.0)
#define psi 90.70*(3.141592654/180.0)
#define beginning -50
#define end 50
#define binsize 3

float chisq()
{

    printf("Initializing likelihood calculator...");
    float chisqval = 0.0;
    int i, j;
    int count = 0;
    int counttest = 1;

    FILE* f;

    float x[nbody];
    float y[nbody];
    float z[nbody];

    float l[nbody];
    float b[nbody];

    float lambda[nbody];
    float beta[nbody];

    bodyptr p;

    // Histogram prep
    int index1, index2;
    int maxindex1 = (int)end / binsize;
    int maxindex2 = (int)abs(beginning) / binsize;

    float histodata1[maxindex1 + 1], histodata2[maxindex2 + 1], histodata1sort[maxindex1 + 1], histodata2sort[maxindex2 + 1];

    // Zero all of the histogram arrays
    for (i = 0; i <= maxindex1; i++)
    {
        histodata1[i] = 0.0;
        histodata1sort[i] = 0.0;
    }
    for (i = 0; i <= maxindex2; i++)
    {
        histodata2[i] = 0.0;
        histodata2sort[i] = 0.0;
    }


    printf("done\n");

    printf("Importing simulation results...");
    for (p = bodytab; p < bodytab + nbody; p++)
    {

        x[counttest-1] = Pos(p)[0];
        y[counttest-1] = Pos(p)[1];
        z[counttest-1] = Pos(p)[2];

        counttest++;
    }
    printf("done\n");

    printf("Transforming simulation results...");
    for (i = 0; i <= counttest - 2; i++)
    {
        count++;

        // Convert to (l,b) (involves convert x to Sun-centered)
        float r;

        x[count-1] += r0;

        r = pow(x[count-1] * x[count-1] + y[count-1] * y[count-1] + z[count-1] * z[count-1], 0.5);

        // Leave in radians to make rotation easier
        b[count-1] = atan2(z[count-1], pow(x[count-1] * x[count-1] + y[count-1] * y[count-1], 0.5));
        l[count-1] = atan2(y[count-1], x[count-1]);

        // Convert to (lambda, beta) (involves a rotation using the Newberg et al (2009) rotation matricies)

        beta[count-1] = (180.0 / 3.141592654) * asin( sin(theta) * sin(phi) * cos(b[count-1]) * cos(l[count-1]) - sin(theta) * cos(phi) * cos(b[count-1]) * sin(l[count-1]) + cos(theta) * sin(b[count-1]) );

        lambda[count-1] = (180.0 / 3.141592654) * atan2( (-sin(psi) * cos(phi) - cos(theta) * sin(phi) * cos(psi)) * cos(b[count-1]) * cos(l[count-1]) + (-sin(psi) * sin(phi) + cos(theta) * cos(phi) * cos(psi)) * cos(b[count-1]) * sin(l[count-1]) + cos(psi) * sin(theta) * sin(b[count-1]), (cos(psi) * cos(phi) - cos(theta) * sin(phi) * sin(psi)) * cos(b[count-1]) * cos(l[count-1]) + (cos(psi) * sin(phi) + cos(theta) * cos(phi) * sin(psi)) * cos(b[count-1]) * sin(l[count-1]) + sin(psi) * sin(theta) * sin(b[count-1]) );

        // Create the histogram
        if (lambda[count-1] > 0 && lambda[count-1] < end)
        {
            index1 = (int)(lambda[count-1] / binsize);
            histodata1[index1]++;
            histodata1sort[index1]++;
        }
        else if (lambda[count-1] > beginning)
        {
            index2 = abs((int)(lambda[count-1] / binsize));
            histodata2[abs(index2)]++;
            histodata2sort[abs(index2)]++;
        }
    }
    printf("done\n");

    // Sort the two histogram arrays to find the largest entry (this is why I used 2 arrays, because qsort overwrites it)

    qsort(histodata1sort, maxindex1, sizeof(int), compare);
    qsort(histodata2sort, maxindex2, sizeof(int), compare);

    // Get the single largest bin so we can normalize over it
    int largestbin;
    int largestbin1 = histodata1sort[maxindex1+1];
    int largestbin2 = histodata2sort[maxindex2+1];

    if (largestbin1 >= largestbin2)
    {
        largestbin = largestbin1;
    }
    else
    {
        largestbin = largestbin2;
    }

    printf("Largest bin: %i\n", largestbin);

    printf("Outputting to disk...");

    f = fopen("hist", "w");

    if (f == NULL)
    {
        printf("There was an error writing to the file\n");
        exit(0);
    }

    printf("...file open...");

    // Print out the histogram
    float foo;
    i = 0;
    for (foo = -binsize; foo >= beginning; foo -= binsize)
    {
        fprintf(f, "%f %f\n", foo + (binsize / 2.0) , histodata2[i] / ((float)largestbin));
        i++;
    }
    i = 0;
    for (foo = 0; foo <= end; foo += binsize)
    {
        fprintf(f, "%f %f\n", foo + (binsize / 2.0) , histodata1[i] / ((float)largestbin));
        i++;
    }
    fclose(f);

    printf("done\n");

    // Calculate the chisq value by reading a file called "histogram" (the real data histogram should already be normalized)

    f = fopen("histogram", "r");

    if (f == NULL)
    {
        printf("histogram file not found...exiting\n");
        exit(0);
    }

    float* fileLambda = NULL, *fileCount = NULL, *fileCountErr = NULL;
    int filecount = 0;

    do
    {
        filecount++;
        fileLambda = (float*) realloc (fileLambda, filecount * sizeof(float));
        fileCount = (float*) realloc (fileCount, filecount * sizeof(float));
        fileCountErr = (float*) realloc (fileCountErr, filecount * sizeof(float));
    }
    while (fscanf(f, "%f %f %f", &fileLambda[filecount-1], &fileCount[filecount-1], &fileCountErr[filecount-1]) != EOF);

    fclose(f);

    // Calculate the chisq
    i = 0;
    for (foo = -binsize; foo >= beginning; foo -= binsize)
    {
        for (j = 0; j <= filecount - 1; j++)
        {
            // Don't include bins with zero, or zero error, this means there were no counts in that bin to begin with
            if (fileLambda[j] == foo + (binsize / 2.0) && histodata2[i] != 0 && fileCountErr[j] != 0)
            {
                chisqval += ((fileCount[j] - (histodata2[i] / largestbin)) / fileCountErr[j]) * ((fileCount[j] - (histodata2[i] / largestbin)) / fileCountErr[j]);
            }
        }
        i++;
    }
    i = 0;
    for (foo = 0; foo <= end; foo += binsize)
    {
        for (j = 0; j <= filecount - 1; j++)
        {
            // Don't include bins with zero, or zero error, this means there were no counts in that bin to begin with
            if (fileLambda[j] == foo + (binsize / 2.0) && histodata1[i] != 0 && fileCountErr[j] != 0)
            {
                chisqval += ((fileCount[j] - (histodata1[i] / largestbin)) / fileCountErr[j]) * ((fileCount[j] - (histodata1[i] / largestbin)) / fileCountErr[j]);
            }
        }
        i++;
    }

    free(fileLambda);
    free(fileCount);
    free(fileCountErr);

    printf("CHISQ = %f\n", chisqval);


    // MAXIMUM likelihood, multiply by -1
    float likelihood = -1.0 * chisqval;
    printf("likelihood = %f\n", likelihood);
    return likelihood;
}

int compare (const void* a, const void* b)
{
    return ( *(int*)a - *(int*)b );
}

