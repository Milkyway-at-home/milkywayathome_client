/*This program tests code for the momentum calculation from a simulation state, 
the momentum likelihood calculation, and the reading of momentum data*/

#include "nbody_mass.h"
#include "nbody_histogram.h"
#include "nbody_types.h"
#include "nbody_checkpoint.h"
#include "nbody.h"
#include <stdio.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

//defined in nbody_checkpoint.c and not the header, so need to define it again here
typedef struct
{
    char header[128];                     /* "mwnbody" */
    uint32_t majorVersion, minorVersion;  /* Version check */
    uint32_t nbody;
    uint32_t step;
    uint32_t realSize;                   /* Does the checkpoint use float or double */
    uint32_t ptrSize;
    uint32_t nOrbitTrace;
    uint32_t nShiftLMC;
    uint32_t treeIncest;
    real rsize;
    NBodyCtx ctx;
} NBodyCheckpointHeader;

//Function to read in bestlikelihood body state from checkpoint, so the same checkpoint can be used even as versions change
int read_bestbodytab_from_checkpoint(const char* checkpointFile, Body** bodytab_out, int* nbody_out) {
    int fd = open(checkpointFile, O_RDONLY);
    if (fd == -1) {
        perror("Failed to open checkpoint file");
        return 1;
    }

    struct stat sb;
    if (fstat(fd, &sb) == -1) {
        perror("Failed to stat checkpoint file");
        close(fd);
        return 1;
    }

    char* mptr = mmap(NULL, sb.st_size, PROT_READ, MAP_SHARED, fd, 0);
    if (mptr == MAP_FAILED) {
        perror("Failed to mmap checkpoint file");
        close(fd);
        return 1;
    }

    // Read header
    NBodyCheckpointHeader cpHdr;
    memcpy(&cpHdr, mptr, 800); //size of header hardcoded in as it will change with updates
    if (strncmp(cpHdr.header, "mwnbody", 7) != 0) {
        fprintf(stderr, "Invalid checkpoint header\n");
        munmap(mptr, sb.st_size);
        close(fd);
        return 1;
    }

    uint32_t nbody = cpHdr.nbody;

    // Find offset to bestLikelihoodBodyTab
    char* p = mptr + 800; //size of header hardcoded in as it will change with updates
    size_t sizeOfData;
    memcpy(&sizeOfData, p, sizeof(size_t));
    p += sizeof(size_t);
    p += sizeOfData; // skip extractedSt

    p += nbody * sizeof(Body); // skip normal bodytab

    // Now p points to best likelihood bodytab
    Body* bestbodytab = (Body*)malloc(nbody * sizeof(Body));
    if (!bestbodytab) {
        fprintf(stderr, "Failed to allocate bodytab\n");
        munmap(mptr, sb.st_size);
        close(fd);
        return 1;
    }
    memcpy(bestbodytab, p, nbody * sizeof(Body));

    *bodytab_out = bestbodytab;
    *nbody_out = nbody;

    munmap(mptr, sb.st_size);
    close(fd);
    return 0;
}

int main()
{
    MainStruct* data; // Input histogram
    MainStruct* histogram; // Histogram calculated from simulation state
    data = nbReadHistogram("./average_bins_test.hist"); // input momentum data for the test is included in this histogram
    histogram = nbReadHistogram("./average_bins_test.hist"); // reading one in so its initialized

    if(histogram == NULL || data == NULL)
    {
        printf("\tFailed to read in histogram for momentum test\n");
        return 1;
    }

    NBodyFlags nbf;
    NBodyCtx ctx;
    NBodyState st;

    /* A known nbody state will be read in from a checkpoint file. Momentum values calculated 
    will be checked against those calculated from the corresponding output file in python*/
    nbf.checkpointFileName = "./momentum_test_checkpoint";
    if(read_bestbodytab_from_checkpoint(nbf.checkpointFileName, &st.bodytab, &st.nbody))
    {
        printf("\tFailed to read best likelihood bodystate from checkpoint file\n");
        return 1;
    }

    //initialize needed context values
    ctx.sunGCDist = 8.0;
    ctx.MomentumSigma = 2.5;
    ctx.IterMax = 6;
    ctx.MomentumCorrect = 1.111;

    /* Set desired parameters for test */
    data->histograms[0]->params.nRange = 0;
    histogram->histograms[0]->params.lambdaStart = -33;
    histogram->histograms[0]->params.lambdaEnd = 48;
    histogram->histograms[0]->params.betaStart = -15;
    histogram->histograms[0]->params.betaEnd = 15;
    histogram->histograms[0]->params.phi=128.790000;
    histogram->histograms[0]->params.theta=54.390000;
    histogram->histograms[0]->params.psi=90.700000;

    nbCalcMomentum(&st, &ctx, data->histograms[0], histogram->histograms[0]);

    mwvector L = histogram->histograms[0]->params.L;
    mwvector LErr = histogram->histograms[0]->params.LErr;

    // Values calculated from the output file in python are:
    // L = [1527.6464608614306, 1597.9471123610872, 3393.8635866082823]
    // LErr = [125.73297168744209, 173.28623918527114, 50.229315701648076]

    if (L.x > 1527.6465 || L.x < 1527.6464 || L.y > 1597.9472 || L.y < 1597.9471 || L.z > 3393.8636 || L.z < 3393.8635)
    {
        printf("\tMomentum L calculation failed\n");
        printf("Calculated L = [%.15f, %.15f, %.15f]\n", L.x, L.y, L.z);
        return 1;
    }
    if (LErr.x > 125.7330 || LErr.x < 125.7329 || LErr.y > 173.2863 || LErr.y < 173.2862 || LErr.z > 50.2294 || LErr.z < 50.2293)
    {
        printf("\tMomentum LErr calculation failed\n");
        printf("Calculated LErr = [%.15f, %.15f, %.15f]\n", LErr.x, LErr.y, LErr.z);
        return 1;
    }

    /* Make sure reading in with EMD range also works*/
    data->histograms[0]->params.nRange = 4;
    data->histograms[0]->params.EMDRange[0] = -33;
    data->histograms[0]->params.EMDRange[1] = 0;
    data->histograms[0]->params.EMDRange[2] = 0;
    data->histograms[0]->params.EMDRange[3] = 48;

    nbCalcMomentum(&st, &ctx, data->histograms[0], histogram->histograms[0]);

    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wfloat-equal"
    if (L.x != histogram->histograms[0]->params.L.x
        || L.y != histogram->histograms[0]->params.L.y
        || L.z != histogram->histograms[0]->params.L.z)
    {
        printf("\tUsing EMD range did not return the same momentum result\n");
        return 1;
    }
    if (LErr.x != histogram->histograms[0]->params.LErr.x
        || LErr.y != histogram->histograms[0]->params.LErr.y
        || LErr.z != histogram->histograms[0]->params.LErr.z)
    {
        printf("\tUsing EMD range did not return the same momentum error\n");
        return 1;
    }
    #pragma GCC diagnostic pop
    //Final likelihood calculated in python is 771983.400390398
    //Note that likelihood calculated at this point is positive, and is made negative later

    real likelihood = nbMomentumLikelihood(data->histograms[0], histogram->histograms[0]);
    if (likelihood > 771983.4004 || likelihood < 771983.4003)
    {
        printf("\tMomentum likelihood calculation failed\n");
        return 1;
    }

    return 0;
}