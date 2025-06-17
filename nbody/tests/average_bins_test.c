/*This program tests the likelihood calculation of the beta dispersion when more than one Lambda bin is considered.
It only applies to an odd number of bins centered on the comparison bin, and will reduce an even input to the next lowest odd number*/

#include "nbody_mass.h"
#include "nbody_histogram.h"
#include "nbody_types.h"
#include <stdio.h>

int main()
{
    MainStruct* histogram; // Input histogram
    MainStruct* histogram2; // Comparison histogram for average bins
    MainStruct* histogram3; // Comparison histogram for average bin comparison with one bad bin excluded
    MainStruct* histogram4; // Comparison histogram for avgBin values that result in default behavior
    MainStruct* histogram5; // Input histogram used against hist3 to test behavior on histogram boundary

    histogram  = nbReadHistogram("./average_bins_test.hist"); // Test that betaDispBins can be read in from a histogram
    histogram2 = nbReadHistogram("./average_bins_test.hist");
    histogram3 = nbReadHistogram("./average_bins_test.hist");
    histogram4 = nbReadHistogram("./average_bins_test.hist");
    histogram5 = nbReadHistogram("./average_bins_test.hist");

    //if getting segFault, ensure that test hist is in /nbody/tests

    histogram2->histograms[1]->lambdaBins = 3;
    histogram2->histograms[1]->betaBins = 1;
    histogram2->histograms[1]->totalNum = 10;
    histogram2->histograms[1]->totalSimulated = 3;

    histogram2->histograms[1]->data[0].useBin= 1;
    histogram2->histograms[1]->data[1].useBin= 1;
    histogram2->histograms[1]->data[2].useBin= 1;
    histogram2->histograms[1]->data[0].variable= 10;  
    histogram2->histograms[1]->data[1].variable= 22;
    histogram2->histograms[1]->data[2].variable= 10;
    histogram2->histograms[1]->data[0].err= 3*mw_sqrt(5.0);      
    histogram2->histograms[1]->data[1].err= 3*mw_sqrt(15.0); 
    histogram2->histograms[1]->data[2].err= 3*mw_sqrt(5.0); 

    histogram3->histograms[1]->lambdaBins = 3;
    histogram3->histograms[1]->betaBins = 1;
    histogram3->histograms[1]->totalNum = 10;
    histogram3->histograms[1]->totalSimulated = 3;

    histogram3->histograms[1]->data[0].useBin= 1;
    histogram3->histograms[1]->data[1].useBin= 1;
    histogram3->histograms[1]->data[2].useBin= 1;
    histogram3->histograms[1]->data[0].variable= 0;  
    histogram3->histograms[1]->data[1].variable= 12;
    histogram3->histograms[1]->data[2].variable= 16;
    histogram3->histograms[1]->data[0].err= 0;      
    histogram3->histograms[1]->data[1].err= 6; 
    histogram3->histograms[1]->data[2].err= 8; 

    histogram4->histograms[1]->lambdaBins = 3;
    histogram4->histograms[1]->betaBins = 1;
    histogram4->histograms[1]->totalNum = 10;
    histogram4->histograms[1]->totalSimulated = 3;

    histogram4->histograms[1]->data[0].useBin= 1;
    histogram4->histograms[1]->data[1].useBin= 1;
    histogram4->histograms[1]->data[2].useBin= 1;
    histogram4->histograms[1]->data[0].variable= 6;  
    histogram4->histograms[1]->data[1].variable= 14;
    histogram4->histograms[1]->data[2].variable= 4;
    histogram4->histograms[1]->data[0].err= 2;      
    histogram4->histograms[1]->data[1].err= 5; 
    histogram4->histograms[1]->data[2].err= 12; 

    histogram5->histograms[1]->lambdaBins = 3;
    histogram5->histograms[1]->betaBins = 1;
    histogram5->histograms[1]->totalNum = 10;
    histogram5->histograms[1]->totalSimulated = 3;
    histogram5->histograms[1]->betaDispBins = 3;

    histogram5->histograms[1]->data[0].useBin= 1;
    histogram5->histograms[1]->data[1].useBin= 1;
    histogram5->histograms[1]->data[2].useBin= 1;
    histogram5->histograms[1]->data[0].variable= 0;  
    histogram5->histograms[1]->data[1].variable= 0;
    histogram5->histograms[1]->data[2].variable= 1;
    histogram5->histograms[1]->data[0].err= -1;      
    histogram5->histograms[1]->data[1].err= -1; 
    histogram5->histograms[1]->data[2].err= 12; 

    real avgScore =  nbLikelihood(histogram->histograms[1], histogram2->histograms[1], histogram->histograms[1]->betaDispBins); //var 1 = 1, var 2 = 14, err 1 = 12, err 2 = 5, should give total difference and error of 13 -> 1 sigma difference -> likelihood 1
    real avgScore2 =  nbLikelihood(histogram->histograms[1], histogram3->histograms[1], histogram->histograms[1]->betaDispBins); //same numbers as above
    real avgScore3 =  nbLikelihood(histogram->histograms[1], histogram2->histograms[1], 4); //same numbers as above
    real score =  nbLikelihood(histogram->histograms[1], histogram4->histograms[1], 0); //only center bin considered, so same numbers as above
    real score1 =  nbLikelihood(histogram->histograms[1], histogram4->histograms[1], 1); //only center bin considered, so same numbers as above
    real score2 =  nbLikelihood(histogram->histograms[1], histogram4->histograms[1], 2); //only center bin considered, so same numbers as above
    real boundScore =  nbLikelihood(histogram5->histograms[1], histogram3->histograms[1], histogram5->histograms[1]->betaDispBins); //same numbers as above, only bottom two bins of hist3 should be considered
    if (avgScore < .499 || avgScore > .501)
    {
        printf("\tFailed standard average bin calculation: likelihood was %2.2f, expected .5", avgScore);
        return 1;
    }
    else if (avgScore2 < .499 || avgScore2 > .501)
    {
        printf("\tFailed average bin calculation with missing data: likelihood was %2.2f, expected .5", avgScore2);
        return 1;
    }
    else if (avgScore2 < .499 || avgScore2 > .501)
    {
        printf("\tFailed average bin calculation with avgBins=4: likelihood was %2.2f, expected .5", avgScore2);
        return 1;
    }
    else if (score < .499 || score > .501)
    {
        printf("\tFailed standard bin calculation with avgBins=0: likelihood was %2.2f, expected .5", score);
        return 1;
    }
    else if (score1 < .499 || score1 > .501)
    {
        printf("\tFailed standard bin calculation with avgBins=1: likelihood was %2.2f, expected .5", score1);
        return 1;
    }
    else if (score2 < .499 || score2 > .501)
    {
        printf("\tFailed standard bin calculation with avgBins=2: likelihood was %2.2f, expected .5", avgScore);
        return 1;
    }
    else if (boundScore < .499 || boundScore > .501)
    {
        printf("\tFailed average bin calculation on boundary of histogram: likelihood was %2.2f, expected .5", avgScore);
        return 1;
    }
    else // all likelihoods found were what we expected, test passed
    {
        return 0;
    }
}    
