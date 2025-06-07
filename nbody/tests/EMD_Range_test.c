/*This test tests the ability to read in lambda ranges from the histogram and use these as the ranges to calculate EMD*/
/*When no range is given, should default to using all given lambda bins*/

#include "nbody_mass.h"
#include "nbody_histogram.h"
#include "nbody_types.h"
#include "nbody_emd.h"
#include <stdio.h>

int main()
{
    MainStruct* histogram_1;
    MainStruct* histogram_2;
    real scoreWithRange;
    real scoreWithoutRange;

    histogram_1 = nbReadHistogram("./EMD_range_1.hist"); // Has a lambda range defined that should be used 
    histogram_2 = nbReadHistogram("./EMD_range_2.hist"); // Has no lambda range defined and should default to using all bins

    scoreWithRange = nbMatchEMD(histogram_1, histogram_2); // With given range, should be an exact match for EMD
    scoreWithoutRange = nbMatchEMD(histogram_2, histogram_1); // Should do a normal comparison and give the likelihood found before EMD range was added

    mw_printf("Score with range: %f\n", scoreWithRange);
    mw_printf("Score without range: %f\n", scoreWithoutRange);

    if(scoreWithRange != 0)
    {
        printf("EMD calculated with given ranges not correct, gave %f should be 0 \n", scoreWithRange);
        return 1;
    }

    if(scoreWithoutRange < 5.449190 || scoreWithoutRange > 5.449192) // should expect score of -5.449191188301354
    {
        printf("EMD calculated with whole histogram is not correct; either the function could not default to normal behavior or the EMD likelihood calculation has been changed. Got %1.15f, expected 5.449191188301354 \n", scoreWithoutRange);
        return 1;
    }
    else //both cases return expected likelihood
    {
        return 0;
    }
}