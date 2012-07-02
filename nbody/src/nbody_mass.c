#include "milkyway_math.h"
#include "nbody_types.h"


/*In order to decrease the size of the numbers
 * computed all these functions are
 * calculated in log space*/

static real factorial(int n) {
    if (n <= 1){
        return 0;
    }
    return mw_log(n) + factorial(n-1);
}
    
static real choose(int n, int c) {
	unsigned int i;
	real result = 0;
	/*This for loop calulates log(n!/(n-c)!)*/
	for(i = n-c+1; i <= (unsigned int)n; i++){
		result += mw_log(i);
	}
	result -= factorial(c);
    return result;
}

real probability_match(int n, int k, real pobs){
    real result;
    result = choose(n, k);
    result += mw_log(pobs) * (real) k;
    result += mw_log(1-pobs) * (real) (n-k);
    return result;
}
