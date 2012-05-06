#include "milkyway_math.h"
#include "nbody_types.h"

int factorial(int n) {
    if (n <= 0){
        return 1;
    }
    return n * factorial(n-1);
}
    
real choose(int n, int c) {
    return (real)factorial(n) / (real)(factorial(c) * (real)factorial(n-c));
}

real propability_match(int n, int k, real pobs){
    real result;
    result = choose(n, k);
    result *= mw_pow(pobs, (real)k);
    result *= mw_pow((1-pobs), (real)(n-k));
    result = mw_log(result);
    return result;
}
