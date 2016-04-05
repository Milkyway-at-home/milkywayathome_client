/*
 * Copyright (c) 2012 Rensselaer Polytechnic Institute
 *
 * This file is part of Milkway@Home.
 *
 * Milkyway@Home is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Milkyway@Home is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Milkyway@Home.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "nbody_mass.h"

#include "milkyway_math.h"
#include "nbody_types.h"


/*In order to decrease the size of the numbers
 * computed all these functions are
 * calculated in log space*/

static real factorial(int n){
     int counter;
     real result = 0.0;

     for (counter = n; counter >= 1; counter--)
       {
          result += mw_log((real) counter);
       }

     return result;
  }


static real choose(int n, int c)
{
    unsigned int i;
    real result = 0.0;
    
    /* This for loop calulates log(n!/(n-c)!) */
    for (i = n - c + 1; i <= (unsigned int) n; ++i)
    {
        result += mw_log(i);
    }
//     mw_printf("i = % 10.10i\n",i);
//     mw_printf("c = % 10.10i\n", c);
//     mw_printf("n - c + 1 = % 10.10i\n", n - c + 1);
//     mw_printf("result = % 10.10f\n", result);
//     mw_printf("factorial(c) = % 10.10f\n",factorial(c));
    result -= factorial(c);
    return result;
}

real probability_match(int n, real ktmp, real pobs)
{
    real result = 0.0;
    real result1, result2, result3;
    // result +=  (real) mw_log(choose(n, k));                                                                                                                                                        
    // result += mw_log(pobs) * (real) k;                                                                                                                                                             
    // result += mw_log(1.0 - pobs) * (real)(n - k);                                                                                                                                                  

    /*
     * Previously, this function took in k as an int. Bad move.
     * This function was called twice, one of which sent a real valued k: (int) k1 and (real) k2
     * That real k2 was then converted to int. Could result in converted (int) k1 != (int) k2 when k1 = k2. 
     * Special result was poor likelihood for some histograms when check against themselves!
     * General results: unknown. But probably not good. (most likely caused different machines to report
     * different likelihood values).
     * 
     */
    int k = (int) mw_round(ktmp);    //patch. See above. 
//     mw_printf("prob n = %i\n", n);
//     mw_printf("prob ktmp = %f\n", ktmp);
//     mw_printf("prob pObs = %f\n", pobs);
//     mw_printf("prob k = %i\n", k);
    //The previous calculation does not return the right values.  Furthermore, we need a zeroed metric.                                                                                              
    result =  (real) choose(n, k);
    result += k * mw_log(pobs); 
    result += (n - k) * mw_log(1.0 - pobs);
    
//     result1 =  (real) choose(n, k);
//     result2 = k * mw_log(pobs); 
//     result3 = (n - k) * mw_log(1.0 - pobs);
//     mw_printf("factorial result = %f\n", result1);
//     mw_printf("p^k result = %f\n", result2);
//     mw_printf("(1-p)^n-k result = %f\n", result3);
//     mw_printf("sum = %f\n", result);
    
//     mw_printf("returned value= %f\n\n", mw_pow(M_E, result));
//     mw_printf("Coeff = %10.10f\n", choose(n,k));
//     mw_printf("Term 1 = %10.10f\n", k * mw_log(pobs));
//     mw_printf("Term 2 = %10.10f\n", (n-k) * mw_log(1 - pobs));
//     mw_printf("result = %10.10f\n",result);
    return mw_exp(result);
}


