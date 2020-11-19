/*********************** DISCLAIMER *************************
* This code is included for its potential usefulness but    *
* has not been fully tested. It is not currently used in    *
* any other milkyway@home files. Use at your own risk.      *
************************************************************/
#ifndef _NBODY_DENSITY_H_
#define _NBODY_DENSITY_H_

#include "nbody_types.h"

#ifdef __cplusplus
extern "C" {
#endif

real nbExtDensity(const Potential* pot, mwvector pos);

#ifdef __cplusplus
}
#endif

#endif /* _NBODY_DENSITY_H_ */

