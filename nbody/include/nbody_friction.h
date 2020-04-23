/*********************** DISCLAIMER *************************
* This code is included for its potential usefulness but    *
* has not been fully tested. It is not currently used in    *
* any other milkyway@home files. Use at your own risk.      *
************************************************************/
#ifndef _NBODY_FRICTION_H_
#define _NBODY_FRICTION_H_


#ifdef __cplusplus
extern "C" {
#endif

mwvector dynamicalFriction(mwvector pos, mwvector vel, real mass, const Potential* pot);

#ifdef __cplusplus
}
#endif

#endif /* _NBODY_FRICTION_H*/

