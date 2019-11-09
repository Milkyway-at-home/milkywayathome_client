


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

