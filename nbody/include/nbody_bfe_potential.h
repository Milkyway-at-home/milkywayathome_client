#ifndef _NBODY_BFE_POTENTIAL_H_
#define _NBODY_BFE_POTENTIAL_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "milkyway_math.h"

/* struct exp_bfe_t is declared but deliberately never defined. */
struct exp_bfe_t;
typedef struct exp_bfe_t exp_bfe_t;

exp_bfe_t *exp_bfe_open(const char *yaml_filename);

/* Note that xyz.w is ignored */
mwvector exp_bfe_get_acceleration(exp_bfe_t *exp_bfe, mwvector xyz, real t);

void exp_bfe_close(exp_bfe_t *exp_bfe);

#ifdef __cplusplus
} // end of extern "C"
#endif

#endif /* _NBODY_BFE_POTENTIAL_H_ */
