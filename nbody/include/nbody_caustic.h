#ifndef _NBODY_CAUSTIC_H_
#define _NBODY_CAUSTIC_H_


#include <complex.h>
#include "nbody_config.h"
#include "nbody_types.h"
#include "nbody_lua.h"
#include "nbody_util.h"
#include "milkyway_math.h"
#include "milkyway_extra.h"

mwvector causticHaloAccel(const Halo* h, mwvector pos, real r);


#endif