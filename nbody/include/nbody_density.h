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

static inline real hernquistSphericalDensity(const Spherical* sph, real r);
static inline real plummerSphericalDensity(const Spherical* sph, real r);

static inline real miyamotoNagaiDiskDensity(const Disk* disk, mwvector pos);
static inline real doubleExponentialDiskDensity(const Disk* disk, mwvector pos);
static inline real sech2ExponentialDiskDensity(const Disk* disk, mwvector pos);

static inline real logarithmicHaloDensity(const Halo* h, mwvector pos);
static inline real triaxialHaloDensity(const Halo* h, mwvector pos);
static inline real hernquistHaloDensity(const Halo* h,  real r);
static inline real plummerHaloDensity(const Halo* h, real r);
static inline real NFWMHaloDensity(const Halo* h,  real r);
static inline real allenSantillanHaloDensity(const Halo* h, real r);
static inline real wilkinsonEvansHaloDensity(const Halo* h, real r);
static inline real ninkovicHaloDensity(const Halo* h, real r);
static inline real KVHalo(const Halo* h, real r);

#ifdef __cplusplus
}
#endif

#endif /* _NBODY_DENSITY_H_ */

