/*
 *  Copyright (c) 2010-2011 Matthew Arsenault
 *  Copyright (c) 2010-2011 Rensselaer Polytechnic Institute
 *
 *  This file is part of Milkway@Home.
 *
 *  Milkway@Home is free software: you may copy, redistribute and/or modify it
 *  under the terms of the GNU General Public License as published by the
 *  Free Software Foundation, either version 3 of the License, or (at your
 *  option) any later version.
 *
 *  This file is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _NBODY_POTENTIAL_TYPES_H_
#define _NBODY_POTENTIAL_TYPES_H_


#include "nbody_config.h"
#include "milkyway_math.h"
#include "milkyway_extra.h"

#define _SPHERICAL 0

typedef enum
{
    InvalidSpherical   = InvalidEnum,
    SphericalPotential = _SPHERICAL
} spherical_t;

/* Spherical potential */
typedef struct MW_ALIGN_TYPE
{
    spherical_t type;
    real mass;
    real scale;
} Spherical;

#define SPHERICAL_TYPE "Spherical"


/* Can't get the enum value in preprocessor, so do this */
#define _MN_DISK 0
#define _EXP_DISK 1


/* Supported disk models */
typedef enum
{
    InvalidDisk       = InvalidEnum,
    MiyamotoNagaiDisk = _MN_DISK,
    ExponentialDisk   = _EXP_DISK
} disk_t;

typedef struct MW_ALIGN_TYPE
{
    disk_t type;
    real mass;         /* disk mass */
    real scaleLength;  /* "a" for M-N, "b" for exp disk */
    real scaleHeight;  /* unused for exponential disk. "b" for Miyamoto-Nagai disk */
} Disk;

#define DISK_TYPE "Disk"

/* Supported halo models */

/* Can't get the enum value in preprocessor, so do this */
#define _LOG_HALO 0
#define _NFW_HALO 1
#define _TRIAXIAL_HALO 2
typedef enum
{
    InvalidHalo     = InvalidEnum,
    LogarithmicHalo = _LOG_HALO,
    NFWHalo         = _NFW_HALO,
    TriaxialHalo    = _TRIAXIAL_HALO
} halo_t;

typedef struct MW_ALIGN_TYPE
{
    halo_t type;
    real vhalo;         /* common to all 3 halos */
    real scaleLength;   /* common to all 3 halos */
    real flattenZ;      /* used by logarithmic and triaxial */
    real flattenY;      /* used by triaxial */
    real flattenX;      /* used by triaxial */
    real triaxAngle;    /* used by triaxial */

    real c1;           /* Constants calculated for triaxial from other params */
    real c2;           /* TODO: Lots more stuff could be cached, but should be done less stupidly */
    real c3;
} Halo;

#define HALO_TYPE "Halo"

typedef struct MW_ALIGN_TYPE
{
    Spherical sphere[1];
    Disk disk;
    Halo halo;
    void* rings;       /* currently unused */
} Potential;

#define POTENTIAL_TYPE "Potential"



#define EMPTY_SPHERICAL { InvalidSpherical, 0.0, 0.0 }
#define EMPTY_DISK { InvalidDisk, 0.0, 0.0, 0.0 }
#define EMPTY_HALO { InvalidHalo, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
#define EMPTY_POTENTIAL { {EMPTY_SPHERICAL}, EMPTY_DISK, EMPTY_HALO, NULL }


#endif /* _NBODY_POTENTIAL_TYPES_H_ */

