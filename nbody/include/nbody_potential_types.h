/*
 *  Copyright (c) 2010-2011 Matthew Arsenault
 *  Copyright (c) 2010-2011 Rensselaer Polytechnic Institute
 *  Copyright (c) 2016 Siddhartha Shelton
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

#define _NO_SPHERICAL 0
#define _HERN_SPHERICAL 1
#define _PLUMMER_SPHERICAL 2

typedef enum
{
    InvalidSpherical   = InvalidEnum,
    NoSpherical = _NO_SPHERICAL,
    HernquistSpherical = _HERN_SPHERICAL,
    PlummerSpherical = _PLUMMER_SPHERICAL
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
#define _NO_DISK 0
#define _MN_DISK 1
#define _FREEMAN_DISK 2
#define _DOUBEXPO_DISK 3
#define _SECHEXPO_DISK 4
#define _BAR 5


/* Supported disk models */
typedef enum
{
    InvalidDisk           = InvalidEnum,
    NoDisk                = _NO_DISK,
    MiyamotoNagaiDisk     = _MN_DISK,
    FreemanDisk           = _FREEMAN_DISK,
    DoubleExponentialDisk = _DOUBEXPO_DISK,
    Sech2ExponentialDisk  = _SECHEXPO_DISK,
    OrbitingBar  = _BAR
} disk_t;

typedef struct MW_ALIGN_TYPE
{
    disk_t type;
    real mass;         /* disk mass */
    real scaleLength;  /* "a" for M-N, "b" for exp disk */
    real scaleHeight;  /* unused for exponential disk. "b" for Miyamoto-Nagai disk */
    real patternSpeed; //for bars only
    real startAngle;   //for bars only
} Disk;

#define DISK_TYPE "Disk"

/* Supported halo models */

/* Can't get the enum value in preprocessor, so do this */
#define _NO_HALO 0
#define _LOG_HALO 1
#define _NFW_HALO 2
#define _TRIAXIAL_HALO 3
#define _CAUSTIC_HALO 4
#define _AS_HALO 5
#define _WE_HALO 6
#define _NFWM_HALO 7
#define _PLUMMER_HALO 8
#define _HERNQUIST_HALO 9
#define _NINKOVIC_HALO 10
typedef enum
{
    InvalidHalo        = InvalidEnum,
    NoHalo             = _NO_HALO,
    LogarithmicHalo    = _LOG_HALO,
    NFWHalo            = _NFW_HALO,
    TriaxialHalo       = _TRIAXIAL_HALO,
    CausticHalo        = _CAUSTIC_HALO,
    AllenSantillanHalo = _AS_HALO,
    WilkinsonEvansHalo = _WE_HALO,
    NFWMassHalo        = _NFWM_HALO,
    PlummerHalo        = _PLUMMER_HALO,
    HernquistHalo      = _HERNQUIST_HALO,
    NinkovicHalo       = _NINKOVIC_HALO
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

    real c1;            /* Constants calculated for triaxial from other params */
    real c2;            /* TODO: Lots more stuff could be cached, but should be done less stupidly */
    real c3;

    real mass;
    real gamma;
    real lambda;

    real rho0;          /*used by Ninkovic Halo*/
} Halo;

#define HALO_TYPE "Halo"

 /* Supported Dwarf Galaxy models */
#define _PLUMMER_DWARF 0
#define _NFW_DWARF 1
#define _GEN_HERN_DWARF 2
#define _EINASTO_DWARF 3
typedef enum
{
    InvalidDwarf       = InvalidEnum,
    Plummer            = _PLUMMER_DWARF,
    NFW                = _NFW_DWARF,
    General_Hernquist  = _GEN_HERN_DWARF,
    Einasto            = _EINASTO_DWARF
} dwarf_t;

typedef struct MW_ALIGN_TYPE
{
    dwarf_t type;
    real mass;        
    real scaleLength;   
    real n; //used by einasto
    real p0; //used by nfw
    real r200; // virial radius
} Dwarf;

#define DWARF_TYPE "Dwarf"


typedef struct MW_ALIGN_TYPE
{
    Spherical sphere[1];
    Disk disk;
    Disk disk2;
    Halo halo;
    void* rings;       /* currently unused */
} Potential;

#define POTENTIAL_TYPE "Potential"


#define EMPTY_SPHERICAL { InvalidSpherical, 0.0, 0.0 }
#define EMPTY_DISK { InvalidDisk, 0.0, 0.0, 0.0, 0.0, 0.0 }
#define EMPTY_DISK2 { InvalidDisk, 0.0, 0.0, 0.0, 0.0, 0.0 }
#define EMPTY_HALO { InvalidHalo, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
#define EMPTY_DWARF { InvalidDwarf, 0.0, 0.0, 0.0, 0.0, 0.0 }
#define EMPTY_POTENTIAL { {EMPTY_SPHERICAL}, EMPTY_DISK, EMPTY_DISK2, EMPTY_HALO, NULL }

#endif /* _NBODY_POTENTIAL_TYPES_H_ */

