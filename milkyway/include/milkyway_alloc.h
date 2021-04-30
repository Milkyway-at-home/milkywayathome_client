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

#ifndef _MILKYWAY_ALLOC_H_
#define _MILKYWAY_ALLOC_H_

#include "milkyway_config.h"
#include <stddef.h>


#ifdef __cplusplus
extern "C" {
#endif

#if HAVE___MINGW_ALIGNED_MALLOC
  #define mwFreeA __mingw_aligned_free
#elif HAVE__ALIGNED_MALLOC
  #define mwFreeA _aligned_free
#elif HAVE_POSIX_MEMALIGN
  #define mwFreeA free
#else
  #define mwFreeA free
#endif

#if HAVE___MINGW_ALIGNED_MALLOC
  #define _aligned_malloc __mingw_aligned_malloc
#endif


/* Allocations with that abort everything on failure */
void* mwMalloc(size_t size);
void* mwCalloc(size_t count, size_t size);
void* mwRealloc(void* ptr, size_t size);

/* Safe allocations aligned to 16 */
void* mwMallocA(size_t size);
void* mwCallocA(size_t count, size_t size);


/* Checks that we can actually align to 16/32 */
int mwAllocA32Safe(void);
int mwAllocA16Safe(void);

#ifdef __cplusplus
}
#endif

#endif /* _MILKYWAY_ALLOC_H_ */

