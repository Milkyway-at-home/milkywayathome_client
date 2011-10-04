/*
 * Copyright (c) 2011 Matthew Arsenault
 * Copyright (c) 2011 Rensselaer Polytechnic Institute
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

#ifndef _REPLACE_AMD_IL_H_
#define _REPLACE_AMD_IL_H_

#include <stddef.h>
#include "milkyway_cl.h"

#ifdef __cplusplus
extern "C" {
#endif

unsigned char* getModifiedAMDBinary(unsigned char* bin,
                                    size_t binSize,
                                    int nStream,
                                    MWCALtargetEnum target,
                                    size_t* newBinLenOut);

#ifdef __cplusplus
}
#endif

#endif /* _REPLACE_AMD_IL_H_ */

