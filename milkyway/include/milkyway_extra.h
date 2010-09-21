/* Copyright 2010 Matthew Arsenault, Travis Desell, Boleslaw
Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail and
Rensselaer Polytechnic Institute.

This file is part of Milkway@Home.

Milkyway@Home is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Milkyway@Home is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Milkyway@Home.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef _MILKYWAY_EXTRA_H_
#define _MILKYWAY_EXTRA_H_

/* The ATI CL compiler is dumb and has a problem with static anything */
#ifndef __ATI_CL__ /* FIXME: What is actually defined? */
  #define _MW_STATIC static
#else
  #define _MW_STATIC
#endif /* __ATI_CL__ */

#ifdef _MSC_VER
  #define strdup _strdup
  #define isnan _isnan
  #define copysign _copysign
  #define access _access
#endif /* _MSC_VER */

#endif /* _MILKYWAY_EXTRA_H_ */

