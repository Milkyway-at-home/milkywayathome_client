/*
 * Copyright (c) 2012 Matthew Arsenault
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

#ifndef _NBODY_GL_INCLUDES_H_
#define _NBODY_GL_INCLUDES_H_

#include "nbody_config.h"

#define GLFW_INCLUDE_GLCOREARB
#if USE_GL3W
  #include <GL3/gl3w.h>
  #define GLFW_NO_GLU 1
  #include <GL/glfw3.h>
#else
  #define GL3_PROTOTYPES
  #define GLFW_INCLUDE_GL3 1
  #define GLFW_NO_GLU 1
  #include <GL/glfw3.h>
#endif

#include "milkyway_gcc_diagnostics.h"

/* Only recentish GCC's support #pragma GCC diagnostic push / pop */
GCC_DIAG_OFF(shadow)
GCC_DIAG_OFF(float-equal)
GCC_DIAG_OFF(type-limits)
GCC_DIAG_OFF(unused-parameter)
GCC_DIAG_OFF(switch-default)
/* GCC_DIAG_OFF(unknown-pragmas) */

#include <glm/glm.hpp>
#include <glm/ext.hpp>
#include <glm/gtc/random.hpp>
#include <glm/gtx/fast_square_root.hpp>
#include <glm/gtx/quaternion.hpp>

#include "MousePoles.h"

GCC_DIAG_ON(shadow)
GCC_DIAG_ON(float-equal)
GCC_DIAG_ON(type-limits)
GCC_DIAG_ON(unused-parameter)
GCC_DIAG_ON(switch-default)
/* GCC_DIAG_ON(unknown-pragmas) */


#endif /* _NBODY_GL_INCLUDES_H_ */



