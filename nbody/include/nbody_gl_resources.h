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

#ifndef _NBODY_GL_RESOURCES_H_
#define _NBODY_GL_RESOURCES_H_

#ifdef __cplusplus
#define EXTERNC extern "C"
#else
#define EXTERNC extern
#endif


EXTERNC unsigned char particle_vertex_glsl[];
EXTERNC size_t particle_vertex_glsl_len;

EXTERNC unsigned char particle_texture_fragment_glsl[];
EXTERNC size_t particle_texture_fragment_glsl_len;

EXTERNC unsigned char particle_point_fragment_glsl[];
EXTERNC size_t particle_point_fragment_glsl_len;


EXTERNC unsigned char axes_vertex_glsl[];
EXTERNC size_t axes_vertex_glsl_len;

EXTERNC unsigned char axes_fragment_glsl[];
EXTERNC size_t axes_fragment_glsl_len;


EXTERNC unsigned char galaxy_vertex_glsl[];
EXTERNC size_t galaxy_vertex_glsl_len;

EXTERNC unsigned char galaxy_fragment_glsl[];
EXTERNC size_t galaxy_fragment_glsl_len;


EXTERNC unsigned char text_vertex_glsl[];
EXTERNC size_t text_vertex_glsl_len;

EXTERNC unsigned char text_fragment_glsl[];
EXTERNC size_t text_fragment_glsl_len;

#if 0
EXTERNC const struct {
  unsigned int 	 width;
  unsigned int 	 height;
  unsigned int 	 bytes_per_pixel; /* 3:RGB, 4:RGBA */
  unsigned char	 pixel_data[1024 * 1024 * 3 + 1];
} milkywayImage;
#else
EXTERNC const struct {
  unsigned int 	 width;
  unsigned int 	 height;
  unsigned int 	 bytes_per_pixel; /* 3:RGB, 4:RGBA */
  unsigned char	 pixel_data[1024 * 1024 * 4 + 1];
} milkywayImage;
#endif


#undef EXTERNC

#endif /* _NBODY_GL_RESOURCES_H_ */

