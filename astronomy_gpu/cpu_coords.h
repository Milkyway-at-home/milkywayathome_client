/*
 * Copyright 2008, 2009 Travis Desell, Dave Przybylo, Nathan Cole,
 * Boleslaw Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail
 * and Rensselaer Polytechnic Institute.
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
 * */

#ifndef CPU_COORDS_H
#define CPU_COORDS_H

#include "../astronomy/parameters.h"


void populate_lb(int sgr_coordinates, int wedge,
		 int mu_steps, double mu_min, double mu_step_size,
		 int nu_steps, double nu_min, double nu_step_size,
		 double *sinb, double *sinl, double *cosb, double *cosl);

void gc_eq_gal_lb(int wedge, double amu_rad, double anu_rad, double *cpu__lb);
void gc_eq_gal(int wedge, double amu_rad, double anu_rad, double *glong, double *glat);
void gc_sgr_gal_lb(int wedge, double amu_rad, double anu_rad, double *cpu__lb);
void gc_sgr_gal(int wedge, double amu_rad, double anu_rad, double *glong, double *glat);
void cpu__gc_to_lb(int wedge, INTEGRAL *integral, double **cpu__lb);
void cpu__gc_eq_gal_lb(int wedge, int mu_steps, double mu_min, double mu_step_size, int nu_steps, double nu_min, double nu_step_size, double **cpu__lb);
void cpu__gc_sgr_gal_lb(int wedge, int mu_steps, double mu_min, double mu_step_size, int nu_steps, double nu_min, double nu_step_size, double **cpu__lb);

#endif
