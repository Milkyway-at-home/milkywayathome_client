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

#ifndef CPU_R_CONSTANTS_H
#define CPU_R_CONSTANTS_H

#include "parameters.h"

void cpu__r_constants(int n_convolve, INTEGRAL* integral, double** cpu__V, double** cpu__r_consts);
void cpu__r_constants(int n_convolve,
                      int r_steps,
                      double r_min,
                      double r_step_size,
                      int mu_steps,
                      double mu_min,
                      double mu_step_size,
                      int nu_steps,
                      double nu_min,
                      double nu_step_size,
                      double** cpu__V,
                      double** cpu__r_consts);

#endif

