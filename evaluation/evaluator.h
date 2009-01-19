/*
Copyright 2008, 2009 Travis Desell, Dave Przybylo, Nathan Cole,
Boleslaw Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail
and Rensselaer Polytechnic Institute.

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

#ifndef FGDO_EVALUATOR_H
#define FGDO_EVALUATOR_H

extern double (*evaluate)(double*);

void evaluator__init(int *number_arguments, char*** arguments, void (*read_data)(int, int));
void evaluator__init_integral(void (*i_f)(double*, double**), int i_p_l, void (*i_c)(double*, int, double**), int i_r_l);
void evaluator__init_likelihood(void (*l_f)(double*, double**), int l_p_l, double (*l_c)(double*, int), int l_r_l);

#endif
