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

#ifndef COORDS_H
#define COORDS_H

#define rmat00 -0.054875539726
#define rmat01 -0.873437108010
#define rmat02 -0.483834985808
#define rmat10  0.494109453312
#define rmat11 -0.444829589425
#define rmat12  0.746982251810
#define rmat20 -0.867666135858
#define rmat21 -0.198076386122
#define rmat22  0.455983795705


#define stripe_separation 2.5
#define survey_center_dec 32.5
#define survey_center_ra 185.0

#define F_A_NODE 	95.0	//survey_center_ra - 90.0
#define F_A_NODE_RAD	(float)1.65806281566619873046875

#define D_A_NODE	95.0
#define D_A_NODE_RAD	D_A_NODE * D_DEG2RAD

double d_get_incl(int wedge);
double d_get_incl_rad(int wedge);

float f_get_incl(int wedge);
float f_get_incl_rad(int wedge);

void d_lbr2xyz(const double *lbr, double *xyz);

#endif
