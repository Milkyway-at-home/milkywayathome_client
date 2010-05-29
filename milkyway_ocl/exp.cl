/*
Copyright 2010 Anthony Waters, Travis Desell,
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

//doesn't appear to work as a macro with nvidia
//cephes dp library
#define EXP(x)								\
  double px = (double) (int)(x * 1.4426950408889634073599246		\
			     + 0.5) - 1; //log2e, estimated floor	\
  int n = px;								\
  x -= px * 6.93145751953125E-1; //C1					\
  x -= px * 1.42860682030941723212E-6; //C2				\
  double xx = x * x;							\
  px = x * (( 1.26177193074810590878E-4 * xx +				\
	      3.02994407707441961300E-2) * xx +				\
	    9.99999999999999999910E-1); //polevl			\
  x = px / ((((3.00198505138664455042E-6 * xx +				\
	       2.52448340349684104192E-3) * xx +			\
	      2.27265548208155028766E-1) * xx +				\
	     2.00000000000000000009E0) - px); //polevl(xx,Q,3)		\
  x = 1.0 + 2.0 * x;							\
  //printf("x:%.15f n:=%d\n", x, n);					\
  //ldexp, x * 2^n							\
  if (n < 0)								\
    {									\
      if (n < -60)							\
	x = -0.0;							\
      else								\
	{								\
	  n = n * -1;							\
	  if (n > 30)							\
	    {								\
	      x = x / (double) (2 << (30-1));				\
	      n -= 30;							\
	    }								\
	  x = x / (double) (2 << (n-1));				\
	}								\
    }									\
  else if (n > 0)							\
    {									\
      x = x * (double) (2 << (n-1));					\
    }
