# Copyright 2010 Matthew Arsenault, Travis Desell, Dave Przybylo,
# Nathan Cole, Boleslaw Szymanski, Heidi Newberg, Carlos Varela, Malik
# Magdon-Ismail and Rensselaer Polytechnic Institute.

# This file is part of Milkway@Home.

# Milkyway@Home is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# Milkyway@Home is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with Milkyway@Home.  If not, see <http://www.gnu.org/licenses/>.
#

#fdlibm seems to not have a single, nice home, so just include the
#entire source

cmake_minimum_required(VERSION 2.8)
project(fdlibm)

set(fdlibm_MAJOR_VERSION 5)
set(fdlibm_MINOR_VERSION 3)
set(fdlibm_VERSION "${fdlibm_MAJOR_VERSION}.${fdlibm_MINOR_VERSION}")

include_directories(${PROJECT_INCLUDE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
link_directories(${LIBRARY_OUTPUT_PATH})

# Note taken from original Makefile:
#
#  @(#)Makefile 1.4 95/01/18
#
#  ====================================================
#  Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
#
#  Developed at SunSoft, a Sun Microsystems, Inc. business.
#  Permission to use, copy, modify, and distribute this
#  software is freely granted, provided that this notice
#  is preserved.
#  ====================================================
#
#
# There are two options in making libm at fdlibm compile time:
# 	_IEEE_LIBM 	--- IEEE libm; smaller, and somewhat faster
#	_MULTI_LIBM	--- Support multi-standard at runtime by
#			    imposing wrapper functions defined in
#			    fdlibm.h:
#				_IEEE_MODE 	-- IEEE
#				_XOPEN_MODE 	-- X/OPEN
#				_POSIX_MODE 	-- POSIX/ANSI
#				_SVID3_MODE 	-- SVID
#
# Here is how to set up CFLAGS to create the desired libm at
# compile time:
#
# 	CFLAGS = -D_IEEE_LIBM		... IEEE libm (recommended)
#	CFLAGS = -D_SVID3_MODE	... Multi-standard supported
#					    libm with SVID as the
#					    default standard
#	CFLAGS = -D_XOPEN_MODE	... Multi-standard supported
#					    libm with XOPEN as the
#					    default standard
#	CFLAGS = -D_POSIX_MODE	... Multi-standard supported
#					    libm with POSIX as the
#					    default standard
#	CFLAGS = 			... Multi-standard supported
#					    libm with IEEE as the
#					    default standard
#
# NOTE: if scalb's second arguement is an int, then one must
# define _SCALB_INT in CFLAGS. The default prototype of scalb
# is double scalb(double, double)
#

option(USE_IEEE_LIBM "Use IEEE libm (recommended)" ON)
#TODO: Do we care about the others? I doubt it.

mark_as_advanced(USE_IEEE_LIBM)

set(fdlibm_lib_src
  e_acos.c     e_lgamma.c     k_tan.c       s_isnan.c        s_tanh.c     w_jn.c
  e_acosh.c    e_lgamma_r.c   s_asinh.c     s_ldexp.c        w_acos.c     w_lgamma.c
  e_asin.c     e_log10.c      s_atan.c      s_lib_version.c  w_acosh.c    w_lgamma_r.c
  e_atan2.c    e_log.c        s_cbrt.c      s_log1p.c        w_asin.c     w_log10.c
  e_atanh.c    e_pow.c        s_ceil.c      s_logb.c         w_atan2.c    w_log.c
  e_cosh.c     e_remainder.c  s_copysign.c  s_matherr.c      w_atanh.c    w_pow.c
  e_exp.c      e_rem_pio2.c   s_cos.c       s_modf.c         w_cosh.c     w_remainder.c
  e_fmod.c     e_scalb.c      s_erf.c       s_nextafter.c    w_exp.c      w_scalb.c
  e_gamma.c    e_sinh.c       s_expm1.c     s_rint.c         w_fmod.c     w_sinh.c
  e_gamma_r.c  e_sqrt.c       s_fabs.c      s_scalbn.c       w_gamma.c    w_sqrt.c
  e_hypot.c    k_cos.c        s_finite.c    s_signgam.c      w_gamma_r.c
  e_j0.c       k_rem_pio2.c   s_floor.c     s_significand.c  w_hypot.c
  e_j1.c       k_sin.c        s_frexp.c     s_sin.c          w_j0.c
  e_jn.c       k_standard.c   s_ilogb.c     s_tan.c          w_j1.c)

#TODO: Other acceptable flags from CFLAGS, and required march on 32 bit
add_library(fdlibm ${fdlibm_lib_src})
set_target_properties(fdlibm
  PROPERTIES
    COMPILE_DEFINITIONS "_IEEE_LIBM"
    COMPILE_FLAGS
     "-mfpmath=sse -msse -msse2 -O3 -std=c99 -fno-strict-aliasing -fomit-frame-pointer -momit-leaf-frame-pointer -w")


