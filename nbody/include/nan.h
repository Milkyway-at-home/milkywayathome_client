/* `NAN' constant for IEEE 754 machines.

Copyright (C) 1992 Free Software Foundation, Inc.
This file is part of the GNU C Library.

The GNU C Library is free software; you can redistribute it and/or
modify it under the terms of the GNU Library General Public License as
published by the Free Software Foundation; either version 2 of the
License, or (at your option) any later version.

The GNU C Library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Library General Public License for more details.

You should have received a copy of the GNU Library General Public
License along with the GNU C Library; see the file COPYING.LIB.  If
not, write to the Free Software Foundation, Inc., 675 Mass Ave,
Cambridge, MA 02139, USA.  */

#ifndef	_NAN_H

#define	_NAN_H	1

/* IEEE Not A Number.  */

#include <endian.h>

#if	__BYTE_ORDER == __BIG_ENDIAN
#define	__nan_bytes		{ 0x7f, 0xf8, 0, 0, 0, 0, 0, 0 }
#endif

#if	__BYTE_ORDER == __LITTLE_ENDIAN
#define	__nan_bytes		{ 0, 0, 0, 0, 0, 0, 0xf8, 0x7f }
#endif

#ifdef __i386__
/* Signal double NaN. */
#define	__snan_bytes		{ 0, 0, 0x0f, 0, 0, 0, 0xf0, 0x7f }
/* Quiet double NaN, It is also indefinite(?). */
#define	__qnan_bytes		{ 0, 0, 0, 0, 0, 0, 0xf8, 0x7f }
/* Signal float NaN. */
#define	__snanf_bytes		{ 0, 0x0f, 0x80, 0x7f }
/* Quiet float NaN. It is also indefinite(?). */
#define	__qnanf_bytes		{ 0, 0, 0xc0, 0x7f }
#endif

#ifdef	__GNUC__
#define	NAN \
  (__extension__ ((union { unsigned char __c[8];		      \
			   double __d; })			      \
		  { __nan_bytes }).__d)

#ifdef __i386__
#define	_SNAN \
  (__extension__ ((union { unsigned char __c[8];		      \
			   double __d; })			      \
		  { __snan_bytes }).__d)
#define	_QNAN \
  (__extension__ ((union { unsigned char __c[8];		      \
			   double __d; })			      \
		  { __qnan_bytes }).__d)
#define	_SNANF \
  (__extension__ ((union { unsigned char __c[4];		      \
			   float __f; })			      \
		  { __snanf_bytes }).__f)
#define	_QNANF \
  (__extension__ ((union { unsigned char __c[4];		      \
			   float __f; })			      \
		  { __qnanf_bytes }).__f)
#endif
#else	/* Not GCC.  */
static __const char __nan[8] = __nan_bytes;
#define	NAN	(*(__const double *) __nan)

#ifdef __i386__
static __const char __snan[8] = __snan_bytes;
#define	_SNAN	(*(__const double *) __snan)
static __const char __qnan[8] = __qnan_bytes;
#define	_QNAN	(*(__const double *) __qnan)
static __const char __snanf[4] = __snanf_bytes;
#define	_SNANF	(*(__const float *) __snanf)
static __const char __qnanf[4] = __qnanf_bytes;
#define	_QNANF	(*(__const float *) __qnanf)
#endif
#endif	/* GCC.  */

#endif	/* nan.h */
