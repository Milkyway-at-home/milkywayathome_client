# Multi precision toolbox for Scilab
# Copyright (C) 2009 - Jonathan Blanchard
#
# This file must be used under the terms of the CeCILL.
# This source file is licensed as described in the file COPYING, which
# you should have received as part of this distribution.  The terms
# are also available at
# http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt

# - Attempts to find GNU MPFR
# The following variables will be defined
#
#  LIBMPFR_FOUND - MPFR was found
#  LIBMPFR_INCLUDE_DIR - The MPFR include directory
#  LIBMPFR_LIB - The MPFR shared library full path
#  LIBMPFR_DIR - MPFR main path, i.e. the one that contains include and lib.

# Check if LIBMPFR_DIR is defined and use that path first.
IF(LIBMPFR_DIR)
	FIND_PATH(LIBMPFR_INCLUDE_DIR mpfr.h PATHS ${LIBMPFR_DIR}/include NO_DEFAULT_PATH)
	FIND_LIBRARY(LIBMPFR_LIB mpfr PATHS ${LIBMPFR_DIR}/lib NO_DEFAULT_PATH)

	IF(NOT LIBMPFR_LIB)
		MESSAGE(STATUS "Warning : MPFR not found in the path specified in LIBMPFR_DIR")
		UNSET(LIBMPFR_DIR)
	ENDIF()
ENDIF()

FIND_PATH(LIBMPFR_INCLUDE_DIR mpfr.h)

FIND_LIBRARY(LIBMPFR_LIB NAMES mpfr)

# Get the root MPFR path and add a cache entry.
GET_FILENAME_COMPONENT(LIBMPFR_DIR ${LIBMPFR_INCLUDE_DIR} PATH)
SET(LIBMPFR_DIR ${LIBMPFR_DIR} CACHE PATH "MPFR root directory.")

INCLUDE(FindPackageHandleStandardArgs)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(LibMPFR DEFAULT_MSG LIBMPFR_LIB LIBMPFR_INCLUDE_DIR)

MARK_AS_ADVANCED(LIBMPFR_INCLUDE_DIR LIBMPFR_LIB)
