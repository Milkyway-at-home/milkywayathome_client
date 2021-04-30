# Multi precision toolbox for Scilab
# Copyright (C) 2009 - Jonathan Blanchard
#
# This file must be used under the terms of the CeCILL.
# This source file is licensed as described in the file COPYING, which
# you should have received as part of this distribution.  The terms
# are also available at
# http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt

# - Attempts to find GMP the GNU Multiple Precision Arithmetic Library
# The following variables will be defined
#
#  LIBGMP_FOUND - GMP was found
#  LIBGMP_INCLUDE_DIR - The GMP include directory
#  LIBGMP_LIB - The GMP shared library full path
#  LIBGMP_DIR - GMP main path, i.e. the one that contains include and lib.

# Check if LIBGMP_DIR is defined and use that path first.
IF(LIBGMP_DIR)
	FIND_PATH(LIBGMP_INCLUDE_DIR gmp.h PATHS ${LIBGMP_DIR}/include NO_DEFAULT_PATH)
	FIND_LIBRARY(LIBGMP_LIB gmp PATHS ${LIBGMP_DIR}/lib NO_DEFAULT_PATH)

	IF(NOT LIBGMP_LIB)
		MESSAGE(STATUS "Warning : GMP not found in the path specified in LIBGMP_DIR")
		UNSET(LIBGMP_DIR)
	ENDIF()
ENDIF()

FIND_PATH(LIBGMP_INCLUDE_DIR gmp.h)

FIND_LIBRARY(LIBGMP_LIB NAMES gmp)

# Get the root GMP path and add a cache entry.
GET_FILENAME_COMPONENT(LIBGMP_DIR ${LIBGMP_INCLUDE_DIR} PATH)
SET(LIBGMP_DIR ${LIBGMP_DIR} CACHE PATH "GMP root directory.")

INCLUDE(FindPackageHandleStandardArgs)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(LibGMP DEFAULT_MSG LIBGMP_LIB LIBGMP_INCLUDE_DIR)

MARK_AS_ADVANCED(LIBGMP_INCLUDE_DIR LIBGMP_LIB)
