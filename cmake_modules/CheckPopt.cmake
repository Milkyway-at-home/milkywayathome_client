# Multi precision toolbox for Scilab
# Copyright (C) 2009 - Jonathan Blanchard
#
# This file must be used under the terms of the CeCILL.
# This source file is licensed as described in the file COPYING, which
# you should have received as part of this distribution.  The terms
# are also available at
# http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt

# Check if popt is available.

FUNCTION(CHECK_POPT)

    IF(DEFINED HAVE_POPT)
        RETURN()
    ENDIF()

    find_path(POPT_INCLUDE_DIR popt.h)
    find_library(POPT_LIBRARY popt)

    if(POPT_INCLUDE_DIR AND POPT_LIBRARY)
        set(POPT_FOUND TRUE)
    endif(POPT_INCLUDE_DIR AND POPT_LIBRARY)



    IF(HAVE_POPT_H)
        MESSAGE(STATUS "Checking if popt is available")

        TRY_COMPILE(POPT_CHECK ${CMAKE_BINARY_DIR} ${CMAKE_MODULE_PATH}/poptcheck.c)
ck.c
    ENDIF()



    IF(POPT_CHECK)
        SET(HAVE_POPT 1 CACHE BOOL "Set if popt is available.")
        MESSAGE(STATUS "Checking if popt is available - yes")
    ELSEIF(HAVE_POPT_H)
        SET(HAVE_POPT 0 CACHE BOOL "Set if popt is available.")
        MESSAGE(STATUS "Checking if popt is available - error")
    ELSE()
        SET(HAVE_POPT 0 CACHE BOOL "Set if popt is available.")
        MESSAGE(STATUS "Checking if popt is available - no")
    ENDIF()

    MARK_AS_ADVANCED(HAVE_POPT)

ENDFUNCTION()

