# Multi precision toolbox for Scilab
# Copyright (C) 2009 - Jonathan Blanchard
#
# This file must be used under the terms of the CeCILL.
# This source file is licensed as described in the file COPYING, which
# you should have received as part of this distribution.  The terms
# are also available at
# http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt

# Check if the compiler supports C99.


FUNCTION(CHECK_C99_SUPPORT)

    IF(DEFINED HAVE_C99_SUPPORT)
        RETURN()
    ENDIF()

    MESSAGE(STATUS "Checking if the C compiler supports C99")

    TRY_COMPILE(C99_CHECK ${CMAKE_BINARY_DIR} ${CMAKE_MODULE_PATH}/c99check.c)

    IF(C99_CHECK)
        SET(HAVE_C99_SUPPORT 1 CACHE STRING "Status of C99 support.")
        MESSAGE(STATUS "Checking if the C compiler supports C99 - yes")
    ELSE()
        SET(HAVE_C99_SUPPORT 0 CACHE STRING "Status of C99 support.")
        MESSAGE(STATUS "Checking if the C compiler supports C99 - no")
    ENDIF()

    MARK_AS_ADVANCED(HAVE_C99_SUPPORT)

ENDFUNCTION()

