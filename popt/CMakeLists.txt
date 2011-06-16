#
# Taken from https://github.com/arsenm/POPT
# fc8836e177f3df560980dd923d565424326cbd1f (Sun Jun 5 14:06:11 2011 -0400)
#
# which was forked off of https://github.com/devzero2000/POPT
# at 30e86359f61b13b6ca0a (Thu May 5 10:16:38 2011 +0000)
#
# With the following changes:
#  - Renamed config.h to popt_config.h
#  - poptReadDefaultConfig #if 0'd out to remove need for configmake.h
#  - Tree reorganization
#


cmake_minimum_required(VERSION 2.6)
project(popt)

include(CheckIncludeFiles)
include(CheckLibraryExists)
include(CheckFunctionExists)
include(IgnoreWarningFlags)
include(CheckSSE2)

# We want this to be small, so use -Os unless we're debugging
if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_BUILD_TYPE MinSizeRel)
endif()

set(POPT_VERSION_MAJOR 1)
set(POPT_VERSION_MINOR 17)
set(POPT_VERSION "${POPT_VERSION_MAJOR}.${POPT_VERSION_MINOR}")

if(CMAKE_COMPILER_IS_GNUCC)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -std=gnu99")
endif()

check_function_exists(dcgettext HAVE_DCGETTEXT)

check_include_files(dlfcn.h HAVE_DLFCN_H)
check_include_files(float.h HAVE_FLOAT_H)
check_include_files(fnmatch.h HAVE_FNMATCH_H)

check_function_exists(gettext HAVE_GETTEXT)
check_include_files(glob.h HAVE_GLOB_H)

# We don't need anything libintl or libiconv provide, and they
# massively increase the final binary size, so just pretend we never have them

#check_function_exists(iconv HAVE_ICONV)

set(HAVE_ICONV FALSE)

check_include_files(langinfo.h HAVE_LANGINFO_H)
#check_include_files(libintl.h HAVE_LIBINTL_H)
set(HAVE_LIBINTL_H FALSE)

check_include_files(mcheck.h HAVE_MCHECK_H)
check_include_files(memory.h HAVE_MEMORY_H)

check_function_exists(mtrace HAVE_MTRACE)
check_function_exists(srandom HAVE_SRANDOM)

check_include_files(stdint.h HAVE_STDINT_H)
check_include_files(stdlib.h HAVE_STDLIB_H)

check_function_exists(stpcpy HAVE_STPCPY)
check_function_exists(strerror HAVE_STRERROR)

check_include_files(unistd.h HAVE_UNISTD_H)

check_function_exists(vasprintf HAVE_VASPRINTF)
check_function_exists(__secure_getenv HAVE___SECURE_GETENV)

set(libpopt_src src/popt.c
                src/poptparse.c
                src/poptconfig.c
                src/popthelp.c
                src/poptint.c)

set(libpopt_headers include/popt.h
                    include/poptint.h
                    include/system.h)
set(libpopt_public_headers include/popt.h)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/popt_config.h.in ${PROJECT_BINARY_DIR}/include/popt_config.h)

include_directories(${PROJECT_BINARY_DIR}/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

add_definitions("-DHAVE_CONFIG_H")

add_library(popt STATIC ${libpopt_src} ${libpopt_headers})


set_target_properties(popt PROPERTIES
                             PUBLIC_HEADER ${libpopt_public_headers}
                             COMPILE_FLAGS "${DISABLE_WARNINGS_FLAG}")

set(POPT_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include" CACHE INTERNAL "popt headers")
set(POPT_LIBRARY popt CACHE INTERNAL "popt library")

########### install files ###############

# install(TARGETS popt
#           ARCHIVE       DESTINATION lib
#           PUBLIC_HEADER DESTINATION include)

add_subdirectory(tests EXCLUDE_FROM_ALL)

maybe_disable_ssen(popt)
