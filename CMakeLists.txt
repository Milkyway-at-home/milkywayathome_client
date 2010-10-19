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
#
# You should have received a copy of the GNU General Public License
# along with Milkyway@Home.  If not, see <http://www.gnu.org/licenses/>.
#

cmake_minimum_required (VERSION 2.6)
project(milkywayathome)

enable_testing()

# Find custom cmake modules
set(CMAKE_MODULE_PATH "${milkywayathome_SOURCE_DIR}/cmake_modules")
set(PROJECT_TOOL_DIR "${milkywayathome_SOURCE_DIR}/tools")

include(MinGWPaths)
include(StaticBOINC)
include(PrintInfo)
include(DownloadCrlibm)
set_mingw_paths()

# Project build selection
option(SEPARATION "Build separation" ON)
option(NBODY "Build nbody" ON)


# Normal options
option(BOINC_APPLICATION
       "Run as BOINC application or independently" ON)

option(BOINC_RELEASE
       "Build with all options for making an actual release" ON)

option(DOUBLEPREC
       "Build with double precision" ON)


# 'Secret' options
option(DISABLE_DENORMALS
       "Disable denormal handling with SSE" OFF)
mark_as_advanced(DISABLE_DENORMALS)

option(STATIC_EVERYTHING
       "Prefer static libraries" ON)
mark_as_advanced(STATIC_EVERYTHING)

option(BOINC_RELEASE_NAME
       "Use BOINC platform binary names" OFF)
mark_as_advanced(BOINC_RELEASE_NAME)

option(DOWNLOAD_CRLIBM
       "Download crlibm and install to thirdparty" OFF)
mark_as_advanced(DOWNLOAD_CRLIBM)


# Find libraries installed to this prefix. This way we can just dump
# dependencies here on windows where everything is painful.
set(MILKYWAY_ROOT ${PROJECT_SOURCE_DIR})
list(INSERT CMAKE_PREFIX_PATH 0 ${MILKYWAY_ROOT})

# Stuff which the build will download and use from there.
set(MILKYWAY_THIRDPARTY "${MILKYWAY_ROOT}/thirdparty")

if(DOWNLOAD_CRLIBM)
  download_build_crlibm()
endif()


# Set output to be in this root
set(EXECUTABLE_OUTPUT_PATH ${MILKYWAY_ROOT}/bin)
set(LIBRARY_OUTPUT_PATH ${MILKYWAY_ROOT}/lib)
set(PROJECT_INCLUDE_DIR ${MILKYWAY_ROOT}/include)

include_directories(${PROJECT_INCLUDE_DIR})
link_directories(${LIBRARY_OUTPUT_PATH})

# Avoid issue I'm too lazy to fix where using the config header, some
# files were getting single and others double
if(DOUBLEPREC)
  add_definitions("-DDOUBLEPREC=1")
endif()

if(DISABLE_DENORMALS)
  add_definitions("-DDISABLE_DENORMALS=1")
endif()

if(BOINC_RELEASE)
  set(STATIC_EVERYTHING ON)
  set(BOINC_RELEASE_NAME ON)
  set(BOINC_APPLICATION ON)
  set(DOUBLEPREC ON)
endif()



#choose release build type by default
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
      "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
endif(NOT CMAKE_BUILD_TYPE)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(MW_ENABLE_DEBUG ON)
endif()


maybe_static()

if(UNIX AND NOT APPLE)
  set(NEED_POPT_DEPS 0)
else()
  set(NEED_POPT_DEPS 1)
endif()

#FIXME: Somehow test for and require the Apple gcc. The macports gccs
#seem to not work when you try to link against frameworks,
#CoreFoundation etc. There's probably a better way to check and handle
#this.
if(APPLE AND BOINC_APPLICATION)
  include(CheckAppleGCC)
  check_apple_gcc()
  if(NOT HAVE_APPLE_GCC)

    message(WARNING
            "You should specify a compiler by setting the CC environment variable or "
            "with -DCMAKE_C_COMPILER=/path/to/compiler. This should be /usr/bin/gcc usually."
            "Latest clang also works")
    message(FATAL_ERROR "Attempting to use a non-Apple GCC. An Apple compiler is required on OS X.")
  endif()

  if(BOINC_RELEASE)
    # We don't want to bother compiling twice until we need the universal binary for releases
    #set(CMAKE_OSX_ARCHITECTURES i386;x86_64)
    #set(CMAKE_OSX_ARCHITECTURES x86_64)
  endif()
endif()

#Some kind of workaround for a popular cmake bug where make test
#doesn't behave as anyone would ever expect.
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND})

option(SEPARATION_OPENCL
        "Build the OpenCL separation" OFF)

option(NBODY_OPENCL
        "Build the OpenCL nbody" OFF)

option(MILKYWAY_OPENCL
        "Build the OpenCL nbody" OFF)
mark_as_advanced(MILKYWAY_OPENCL)
#FIXME: Be able to build one with / without openl
if(NBODY_OPENCL OR SEPARATION_OPENCL)
  set(MILKYWAY_OPENCL ON)
  add_definitions("-DMILKYWAY_OPENCL=1")
endif()

option(NBODY
    "Build the n-body code" ON)

option(BOINC_APP_GRAPHICS
       "BOINC app graphics" OFF)


set(gcc_extra_warnings "-Wshadow -Wpointer-arith -Wcast-align -common-Wwrite-strings -compare -fno-common -Wstrict-aliasing -Wstrict-aliasing=3 -Wswitch-enum -Wswitch-default -Wmissing-declarations -Wfloat-equal -Wcomment")

#TODO: Other compilers
if(CMAKE_COMPILER_IS_GNUCC) # I think this could maybe be IS_UNIX
  set(CMAKE_C_FLAGS "-Wall -Wextra ${CMAKE_C_FLAGS} -funroll-loops -maccumulate-outgoing-args")

  set(CMAKE_C_FLAGS_DEBUG "-Wall -Wextra ${gcc_extra_warnings} -O0 -g ${CMAKE_C_FLAGS_DEBUG}")
endif(CMAKE_COMPILER_IS_GNUCC)


#Set math assumptions that we can make. We don't ever use the error
#checking, or the other rounding modes.
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mhard-float -fno-rounding-math -fno-math-errno")


if(NOT DOUBLEPREC)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsingle-precision-constant")
endif()


print_separator()
print_build_info()
print_build_options()
print_separator()


if(NOT SEPARATION AND NOT NBODY)
  message(FATAL_ERROR "Nothing will be built!")
endif()


add_subdirectory(milkyway)
include_directories("${PROJECT_SOURCE_DIR}/milkyway/include")

if(SEPARATION)
  add_subdirectory(separation)
endif()

if(NBODY)
  add_subdirectory(nbody)
endif()


