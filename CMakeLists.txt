# Copyright 2010 Matthew Arsenault, Travis Desell, Dave Przybylo,
# Nathan Cole, Boleslaw Szymanski, Heidi Newberg, Carlos Varela, Malik
# Magdon-Ismail and Rensselaer Polytechnic Institute.
#
# This file is part of Milkway@Home.
#
# Milkyway@Home is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Milkyway@Home is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Milkyway@Home.  If not, see <http://www.gnu.org/licenses/>.
#

cmake_minimum_required(VERSION 2.6)
project(milkywayathome_client)
enable_testing()

if(NOT MSVC)
  # Use the separate generator for MSVC
  # FIXME: This doesn't really work. The tests that cmake runs before anything else will detect 64bitness
  option(BUILD_32 "32 bit build" OFF)
endif()

if(BUILD_32)
  if(APPLE)
    # TODO: If using this, it seems only the Apple gcc/llvm-gcc works
    set(CMAKE_OSX_ARCHITECTURES i386)
  endif()

  set(CMAKE_C_FLAGS "-m32 ${CMAKE_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "-m32 ${CMAKE_C_FLAGS}")
endif()



# Find custom cmake modules
set(MILKYWAYATHOME_CLIENT_CMAKE_MODULES "${milkywayathome_client_SOURCE_DIR}/cmake_modules")
list(APPEND CMAKE_MODULE_PATH ${MILKYWAYATHOME_CLIENT_CMAKE_MODULES})
set(PROJECT_TOOL_DIR "${milkywayathome_client_SOURCE_DIR}/tools")

include(CMakeDependentOption)
include(CheckLibraryExists)


include(CompilerID)

if(WIN32)
  set(DISABLE_WIN32_BULLSHIT_DEFINES "-D_CRT_SECURE_NO_WARNINGS=1 -DWIN32_LEAN_AND_MEAN=1 -DVC_EXTRALEAN=1 -DNOMINMAX=1")
  if(MSVC)
    string(REPLACE "-D" "/D " DISABLE_WIN32_BULLSHIT_DEFINES ${DISABLE_WIN32_BULLSHIT_DEFINES})
  endif()
  add_definitions(${DISABLE_WIN32_BULLSHIT_DEFINES})
elseif(APPLE)
  add_definitions(-DOS_OBJECT_USE_OBJC=0 -DOS_OBJECT_HAVE_OBJC_SUPPORT=0)
endif()
if(CXX_COMPILER_IS_SUN)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -features=zla")
endif()


include(MinGWPaths)
include(MWLink)
include(PrintInfo)
include(CheckCFlag)
include(MSVCFlags)
include(PreferStatic)
include(CMakeDependentOption)
include(CheckFunctionExists)
include(Profiling)
include(BOINCPaths)

#These must be after CheckCFlags is included

#TODO: Make sure to use Apple GCC for releases
if(HAVE_FLAG_STATIC_LIBGCC AND NOT APPLE)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -static-libgcc")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static-libgcc")
endif()

if(HAVE_FLAG_PTHREAD)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pthread")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
endif()


check_function_exists(stpcpy HAVE_STPCPY)
check_library_exists(rt shm_open "" HAVE_LIBRT)
check_library_exists(dl dlopen "" HAVE_LIBDL)


option(MSVC_USE_STATIC_CRT "Use static CRT" ON)
mark_as_advanced(MSVC_USE_STATIC_CRT)

if(MSVC)
  set_msvc_extra_link_flags()
  set_msvc_find_paths()
  if(MSVC_USE_STATIC_CRT)
    # http://msdn.microsoft.com/en-us/library/2kzt1wy3(v=VS.100).aspx
    set_msvc_mt()
  endif()
endif()

# Project build selection
option(SEPARATION "Build separation" ON)
option(NBODY "Build nbody" ON)

# Normal options
option(BOINC_APPLICATION
       "Run as BOINC application or independently" ON)

cmake_dependent_option(INSTALL_BOINC "make install should go to BOINC project directory (e.g. /var/lib/boinc-client/projects/milkyway.cs.rpi.edu/" ON
  "BOINC_APPLICATION" OFF)
mark_as_advanced(INSTALL_BOINC)



option(DOUBLEPREC
       "Build with double precision" ON)

option(NVIDIA_OPENCL "Using Nvidia OpenCL" OFF)
mark_as_advanced(NVIDIA_OPENCL)

option(AMD_OPENCL "Using AMD OpenCL" OFF)
mark_as_advanced(AMD_OPENCL)



option(DISABLE_DENORMALS
       "Disable denormal handling with SSE" OFF)
mark_as_advanced(DISABLE_DENORMALS)

option(BOINC_RELEASE_NAMES
       "Use BOINC platform binary names" OFF)
mark_as_advanced(BOINC_RELEASE_NAMES)

# -fno-errno-math

# Find libraries installed to this prefix. This way we can just dump
# dependencies here on windows where everything is painful.
set(MILKYWAY_ROOT ${PROJECT_BINARY_DIR})
list(INSERT CMAKE_PREFIX_PATH 0 ${MILKYWAY_ROOT})


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

#choose release build type by default
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
      "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel Profile."
      FORCE)
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(MW_ENABLE_DEBUG ON)
endif()

option(NBODY
        "Build the n-body code" ON)

option(NBODY_BLENDER_OUTPUT
        "Builds with the ability to output for use with Blender Visualizer." OFF)
mark_as_advanced(NBODY_BLENDER_OUTPUT)

cmake_dependent_option(NBODY_CRLIBM "Use crlibm for maths" ON
                         "DOUBLEPREC" ON)
mark_as_advanced(NBODY_CRLIBM)


option(MILKYWAY_IPHONE_APP "April fools?" OFF)
mark_as_advanced(MILKYWAY_IPHONE_APP)
if(MILKYWAY_IPHONE_APP)
  add_definitions("-DMILKYWAY_IPHONE_APP")
endif()


################################################################################



set(gcc_warnings "-Wall"
                 "-Wextra"
                 "-Wshadow"
                 "-Wredundant-decls"
                 "-Winline"
                 "-Wdisabled-optimization"
                 "-Wpointer-arith"
                 "-Wcast-align"
                 "-Wstrict-aliasing"
                 "-Wstrict-aliasing=3"
                 "-Wswitch-enum"
                 "-Wswitch-default"
                 "-Wfloat-equal"
                 "-Wwrite-strings"
                 "-Wcomment"

                 # pragma to disable -Wunknown-pragmas seems to not work
                 "-Wno-unknown-pragmas"
               #  "-Wsuggest-attribute=pure"
               #  "-Wsuggest-attribute=const"
               #  "-Wsuggest-attribute=noreturn"

               # "-Wstrict-overflow=5"
               # "-Wmissing-declarations"
               # "-Wunreachable-code"
               )
add_flag_if_supported("-Wno-psabi")
add_flag_if_supported("-std=gnu99")
add_flag_if_supported("-xc99")
# set(clang_analyzer "--analyze"
#                    -Xanalyzer -analyzer-eagerly-assume
#                    -Xanalyzer -analyzer-check-dead-stores
#                    -Xanalyzer -analyzer-check-llvm-conventions
#                    -Xanalyzer -analyzer-check-security-syntactic
#                    -Xanalyzer -analyzer-opt-analyze-nested-blocks
#                    -Xanalyzer -warn-uninit-values
#                    -Xanalyzer -warn-sizeof-pointer
#                    #-Xanalyzer -analyzer-check-buffer-overflows
#                    -Xanalyzer -analyzer-opt-analyze-headers
#                    -Xanalyzer -analyzer-experimental-checks
#                    #-Xanalyzer -analyzer-display-progress
#                    -Xanalyzer -analyzer-no-purge-dead)


set(gcc_extra_flags "-fno-common -funswitch-loops ")
set(gcc_only_extra_flags "-freg-struct-return ")



add_flag_if_supported("-maccumulate-outgoing-args")


# ICC does bad things by default
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${SAFE_MATH_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SAFE_MATH_FLAGS}")

#Set math assumptions that we can make. We don't ever use the error
#checking, or the other rounding modes.
set(extra_gcc_math_flags "-fno-rounding-math"
                         "-fno-math-errno")

if(NOT DOUBLEPREC)
  list(APPEND gcc_extra_flags "-fsingle-precision-constant")
endif()

if(CMAKE_COMPILER_IS_GNUCC OR C_COMPILER_IS_CLANG)
  foreach(flag ${gcc_extra_flags} ${gcc_warnings} ${extra_gcc_math_flags})
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${flag}")
  endforeach()
endif()

if(CMAKE_COMPILER_IS_GNUCC AND NOT C_COMPILER_IS_CLANG)
    foreach(flag ${gcc_only_extra_flags})
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${flag}")
  endforeach()
endif()

if(CMAKE_COMPILER_IS_GNUCC OR CXX_COMPILER_IS_CLANG)
  foreach(flag ${gcc_extra_flags} ${gcc_warnings} ${extra_gcc_math_flags})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${flag}")
  endforeach()
endif()


if(C_COMPILER_IS_CLANG)
  foreach(flag ${clang_warnings})
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${flag}")
  endforeach()

  if(MW_ENABLE_DEBUG)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ftrapv -fsanitize=undefined-trap -fsanitize-undefined-trap-on-error")
  endif()
endif()



# Clang currently missing many of the gcc warnings and options, so
# shut up the warnings about them
set(shut_up_unknown_flags "-Wno-unknown-warning-option -Qunused-arguments")

if(C_COMPILER_IS_CLANG)
  set(CMAKE_C_FLAGS "${shut_up_unknown_flags} ${CMAKE_C_FLAGS}")
endif()

if(CXX_COMPILER_IS_CLANG)
  set(CMAKE_CXX_FLAGS "${shut_up_unknown_flags} ${CMAKE_CXX_FLAGS}")
endif()

if(C_COMPILER_IS_PATHCC)
  # C99 inline apparently not supported and without this there are billions of warnings about it
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fgnu89-inline")
endif()

if(NOT SEPARATION AND NOT NBODY)
  message(FATAL_ERROR "Nothing will be built!")
endif()

if(WIN32)
  set(windows_version 0x0500) # Windows 2000
  add_definitions(-D_WINDOWS)
  add_definitions(-D_WIN32_WINDOWS=${windows_version})
  add_definitions(-D_WIN32_WINNT=${windows_version})
  add_definitions(-DWIN32_WINNT=${windows_version})
  add_definitions(-DWINVER=${windows_version})
endif()

# locally installed headers, and generated *_config.h
set(MILKYWAY_INSTALL_INCLUDE_DIR "${MILKYWAY_ROOT}/include" CACHE INTERNAL "Locally installed headers")



find_package(OpenCL)
find_package(OpenMP)
find_package(OpenGL)
find_package(Curses)
find_package(OpenSSL)

if(APPLE)
  cmake_dependent_option(SEPARATION_OPENCL "Build separation with OpenCL" OFF
                                            "OPENCL_FOUND" OFF)
else()
    cmake_dependent_option(SEPARATION_OPENCL "Build separation with OpenCL" ON
                                             "OPENCL_FOUND" OFF)
endif()


cmake_dependent_option(NBODY_OPENCL "Build nbody with OpenCL" OFF
                                     "OPENCL_FOUND" OFF)

if(NBODY_OPENCL OR SEPARATION_OPENCL)
  set(MILKYWAY_OPENCL ON CACHE INTERNAL "Need to build milkyway_cl")
endif()


cmake_dependent_option(NBODY_OPENMP "Use OpenMP for nbody" ON
                                    "OPENMP_FOUND" OFF)

cmake_dependent_option(NBODY_GL "Build nbody visualizer" OFF
                                "OPENGL_FOUND;OPENGL_GLU_FOUND" OFF)

# Required on Windows, optional on others
# OS X GL 3.2 requires Lion anyway so there's no sense in using it by default
if(APPLE)
  set(USE_GL3W_DEFAULT OFF)
else()
  set(USE_GL3W_DEFAULT ON)
endif()
cmake_dependent_option(USE_GL3W "Use GL3W for OpenGL 3 loading" ${USE_GL3W_DEFAULT}
                       "NOT WIN32" ON)


print_separator()
print_build_info()
print_build_options()
print_separator()


include(SpecialLibs) # Do this after options


add_subdirectory(tools)
add_subdirectory(popt)
add_subdirectory(lua)

add_subdirectory(dSFMT)
add_definitions("${DSFMT_FLAGS}")

if(SEPARATION_OPENCL)
  add_subdirectory(libelf)
endif()

add_subdirectory(openpa)

if(NBODY_GL)
  if(NOT APPLE)
    set(OPENGL3_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/third_party/include)
    include_directories(${OPENGL3_INCLUDE_DIR})
  endif()

  if(USE_GL3W)
    # Header and source generated to load all OpenGL3 functions via
    # function pointer to work around Windows
    set(GL3W_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/third_party/include)
    include_directories(${GL3W_INCLUDE_DIR})
    add_library(gl3w STATIC ${CMAKE_CURRENT_SOURCE_DIR}/third_party/src/gl3w.c)
    set(GL3W_LIBRARY gl3w)
  endif()
  
  if(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/glfw/CMakeLists.txt)
    message(FATAL_ERROR "GLFW submodule is missing. Did you forget to git submodule init; git submodule update?")
  endif()

  set(BUILD_SHARED_LIBS OFF)
  set(GLFW_BUILD_EXAMPLES OFF)
  set(GLFW_BUILD_TESTS OFF)
  add_subdirectory(glfw)
  set(GLFW_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/glfw/include")
  set(GLFW_LIBRARY glfw)
  # TODO: We require 10.7 on OS X to use OpenGL 3.2
  # Also The 10.7 Cocoa headers also require clang or Apple llvm-gcc If we're using
  # macports gcc we can build glfw as an external project with the
  # other compiler


  if(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/glm/CMakeLists.txt)
    message(FATAL_ERROR "GLM submodule is missing. Did you forget to git submodule init; git submodule update?")
  endif()
  set(GLM_PACKAGE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/glm)
  list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/glm/util)
  find_package(GLM)
endif()


if(BOINC_APPLICATION)
  if(EXISTS ${PROJECT_SOURCE_DIR}/boinc/boinc/api/boinc_api.cpp)  # The submodule is there
    message(STATUS "Using BOINC from submodule")
    add_subdirectory(boinc)
  else()
    message(STATUS "Trying to find BOINC")
    find_package_static(BOINC REQUIRED)
  endif()

  include_directories(${BOINC_INCLUDE_DIR})
endif()

include_directories(${POPT_INCLUDE_DIR})
include_directories(${LUA_INCLUDE_DIR})
add_subdirectory(milkyway)

if(SEPARATION)
  add_subdirectory(separation)
endif()

if(NBODY)
  if(NBODY_CRLIBM)
    set(CRLIBM_QUIET_BUILD TRUE)
    add_subdirectory(crlibm)
    set(CRLIBM_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/crlibm)
    set(CRLIBM_LIBRARY crlibm)
  endif()
  add_subdirectory(nbody)
endif()

if(${CMAKE_VERSION} VERSION_GREATER "2.8.0" OR ${CMAKE_VERSION} VERSION_EQUAL "2.8.0")
  find_package(Git)
endif()

find_program(XZ_COMMAND xz)

if(GIT_FOUND AND XZ_COMMAND)
  macro(make_release_tarball targetname releasename)
    add_custom_target(${targetname} COMMAND ${GIT_EXECUTABLE} archive --prefix="${releasename}/" --format=tar HEAD | ${XZ_COMMAND} -9 > "${releasename}.tar.xz")
  endmacro()

  make_release_tarball("nbody_release" "milkyway_nbody_${nbody_VERSION}")
  make_release_tarball("separation_release" "milkyway_separation_${separation_VERSION}")
endif()

if(0)
  add_subdirectory(LModL EXCLUDE_FROM_ALL)
  add_subdirectory(simulation EXCLUDE_FROM_ALL)
endif()

add_custom_command(OUTPUT "${PROJECT_BINARY_DIR}/tools/app_info.xml"
                    WORKING_DIRECTORY "${PROJECT_BINARY_DIR}/tools"
                    COMMAND ${LUA_INTERPRETER} "${PROJECT_SOURCE_DIR}/tools/GenerateAppInfo.lua"
                                               "$SEPARATION_BIN_NAME"
                                               "$NBODY_BIN_NAME"
                                               "$NBODY_GL_BIN_NAME")



add_custom_target(app_info DEPENDS "${PROJECT_BINARY_DIR}/tools/app_info.xml")

if(INSTALL_BOINC)
  # FIXME:
  # Need to rename and install libraries, e.g. libgomp
  # Need to chown to boinc user
  install(FILES "${PROJECT_BINARY_DIR}/tools/app_info.xml" DESTINATION ${MILKYWAY_PROJECT_DIRECTORY})
endif()

option(TAO "Build TAO search package" OFF)
if(TAO)
  add_subdirectory(mwtao)
endif()



# Workaround for MSVS10 to avoid the Dialog Hell
# FIXME: This could be removed with future version of CMake.
if(MSVC_VERSION EQUAL 1600)
  set(MILKYWAYATHOME_CLIENT_SLN_FILENAME "${CMAKE_CURRENT_BINARY_DIR}/milkywayathome_client.sln")
  if(EXISTS "${MILKYWAYATHOME_CLIENT_SLN_FILENAME}")
    file(APPEND "${MILKYWAYATHOME_CLIENT_SLN_FILENAME}" "\n# This should be regenerated!\n")
  endif()
endif()
