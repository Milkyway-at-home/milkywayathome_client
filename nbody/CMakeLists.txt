# Copyright 2010-2011 Matthew Arsenault, Travis Desell, Dave Przybylo,
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
project(milkyway_nbody)

set(nbody_VERSION_MAJOR 0
          CACHE INTERNAL "N-body version number")
set(nbody_VERSION_MINOR 51
          CACHE INTERNAL "N-body version number")
set(nbody_VERSION "${nbody_VERSION_MAJOR}.${nbody_VERSION_MINOR}"
          CACHE INTERNAL "N-body version number")

include(MaybeDlCheck)
include(StaticBOINC)
include(SpecialLibs)
include(PrintInfo)
include(CheckSSE2)
include(BOINCPaths)
include(CMakeDependentOption)

cmake_dependent_option(NBODY_STATIC "Build Nbody as fully static binary" ON
                                      "NOT NBODY_OPENCL" OFF)
mark_as_advanced(NBODY_STATIC)
maybe_static(${NBODY_STATIC})

if(OPENGL_FOUND AND OPENGL_GLU_FOUND AND GLUT_FOUND)
  set(NBODY_GL ON CACHE INTERNAL "Build N-body graphics")
endif()

include_directories(${MILKYWAY_INCLUDE_DIR})
include_directories(${MILKYWAY_INSTALL_INCLUDE_DIR})
include_directories(${DSFMT_INCLUDE_DIR})


set(NBODY_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include")
include_directories(${NBODY_INCLUDE_DIR})


if(OPENMP_FOUND AND NBODY_OPENMP)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
elseif(NOT OPENMP_FOUND AND NBODY_OPENMP)
  message(WARNING "Did not find OpenMP support, but was enabled. Continuing without OpenMP")
endif()

if(BOINC_RELEASE)
  set(NBODY_CRLIBM ON)
endif()

if(NOT DOUBLEPREC AND NBODY_CRLIBM)
  message(WARNING "Attempting to use single precision with crlibm. Disabling use of crlibm")
  set(NBODY_CRLIBM OFF)
endif()

if(NBODY_OPENCL)
  if(NOT OPENCL_FOUND)
    message(WARNING "N-body OpenCL was requested, but OpenCL not found. Disabling OpenCL.")
    set(NBODY_OPENCL OFF)
  else()
    include_directories(${OPENCL_DIRECTORIES})
  endif()
endif()

# Handle downloading, extracting and building dSFMT

print_separator()
message("   Building:          ${PROJECT_NAME}")
message("   BOINC Application: ${BOINC_APPLICATION}")
message("   Double precision:  ${DOUBLEPREC}")
message("   Enable crlibm:     ${NBODY_CRLIBM}")
message("  ----")
print_libs()
print_separator()

if(BOINC_APPLICATION AND NOT BOINC_LIBRARIES)
  message(FATAL_ERROR "Trying to build as a BOINC application, but the BOINC library wasn't found")
endif()



#CHECKME: We don't actually need special number math
#except for assignment in a couple places. Not sure how useful this is
#anyway.
#if(NOT ENABLE_FDLIBM)
#  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ffinite-math-only")
#endif()

set(NBODY_SRC_DIR "${PROJECT_SOURCE_DIR}/src/")
set(NBODY_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include/")
set(nbody_lib_src ${NBODY_SRC_DIR}/nbody_chisq.c
                  ${NBODY_SRC_DIR}/nbody_grav.c
                  ${NBODY_SRC_DIR}/nbody_io.c
                  ${NBODY_SRC_DIR}/nbody_types.c
                  ${NBODY_SRC_DIR}/nbody_tree.c
                  ${NBODY_SRC_DIR}/nbody_orbit_integrator.c
                  ${NBODY_SRC_DIR}/nbody_potential.c
                  ${NBODY_SRC_DIR}/nbody_simulation.c
                  ${NBODY_SRC_DIR}/nbody_check_params.c
                  ${NBODY_SRC_DIR}/nbody_plummer.c
                  ${NBODY_SRC_DIR}/nbody_step.c
                  ${NBODY_SRC_DIR}/nbody_show.c
                  ${NBODY_SRC_DIR}/nbody_checkpoint.c
                  ${NBODY_SRC_DIR}/nbody_defaults.c
                  ${NBODY_SRC_DIR}/nbody_coordinates.c
                  ${NBODY_SRC_DIR}/nbody_shmem.c
                  ${NBODY_SRC_DIR}/nbody_util.c)


set(nbody_lib_headers ${NBODY_INCLUDE_DIR}/nbody_chisq.h
                      ${NBODY_INCLUDE_DIR}/nbody_grav.h
                      ${NBODY_INCLUDE_DIR}/nbody_config.h.in
                      ${NBODY_INCLUDE_DIR}/nbody_io.h
                      ${NBODY_INCLUDE_DIR}/nbody_tree.h
                      ${NBODY_INCLUDE_DIR}/nbody_orbit_integrator.h
                      ${NBODY_INCLUDE_DIR}/nbody_potential.h
                      ${NBODY_INCLUDE_DIR}/nbody_check_params.h
                      ${NBODY_INCLUDE_DIR}/nbody_plummer.h
                      ${NBODY_INCLUDE_DIR}/nbody.h
                      ${NBODY_INCLUDE_DIR}/nbody_show.h
                      ${NBODY_INCLUDE_DIR}/nbody_priv.h
                      ${NBODY_INCLUDE_DIR}/nbody_types.h
                      ${NBODY_INCLUDE_DIR}/nbody_checkpoint.h
                      ${NBODY_INCLUDE_DIR}/nbody_defaults.h
                      ${NBODY_INCLUDE_DIR}/nbody_coordinates.h
                      ${NBODY_INCLUDE_DIR}/nbody_shmem.h
                      ${NBODY_INCLUDE_DIR}/nbody_util.h
                      ${NBODY_INCLUDE_DIR}/nbody_graphics.h)


set(nbody_lua_src ${NBODY_SRC_DIR}/nbody_lua.c
                  ${NBODY_SRC_DIR}/nbody_lua_models.c

                  ${NBODY_SRC_DIR}/nbody_lua_types/nbody_lua_types.c
                  ${NBODY_SRC_DIR}/nbody_lua_types/nbody_lua_type_marshal.c

                  ${NBODY_SRC_DIR}/nbody_lua_types/nbody_lua_nbodystate.c

                  ${NBODY_SRC_DIR}/nbody_lua_types/nbody_lua_nbodyctx.c
                  ${NBODY_SRC_DIR}/nbody_lua_types/nbody_lua_body.c
                  ${NBODY_SRC_DIR}/nbody_lua_types/nbody_lua_halo.c
                  ${NBODY_SRC_DIR}/nbody_lua_types/nbody_lua_disk.c
                  ${NBODY_SRC_DIR}/nbody_lua_types/nbody_lua_spherical.c
                  ${NBODY_SRC_DIR}/nbody_lua_types/nbody_lua_potential.c
                  ${NBODY_SRC_DIR}/nbody_lua_types/nbody_lua_histogram_params.c)



set(nbody_lua_headers ${NBODY_INCLUDE_DIR}/nbody_lua.h
                      ${NBODY_INCLUDE_DIR}/nbody_lua_models.h
                      ${NBODY_INCLUDE_DIR}/nbody_lua_types.h
                      ${NBODY_INCLUDE_DIR}/nbody_lua_type_marshal.h

                      ${NBODY_INCLUDE_DIR}/nbody_lua_nbodystate.h

                      ${NBODY_INCLUDE_DIR}/nbody_lua_nbodyctx.h
                      ${NBODY_INCLUDE_DIR}/nbody_lua_body.h
                      ${NBODY_INCLUDE_DIR}/nbody_lua_halo.h
                      ${NBODY_INCLUDE_DIR}/nbody_lua_disk.h
                      ${NBODY_INCLUDE_DIR}/nbody_lua_spherical.h
                      ${NBODY_INCLUDE_DIR}/nbody_lua_potential.h
                      ${NBODY_INCLUDE_DIR}/nbody_lua_histogram_params.h)


add_library(nbody_lua STATIC ${nbody_lua_src} ${nbody_lua_headers})
set_target_properties(nbody_lua PROPERTIES
                        COMPILE_DEFINITIONS "NBODY_LUA_TYPES_COMPILATION")



#TODO: The tests dSFMT has

# Require SSE and SSE2 (needed for doubles) This avoids a host of
# issues with x87, ranging from different intermediate precisions on
# different places, and removes many of the places where different
# compiler optimizations can change the precision.
#
# Some of these issues are summarized here: There are others to be
# aware of.
# http://www.yosefk.com/blog/consistency-how-to-defeat-the-purpose-of-ieee-floating-point.html
#
# Also as a nice advantage, SSE for a quick benchmark I ran is ~20%
# faster than using the antique x87 junk. There is also a later
# opportunity to do hand vectorization for an even bigger improvement.
#
# The GPU floating point stuff is
# "similar" to SSE, but not SSE.  OpenCL demands a stricter
# conformance to standards for doubles, so hopefully when we start
# using GPUs the results will be the same as using SSE/SSE2 and
# doubles. No promises however.
#
# Summary: we abandon old, old processors for the sake of convenience
# and consistency, and a nice speed boost on modern x86/x86_64
# hardware.
#

set(CMAKE_C_FLAGS "-std=gnu99 ${SSE2_FLAGS} ${CMAKE_C_FLAGS}")

set(nbody_link_libs )

if(NBODY_OPENCL)
  # set(c_kernel_dir "${CMAKE_CURRENT_SOURCE_DIR}/ckernels")
  # list(APPEND nbody_lib_src "${PROJECT_SOURCE_DIR}/src/gravmap_opencl.c")

  # list(APPEND nbody_lib_src "${c_kernel_dir}/cl_gravmap.c")
  # set_source_files_properties("${c_kernel_dir}/cl_gravmap.c"
  #                             PROPERTIES GENERATED TRUE)

  # list(APPEND nbody_link_libs ${OPENCL_LIBRARIES})
endif()

add_library(nbody_simulation STATIC ${nbody_lib_src} ${nbody_lib_headers})


if(DOWNLOAD_CRLIBM)
  add_dependencies(nbody_simulation crlibm_build)
  set(CRLIBM_INCLUDE_DIR "${MILKYWAY_ROOT}/include"})
  set(CRLIBM_LIBRARY "${LIBRARY_OUTPUT_PATH}/libcrlibm${CMAKE_STATIC_LIBRARY_SUFFIX}")
else()
  if(BOINC_RELEASE AND NOT CRLIBM_FOUND)
    message(FATAL "Trying to release without crlibm")
  endif()
endif()

set(nbody_exe_link_libs nbody_simulation
                        nbody_lua
                        milkyway_lua
                        nbody_simulation
                        milkyway

                        ${POPT_LIBRARY}
                        ${LUA_LIBRARIES}
                        m
                        ${OS_SPECIFIC_LIBS}
                        ${nbody_link_libs})

if(BOINC_APPLICATION)
  list(INSERT nbody_exe_link_libs 5 ${BOINC_LIBRARIES})
endif()


if(NBODY_CRLIBM)
  include_directories(${CRLIBM_INCLUDE_DIR})
  list(APPEND nbody_exe_link_libs ${CRLIBM_LIBRARY})
endif()

if(NOT HAVE_SSE2 AND SYSTEM_IS_x86)
  message(FATAL_ERROR "SSE2 is required for x86 systems")
endif()

add_executable(milkyway_nbody ${NBODY_SRC_DIR}/main.c)
milkyway_link(milkyway_nbody ${BOINC_APPLICATION} ${NBODY_STATIC} "${nbody_exe_link_libs}")

if(NBODY_GL AND WIN32 AND NOT BOINC_APPLICATION)
  message(FATAL "N-body visualizer unsupported on Windows without BOINC")
endif()
if(NBODY_GL AND BOINC_APPLICATION AND NOT BOINC_GRAPHICS_FOUND)
  message(FATAL "BOINC graphics library not found")
endif()

if(NBODY_GL)
  include_directories(${OPENGL_INCLUDE_DIR})
  include_directories(${GLUT_INCLUDE_PATH})

  set(nbody_screensaver_src ${NBODY_SRC_DIR}/nbody_gl.c
                            ${NBODY_SRC_DIR}/nbody_gl_util.c
                            ${NBODY_SRC_DIR}/graphics_main.c)
  set(nbody_screensaver_headers ${NBODY_INCLUDE_DIR}/nbody_gl.h
                                ${NBODY_INCLUDE_DIR}/nbody_gl_util.h
                                ${NBODY_INCLUDE_DIR}/nbody_graphics.h)

  set(vis_link_libs milkyway ${POPT_LIBRARY} ${GLUT_LIBRARY} ${OPENGL_LIBRARIES})
  if(BOINC_APPLICATION)
    list(APPEND vis_link_libs ${BOINC_GRAPHICS_LIBRARY} ${BOINC_LIBRARIES})
  endif()
  list(APPEND vis_link_libs ${OS_SPECIFIC_LIBS} m)

  add_executable(milkyway_nbody_graphics ${nbody_screensaver_src} ${nbody_screensaver_headers})
  target_link_libraries(milkyway_nbody_graphics ${vis_link_libs})
endif()


if(BOINC_RELEASE_NAMES)
  if(NBODY_OPENMP)
    set(NBODY_PLAN "mt")
  else()
    if(HAVE_SSE2)
      set(NBODY_PLAN "sse2")
    endif()
  endif()

  get_boinc_bin_name("milkyway_nbody" "${nbody_VERSION}" "${NBODY_PLAN}")
  set_target_properties(milkyway_nbody
                          PROPERTIES
						    OUTPUT_NAME ${BOINC_BIN_NAME})

  if(NBODY_GL)
    get_boinc_bin_name("milkyway_nbody_graphics" "${nbody_VERSION}" "")
    set_target_properties(milkyway_nbody_graphics
                            PROPERTIES
                              OUTPUT_NAME ${BOINC_BIN_NAME})
  endif()
endif()

install(TARGETS milkyway_nbody
          RUNTIME DESTINATION ${MILKYWAY_PROJECT_DIRECTORY})

if(NBODY_GL)
  install(TARGETS milkyway_nbody_graphics
          RUNTIME DESTINATION ${MILKYWAY_PROJECT_DIRECTORY})
  get_target_property(NBODY_GRAPHICS_BIN_NAME milkyway_nbody_graphics OUTPUT_NAME)
  set(NBODY_GRAPHICS_BIN_NAME ${NBODY_GRAPHICS_BIN_NAME} PARENT_SCOPE)
endif()


configure_file (
  "${NBODY_INCLUDE_DIR}/nbody_config.h.in"
  "${PROJECT_INCLUDE_DIR}/nbody_config.h")


find_package(OpenSSL)
message("   OpenSSL Libraries  ${OPENSSL_LIBRARIES}")

if(NOT OPENSSL_FOUND)
  message(WARNING "OpenSSL not found, not using N-body tests")
else()
  add_subdirectory(tests EXCLUDE_FROM_ALL)
endif()


