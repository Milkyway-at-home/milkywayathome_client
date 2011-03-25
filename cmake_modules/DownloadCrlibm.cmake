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

include(MaybeDlCheck)
include(CheckSSE2)

macro(download_build_crlibm)
  set(crlibmVer "1.0beta4")
  set(crlibmMD5 "8ecabd55d4a2d08030eb770ae2a5670a")
  set(crlibmTar "crlibm-${crlibmVer}.tar.gz")
  set(crlibmSrcPath "${MILKYWAY_THIRDPARTY}")

  maybe_dl_check("crlibm"
                 "${crlibmVer}"
                 "${crlibmMD5}"
                 "http://lipforge.ens-lyon.fr/frs/download.php/162/${crlibmTar}"
                 "${MILKYWAY_THIRDPARTY}"
                 "${crlibmTar}")

 #TODO: Set march, mtune etc. from user flags, but avoid some others
  set(crlibm_flags "-std=c99 -O3 ${SSE2_FLAGS}")

  if(HAVE_SSE2) # x86/x86_64 systems
    set(crlibm_configure_sse2 "--enable-sse2")
    # doesn't apply and breaks the build on PPC
  endif()

  set(crlibm_file "${LIBRARY_OUTPUT_PATH}/libcrlibm${CMAKE_STATIC_LIBRARY_SUFFIX}")
  set(crlibm_header "${PROJECT_INCLUDE_DIR}/crlibm.h")
  set(scs_file "${LIBRARY_OUTPUT_PATH}/libscs${CMAKE_STATIC_LIBRARY_SUFFIX}")
  add_custom_command(
         OUTPUT "${crlibm_file}" "${crlibm_header}"
         COMMAND "${crlibmSrcPath}/crlibm-${crlibmVer}/configure"
                 "${crlibm_configure_sse2}"
                 "--prefix=${MILKYWAY_ROOT}"
                 "CC=${CMAKE_C_COMPILER}"
                 "MAKE=${CMAKE_BUILD_TOOL}"
                 "CFLAGS=${crlibm_flags}"
         COMMAND "${CMAKE_BUILD_TOOL}"
         COMMAND "${CMAKE_BUILD_TOOL}" "install"
         WORKING_DIRECTORY "${crlibmSrcPath}/crlibm-${crlibmVer}"
         COMMENT "Building crlibm")

  add_custom_target(crlibm_build DEPENDS "${crlibm_file}")
endmacro()

