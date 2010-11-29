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

macro(install_dsfmt)
  set(dsfmtVer "2.1")
  set(randSrcMD5 "b3a38dac7fd8996a70d02edc4432dd75")
  set(randSrcPath "${MILKYWAY_THIRDPARTY}")
  set(dsfmtTar "dSFMT-src-${dsfmtVer}.tar.gz")
  maybe_dl_check("dSFMT"
    "${dsfmtVer}"
    "${randSrcMD5}"
    "http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/SFMT/${dsfmtTar}"
    "${randSrcPath}"
    "${dsfmtTar}")

  set(dsfmt_mexp "19937") #Value it assumes if you don't specify it; stops warning
  add_definitions("-DDSFMT_MEXP=${dsfmt_mexp}")

  if(NOT MSVC)
    set(dsmft_flags "-O3 -DNDEBUG -finline-functions -fomit-frame-pointer -fno-strict-aliasing --param max-inline-insns-single=1800 -std=c99 ${SSE2_FLAGS}")
  else()
      set(dsmft_flags "/O2 -DNDEBUG ${SSE2_FLAGS}")
  endif()

  set(dsmft_src "${MILKYWAY_THIRDPARTY}/dSFMT-src-${dsfmtVer}/dSFMT.c")
  set(DSMFT_INCLUDE_DIR "${MILKYWAY_THIRDPARTY}/dSFMT-src-${dsfmtVer}/")
  set(dsmft_hdr "${DSMFT_INCLUDE_DIR}/dSFMT.h")

  include_directories("${MILKYWAY_THIRDPARTY}/dSFMT-src-${dsfmtVer}/")

  add_library(dsfmt STATIC ${dsmft_src})

  set_property(SOURCE ${dsmft_src}
               PROPERTY COMPILE_DEFINITIONS "DSFMT_MEXP=${dsfmt_mexp}")

  set_property(SOURCE ${dsmft_src}
               PROPERTY COMPILE_FLAGS "${dsmft_flags}")



  #file(INSTALL ${dsmft_hdr} DESTINATION include)
  #file(INSTALL dsfmt DESTINATION lib)
endmacro()

