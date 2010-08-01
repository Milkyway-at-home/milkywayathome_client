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

function(maybe_dl_check name version md5Hash url srcPath srcTar)
  if(NOT EXISTS "${srcPath}/${srcTar}")
    message(STATUS "Downloading ${name}")
    file(DOWNLOAD ${url} "${srcPath}/${srcTar}"
         TIMEOUT 60
         EXPECTED_MD5 ${md5Hash}
         LOG "Downloading ${name}"
         SHOW_PROGRESS)
  else()
    message(STATUS "Already have ${name}")
  endif()

  message(STATUS "Extracting ${name} source")
  execute_process(
    COMMAND "${CMAKE_COMMAND}" -E tar xzf "${srcPath}/${srcTar}" "${name}"
    WORKING_DIRECTORY "${srcPath}")
endfunction()


