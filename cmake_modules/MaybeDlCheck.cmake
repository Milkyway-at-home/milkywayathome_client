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

function(maybe_dl_check name md5Hash url tarName)
  if(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/${name}")
    message(STATUS "Downloading ${name}")
    # Download somewhere you don't see
    file(DOWNLOAD ${url} "${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_FILES_DIRECTORY}/${tarName}"
         TIMEOUT 300
         EXPECTED_MD5 ${md5Hash}
         LOG "Downloading ${name}"
         SHOW_PROGRESS)
     message(STATUS "Extracting ${name}")
     execute_process(COMMAND
                      "${CMAKE_COMMAND}" -E tar xvjf "${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_FILES_DIRECTORY}/${tarName}"
                      WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
  else()
    message(STATUS "Already have ${name}")
  endif()
endfunction()


