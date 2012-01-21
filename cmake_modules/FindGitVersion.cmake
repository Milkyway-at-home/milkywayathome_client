# Copyright 2011 Matthew Arsenault, Travis Desell, Dave Przybylo,
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


find_package(Git)

if(IS_DIRECTORY ".svn")
  execute_process(COMMAND svnversion ${repodir}
                    OUTPUT_VARIABLE SVN_REVISION
                    OUTPUT_STRIP_TRAILING_WHITESPACE
                    ERROR_STRIP_TRAILING_WHITESPACE)

  if(SVN_REVISION STREQUAL "exported")
    message(FATAL_ERROR "Did not find subversion version")
  endif()
  set(REPOSITORY_REVISION ${SVN_REVISION})

elseif(IS_DIRECTORY ".git")
  execute_process(COMMAND ${GIT_EXECUTABLE} rev-list HEAD
                    OUTPUT_VARIABLE GIT_REVISION
                    OUTPUT_STRIP_TRAILING_WHITESPACE
                    ERROR_STRIP_TRAILING_WHITESPACE)
  set(REPOSITORY_REVISION ${GIT_REVISION})
else()
  message(FATAL_ERROR "Directory does not contain .svn or .git")
endif()

string(TOUPPER HEADER_BASENAME uppername)

file(WRITE ${HEADER_BASENAME}.h
  "
   #ifndef _${uppername}_H_
   #define _${uppername}_H_

   #define SVN_REVISION \"${REPOSITORY_REVISION}\"

   #endif /* _${uppername}_H_ */
"
)


