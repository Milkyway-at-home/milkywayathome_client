#
# Copyright 2011 Rensselaer Polytechnic Institute.
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


find_program(TEXINFO_COMMAND texi2pdf)

if(TEXINFO_COMMAND)
  add_custom_target(doc
                    COMMAND ${TEXINFO_COMMAND}
                    "${CMAKE_CURRENT_SOURCE_DIR}/milkyway_nbody.texinfo"
                    --output="${CMAKE_CURRENT_BINARY_DIR}/milkyway_nbody.pdf")
else()
  message(WARNING "Could not find texi2pdf. Documentation cannot be built")
endif()

