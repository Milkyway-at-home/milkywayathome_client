#
# Copyright (C) 2011 Matthew Arsenault
# Copyright (C) 2011 Rensselaer Polytechnic Institute.
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

set(MILKYWAY_PROJECT_HOST "milkyway.cs.rpi.edu")

if(WIN32)
  set(BOINC_PROJECT_DIRECTORY "$ENV{SYSTEMROOT}/ProgramData/BOINC Data/projects/")
elseif(APPLE)
  set(BOINC_PROJECT_DIRECTORY "/Library/Application Support/BOINC Data/projects/")
else()
  set(BOINC_PROJECT_DIRECTORY "/var/lib/boinc/projects/")
endif()

set(BOINC_USER "boinc")
set(BOINC_GROUP "boinc")

set(MILKYWAY_PROJECT_DIRECTORY "${BOINC_PROJECT_DIRECTORY}/${MILKYWAY_PROJECT_HOST}_milkyway")




