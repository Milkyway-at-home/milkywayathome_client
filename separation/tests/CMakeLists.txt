#
# Copyright (C) 2011 Matthew Arsenault
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

include(MaybeDLCheck)

maybe_dl_check("stars" "9832cf1dd6e1e11083f9594fd4bf657b"
                       "http://milkyway.cs.rpi.edu/milkyway/download/stars.tar.bz2"
                       "stars.tar.bz2")

# TODO: Test different paths
add_test(NAME tests
           COMMAND ${LUA_INTERPRETER} "SeparationTests.lua"
                                      "${EXECUTABLE_OUTPUT_PATH}/${separation_base_bin_name}"
                                      "")


