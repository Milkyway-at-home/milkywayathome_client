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

add_custom_command(OUTPUT "stars.tar.bz2"
                     COMMAND ${CMAKE_COMMAND}
                       -DCMAKE_MODULE_PATH=${CMAKE_MODULE_PATH}
                       -DCMAKE_FILES_DIRECTORY=${CMAKE_FILES_DIRECTORY}
                       -DCMAKE_CURRENT_BINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}
                       -DNAME="stars"
                       -DTARNAME="stars.tar.bz2"
                       -DMD5_FILE="9832cf1dd6e1e11083f9594fd4bf657b"
                       -DURL="http://milkyway.cs.rpi.edu/milkyway/download/stars.tar.bz2"
                       -P "${MILKYWAYATHOME_CLIENT_CMAKE_MODULES}/DownloadStars.cmake")


# TODO: Test different paths, such as no intrinsics etc.
add_test(NAME tests
           WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/tests"
           COMMAND $<TARGET_FILE:lua> "${PROJECT_SOURCE_DIR}/tests/SeparationTests.lua"
                                       $<TARGET_FILE:milkyway_separation>
                                       "${PROJECT_SOURCE_DIR}/tests"
                                       "")

add_custom_target(test_data DEPENDS "stars.tar.bz2")
# FIXME: How to add dependency on tests of test_data?

