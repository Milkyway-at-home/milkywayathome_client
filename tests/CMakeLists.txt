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

cmake_minimum_required (VERSION 2.6)

set(test_app ${PROJECT_BINARY_DIR}/bin/milkyway_client)


#TODO: Fix milkyway_client to use filename argument instead of stupid preprocessor / hardcoded filename
set(small_test_args1 -np 20 -p 0.571713 12.312119 -3.305187 148.010257 22.453902 0.42035 -0.468858 0.760579 -1.361644 177.884238 23.882892 1.210639 -1.611974 8.534378 -1.361644 177.884238 10.882892 1.210639 -1.611974 8.534378 )

set(small_test_args2 -np 8 -p 0.34217373320392042 25.9517910846623 -2.1709414738826602 38.272511356953906 30.225190442596112 2.2149060013372885 0.32316169064291655 2.7740244716285285)

set(small_test_args3 -np 8 -p 0.40587961154742185 17.529961843393409 -1.8575145272144837 29.360893891378243 31.228263575178566 -1.551741065334 .064096152599308373 2.55428209912781)

foreach(i RANGE 1 3)
  add_test(small_test_${i} ${test_app} ${small_test_args${i}})
endforeach()


