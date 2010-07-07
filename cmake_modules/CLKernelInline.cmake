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

function(inline_kernel name file c_kernel_dir)
  file(READ "${dir}/${file}" str)

  set(name "cl_${name}")

  #Remove the comments. Otherwise, C++ comments will comment out the
  #rest of the program unless we do more work.
  string(REGEX REPLACE
             "(/\\*([^*]|[\r\n]|(\\*+([^*/]|[\r\n])))*\\*+/)|(//[^\r\n]*)"
             "" # No comment
             str_stripped

             #Quote this to prevent ;'s being interpreted as list separators
             "${str}")

  #Escape special characters
  #TODO: Other things that need escaping
  string(REGEX REPLACE
             "(\n)|(\r\n)"
             "\\\\n"
             str_escaped
             "${str_stripped}")

  string(REGEX REPLACE
             "(\")"
             "\\\\\""
             str_escaped
             "${str_escaped}")

  set(cfile "\n#include \"${name}.h\"\n
           \nconst char* ${name}_src = \"${str_escaped}\";
           \n\n")

  string(TOUPPER "${name}" upname)
  set(hfile "\n#ifndef _${upname}_H_\n#define _${upname}_H_
            \n#ifdef __cplusplus\nextern \"C\" {\n#endif
            \nextern const char* ${name}_src;
            \n#ifdef __cplusplus\n} \n#endif
            \n#endif /* _${upname}_H_ */\n\n")

  file(WRITE "${c_kernel_dir}/${name}.c" "${cfile}")
  file(WRITE "${c_kernel_dir}/${name}.h" "${hfile}")
endfunction()

