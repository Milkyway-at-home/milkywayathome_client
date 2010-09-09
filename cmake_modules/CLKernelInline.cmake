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

# Run the C preprocessor to dump everything into 1 file, which we then
# can embed in the source and ship that without needing to send source
# files.

#FIXME: Passing a list of include arguments wasn't working and I don't
#feel like fighting with cmake now.
function(inline_kernel kernel_file name c_kernel_dir include_dir1 include_dir2 include_dir3)

  #FIXME: Less hacky than definining __OPENCL_VERSION__
  #This is also somewhat sketchy in the first place.
  #Run preprocessor
  execute_process(COMMAND ${CMAKE_C_COMPILER} -I${include_dir1} -I${include_dir2} -I${include_dir3} -std=c99 -D__OPENCL_VERSION__=100 -DSEPARATION_INLINE_KERNEL=1 -E -x c ${kernel_file}
                   OUTPUT_VARIABLE kernel_cpp OUTPUT_STRIP_TRAILING_WHITESPACE)

  #file(WRITE "${c_kernel_dir}/${name}.cl" "${kernel_cpp}")

  #Escape special characters
  #TODO: Other things that need escaping
  string(REGEX REPLACE
             "(\n)|(\r\n)"
             "\\\\n"
             str_escaped
             "${kernel_cpp}")

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

