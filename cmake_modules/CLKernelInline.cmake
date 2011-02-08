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


# Visual Studio has some kind of idiotic string literal length
# limitation of 16380 characters which you work around by breaking it
# up into multiple literals which get concatenated
# Cmake also makes this painful to do. This function makes me want to cry.
function(visual_studio_is_shit str)
  # FIXME: if you get unlucky, this will break on special characters and escaping will get all broken.
  # I gave up on splitting on newline with cmake regex.
  set(chunk_size 12000)
  string(LENGTH "${str}" str_len)
  math(EXPR num_str_chunks "${str_len} / ${chunk_size}")
  math(EXPR len_chunks "${num_str_chunks} * ${chunk_size}")

  # Horrible work around lack of ceil()
  if (${len_chunks} LESS ${str_len})
    math(EXPR num_str_chunks "${num_str_chunks} + 1")
  endif()

  set(new_str "")
  set(chunk_len ${chunk_size})
  foreach(i RANGE 1 ${num_str_chunks})
    math(EXPR lower "(${i} - 1) * ${chunk_size}")
    math(EXPR upper "${i} * ${chunk_size}")

    if(${upper} GREATER ${str_len})
      math(EXPR chunk_len "${str_len} - ${lower}")
    endif()

    string(SUBSTRING "${str}" "${lower}" "${chunk_len}" sub_str)
    set(new_str "${new_str}\"${sub_str}\"\n")
  endforeach()
  set(BROKEN_UP_LITERAL "${new_str}" PARENT_SCOPE)
endfunction()

function(inline_kernel name full_src outfile)
  #Escape special characters
  #TODO: Other things that need escaping

  # Escape escaped newlines, e.g. those in macros
  string(REGEX REPLACE
             "\\\\(\n)|(\r\n)"
             "\\\\\\\\\\n"
             str_escaped
             "${full_src}")

  string(REGEX REPLACE
             "(\n)|(\r\n)"
             "\\\\n"
             str_escaped
             "${str_escaped}")

  string(REGEX REPLACE
             "(\")"
             "\\\\\""
             str_escaped
             "${str_escaped}")

  visual_studio_is_shit("${str_escaped}")
  set(cfile "\nconst char* ${name} = ${BROKEN_UP_LITERAL};\n\n")

  file(WRITE "${outfile}" "${cfile}")
endfunction()


macro(strip_includes)
  string(REGEX REPLACE "[ \t]*#include[ \t]*[\"<][^\">]*[\">]" "" str "${str}")
endmacro()

# Warning: This won't work if there are comments in strings
macro(strip_comments)
  #Remove the comments. Otherwise, C++ comments will comment out the
  #rest of the program unless we do more work.
  string(REGEX REPLACE
             "(/\\*([^*]|[\r\n]|(\\*+([^*/]|[\r\n])))*\\*+/)|(//[^\r\n]*)"
             "" # No comment
             str
             #Quote this to prevent ;'s being interpreted as list separators
             "${str}")
  string(REGEX REPLACE "\n\n" "\n" str "${str}")
endmacro()


# It's a Webkit joke

# Need to pull in all header file dependencies for the kernel, and
# strip out the #include's. Then we can pack it into a string, and
# then we don't need to deal with distributing extra files. Ideally we
# could also strip comments, extra white space, rename things to
# save space etc. Hopefully someday someone will write an OpenCL minifier
# that will do this for you. There seem to be GLSL ones that don't
# seem to work.

# We can't ship headers anyway. Because of the stupid way BOINC
# handles files, they won't persist in the project directory if we
# just put them there. They will either be deleted, or will have a
# stupid "XML link" thing which it uses, which obviously won't work
# with the CL compiler.
function(all_in_one_file src_list)
  set(all_in_one_string "")
  foreach(i ${src_list})
    file(READ ${i} str)
    strip_includes()
    strip_comments()
    set(all_in_one_string "${all_in_one_string}${str}")
  endforeach()
  set(ALL_IN_ONE_STRING "${all_in_one_string}" PARENT_SCOPE)
endfunction()

