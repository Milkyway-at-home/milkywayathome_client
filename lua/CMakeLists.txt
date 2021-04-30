# Copyright (C) 2007-2011 LuaDist.
# Created by Peter Drahoš, Peter Kapec
# Redistribution and use of this file is allowed according to the terms of the MIT license.
# For details see the COPYRIGHT file distributed with LuaDist.
# Please note that the package source code is licensed under its own license.

project( lua C )
cmake_minimum_required( VERSION 2.6 )

include(PreferStatic)
include(IgnoreWarningFlags)
include(CPUNameTest)

if(APPLE)
    # Apple GCC extension that breaks with macports GCC on PowerPC
    # It is the same as __attribute__((visibility("hidden")))
    try_compile(HAVE_PRIVATE_EXTERN ${PROJECT_SOURCE_DIR}/cmake_modules ${PROJECT_SOURCE_DIR}/cmake_modules/test_private_extern.c)
    if(NOT HAVE_PRIVATE_EXTERN)
      add_definitions("-D__private_extern__=__attribute__((visibility(\"hidden\")))")
    endif()
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  add_definitions("-DLUA_USE_APICHECK")
endif()


## CONFIGURATION
# Default configuration (we assume POSIX by default)
set( LUA_PATH "LUA_PATH" CACHE INTERNAL "Environment variable to use as package.path." )
set( LUA_CPATH "LUA_CPATH" CACHE INTERNAL "Environment variable to use as package.cpath." )
set( LUA_INIT "LUA_INIT" CACHE INTERNAL "Environment variable for initial script." )

#option( LUA_ANSI "Use only ansi features." OFF )
option( LUA_USE_RELATIVE_LOADLIB "Use modified loadlib.c with support for relative paths on posix systems." ON)
mark_as_advanced(LUA_USE_RELATIVE_LOADLIB)

set( LUA_IDSIZE 60 CACHE INTERNAL "gives the maximum size for the description of the source." )
set( LUA_PROMPT "> " CACHE INTERNAL "Is the default prompt used by stand-alone Lua." )
set( LUA_PROMPT2 ">> " CACHE INTERNAL "Is the default continuation prompt used by stand-alone Lua." )
set( LUA_MAXINPUT 512 CACHE INTERNAL "Is the maximum length for an input line in the stand-alone interpreter.")

#2DO: LUAI_* and LUAL_* settings, for now defaults are used.

set( LUA_DIRSEP "/" )
set( LUA_MODULE_SUFFIX ${CMAKE_SHARED_MODULE_SUFFIX} )
set( LUA_LDIR ${INSTALL_LMOD} )
set( LUA_CDIR ${INSTALL_CMOD} )

if( LUA_USE_RELATIVE_LOADLIB )
  # This will set up relative paths to lib
  string( REGEX REPLACE "[^!/]+" ".." LUA_DIR "!/${INSTALL_BIN}/" )
else()
  # Direct path to installation
  set( LUA_DIR ${CMAKE_INSTALL_PREFIX} CACHE STRING "Destination from which modules will be resolved. See INSTALL_LMOD and INSTALL_CMOD.")
endif()

set( LUA_PATH_DEFAULT "./?.lua;${LUA_DIR}${LUA_LDIR}/?.lua;${LUA_DIR}${LUA_LDIR}/?/init.lua" )
set( LUA_CPATH_DEFAULT "./?${LUA_MODULE_SUFFIX};${LUA_DIR}${LUA_CDIR}/?${LUA_MODULE_SUFFIX};${LUA_DIR}${LUA_CDIR}/loadall${LUA_MODULE_SUFFIX}" )

if( WIN32 AND NOT CYGWIN )
  # Windows systems
  option( LUA_WIN "Windows specific build." ON )
  option( LUA_BUILD_AS_DLL "Build Lua library as Dll." ON )
  # Paths (Double escapes needed)
  set( LUA_DIRSEP "\\\\" )
  string( REPLACE "/" ${LUA_DIRSEP} LUA_DIR "${LUA_DIR}" )
  string( REPLACE "/" ${LUA_DIRSEP} LUA_LDIR "${LUA_LDIR}" )
  string( REPLACE "/" ${LUA_DIRSEP} LUA_CDIR "${LUA_CDIR}" )
  string( REPLACE "/" ${LUA_DIRSEP} LUA_PATH_DEFAULT "${LUA_PATH_DEFAULT}" )
  string( REPLACE "/" ${LUA_DIRSEP} LUA_CPATH_DEFAULT "${LUA_CPATH_DEFAULT}" )
else()
  # Posix systems (incl. Cygwin)
  option( LUA_USE_POSIX "Use POSIX functionality." ON )
  option( LUA_USE_DLOPEN "Use dynamic linker to load modules." OFF)
  option( LUA_USE_MKSTEMP "Use mkstep." ON )
  option( LUA_USE_ISATTY "Use tty." ON )
  option( LUA_USE_POPEN "Use popen." ON )
  option( LUA_USE_ULONGJMP "Use ulongjmp" ON)

  mark_as_advanced(LUA_USE_POSIX LUA_USE_DLOPEN LUA_USE_MKSTEMP LUA_USE_ISATTY LUA_USE_POPEN LUA_USE_ULONGJMP)
endif()

## SETUP
# Optional libraries
find_package_static( Readline )

if( READLINE_FOUND )
  option(LUA_USE_READLINE "Use readline in the Lua CLI." OFF)
  mark_as_advanced(LUA_USE_READLINE)
endif()

find_package(Curses)
if( CURSES_FOUND )
  option( LUA_USE_CURSES "Use curses in the Lua CLI." OFF )
  mark_as_advanced(LUA_USE_CURSES)
endif()

# Setup needed variables and libraries
if( LUA_USE_POSIX )
  # On POSIX Lua links to standard math library "m"
  list( APPEND LIBS m )
endif()

find_library(HAVE_DL dl)

if( LUA_USE_DLOPEN AND HAVE_DL)
  # Link to dynamic linker library "dl"
  list( APPEND LIBS dl )
endif()

if( LUA_USE_READLINE )
  # Add readline
  include_directories( ${READLINE_INCLUDE_DIR} )
  list( APPEND LIBS ${READLINE_LIBRARY} )
endif()

if( LUA_USE_CURSES )
  # Add curses
  include_directories( ${CURSES_INCLUDE_DIR} )
  list( APPEND LIBS ${CURSES_LIBRARY} )
endif()

## SOURCES
# Generate luaconf.h
configure_file( include/luaconf.h.in ${PROJECT_INCLUDE_DIR}/luaconf.h )

# Sources and headers
include_directories( include ${CMAKE_CURRENT_BINARY_DIR} )
set( SRC_LIB src/lapi.c src/lcode.c src/ldebug.c src/ldo.c src/ldump.c src/lfunc.c src/lgc.c src/llex.c src/lmem.c src/lobject.c src/lopcodes.c src/lparser.c src/lstate.c src/lstring.c src/ltable.c src/ltm.c src/lundump.c src/lvm.c src/lzio.c src/lauxlib.c src/lbaselib.c src/ldblib.c src/liolib.c src/lmathlib.c src/loslib.c src/ltablib.c src/lstrlib.c src/linit.c )
set( SRC_LUA src/lua.c )
set( SRC_LUAC src/luac.c src/print.c )

if( LUA_USE_RELATIVE_LOADLIB )
  # Use modified loadlib
  list( APPEND SRC_LIB src/loadlib_rel.c )
else()
  list( APPEND SRC_LIB src/loadlib.c )
endif()

## BUILD
add_library( lua51 STATIC ${SRC_LIB} )

add_executable( lua ${SRC_LUA} )
target_link_libraries( lua lua51 ${LIBS})

set_target_properties(lua51 PROPERTIES COMPILE_FLAGS "${DISABLE_WARNINGS_FLAG}")
set_target_properties(lua PROPERTIES COMPILE_FLAGS "${DISABLE_WARNINGS_FLAG}")

if(NOT ALWAYS_HAVE_SSE2)
  disable_sse2(lua51)
endif()


set(LUA_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include" CACHE INTERNAL "lua headers")
set(LUA_INTERPRETER $<TARGET_FILE:lua> CACHE INTERNAL "Lua interpreter")
set(LUA_LIBRARIES lua51 ${LIBS} CACHE INTERNAL "Lua libraries")

# install(TARGETS lua lua51
#           RUNTIME       DESTINATION bin
#           ARCHIVE       DESTINATION lib
#           PUBLIC_HEADER DESTINATION include)

