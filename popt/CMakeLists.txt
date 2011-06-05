
cmake_minimum_required(VERSION 2.6)
project(popt)

include(CheckIncludeFiles)
include(CheckLibraryExists)
include(CheckFunctionExists)
include(IgnoreWarningFlags)

# We want this to be small, so use -Os unless we're debugging
if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_BUILD_TYPE MinSizeRel)
endif()


set(POPT_VERSION_MAJOR 1)
set(POPT_VERSION_MINOR 17)
set(POPT_VERSION "${POPT_VERSION_MAJOR}.${POPT_VERSION_MINOR}")

if(CMAKE_COMPILER_IS_GNUCC)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -std=gnu99")
endif()

check_function_exists(CFPreferencesCopyAppValue HAVE_CFPREFERENCESCOPYAPPVALUE)
check_function_exists(dcgettext HAVE_DCGETTEXT)

check_include_files(dlfcn.h HAVE_DLFCN_H)
check_include_files(float.h HAVE_FLOAT_H)
check_include_files(fnmatch.h HAVE_FNMATCH_H)

check_function_exists(geteuid HAVE_GETEUID)
check_function_exists(gettext HAVE_GETTEXT)
check_function_exists(getuid HAVE_GETUID)

check_include_files(glob.h HAVE_GLOB_H)

# We don't need anything libintl or libiconv provide, and they
# massively increase the final binary size, so just pretend we never have them

#check_function_exists(iconv HAVE_ICONV)

set(HAVE_ICONV FALSE)

check_include_files(inttypes.h HAVE_INTTYPES_H)
check_include_files(langinfo.h HAVE_LANGINFO_H)
#check_include_files(libintl.h HAVE_LIBINTL_H)
set(HAVE_LIBINTL_H FALSE)

check_include_files(mcheck.h HAVE_MCHECK_H)
check_include_files(memory.h HAVE_MEMORY_H)

check_function_exists(mtrace HAVE_MTRACE)
check_function_exists(setregid HAVE_SETREGID)
check_function_exists(srandom HAVE_SRANDOM)

check_include_files(stdint.h HAVE_STDINT_H)
check_include_files(stdlib.h HAVE_STDLIB_H)

check_function_exists(stpcpy HAVE_STPCPY)
check_function_exists(strerror HAVE_STRERROR)

check_include_files(strings.h HAVE_STRINGS_H)
check_include_files(string.h HAVE_STRING_H)
check_include_files(sys/stat.h HAVE_SYS_STAT_H)
check_include_files(sys/types.h HAVE_SYS_TYPES_H)
check_include_files(unistd.h HAVE_UNISTD_H)

check_function_exists(vasprintf HAVE_VASPRINTF)
check_function_exists(__secure_getenv HAVE___SECURE_GETENV)

set(libpopt_src popt/popt.c
                popt/poptparse.c
                popt/poptconfig.c
                popt/popthelp.c
                popt/poptint.c)

set(libpopt_headers popt/popt.h
                    popt/poptint.h
                    popt/system.h)
set(libpopt_public_headers popt/popt.h)

set(LINK_FLAGS "-no-undefined ${LINK_FLAGS}")

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake ${PROJECT_BINARY_DIR}/config.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/configmake.h.cmake ${PROJECT_BINARY_DIR}/configmake.h)

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/popt)

add_definitions("-DHAVE_CONFIG_H")

add_library(popt STATIC ${libpopt_src} ${libpopt_headers})


set_target_properties(popt PROPERTIES
                             PUBLIC_HEADER ${libpopt_public_headers}
                             COMPILE_FLAGS "${DISABLE_WARNINGS_FLAG}")

set(POPT_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/popt" CACHE INTERNAL "popt headers")
set(POPT_LIBRARY popt CACHE INTERNAL "popt library")

########### tests ###############

add_executable(test1 popt/test1.c)
target_link_libraries(test1 popt)
set_target_properties(test1 PROPERTIES EXCLUDE_FROM_ALL TRUE)

add_executable(test2 popt/test2.c)
target_link_libraries(test2 popt)
set_target_properties(test2 PROPERTIES EXCLUDE_FROM_ALL TRUE)

add_executable(test3 popt/test3.c)
target_link_libraries(test3 popt)
set_target_properties(test3 PROPERTIES EXCLUDE_FROM_ALL TRUE)

add_executable(tdict popt/tdict.c)
target_link_libraries(tdict popt)
set_target_properties(tdict PROPERTIES EXCLUDE_FROM_ALL TRUE)

########### install files ###############

install(TARGETS popt
          ARCHIVE       DESTINATION lib
          PUBLIC_HEADER DESTINATION include)

