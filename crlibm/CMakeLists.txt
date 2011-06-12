#
# Taken from https://github.com/Milkyway-at-home/crlibm
#
# 9c4c1c7e5f87f29535488b605ce38777e1a8b004 (Mon Jun 6 22:30:32 2011 -0400)
#
# Stripped out documentation and some extra stuff
# Go there for actually working on it or documentation
#


cmake_minimum_required(VERSION 2.6)
project(crlibm C)

include(CheckIncludeFile)
include(CheckFunctionExists)
include(CPUNameTest)

option(CRLIBM_QUIET_BUILD "Quietly build" ON)

if(NOT CRLIBM_QUIET_BUILD)
  message("Trying to build crlibm on an ${CMAKE_SYSTEM_PROCESSOR} under ${CMAKE_SYSTEM_NAME}, compiler is ${CMAKE_C_COMPILER}")
endif()

  #choose release build type by default
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
      "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
endif()

# What was in the autotools
set(crlibm_headers crlibm.h
                   crlibm_private.h
                   crlibm_msvc_hacks.h
                   triple-double.h
	               exp-td.h
	               expm1.h
                   log.h
                   rem_pio2_accurate.h
                   trigo_fast.h
                   trigo_accurate.h
                   trigpi.h
	               asincos.h
	               pow.h
                   atan_fast.h
                   atan_accurate.h
                   csh_fast.h)

set(crlibm_src crlibm_private.c
               triple-double.c
               exp-td.c
               exp-td-standalone.c
               expm1-standalone.c
               expm1.c
	           log.c
	           log1p.c
	           rem_pio2_accurate.c
	           trigo_fast.c
               trigo_accurate.c
	           trigpi.c
               asincos.c
               pow.c
	           atan_fast.c
               atan_accurate.c
	           csh_fast.c)


# All of the stuff that's here
# set(crlibm_headers crlibm.h
#                    crlibm_private.h
#                    triple-double.h
#                    expm1.h
# 	                 exp-td.h
#                    log-td.h
# 	                 log10-td.h
#                    log2-td.h
# 	                 rem_pio2_accurate.h
#                    trigo_fast.h
#                    trigo_accurate.h
# 	                 asin-td.h
#                    acos-td.h
#                    atan_fast.h
#                    atan_accurate.h
#                    csh_fast.h
#
#                    # headers missing from cmake one
#                    log.h
#                    trigpi.h
#  	                 asincos.h
#                    pow.h)

# set(crlibm_src crlibm_private.c
#                exp-td.c
#                exp-td-standalone.c
# 	             expm1-standalone.c
#                expm1.c
# 	             log-td.c
# 	             log1p.c
#                log10-td.c
#                log2-td.c
#                rem_pio2_accurate.c
# 	             trigo_fast.c
#                trigo_accurate.c
#                asin-td.c
#                # acos-td.c
# 	             atan_fast.c
#                atan_accurate.c
#                csh_fast.c
#
#                # source missing from cmake one
#                triple-double.c
#  	             log.c
#    	         trigpi.c
#                asincos.c
#                pow.c)

set(scs_lib_src scs_lib/scs_private.c
                scs_lib/addition_scs.c
                scs_lib/division_scs.c
                scs_lib/print_scs.c
	            scs_lib/double2scs.c
                scs_lib/zero_scs.c
	            scs_lib/multiplication_scs.c
                scs_lib/scs2double.c)


set(scs_lib_headers scs_lib/tests/tbx_timing.h)

list(APPEND crlibm_src ${scs_lib_src})
list(APPEND crlibm_headers ${scs_lib_headers})

check_include_file(fenv.h HAVE_FENV_H)
check_include_file(float.h HAVE_FLOAT_H)
check_include_file(fpu_control.h HAVE_FPU_CONTROL_H)
check_include_file(inttypes.h HAVE_INTTYPES_H)

add_definitions(-DHAVE_CONFIG_H)


# IBM's Libultim

find_path(IBM_H
	  MathLib.h
	  /usr/include/ /usr/local/include/
	  DOC "Path of Mathlib.h, the include file for IBM's libultim, for comparison purpose"
)

find_library(IBM_LIB
             NAMES ultim
             PATH /usr/lib/ /usr/local/lib/
             DOC "Directory of IBM's libultim, for comparison purpose")

# for compatibility with configure
if(IBM_H AND IBM_LIB)
  add_definitions(-DHAVE_MATHLIB_H)
  include_directories(${IBM_H})   # add include files
  link_directories(${IBM_LIB})    # add link dir and link flag
  set(LIBIBM ultim)
endif()


# Sun's libmcr

find_path(MCR_H
	  mcr.h
	  /usr/include/ /usr/local/include/
	  DOC "Path of mcr.h, the include file for Sun's libmcr, for comparison purpose"
)

find_library(MCR_LIB
             NAMES mcr
             PATH /usr/lib/ /usr/local/lib/
             DOC "Directory of Sun's libmcr, for comparison purpose")

# for compatibility with configure
if(MCR_H AND MCR_LIB)
  add_definitions(-DHAVE_MCR_H)
  include_directories(${MCR_H})   # add include files
  link_directories(${MCR_LIB})    # add link dir and link flag
  set(LIBMCR mcr)
endif()




option(ENABLE_SSE2 "Enable SSE2" ON)

if(SYSTEM_IS_X86)
  if(NOT MSVC)
    set(SSE2_FLAGS "-mfpmath=sse -msse -msse2")
  else()
    set(SSE2_FLAGS "${CMAKE_C_FLAGS} /arch:SSE2")
  endif()
endif()


option(ENABLE_DOUBLE_EXTENDED "Enable double-extended optimised versions of functions if hardware supports it" ON)


add_definitions(-DSCS_NB_WORDS=8 -DSCS_NB_BITS=30)

include_directories(${PROJECT_SOURCE_DIR})

enable_testing()
add_subdirectory(tests EXCLUDE_FROM_ALL)




set(TEST_DIR ${PROJECT_SOURCE_DIR}/tests)

add_test(exp ${TEST_DIR}/blind_test ${TEST_DIR}/exp.testdata)
add_test(log ${TEST_DIR}/blind_test ${TEST_DIR}/log.testdata)
add_test(sin ${TEST_DIR}/blind_test ${TEST_DIR}/sin.testdata)

add_test(cos ${TEST_DIR}/blind_test ${TEST_DIR}/cos.testdata)
add_test(tan ${TEST_DIR}/blind_test ${TEST_DIR}/tan.testdata)
add_test(atan ${TEST_DIR}/blind_test ${TEST_DIR}/atan.testdata)
add_test(sinh ${TEST_DIR}/blind_test ${TEST_DIR}/sinh.testdata)
add_test(cosh ${TEST_DIR}/blind_test ${TEST_DIR}/cosh.testdata)



# First check the host processor

if(CMAKE_SYSTEM_PROCESSOR MATCHES i*86|k6|k7|pentium*|athlon*)
    # CMAKE_SYSTEM_PROCESSOR is unfortunately still i386 on 64 bit OS X
  if(CMAKE_SIZEOF_VOID_P EQUAL 8) # 64-bit
    add_definitions(-DCRLIBM_TYPECPU_AMD64)
  else()
    add_definitions(-DCRLIBM_TYPECPU_X86)
  endif()
endif()

if(CMAKE_SYSTEM_PROCESSOR MATCHES alpha)
  add_definitions(-DCRLIBM_TYPECPU_ALPHA)
endif()

if(CMAKE_SYSTEM_PROCESSOR MATCHES powerpc*)
  add_definitions(-DCRLIBM_TYPECPU_POWERPC)
endif()

if(CMAKE_SYSTEM_PROCESSOR MATCHES sparc)
  add_definitions(-DCRLIBM_TYPECPU_SPARC)
endif()

if(CMAKE_SYSTEM_PROCESSOR MATCHES x86_64|amd64)
  add_definitions(-DCRLIBM_TYPECPU_AMD64)
endif()

if(CMAKE_SYSTEM_PROCESSOR MATCHES ia64)
  add_definitions(-DCRLIBM_TYPECPU_ITANIUM)
endif()


if(ENABLE_SSE2)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${SSE2_FLAGS}")
endif()


option(BUILD_INTERVAL_FUNCTIONS "Something" OFF)


configure_file (
  "${PROJECT_SOURCE_DIR}/crlibm_config.h.in"
  "${PROJECT_SOURCE_DIR}/crlibm_config.h")

if(MSVC)
  set_property(SOURCE ${crlibm_src} PROPERTY LANGUAGE CXX)
endif()

add_library(crlibm STATIC ${crlibm_src} ${crlibm_headers})

set_target_properties(crlibm PROPERTIES PUBLIC_HEADER crlibm.h)

if(CRLIBM_QUIET_BUILD)
  set_target_properties(crlibm PROPERTIES COMPILE_FLAGS "-w")
endif()

# install(TARGETS crlibm
#         ARCHIVE       DESTINATION lib
#         PUBLIC_HEADER DESTINATION include)

