project(libelf)

cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules)

include(CheckTypeSize)
include(CheckIncludeFiles)
include(CheckFunctionExists)
include(CheckSymbolExists)
include(CheckHeaderSTDC)

if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_BUILD_TYPE MinSizeRel)
endif()


set(LIBELF_SRC_DIR libelf-0.8.13)

set(LIBELF_VERSION "LIBELF")


option(ENABLE_SANITY_CHECKS "enable sanity checks by default" ON)
option(ENABLE_NLS "use Native Language Support" OFF)
option(ENABLE_ELF64 "compile in 64-bit support" ON)
option(EXTENDED_FORMAT "support extended file formats" OFF)

mark_as_advanced(ENABLE_ELF64)
mark_as_advanced(ENABLE_EXTENDED_FORMAT)

set(elf_headers lib/elf_repl.h
                lib/byteswap.h
                lib/errors.h
                lib/ext_types.h
                lib/gelf.h
                lib/libelf.h
                lib/nlist.h
                lib/private.h
                lib/verdef.h
                lib/verneed.h)

set(elf_src lib/32.fsize.c
            lib/32.getehdr.c
            lib/32.getphdr.c
            lib/32.getshdr.c
            lib/32.newehdr.c
            lib/32.newphdr.c
            lib/32.xlatetof.c
            lib/64.xlatetof.c
            lib/assert.c
            lib/begin.c
            lib/checksum.c
            lib/cntl.c
            lib/cook.c
            lib/data.c
            lib/end.c
            lib/errmsg.c
            lib/errno.c
            lib/fill.c
            lib/flag.c
            lib/gelfehdr.c
            lib/gelfphdr.c
            lib/gelfshdr.c
            lib/gelftrans.c
            lib/getarhdr.c
            lib/getaroff.c
            lib/getarsym.c
            lib/getbase.c
            lib/getdata.c
            lib/getident.c
            lib/getscn.c
            lib/hash.c
            lib/input.c
            lib/kind.c
            lib/memset.c
            lib/ndxscn.c
            lib/newdata.c
            lib/newscn.c
            lib/next.c
            lib/nextscn.c
            lib/nlist.c
            lib/opt.delscn.c
            lib/rand.c
            lib/rawdata.c
            lib/rawfile.c
            lib/strptr.c
            lib/swap64.c
            lib/update.c
            lib/verdef_32_tof.c
            lib/verdef_32_tom.c
            lib/verdef_64_tof.c
            lib/verdef_64_tom.c
            lib/version.c
            lib/x.elfext.c
            lib/x.movscn.c
            lib/x.remscn.c
            )

set(CMAKE_EXTRA_INCLUDE_FILES sys/types.h)
check_type_size("size_t" SIZE_T)
check_type_size("off_t" OFF_T)
set(CMAKE_EXTRA_INCLUDE_FILES)


#Find 64 bit integer
check_type_size("int64_t" INT64_T)
check_type_size("uint64_t" UINT64_T)
if(INT64_T EQUAL 8)
  set(LIBELF_CV_INT64 "int64_t")
endif()
if(UINT64_T EQUAL 8)
  set(LIBELF_CV_UINT64 "uint64_t")
endif()
if(NOT LIBELF_CV_INT64 OR NOT LIBELF_CV_UINT64)
  check_type_size("long" LONG)
  check_type_size("__int64" __INT64)
  check_type_size("long long" LONG_LONG)

  if(__INT64 EQUAL 8)
    set(LIBELF_CV_INT64 "__int64")
  elseif(LONG_LONG EQUAL 8)
    set(LIBELF_CV_INT64 "long long")
  elseif(LONG EQUAL 8)
    set(LIBELF_CV_INT64 "long")
  else()
    message(FATAL "Did not find a 64 bit int type")
  endif()

  set(LIBELF_CV_UINT64 "unsigned ${LIBELF_CV_INT64}")
endif()

#Find 32 bit integer
check_type_size("int32_t" INT32_T)
check_type_size("uint32_t" UINT32_T)
if(INT32_T EQUAL 4)
  set(LIBELF_CV_INT32 "int32_t")
endif()
if(UINT32_T EQUAL 4)
  set(LIBELF_CV_UINT32 "uint32_t")
endif()
if(NOT LIBELF_CV_INT32 OR NOT LIBELF_CV_UINT32)
  check_type_size("int" LONG)
  check_type_size("__int32" __INT32)
  check_type_size("long" LONG)

  if(__INT32 EQUAL 4)
    set(LIBELF_CV_INT32 "__int32")
  elseif(INT EQUAL 4)
    set(LIBELF_CV_INT32 "__int32")
  elseif(LONG EQUAL 4)
    set(LIBELF_CV_INT32 "long")
  else()
    message(FATAL "Did not find a 32 bit int type")
  endif()

  set(LIBELF_CV_UINT32 "unsigned ${LIBELF_CV_INT32}")
endif()

#Find 16 bit integer
check_type_size("int16_t" INT16_T)
check_type_size("uint16_t" UINT16_T)
if(INT16_T EQUAL 2)
  set(LIBELF_CV_INT16 "int16_t")
endif()
if(UINT16_T EQUAL 2)
  set(LIBELF_CV_UINT16 "uint16_t")
endif()
if(NOT LIBELF_CV_INT16 OR NOT LIBELF_CV_UINT16)
  check_type_size("__int16" __INT16)
  check_type_size("short" LONG)
  check_type_size("int" LONG)

  if(__INT16 EQUAL 2)
    set(LIBELF_CV_INT16 "__int16")
  elseif(SHORT EQUAL 2)
    set(LIBELF_CV_INT16 "short")
  elseif(INT EQUAL 2)
    set(LIBELF_CV_INT16 "int")
  else()
    message(FATAL "Did not find a 16 bit int type")
  endif()

  set(LIBELF_CV_UINT16 "unsigned ${LIBELF_CV_INT16}")
endif()




check_function_exists(mmap HAVE_MMAP)
check_function_exists(catgets HAVE_CATGETS)
check_function_exists(dgettext HAVE_DGETTEXT)
check_function_exists(ftruncate HAVE_FTRUNCATE)
check_function_exists(getpagesize HAVE_GETPAGESIZE)

check_symbol_exists(memcmp string.h HAVE_MEMCMP)
check_symbol_exists(memcpy string.h HAVE_MEMCPY)
check_symbol_exists(memmove string.h HAVE_MEMMOVE)
check_symbol_exists(memset string.h HAVE_MEMSET)



check_include_files(ar.h HAVE_AR_H)
check_include_files(elf.h HAVE_ELF_H)
check_include_files(fcntl.h HAVE_FCNTL_H)
#check_include_files(gelf.h HAVE_GELF_H)

check_include_files(link.h HAVE_LINK_H)
check_include_files(nlist.h HAVE_NLIST_H)
check_include_files(stdint.h HAVE_STDINT_H)

check_include_files(sys/elf.h HAVE_SYS_ELF_H)
check_include_files(sys/link.h HAVE_SYS_LINK_H)
check_include_files(unistd.h HAVE_UNISTD_H)

set(CMAKE_EXTRA_INCLUDE_FILES link.h)
check_type_size("Elf32_Dyn" ELF32_DYN)
if(ELF32_DYN)
  set(__LIBELF_NEED_LINK_H)
endif()

set(CMAKE_EXTRA_INCLUDE_FILES sys/link.h)
check_type_size("Elf32_Dyn" ELF32_DYN)
if(ELF32_DYN)
  set(__LIBELF_NEED_SYS_LINK_H)
endif()
set(CMAKE_EXTRA_INCLUDE_FILES)



# Missing some nonsense I don't care about for Irix
if(HAVE_ELF_H)
  set(CMAKE_EXTRA_INCLUDE_FILES elf.h)
elseif(HAVE_SYS_ELF_H)
  set(CMAKE_EXTRA_INCLUDE_FILES sys/elf.h)
endif()
check_type_size("Elf64_Addr" ELF64_ADDR)
check_type_size("Elf64_Rel" ELF64_Rel)
set(CMAKE_EXTRA_INCLUDE_FILES)


if(ENABLE_ELF64 AND (HAVE_ELF64_ADDR OR APPLE OR WIN32))
  set(__LIBELF64 TRUE)
  if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    set(__LIBELF64_LINUX TRUE)
  endif()
endif()

check_header_stdc()

add_definitions(-DHAVE_CONFIG_H=1)

# Some stuff is defined in config.h and then again in sys_elf.h, and I don't know why
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.h.in ${PROJECT_INCLUDE_DIR}/libelf/config.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/lib/sys_elf.h.in ${PROJECT_INCLUDE_DIR}/libelf/sys_elf.h)

include_directories(${PROJECT_INCLUDE_DIR})
include_directories(${PROJECT_INCLUDE_DIR}/libelf)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib)
add_library(elf STATIC ${elf_src} ${elf_headers})


set(LIBELF_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/lib;${PROJECT_INCLUDE_DIR}" CACHE INTERNAL "libelf headers")
set(LIBELF_LIBRARY elf CACHE INTERNAL "libelf library")




