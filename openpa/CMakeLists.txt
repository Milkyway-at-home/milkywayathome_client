
cmake_minimum_required(VERSION 2.6)
project(openpa)
enable_testing()

set(OPENPA_VERSION_MAJOR 1)
set(OPENPA_VERSION_MINOR 0)
set(OPENPA_VERSION_PATCHLEVEL 3)
set(OPENPA_VERSION "${OPENPA_VERSION_MAJOR}.${OPENPA_VERSION_MINOR}.${OPENPA_VERSION_PATCHLEVEL}")


include(CheckIncludeFiles)
include(CheckLibraryExists)
include(CheckFunctionExists)
include(CheckTypeSize)
include(CheckCSourceCompiles)


if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_BUILD_TYPE MinSizeRel)
endif()





check_include_files(atomic.h OPA_HAVE_ATOMIC_H)
check_include_files(intrin.h OPA_HAVE_INTRIN_H)
check_include_files(pthread.h OPA_HAVE_PTHREAD_H)
check_include_files(stddef.h OPA_HAVE_STDDEF_H)

check_library_exists(pthread pthread_create "" OPA_HAVE_LIBPTHREAD)

set(CMAKE_EXTRA_INCLUDE_FILES pthread.h)

set(CMAKE_REQUIRED_LIBRARIES pthread)
check_function_exists(pthread_yield OPA_HAVE_PTHREAD_YIELD)
set(CMAKE_REQUIRED_LIBRARIES )

set(CMAKE_EXTRA_INCLUDE_FILES sched.h)
check_function_exists(sched_yield OPA_HAVE_SCHED_YIELD)
set(CMAKE_EXTRA_INCLUDES )


check_type_size("int" SIZEOF_INT)
check_type_size("void*" SIZEOF_VOID_P)

#FIXME: These options don't actually work
set(MAX_NTHREADS 100 CACHE STRING "Maximum number of threads for tests")
option(EXPLICIT_EMULATION "Use lock-based emulation" FALSE)
option(USE_LOCK_BASED_PRIMITIVES "Force using lock-based atomic primitives" FALSE)
option(USE_UNSAFE_PRIMITIVES "Unsafe (non-atomic) primitives should be used" FALSE)
option(OPA_HAVE_STRICT_FAIRNESS_CHECKS "If strict checking of atomic operation fairness is desired" FALSE)
mark_as_advanced(MAX_NTHREADS)
mark_as_advanced(EXPLICIT_EMULATION)
mark_as_advanced(USE_LOCK_BASED_PRIMITIVES)
mark_as_advanced(USE_UNSAFE_PRIMITIVES)
mark_as_advanced(HAVE_STRICT_FAIRNESS_CHECKS)

check_c_source_compiles(
  "
  int foo(int) __attribute__ ((pure));

  int a;

  int main()
  {
      return 0;
  }
  "
OPA_HAVE_ATTRIBUTE_PURE)

check_c_source_compiles(
  "
    int foo(int) __attribute__ ((pure));

    int a;

    int main()
    {
        return 0;
    }
  "
  OPA_HAVE_ATTRIBUTE_PURE)

check_c_source_compiles(
  "
    int foo(char*, ...) __attribute__ ((format(printf,1,2)));

    int a;

    int main(int argc, const char* arv[])
    {
        return 0;
    }
  "
  OPA_HAVE_ATTRIBUTE_FORMAT)

if(OPA_HAVE_ATTRIBUTE_FORMAT AND OPA_HAVE_ATTRIBUTE_PURE)
  set(OPA_HAVE_GCC_ATTRIBUTE 1)
endif()

check_c_source_compiles(
  "
    int main(int argc, const char* arv[])
    {
        __asm__(\"ILLEGAL_ASM_STMT\");
        return 0;
    }
  "
  ACCEPTS_BOGUS_ASM)

macro(try_opa_header name header)
  set(CMAKE_REQUIRED_INCLUDES ${PROJECT_SOURCE_DIR}/include;${PROJECT_SOURCE_DIR}/include/primitives)
  set(CMAKE_REQUIRED_DEFINITIONS "-DPRIM_HEADER=\"${header}\""
                                 "-DSIZEOF_INT=${SIZEOF_INT}"
                                 "-DSIZEOF_VOID_P=${SIZEOF_VOID_P}"
                                 "-DHAVE_GCC_ATTRIBUTE=${HAVE_GCC_ATTRIBUTE}")
  check_c_source_compiles(
  "
    #define OPA_SIZEOF_INT SIZEOF_INT
    #define OPA_SIZEOF_VOID_P SIZEOF_VOID_P
    #ifndef _opa_inline
    #define _opa_inline inline
    #endif
    #ifndef _opa_restrict
    #define _opa_restrict restrict
    #endif
    #ifndef _opa_const
    #define _opa_const const
    #endif
    #ifdef HAVE_GCC_ATTRIBUTE
    #define OPA_HAVE_GCC_ATTRIBUTE 1
    #endif
    #include \"opa_util.h\" /* for OPA_ATTRIBUTE and friends */
    #include PRIM_HEADER

    int main(int argc, const char* argv[])
    {
        OPA_int_t a, b;
        int c;

        OPA_store_int(&a, 0);
        OPA_store_int(&b, 1);
        c = OPA_load_int(&a);

        OPA_add_int(&a, 10);
        OPA_incr_int(&a);
        OPA_decr_int(&a);

        c = OPA_decr_and_test_int(&a);
        c = OPA_fetch_and_add_int(&a, 10);
        c = OPA_fetch_and_incr_int(&a);
        c = OPA_fetch_and_decr_int(&a);

        c = OPA_cas_int(&a, 10, 11);
        c = OPA_swap_int(&a, OPA_load_int(&b));

        OPA_write_barrier();
        OPA_read_barrier();
        OPA_read_write_barrier();

        return 0;
    }
  "
  ${name})

  set(CMAKE_REQUIRED_INCLUDES )
  set(CMAKE_REQUIRED_DEFINITIONS )

#  set(CMAKE_REQUIRED_INCLUDES )
#  message("Check out my monster dong ${name} ${MONSTER_DONG}")
endmacro()

if(NOT ACCEPTS_BOGUS_ASM)
  if(NOT CMAKE_CROSSCOMPILING)
    try_opa_header(OPA_HAVE_GCC_X86_32_64 "opa_gcc_intel_32_64.h")
  endif()

  try_opa_header(OPA_HAVE_GCC_X86_32_64_P3 "opa_gcc_intel_32_64_p3.h")
  try_opa_header(OPA_HAVE_GCC_AND_IA64_ASM "opa_gcc_ia64.h")
  try_opa_header(OPA_HAVE_GCC_AND_POWERPC_ASM "opa_gcc_ppc.h")
  try_opa_header(OPA_HAVE_GCC_AND_SICORTEX_ASM "opa_gcc_sicortex.h")
endif()

try_opa_header(OPA_HAVE_GCC_INTRINSIC_ATOMICS "opa_gcc_intrinsics.h")
try_opa_header(OPA_HAVE_NT_INTRINSICS "opa_nt_intrinsics.h")
try_opa_header(OPA_HAVE_SUN_ATOMIC_OPS "opa_sun_atomic_ops.h")



set(opa_src src/opa_primitives.c
            src/opa_queue.c)
set(opa_headers include/opa_primitives.h
                include/opa_queue.h
                include/opa_util.h

                include/primitives/opa_by_lock.h
                include/primitives/opa_gcc_ia64.h
                include/primitives/opa_gcc_intel_32_64.h
                include/primitives/opa_gcc_intel_32_64_p3barrier.h
                include/primitives/opa_gcc_intrinsics.h
                include/primitives/opa_gcc_sicortex.h
                include/primitives/opa_sun_atomic_ops.h
                include/primitives/opa_emulated.h
                include/primitives/opa_gcc_intel_32_64_barrier.h
                include/primitives/opa_gcc_intel_32_64_ops.h
                include/primitives/opa_gcc_intel_32_64_p3.h
                include/primitives/opa_gcc_ppc.h
                include/primitives/opa_nt_intrinsics.h
                include/primitives/opa_unsafe.h

                include/opa_config.h.in)

set(OPA_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include" CACHE INTERNAL "OPA headers")
set(OPA_LIBRARY opa CACHE INTERNAL "OPA library")


message("Dicks ${PROJECT_BINARY_DIR}/include")
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/opa_config.h.in ${PROJECT_INCLUDE_DIR}/opa_config.h)

include_directories(${OPA_INCLUDE_DIR})
include_directories(${PROJECT_INCLUDE_DIR})


add_library(opa ${opa_src} ${opa_headers})

add_subdirectory(test EXCLUDE_FROM_ALL)

