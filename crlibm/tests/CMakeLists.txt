
find_package(GMP)
find_package(MPFR)

if(LIBGMP_FOUND)
  include_directories(${LIBGMP_INCLUDE_DIR})
  add_definitions(-DHAVE_GMP_H)
endif()

if(LIBMPFR_FOUND)
  include_directories(${LIBMPFR_INCLUDE_DIR})
  add_definitions(-DHAVE_MPFR_H)
endif()



add_executable(blind_test blind_test.c test_common.c)

if(LIBGMP_FOUND)
  target_link_libraries(blind_test ${LIBGMP_LIB})
endif()

if(LIBMPFR_FOUND)
  target_link_libraries(blind_test ${LIBMPFR_LIB})
endif()


include_directories(${PROJECT_SOURCE_DIR})

link_directories(${PROJECT_SOURCE_DIR})
target_link_libraries(blind_test crlibm m)

set(test_data acos
              acospi
              asin
              asinpi
              atan
              atanpi
              cos
              cosh
              cospi
              exp
              expm1
              log
              log10
              log1p
              log2
              pow
              sin
              sinh
              sinpi
              tan
              tanpi)

foreach(test ${test_data})
  add_test(NAME ${test} COMMAND $<TARGET_FILE:blind_test> -v ${CRLIBM_TEST_DIR}/${test}.testdata)
endforeach()


