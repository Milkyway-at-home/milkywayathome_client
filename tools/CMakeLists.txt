
cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH "../cmake_modules")

find_package(OpenCL)

set(CMAKE_C_FLAGS "-O3 -std=c99 -mfpmath=sse -msse -msse2 -mno-3dnow ${CMAKE_C_FLAGS}")


add_executable(precision_test PrecisionTest.c)
target_link_libraries(precision_test m)

add_executable(precision_test_ocl PrecisionTest.c)
target_link_libraries(precision_test_ocl ${OPENCL_LIBRARIES})
set_target_properties(precision_test_ocl
                         PROPERTIES
                         COMPILE_DEFINITIONS "TEST_OPENCL=1")


include_directories("../nbody/thirdparty/fdlibm53/")
include_directories("../nbody/include/")
link_directories("../nbody/lib/")

#CHECKME: Why do we need to link to libm? Fails with undefined
#reference to sincos. I don't trust that. It doesn't appear to be used
#in fdlibm though.
add_executable(precision_test_fdlibm PrecisionTest.c)
target_link_libraries(precision_test_fdlibm fdlibm m)
set_target_properties(precision_test_fdlibm
                         PROPERTIES
                         COMPILE_DEFINITIONS "USE_FDLIBM=1")



add_executable(math_bench MathBench.c)
target_link_libraries(math_bench fdlibm crlibm)


