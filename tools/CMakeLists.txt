
cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH "../cmake_modules")

find_package(OpenCL)

add_executable(precision_test PrecisionTest.c)
target_link_libraries(precision_test m)

add_executable(precision_test_ocl PrecisionTest.c)
target_link_libraries(precision_test_ocl ${OPENCL_LIBRARIES})
set_target_properties(precision_test_ocl
                         PROPERTIES
                         COMPILE_DEFINITIONS "TEST_OPENCL=1")

