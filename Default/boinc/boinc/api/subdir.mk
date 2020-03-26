################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../boinc/boinc/api/boinc_api.cpp \
../boinc/boinc/api/boinc_api_fortran.cpp \
../boinc/boinc/api/boinc_opencl.cpp \
../boinc/boinc/api/graphics2.cpp \
../boinc/boinc/api/graphics2_unix.cpp \
../boinc/boinc/api/graphics2_util.cpp \
../boinc/boinc/api/graphics2_win.cpp \
../boinc/boinc/api/graphics_api.cpp \
../boinc/boinc/api/graphics_data.cpp \
../boinc/boinc/api/graphics_impl.cpp \
../boinc/boinc/api/graphics_impl_lib.cpp \
../boinc/boinc/api/graphics_lib.cpp \
../boinc/boinc/api/gutil.cpp \
../boinc/boinc/api/gutil_text.cpp \
../boinc/boinc/api/mac_icon.cpp \
../boinc/boinc/api/make_app_icon_h.cpp \
../boinc/boinc/api/reduce_lib.cpp \
../boinc/boinc/api/reduce_main.cpp \
../boinc/boinc/api/static_graphics.cpp \
../boinc/boinc/api/ttfont.cpp \
../boinc/boinc/api/windows_opengl.cpp \
../boinc/boinc/api/x_opengl.cpp 

OBJS += \
./boinc/boinc/api/boinc_api.o \
./boinc/boinc/api/boinc_api_fortran.o \
./boinc/boinc/api/boinc_opencl.o \
./boinc/boinc/api/graphics2.o \
./boinc/boinc/api/graphics2_unix.o \
./boinc/boinc/api/graphics2_util.o \
./boinc/boinc/api/graphics2_win.o \
./boinc/boinc/api/graphics_api.o \
./boinc/boinc/api/graphics_data.o \
./boinc/boinc/api/graphics_impl.o \
./boinc/boinc/api/graphics_impl_lib.o \
./boinc/boinc/api/graphics_lib.o \
./boinc/boinc/api/gutil.o \
./boinc/boinc/api/gutil_text.o \
./boinc/boinc/api/mac_icon.o \
./boinc/boinc/api/make_app_icon_h.o \
./boinc/boinc/api/reduce_lib.o \
./boinc/boinc/api/reduce_main.o \
./boinc/boinc/api/static_graphics.o \
./boinc/boinc/api/ttfont.o \
./boinc/boinc/api/windows_opengl.o \
./boinc/boinc/api/x_opengl.o 

CPP_DEPS += \
./boinc/boinc/api/boinc_api.d \
./boinc/boinc/api/boinc_api_fortran.d \
./boinc/boinc/api/boinc_opencl.d \
./boinc/boinc/api/graphics2.d \
./boinc/boinc/api/graphics2_unix.d \
./boinc/boinc/api/graphics2_util.d \
./boinc/boinc/api/graphics2_win.d \
./boinc/boinc/api/graphics_api.d \
./boinc/boinc/api/graphics_data.d \
./boinc/boinc/api/graphics_impl.d \
./boinc/boinc/api/graphics_impl_lib.d \
./boinc/boinc/api/graphics_lib.d \
./boinc/boinc/api/gutil.d \
./boinc/boinc/api/gutil_text.d \
./boinc/boinc/api/mac_icon.d \
./boinc/boinc/api/make_app_icon_h.d \
./boinc/boinc/api/reduce_lib.d \
./boinc/boinc/api/reduce_main.d \
./boinc/boinc/api/static_graphics.d \
./boinc/boinc/api/ttfont.d \
./boinc/boinc/api/windows_opengl.d \
./boinc/boinc/api/x_opengl.d 


# Each subdirectory must supply rules for building sources it contributes
boinc/boinc/api/%.o: ../boinc/boinc/api/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


