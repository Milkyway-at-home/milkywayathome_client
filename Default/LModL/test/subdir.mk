################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../LModL/test/coltest.cpp \
../LModL/test/cubetest.cpp \
../LModL/test/cubetest_fullscreen.cpp \
../LModL/test/fpstest.cpp \
../LModL/test/gfxinfo.cpp \
../LModL/test/hsltest.cpp \
../LModL/test/imgrender.cpp \
../LModL/test/quattest.cpp \
../LModL/test/trigtblm.cpp 

C_SRCS += \
../LModL/test/text2bin.c 

OBJS += \
./LModL/test/coltest.o \
./LModL/test/cubetest.o \
./LModL/test/cubetest_fullscreen.o \
./LModL/test/fpstest.o \
./LModL/test/gfxinfo.o \
./LModL/test/hsltest.o \
./LModL/test/imgrender.o \
./LModL/test/quattest.o \
./LModL/test/text2bin.o \
./LModL/test/trigtblm.o 

CPP_DEPS += \
./LModL/test/coltest.d \
./LModL/test/cubetest.d \
./LModL/test/cubetest_fullscreen.d \
./LModL/test/fpstest.d \
./LModL/test/gfxinfo.d \
./LModL/test/hsltest.d \
./LModL/test/imgrender.d \
./LModL/test/quattest.d \
./LModL/test/trigtblm.d 

C_DEPS += \
./LModL/test/text2bin.d 


# Each subdirectory must supply rules for building sources it contributes
LModL/test/%.o: ../LModL/test/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

LModL/test/%.o: ../LModL/test/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


