################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../LModL/src/3dproj.cpp \
../LModL/src/binfile.cpp \
../LModL/src/draw.cpp \
../LModL/src/draw32.cpp \
../LModL/src/draw32p.cpp \
../LModL/src/draw8.cpp \
../LModL/src/drawcore.cpp \
../LModL/src/drawhalo.cpp \
../LModL/src/imgplot.cpp \
../LModL/src/rgbconv.cpp \
../LModL/src/trigtbl.cpp 

OBJS += \
./LModL/src/3dproj.o \
./LModL/src/binfile.o \
./LModL/src/draw.o \
./LModL/src/draw32.o \
./LModL/src/draw32p.o \
./LModL/src/draw8.o \
./LModL/src/drawcore.o \
./LModL/src/drawhalo.o \
./LModL/src/imgplot.o \
./LModL/src/rgbconv.o \
./LModL/src/trigtbl.o 

CPP_DEPS += \
./LModL/src/3dproj.d \
./LModL/src/binfile.d \
./LModL/src/draw.d \
./LModL/src/draw32.d \
./LModL/src/draw32p.d \
./LModL/src/draw8.d \
./LModL/src/drawcore.d \
./LModL/src/drawhalo.d \
./LModL/src/imgplot.d \
./LModL/src/rgbconv.d \
./LModL/src/trigtbl.d 


# Each subdirectory must supply rules for building sources it contributes
LModL/src/%.o: ../LModL/src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


