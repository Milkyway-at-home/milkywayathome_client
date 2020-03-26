################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../simulation/src/lbatoxyz.cpp \
../simulation/src/mwdemo.cpp \
../simulation/src/mwdemo_fullscreen.cpp \
../simulation/src/nbody.cpp \
../simulation/src/readcel.cpp \
../simulation/src/render.cpp \
../simulation/src/rmvdup.cpp \
../simulation/src/showstar.cpp \
../simulation/src/ssdemo.cpp 

OBJS += \
./simulation/src/lbatoxyz.o \
./simulation/src/mwdemo.o \
./simulation/src/mwdemo_fullscreen.o \
./simulation/src/nbody.o \
./simulation/src/readcel.o \
./simulation/src/render.o \
./simulation/src/rmvdup.o \
./simulation/src/showstar.o \
./simulation/src/ssdemo.o 

CPP_DEPS += \
./simulation/src/lbatoxyz.d \
./simulation/src/mwdemo.d \
./simulation/src/mwdemo_fullscreen.d \
./simulation/src/nbody.d \
./simulation/src/readcel.d \
./simulation/src/render.d \
./simulation/src/rmvdup.d \
./simulation/src/showstar.d \
./simulation/src/ssdemo.d 


# Each subdirectory must supply rules for building sources it contributes
simulation/src/%.o: ../simulation/src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


