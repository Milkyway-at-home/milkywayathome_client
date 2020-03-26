################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../glfw/examples/boing.c \
../glfw/examples/gears.c \
../glfw/examples/heightmap.c \
../glfw/examples/simple.c \
../glfw/examples/splitview.c \
../glfw/examples/wave.c 

OBJS += \
./glfw/examples/boing.o \
./glfw/examples/gears.o \
./glfw/examples/heightmap.o \
./glfw/examples/simple.o \
./glfw/examples/splitview.o \
./glfw/examples/wave.o 

C_DEPS += \
./glfw/examples/boing.d \
./glfw/examples/gears.d \
./glfw/examples/heightmap.d \
./glfw/examples/simple.d \
./glfw/examples/splitview.d \
./glfw/examples/wave.d 


# Each subdirectory must supply rules for building sources it contributes
glfw/examples/%.o: ../glfw/examples/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


