################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../glfw/deps/getopt.c \
../glfw/deps/tinycthread.c 

OBJS += \
./glfw/deps/getopt.o \
./glfw/deps/tinycthread.o 

C_DEPS += \
./glfw/deps/getopt.d \
./glfw/deps/tinycthread.d 


# Each subdirectory must supply rules for building sources it contributes
glfw/deps/%.o: ../glfw/deps/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


