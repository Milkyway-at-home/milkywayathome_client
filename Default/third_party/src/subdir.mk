################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../third_party/src/gl3w.c 

OBJS += \
./third_party/src/gl3w.o 

C_DEPS += \
./third_party/src/gl3w.d 


# Each subdirectory must supply rules for building sources it contributes
third_party/src/%.o: ../third_party/src/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


