################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../tools/xxd/xxd.c 

OBJS += \
./tools/xxd/xxd.o 

C_DEPS += \
./tools/xxd/xxd.d 


# Each subdirectory must supply rules for building sources it contributes
tools/xxd/%.o: ../tools/xxd/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


