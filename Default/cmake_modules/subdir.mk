################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../cmake_modules/test_avx.c 

OBJS += \
./cmake_modules/test_avx.o 

C_DEPS += \
./cmake_modules/test_avx.d 


# Each subdirectory must supply rules for building sources it contributes
cmake_modules/%.o: ../cmake_modules/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


