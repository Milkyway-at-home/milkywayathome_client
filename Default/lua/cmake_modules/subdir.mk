################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../lua/cmake_modules/test_private_extern.c 

OBJS += \
./lua/cmake_modules/test_private_extern.o 

C_DEPS += \
./lua/cmake_modules/test_private_extern.d 


# Each subdirectory must supply rules for building sources it contributes
lua/cmake_modules/%.o: ../lua/cmake_modules/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


