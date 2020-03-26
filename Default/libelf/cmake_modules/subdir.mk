################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../libelf/cmake_modules/CheckHeaderSTDC.c 

OBJS += \
./libelf/cmake_modules/CheckHeaderSTDC.o 

C_DEPS += \
./libelf/cmake_modules/CheckHeaderSTDC.d 


# Each subdirectory must supply rules for building sources it contributes
libelf/cmake_modules/%.o: ../libelf/cmake_modules/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


