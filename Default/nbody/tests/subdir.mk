################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../nbody/tests/bessel_test.c \
../nbody/tests/emd_test.c \
../nbody/tests/nbody_test_driver.c 

OBJS += \
./nbody/tests/bessel_test.o \
./nbody/tests/emd_test.o \
./nbody/tests/nbody_test_driver.o 

C_DEPS += \
./nbody/tests/bessel_test.d \
./nbody/tests/emd_test.d \
./nbody/tests/nbody_test_driver.d 


# Each subdirectory must supply rules for building sources it contributes
nbody/tests/%.o: ../nbody/tests/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


