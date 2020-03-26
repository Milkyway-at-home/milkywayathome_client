################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../crlibm/scs_lib/tests/test_accuracy.c \
../crlibm/scs_lib/tests/test_log.c \
../crlibm/scs_lib/tests/test_timing.c 

OBJS += \
./crlibm/scs_lib/tests/test_accuracy.o \
./crlibm/scs_lib/tests/test_log.o \
./crlibm/scs_lib/tests/test_timing.o 

C_DEPS += \
./crlibm/scs_lib/tests/test_accuracy.d \
./crlibm/scs_lib/tests/test_log.d \
./crlibm/scs_lib/tests/test_timing.d 


# Each subdirectory must supply rules for building sources it contributes
crlibm/scs_lib/tests/%.o: ../crlibm/scs_lib/tests/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


