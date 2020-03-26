################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../popt/tests/tdict.c \
../popt/tests/test1.c \
../popt/tests/test2.c \
../popt/tests/test3.c 

OBJS += \
./popt/tests/tdict.o \
./popt/tests/test1.o \
./popt/tests/test2.o \
./popt/tests/test3.o 

C_DEPS += \
./popt/tests/tdict.d \
./popt/tests/test1.d \
./popt/tests/test2.d \
./popt/tests/test3.d 


# Each subdirectory must supply rules for building sources it contributes
popt/tests/%.o: ../popt/tests/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


