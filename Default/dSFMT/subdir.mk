################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../dSFMT/dSFMT.c \
../dSFMT/test.c 

OBJS += \
./dSFMT/dSFMT.o \
./dSFMT/test.o 

C_DEPS += \
./dSFMT/dSFMT.d \
./dSFMT/test.d 


# Each subdirectory must supply rules for building sources it contributes
dSFMT/%.o: ../dSFMT/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


