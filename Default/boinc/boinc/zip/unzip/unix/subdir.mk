################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../boinc/boinc/zip/unzip/unix/unix.c 

OBJS += \
./boinc/boinc/zip/unzip/unix/unix.o 

C_DEPS += \
./boinc/boinc/zip/unzip/unix/unix.d 


# Each subdirectory must supply rules for building sources it contributes
boinc/boinc/zip/unzip/unix/%.o: ../boinc/boinc/zip/unzip/unix/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


