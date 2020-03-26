################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../boinc/boinc/zip/zip/unix/z_unix.c 

OBJS += \
./boinc/boinc/zip/zip/unix/z_unix.o 

C_DEPS += \
./boinc/boinc/zip/zip/unix/z_unix.d 


# Each subdirectory must supply rules for building sources it contributes
boinc/boinc/zip/zip/unix/%.o: ../boinc/boinc/zip/zip/unix/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


