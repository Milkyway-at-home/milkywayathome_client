################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../boinc/boinc/samples/cygwin_fstab/fstab.c 

OBJS += \
./boinc/boinc/samples/cygwin_fstab/fstab.o 

C_DEPS += \
./boinc/boinc/samples/cygwin_fstab/fstab.d 


# Each subdirectory must supply rules for building sources it contributes
boinc/boinc/samples/cygwin_fstab/%.o: ../boinc/boinc/samples/cygwin_fstab/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


