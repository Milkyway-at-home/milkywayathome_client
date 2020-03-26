################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../boinc/boinc/clienttray/tray_win.cpp 

OBJS += \
./boinc/boinc/clienttray/tray_win.o 

CPP_DEPS += \
./boinc/boinc/clienttray/tray_win.d 


# Each subdirectory must supply rules for building sources it contributes
boinc/boinc/clienttray/%.o: ../boinc/boinc/clienttray/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


