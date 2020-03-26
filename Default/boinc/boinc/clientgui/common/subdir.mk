################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../boinc/boinc/clientgui/common/wxPieCtrl.cpp 

OBJS += \
./boinc/boinc/clientgui/common/wxPieCtrl.o 

CPP_DEPS += \
./boinc/boinc/clientgui/common/wxPieCtrl.d 


# Each subdirectory must supply rules for building sources it contributes
boinc/boinc/clientgui/common/%.o: ../boinc/boinc/clientgui/common/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


