################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../boinc/boinc/clientctrl/boincsvcctrl.cpp 

OBJS += \
./boinc/boinc/clientctrl/boincsvcctrl.o 

CPP_DEPS += \
./boinc/boinc/clientctrl/boincsvcctrl.d 


# Each subdirectory must supply rules for building sources it contributes
boinc/boinc/clientctrl/%.o: ../boinc/boinc/clientctrl/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


