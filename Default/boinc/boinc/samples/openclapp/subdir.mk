################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../boinc/boinc/samples/openclapp/openclapp.cpp 

OBJS += \
./boinc/boinc/samples/openclapp/openclapp.o 

CPP_DEPS += \
./boinc/boinc/samples/openclapp/openclapp.d 


# Each subdirectory must supply rules for building sources it contributes
boinc/boinc/samples/openclapp/%.o: ../boinc/boinc/samples/openclapp/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


