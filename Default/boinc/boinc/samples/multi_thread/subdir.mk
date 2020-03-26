################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../boinc/boinc/samples/multi_thread/multi_thread.cpp 

OBJS += \
./boinc/boinc/samples/multi_thread/multi_thread.o 

CPP_DEPS += \
./boinc/boinc/samples/multi_thread/multi_thread.d 


# Each subdirectory must supply rules for building sources it contributes
boinc/boinc/samples/multi_thread/%.o: ../boinc/boinc/samples/multi_thread/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


