################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../boinc/boinc/samples/worker/worker.cpp 

OBJS += \
./boinc/boinc/samples/worker/worker.o 

CPP_DEPS += \
./boinc/boinc/samples/worker/worker.d 


# Each subdirectory must supply rules for building sources it contributes
boinc/boinc/samples/worker/%.o: ../boinc/boinc/samples/worker/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


