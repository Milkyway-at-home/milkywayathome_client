################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../boinc/boinc/samples/sleeper/sleeper.cpp 

OBJS += \
./boinc/boinc/samples/sleeper/sleeper.o 

CPP_DEPS += \
./boinc/boinc/samples/sleeper/sleeper.d 


# Each subdirectory must supply rules for building sources it contributes
boinc/boinc/samples/sleeper/%.o: ../boinc/boinc/samples/sleeper/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


