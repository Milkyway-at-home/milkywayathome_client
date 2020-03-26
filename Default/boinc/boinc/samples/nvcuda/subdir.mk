################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../boinc/boinc/samples/nvcuda/cuda.cpp 

OBJS += \
./boinc/boinc/samples/nvcuda/cuda.o 

CPP_DEPS += \
./boinc/boinc/samples/nvcuda/cuda.d 


# Each subdirectory must supply rules for building sources it contributes
boinc/boinc/samples/nvcuda/%.o: ../boinc/boinc/samples/nvcuda/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


