################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../boinc/boinc/samples/wrappture/fermi.cpp \
../boinc/boinc/samples/wrappture/wrappture.cpp \
../boinc/boinc/samples/wrappture/wrappture_example.cpp 

OBJS += \
./boinc/boinc/samples/wrappture/fermi.o \
./boinc/boinc/samples/wrappture/wrappture.o \
./boinc/boinc/samples/wrappture/wrappture_example.o 

CPP_DEPS += \
./boinc/boinc/samples/wrappture/fermi.d \
./boinc/boinc/samples/wrappture/wrappture.d \
./boinc/boinc/samples/wrappture/wrappture_example.d 


# Each subdirectory must supply rules for building sources it contributes
boinc/boinc/samples/wrappture/%.o: ../boinc/boinc/samples/wrappture/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


