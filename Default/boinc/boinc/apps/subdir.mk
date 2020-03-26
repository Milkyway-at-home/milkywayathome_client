################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../boinc/boinc/apps/1sec.cpp \
../boinc/boinc/apps/concat.cpp \
../boinc/boinc/apps/error.cpp \
../boinc/boinc/apps/upper_case.cpp 

OBJS += \
./boinc/boinc/apps/1sec.o \
./boinc/boinc/apps/concat.o \
./boinc/boinc/apps/error.o \
./boinc/boinc/apps/upper_case.o 

CPP_DEPS += \
./boinc/boinc/apps/1sec.d \
./boinc/boinc/apps/concat.d \
./boinc/boinc/apps/error.d \
./boinc/boinc/apps/upper_case.d 


# Each subdirectory must supply rules for building sources it contributes
boinc/boinc/apps/%.o: ../boinc/boinc/apps/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


