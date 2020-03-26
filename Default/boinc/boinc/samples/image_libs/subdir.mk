################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../boinc/boinc/samples/image_libs/bmplib.cpp \
../boinc/boinc/samples/image_libs/tgalib.cpp 

OBJS += \
./boinc/boinc/samples/image_libs/bmplib.o \
./boinc/boinc/samples/image_libs/tgalib.o 

CPP_DEPS += \
./boinc/boinc/samples/image_libs/bmplib.d \
./boinc/boinc/samples/image_libs/tgalib.d 


# Each subdirectory must supply rules for building sources it contributes
boinc/boinc/samples/image_libs/%.o: ../boinc/boinc/samples/image_libs/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


