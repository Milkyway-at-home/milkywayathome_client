################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../boinc/boinc/samples/example_app/slide_show.cpp \
../boinc/boinc/samples/example_app/uc2.cpp \
../boinc/boinc/samples/example_app/uc2_dll.cpp \
../boinc/boinc/samples/example_app/uc2_graphics.cpp \
../boinc/boinc/samples/example_app/ucn.cpp 

OBJS += \
./boinc/boinc/samples/example_app/slide_show.o \
./boinc/boinc/samples/example_app/uc2.o \
./boinc/boinc/samples/example_app/uc2_dll.o \
./boinc/boinc/samples/example_app/uc2_graphics.o \
./boinc/boinc/samples/example_app/ucn.o 

CPP_DEPS += \
./boinc/boinc/samples/example_app/slide_show.d \
./boinc/boinc/samples/example_app/uc2.d \
./boinc/boinc/samples/example_app/uc2_dll.d \
./boinc/boinc/samples/example_app/uc2_graphics.d \
./boinc/boinc/samples/example_app/ucn.d 


# Each subdirectory must supply rules for building sources it contributes
boinc/boinc/samples/example_app/%.o: ../boinc/boinc/samples/example_app/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


