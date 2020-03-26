################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../boinc/boinc/samples/wrapper/wrapper.cpp 

C_SRCS += \
../boinc/boinc/samples/wrapper/regerror.c \
../boinc/boinc/samples/wrapper/regexp.c \
../boinc/boinc/samples/wrapper/regexp_memory.c \
../boinc/boinc/samples/wrapper/regexp_report.c \
../boinc/boinc/samples/wrapper/regsub.c 

OBJS += \
./boinc/boinc/samples/wrapper/regerror.o \
./boinc/boinc/samples/wrapper/regexp.o \
./boinc/boinc/samples/wrapper/regexp_memory.o \
./boinc/boinc/samples/wrapper/regexp_report.o \
./boinc/boinc/samples/wrapper/regsub.o \
./boinc/boinc/samples/wrapper/wrapper.o 

CPP_DEPS += \
./boinc/boinc/samples/wrapper/wrapper.d 

C_DEPS += \
./boinc/boinc/samples/wrapper/regerror.d \
./boinc/boinc/samples/wrapper/regexp.d \
./boinc/boinc/samples/wrapper/regexp_memory.d \
./boinc/boinc/samples/wrapper/regexp_report.d \
./boinc/boinc/samples/wrapper/regsub.d 


# Each subdirectory must supply rules for building sources it contributes
boinc/boinc/samples/wrapper/%.o: ../boinc/boinc/samples/wrapper/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

boinc/boinc/samples/wrapper/%.o: ../boinc/boinc/samples/wrapper/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


