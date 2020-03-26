################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../boinc/boinc/zip/boinc_zip.cpp \
../boinc/boinc/zip/test.cpp \
../boinc/boinc/zip/testzlibconflict.cpp 

OBJS += \
./boinc/boinc/zip/boinc_zip.o \
./boinc/boinc/zip/test.o \
./boinc/boinc/zip/testzlibconflict.o 

CPP_DEPS += \
./boinc/boinc/zip/boinc_zip.d \
./boinc/boinc/zip/test.d \
./boinc/boinc/zip/testzlibconflict.d 


# Each subdirectory must supply rules for building sources it contributes
boinc/boinc/zip/%.o: ../boinc/boinc/zip/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


