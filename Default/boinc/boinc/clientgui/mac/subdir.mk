################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../boinc/boinc/clientgui/mac/MacBitmapComboBox.cpp \
../boinc/boinc/clientgui/mac/Mac_GUI.cpp \
../boinc/boinc/clientgui/mac/SecurityUtility.cpp \
../boinc/boinc/clientgui/mac/SetVersion.cpp \
../boinc/boinc/clientgui/mac/SetupSecurity.cpp 

OBJS += \
./boinc/boinc/clientgui/mac/MacBitmapComboBox.o \
./boinc/boinc/clientgui/mac/Mac_GUI.o \
./boinc/boinc/clientgui/mac/SecurityUtility.o \
./boinc/boinc/clientgui/mac/SetVersion.o \
./boinc/boinc/clientgui/mac/SetupSecurity.o 

CPP_DEPS += \
./boinc/boinc/clientgui/mac/MacBitmapComboBox.d \
./boinc/boinc/clientgui/mac/Mac_GUI.d \
./boinc/boinc/clientgui/mac/SecurityUtility.d \
./boinc/boinc/clientgui/mac/SetVersion.d \
./boinc/boinc/clientgui/mac/SetupSecurity.d 


# Each subdirectory must supply rules for building sources it contributes
boinc/boinc/clientgui/mac/%.o: ../boinc/boinc/clientgui/mac/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


