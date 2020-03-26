################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../boinc/boinc/mac_installer/AddRemoveUser.cpp \
../boinc/boinc/mac_installer/CustomInstall.cpp \
../boinc/boinc/mac_installer/Installer.cpp \
../boinc/boinc/mac_installer/PostInstall.cpp \
../boinc/boinc/mac_installer/WaitPermissions.cpp \
../boinc/boinc/mac_installer/uninstall.cpp 

OBJS += \
./boinc/boinc/mac_installer/AddRemoveUser.o \
./boinc/boinc/mac_installer/CustomInstall.o \
./boinc/boinc/mac_installer/Installer.o \
./boinc/boinc/mac_installer/PostInstall.o \
./boinc/boinc/mac_installer/WaitPermissions.o \
./boinc/boinc/mac_installer/uninstall.o 

CPP_DEPS += \
./boinc/boinc/mac_installer/AddRemoveUser.d \
./boinc/boinc/mac_installer/CustomInstall.d \
./boinc/boinc/mac_installer/Installer.d \
./boinc/boinc/mac_installer/PostInstall.d \
./boinc/boinc/mac_installer/WaitPermissions.d \
./boinc/boinc/mac_installer/uninstall.d 


# Each subdirectory must supply rules for building sources it contributes
boinc/boinc/mac_installer/%.o: ../boinc/boinc/mac_installer/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


