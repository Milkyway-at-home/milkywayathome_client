################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../boinc/boinc/clientscr/gfx_switcher.cpp \
../boinc/boinc/clientscr/mac_saver_module.cpp \
../boinc/boinc/clientscr/screensaver.cpp \
../boinc/boinc/clientscr/screensaver_win.cpp \
../boinc/boinc/clientscr/screensaver_x11.cpp \
../boinc/boinc/clientscr/ss_app.cpp 

OBJS += \
./boinc/boinc/clientscr/gfx_switcher.o \
./boinc/boinc/clientscr/mac_saver_module.o \
./boinc/boinc/clientscr/screensaver.o \
./boinc/boinc/clientscr/screensaver_win.o \
./boinc/boinc/clientscr/screensaver_x11.o \
./boinc/boinc/clientscr/ss_app.o 

CPP_DEPS += \
./boinc/boinc/clientscr/gfx_switcher.d \
./boinc/boinc/clientscr/mac_saver_module.d \
./boinc/boinc/clientscr/screensaver.d \
./boinc/boinc/clientscr/screensaver_win.d \
./boinc/boinc/clientscr/screensaver_x11.d \
./boinc/boinc/clientscr/ss_app.d 


# Each subdirectory must supply rules for building sources it contributes
boinc/boinc/clientscr/%.o: ../boinc/boinc/clientscr/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


