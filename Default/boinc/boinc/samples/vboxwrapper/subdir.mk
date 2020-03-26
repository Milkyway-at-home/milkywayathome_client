################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../boinc/boinc/samples/vboxwrapper/floppyio.cpp \
../boinc/boinc/samples/vboxwrapper/gbac.cpp \
../boinc/boinc/samples/vboxwrapper/vbox_common.cpp \
../boinc/boinc/samples/vboxwrapper/vbox_mscom42.cpp \
../boinc/boinc/samples/vboxwrapper/vbox_mscom43.cpp \
../boinc/boinc/samples/vboxwrapper/vbox_mscom50.cpp \
../boinc/boinc/samples/vboxwrapper/vbox_mscom51.cpp \
../boinc/boinc/samples/vboxwrapper/vbox_mscom_impl.cpp \
../boinc/boinc/samples/vboxwrapper/vbox_vboxmanage.cpp \
../boinc/boinc/samples/vboxwrapper/vboxcheckpoint.cpp \
../boinc/boinc/samples/vboxwrapper/vboxjob.cpp \
../boinc/boinc/samples/vboxwrapper/vboxlogging.cpp \
../boinc/boinc/samples/vboxwrapper/vboxwrapper.cpp 

OBJS += \
./boinc/boinc/samples/vboxwrapper/floppyio.o \
./boinc/boinc/samples/vboxwrapper/gbac.o \
./boinc/boinc/samples/vboxwrapper/vbox_common.o \
./boinc/boinc/samples/vboxwrapper/vbox_mscom42.o \
./boinc/boinc/samples/vboxwrapper/vbox_mscom43.o \
./boinc/boinc/samples/vboxwrapper/vbox_mscom50.o \
./boinc/boinc/samples/vboxwrapper/vbox_mscom51.o \
./boinc/boinc/samples/vboxwrapper/vbox_mscom_impl.o \
./boinc/boinc/samples/vboxwrapper/vbox_vboxmanage.o \
./boinc/boinc/samples/vboxwrapper/vboxcheckpoint.o \
./boinc/boinc/samples/vboxwrapper/vboxjob.o \
./boinc/boinc/samples/vboxwrapper/vboxlogging.o \
./boinc/boinc/samples/vboxwrapper/vboxwrapper.o 

CPP_DEPS += \
./boinc/boinc/samples/vboxwrapper/floppyio.d \
./boinc/boinc/samples/vboxwrapper/gbac.d \
./boinc/boinc/samples/vboxwrapper/vbox_common.d \
./boinc/boinc/samples/vboxwrapper/vbox_mscom42.d \
./boinc/boinc/samples/vboxwrapper/vbox_mscom43.d \
./boinc/boinc/samples/vboxwrapper/vbox_mscom50.d \
./boinc/boinc/samples/vboxwrapper/vbox_mscom51.d \
./boinc/boinc/samples/vboxwrapper/vbox_mscom_impl.d \
./boinc/boinc/samples/vboxwrapper/vbox_vboxmanage.d \
./boinc/boinc/samples/vboxwrapper/vboxcheckpoint.d \
./boinc/boinc/samples/vboxwrapper/vboxjob.d \
./boinc/boinc/samples/vboxwrapper/vboxlogging.d \
./boinc/boinc/samples/vboxwrapper/vboxwrapper.d 


# Each subdirectory must supply rules for building sources it contributes
boinc/boinc/samples/vboxwrapper/%.o: ../boinc/boinc/samples/vboxwrapper/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


