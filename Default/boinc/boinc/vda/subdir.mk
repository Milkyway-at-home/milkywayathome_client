################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../boinc/boinc/vda/sched_vda.cpp \
../boinc/boinc/vda/ssim.cpp \
../boinc/boinc/vda/stats.cpp \
../boinc/boinc/vda/vda.cpp \
../boinc/boinc/vda/vda_lib.cpp \
../boinc/boinc/vda/vda_lib2.cpp \
../boinc/boinc/vda/vda_policy.cpp \
../boinc/boinc/vda/vdad.cpp 

OBJS += \
./boinc/boinc/vda/sched_vda.o \
./boinc/boinc/vda/ssim.o \
./boinc/boinc/vda/stats.o \
./boinc/boinc/vda/vda.o \
./boinc/boinc/vda/vda_lib.o \
./boinc/boinc/vda/vda_lib2.o \
./boinc/boinc/vda/vda_policy.o \
./boinc/boinc/vda/vdad.o 

CPP_DEPS += \
./boinc/boinc/vda/sched_vda.d \
./boinc/boinc/vda/ssim.d \
./boinc/boinc/vda/stats.d \
./boinc/boinc/vda/vda.d \
./boinc/boinc/vda/vda_lib.d \
./boinc/boinc/vda/vda_lib2.d \
./boinc/boinc/vda/vda_policy.d \
./boinc/boinc/vda/vdad.d 


# Each subdirectory must supply rules for building sources it contributes
boinc/boinc/vda/%.o: ../boinc/boinc/vda/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


