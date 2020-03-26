################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../boinc/boinc/db/boinc_db.cpp \
../boinc/boinc/db/db_base.cpp 

OBJS += \
./boinc/boinc/db/boinc_db.o \
./boinc/boinc/db/db_base.o 

CPP_DEPS += \
./boinc/boinc/db/boinc_db.d \
./boinc/boinc/db/db_base.d 


# Each subdirectory must supply rules for building sources it contributes
boinc/boinc/db/%.o: ../boinc/boinc/db/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


