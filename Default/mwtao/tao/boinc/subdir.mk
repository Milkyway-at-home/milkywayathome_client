################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CXX_SRCS += \
../mwtao/tao/boinc/messages.cxx \
../mwtao/tao/boinc/tao_assimilation_policy.cxx \
../mwtao/tao/boinc/tao_search_status.cxx \
../mwtao/tao/boinc/tao_stop_search.cxx \
../mwtao/tao/boinc/tao_validation_policy.cxx \
../mwtao/tao/boinc/tao_work_generator.cxx \
../mwtao/tao/boinc/workunit_information.cxx 

CXX_DEPS += \
./mwtao/tao/boinc/messages.d \
./mwtao/tao/boinc/tao_assimilation_policy.d \
./mwtao/tao/boinc/tao_search_status.d \
./mwtao/tao/boinc/tao_stop_search.d \
./mwtao/tao/boinc/tao_validation_policy.d \
./mwtao/tao/boinc/tao_work_generator.d \
./mwtao/tao/boinc/workunit_information.d 

OBJS += \
./mwtao/tao/boinc/messages.o \
./mwtao/tao/boinc/tao_assimilation_policy.o \
./mwtao/tao/boinc/tao_search_status.o \
./mwtao/tao/boinc/tao_stop_search.o \
./mwtao/tao/boinc/tao_validation_policy.o \
./mwtao/tao/boinc/tao_work_generator.o \
./mwtao/tao/boinc/workunit_information.o 


# Each subdirectory must supply rules for building sources it contributes
mwtao/tao/boinc/%.o: ../mwtao/tao/boinc/%.cxx
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


