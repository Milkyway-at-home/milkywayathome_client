################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../boinc/boinc/tools/backend_lib.cpp \
../boinc/boinc/tools/cancel_jobs.cpp \
../boinc/boinc/tools/create_work.cpp \
../boinc/boinc/tools/dir_hier_move.cpp \
../boinc/boinc/tools/dir_hier_path.cpp \
../boinc/boinc/tools/hr_db_convert.cpp \
../boinc/boinc/tools/kill_wu.cpp \
../boinc/boinc/tools/poll_wu.cpp \
../boinc/boinc/tools/process_input_template.cpp \
../boinc/boinc/tools/process_result_template.cpp \
../boinc/boinc/tools/remote_submit_test.cpp \
../boinc/boinc/tools/sign_executable.cpp \
../boinc/boinc/tools/updater.cpp 

OBJS += \
./boinc/boinc/tools/backend_lib.o \
./boinc/boinc/tools/cancel_jobs.o \
./boinc/boinc/tools/create_work.o \
./boinc/boinc/tools/dir_hier_move.o \
./boinc/boinc/tools/dir_hier_path.o \
./boinc/boinc/tools/hr_db_convert.o \
./boinc/boinc/tools/kill_wu.o \
./boinc/boinc/tools/poll_wu.o \
./boinc/boinc/tools/process_input_template.o \
./boinc/boinc/tools/process_result_template.o \
./boinc/boinc/tools/remote_submit_test.o \
./boinc/boinc/tools/sign_executable.o \
./boinc/boinc/tools/updater.o 

CPP_DEPS += \
./boinc/boinc/tools/backend_lib.d \
./boinc/boinc/tools/cancel_jobs.d \
./boinc/boinc/tools/create_work.d \
./boinc/boinc/tools/dir_hier_move.d \
./boinc/boinc/tools/dir_hier_path.d \
./boinc/boinc/tools/hr_db_convert.d \
./boinc/boinc/tools/kill_wu.d \
./boinc/boinc/tools/poll_wu.d \
./boinc/boinc/tools/process_input_template.d \
./boinc/boinc/tools/process_result_template.d \
./boinc/boinc/tools/remote_submit_test.d \
./boinc/boinc/tools/sign_executable.d \
./boinc/boinc/tools/updater.d 


# Each subdirectory must supply rules for building sources it contributes
boinc/boinc/tools/%.o: ../boinc/boinc/tools/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


