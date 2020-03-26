################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../boinc/boinc/sched/adjust_user_priority.cpp \
../boinc/boinc/sched/antique_file_deleter.cpp \
../boinc/boinc/sched/assimilator.cpp \
../boinc/boinc/sched/census.cpp \
../boinc/boinc/sched/credit.cpp \
../boinc/boinc/sched/credit_test.cpp \
../boinc/boinc/sched/db_dump.cpp \
../boinc/boinc/sched/db_purge.cpp \
../boinc/boinc/sched/delete_file.cpp \
../boinc/boinc/sched/edf_sim.cpp \
../boinc/boinc/sched/feeder.cpp \
../boinc/boinc/sched/file_deleter.cpp \
../boinc/boinc/sched/file_upload_handler.cpp \
../boinc/boinc/sched/get_file.cpp \
../boinc/boinc/sched/handle_request.cpp \
../boinc/boinc/sched/hr.cpp \
../boinc/boinc/sched/hr_info.cpp \
../boinc/boinc/sched/make_work.cpp \
../boinc/boinc/sched/message_handler.cpp \
../boinc/boinc/sched/plan_class_spec.cpp \
../boinc/boinc/sched/put_file.cpp \
../boinc/boinc/sched/sample_assimilator.cpp \
../boinc/boinc/sched/sample_bitwise_validator.cpp \
../boinc/boinc/sched/sample_dummy_assimilator.cpp \
../boinc/boinc/sched/sample_substr_validator.cpp \
../boinc/boinc/sched/sample_trivial_validator.cpp \
../boinc/boinc/sched/sample_work_generator.cpp \
../boinc/boinc/sched/sched_array.cpp \
../boinc/boinc/sched/sched_assign.cpp \
../boinc/boinc/sched/sched_check.cpp \
../boinc/boinc/sched/sched_config.cpp \
../boinc/boinc/sched/sched_customize.cpp \
../boinc/boinc/sched/sched_driver.cpp \
../boinc/boinc/sched/sched_files.cpp \
../boinc/boinc/sched/sched_hr.cpp \
../boinc/boinc/sched/sched_limit.cpp \
../boinc/boinc/sched/sched_locality.cpp \
../boinc/boinc/sched/sched_main.cpp \
../boinc/boinc/sched/sched_nci.cpp \
../boinc/boinc/sched/sched_resend.cpp \
../boinc/boinc/sched/sched_result.cpp \
../boinc/boinc/sched/sched_score.cpp \
../boinc/boinc/sched/sched_send.cpp \
../boinc/boinc/sched/sched_shmem.cpp \
../boinc/boinc/sched/sched_timezone.cpp \
../boinc/boinc/sched/sched_types.cpp \
../boinc/boinc/sched/sched_util.cpp \
../boinc/boinc/sched/sched_util_basic.cpp \
../boinc/boinc/sched/sched_version.cpp \
../boinc/boinc/sched/script_assimilator.cpp \
../boinc/boinc/sched/script_validator.cpp \
../boinc/boinc/sched/show_shmem.cpp \
../boinc/boinc/sched/single_job_assimilator.cpp \
../boinc/boinc/sched/size_regulator.cpp \
../boinc/boinc/sched/target_batch.cpp \
../boinc/boinc/sched/time_stats_log.cpp \
../boinc/boinc/sched/transitioner.cpp \
../boinc/boinc/sched/trickle_credit.cpp \
../boinc/boinc/sched/trickle_deadline.cpp \
../boinc/boinc/sched/trickle_echo.cpp \
../boinc/boinc/sched/trickle_handler.cpp \
../boinc/boinc/sched/update_stats.cpp \
../boinc/boinc/sched/validate_util.cpp \
../boinc/boinc/sched/validate_util2.cpp \
../boinc/boinc/sched/validator.cpp \
../boinc/boinc/sched/validator_test.cpp \
../boinc/boinc/sched/wu_check.cpp 

OBJS += \
./boinc/boinc/sched/adjust_user_priority.o \
./boinc/boinc/sched/antique_file_deleter.o \
./boinc/boinc/sched/assimilator.o \
./boinc/boinc/sched/census.o \
./boinc/boinc/sched/credit.o \
./boinc/boinc/sched/credit_test.o \
./boinc/boinc/sched/db_dump.o \
./boinc/boinc/sched/db_purge.o \
./boinc/boinc/sched/delete_file.o \
./boinc/boinc/sched/edf_sim.o \
./boinc/boinc/sched/feeder.o \
./boinc/boinc/sched/file_deleter.o \
./boinc/boinc/sched/file_upload_handler.o \
./boinc/boinc/sched/get_file.o \
./boinc/boinc/sched/handle_request.o \
./boinc/boinc/sched/hr.o \
./boinc/boinc/sched/hr_info.o \
./boinc/boinc/sched/make_work.o \
./boinc/boinc/sched/message_handler.o \
./boinc/boinc/sched/plan_class_spec.o \
./boinc/boinc/sched/put_file.o \
./boinc/boinc/sched/sample_assimilator.o \
./boinc/boinc/sched/sample_bitwise_validator.o \
./boinc/boinc/sched/sample_dummy_assimilator.o \
./boinc/boinc/sched/sample_substr_validator.o \
./boinc/boinc/sched/sample_trivial_validator.o \
./boinc/boinc/sched/sample_work_generator.o \
./boinc/boinc/sched/sched_array.o \
./boinc/boinc/sched/sched_assign.o \
./boinc/boinc/sched/sched_check.o \
./boinc/boinc/sched/sched_config.o \
./boinc/boinc/sched/sched_customize.o \
./boinc/boinc/sched/sched_driver.o \
./boinc/boinc/sched/sched_files.o \
./boinc/boinc/sched/sched_hr.o \
./boinc/boinc/sched/sched_limit.o \
./boinc/boinc/sched/sched_locality.o \
./boinc/boinc/sched/sched_main.o \
./boinc/boinc/sched/sched_nci.o \
./boinc/boinc/sched/sched_resend.o \
./boinc/boinc/sched/sched_result.o \
./boinc/boinc/sched/sched_score.o \
./boinc/boinc/sched/sched_send.o \
./boinc/boinc/sched/sched_shmem.o \
./boinc/boinc/sched/sched_timezone.o \
./boinc/boinc/sched/sched_types.o \
./boinc/boinc/sched/sched_util.o \
./boinc/boinc/sched/sched_util_basic.o \
./boinc/boinc/sched/sched_version.o \
./boinc/boinc/sched/script_assimilator.o \
./boinc/boinc/sched/script_validator.o \
./boinc/boinc/sched/show_shmem.o \
./boinc/boinc/sched/single_job_assimilator.o \
./boinc/boinc/sched/size_regulator.o \
./boinc/boinc/sched/target_batch.o \
./boinc/boinc/sched/time_stats_log.o \
./boinc/boinc/sched/transitioner.o \
./boinc/boinc/sched/trickle_credit.o \
./boinc/boinc/sched/trickle_deadline.o \
./boinc/boinc/sched/trickle_echo.o \
./boinc/boinc/sched/trickle_handler.o \
./boinc/boinc/sched/update_stats.o \
./boinc/boinc/sched/validate_util.o \
./boinc/boinc/sched/validate_util2.o \
./boinc/boinc/sched/validator.o \
./boinc/boinc/sched/validator_test.o \
./boinc/boinc/sched/wu_check.o 

CPP_DEPS += \
./boinc/boinc/sched/adjust_user_priority.d \
./boinc/boinc/sched/antique_file_deleter.d \
./boinc/boinc/sched/assimilator.d \
./boinc/boinc/sched/census.d \
./boinc/boinc/sched/credit.d \
./boinc/boinc/sched/credit_test.d \
./boinc/boinc/sched/db_dump.d \
./boinc/boinc/sched/db_purge.d \
./boinc/boinc/sched/delete_file.d \
./boinc/boinc/sched/edf_sim.d \
./boinc/boinc/sched/feeder.d \
./boinc/boinc/sched/file_deleter.d \
./boinc/boinc/sched/file_upload_handler.d \
./boinc/boinc/sched/get_file.d \
./boinc/boinc/sched/handle_request.d \
./boinc/boinc/sched/hr.d \
./boinc/boinc/sched/hr_info.d \
./boinc/boinc/sched/make_work.d \
./boinc/boinc/sched/message_handler.d \
./boinc/boinc/sched/plan_class_spec.d \
./boinc/boinc/sched/put_file.d \
./boinc/boinc/sched/sample_assimilator.d \
./boinc/boinc/sched/sample_bitwise_validator.d \
./boinc/boinc/sched/sample_dummy_assimilator.d \
./boinc/boinc/sched/sample_substr_validator.d \
./boinc/boinc/sched/sample_trivial_validator.d \
./boinc/boinc/sched/sample_work_generator.d \
./boinc/boinc/sched/sched_array.d \
./boinc/boinc/sched/sched_assign.d \
./boinc/boinc/sched/sched_check.d \
./boinc/boinc/sched/sched_config.d \
./boinc/boinc/sched/sched_customize.d \
./boinc/boinc/sched/sched_driver.d \
./boinc/boinc/sched/sched_files.d \
./boinc/boinc/sched/sched_hr.d \
./boinc/boinc/sched/sched_limit.d \
./boinc/boinc/sched/sched_locality.d \
./boinc/boinc/sched/sched_main.d \
./boinc/boinc/sched/sched_nci.d \
./boinc/boinc/sched/sched_resend.d \
./boinc/boinc/sched/sched_result.d \
./boinc/boinc/sched/sched_score.d \
./boinc/boinc/sched/sched_send.d \
./boinc/boinc/sched/sched_shmem.d \
./boinc/boinc/sched/sched_timezone.d \
./boinc/boinc/sched/sched_types.d \
./boinc/boinc/sched/sched_util.d \
./boinc/boinc/sched/sched_util_basic.d \
./boinc/boinc/sched/sched_version.d \
./boinc/boinc/sched/script_assimilator.d \
./boinc/boinc/sched/script_validator.d \
./boinc/boinc/sched/show_shmem.d \
./boinc/boinc/sched/single_job_assimilator.d \
./boinc/boinc/sched/size_regulator.d \
./boinc/boinc/sched/target_batch.d \
./boinc/boinc/sched/time_stats_log.d \
./boinc/boinc/sched/transitioner.d \
./boinc/boinc/sched/trickle_credit.d \
./boinc/boinc/sched/trickle_deadline.d \
./boinc/boinc/sched/trickle_echo.d \
./boinc/boinc/sched/trickle_handler.d \
./boinc/boinc/sched/update_stats.d \
./boinc/boinc/sched/validate_util.d \
./boinc/boinc/sched/validate_util2.d \
./boinc/boinc/sched/validator.d \
./boinc/boinc/sched/validator_test.d \
./boinc/boinc/sched/wu_check.d 


# Each subdirectory must supply rules for building sources it contributes
boinc/boinc/sched/%.o: ../boinc/boinc/sched/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


