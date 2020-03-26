################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../boinc/boinc/client/acct_mgr.cpp \
../boinc/boinc/client/acct_setup.cpp \
../boinc/boinc/client/app.cpp \
../boinc/boinc/client/app_config.cpp \
../boinc/boinc/client/app_control.cpp \
../boinc/boinc/client/app_graphics.cpp \
../boinc/boinc/client/app_start.cpp \
../boinc/boinc/client/async_file.cpp \
../boinc/boinc/client/auto_update.cpp \
../boinc/boinc/client/boinc_cmd.cpp \
../boinc/boinc/client/boinc_log.cpp \
../boinc/boinc/client/check_security.cpp \
../boinc/boinc/client/check_state.cpp \
../boinc/boinc/client/client_msgs.cpp \
../boinc/boinc/client/client_state.cpp \
../boinc/boinc/client/client_types.cpp \
../boinc/boinc/client/coproc_sched.cpp \
../boinc/boinc/client/cpu_sched.cpp \
../boinc/boinc/client/cs_account.cpp \
../boinc/boinc/client/cs_apps.cpp \
../boinc/boinc/client/cs_benchmark.cpp \
../boinc/boinc/client/cs_cmdline.cpp \
../boinc/boinc/client/cs_files.cpp \
../boinc/boinc/client/cs_notice.cpp \
../boinc/boinc/client/cs_platforms.cpp \
../boinc/boinc/client/cs_prefs.cpp \
../boinc/boinc/client/cs_proxy.cpp \
../boinc/boinc/client/cs_scheduler.cpp \
../boinc/boinc/client/cs_statefile.cpp \
../boinc/boinc/client/cs_trickle.cpp \
../boinc/boinc/client/current_version.cpp \
../boinc/boinc/client/dhrystone.cpp \
../boinc/boinc/client/dhrystone2.cpp \
../boinc/boinc/client/file_names.cpp \
../boinc/boinc/client/file_xfer.cpp \
../boinc/boinc/client/gpu_amd.cpp \
../boinc/boinc/client/gpu_detect.cpp \
../boinc/boinc/client/gpu_intel.cpp \
../boinc/boinc/client/gpu_nvidia.cpp \
../boinc/boinc/client/gpu_opencl.cpp \
../boinc/boinc/client/gui_http.cpp \
../boinc/boinc/client/gui_rpc_server.cpp \
../boinc/boinc/client/gui_rpc_server_ops.cpp \
../boinc/boinc/client/hostinfo_network.cpp \
../boinc/boinc/client/hostinfo_unix.cpp \
../boinc/boinc/client/hostinfo_unix_test.cpp \
../boinc/boinc/client/hostinfo_win.cpp \
../boinc/boinc/client/http_curl.cpp \
../boinc/boinc/client/log_flags.cpp \
../boinc/boinc/client/mac_address.cpp \
../boinc/boinc/client/main.cpp \
../boinc/boinc/client/net_stats.cpp \
../boinc/boinc/client/pers_file_xfer.cpp \
../boinc/boinc/client/project.cpp \
../boinc/boinc/client/result.cpp \
../boinc/boinc/client/rr_sim.cpp \
../boinc/boinc/client/rrsim_test.cpp \
../boinc/boinc/client/sandbox.cpp \
../boinc/boinc/client/scheduler_op.cpp \
../boinc/boinc/client/setprojectgrp.cpp \
../boinc/boinc/client/sim.cpp \
../boinc/boinc/client/sim_util.cpp \
../boinc/boinc/client/switcher.cpp \
../boinc/boinc/client/sysmon_win.cpp \
../boinc/boinc/client/thread.cpp \
../boinc/boinc/client/time_stats.cpp \
../boinc/boinc/client/whetstone.cpp \
../boinc/boinc/client/work_fetch.cpp 

OBJS += \
./boinc/boinc/client/acct_mgr.o \
./boinc/boinc/client/acct_setup.o \
./boinc/boinc/client/app.o \
./boinc/boinc/client/app_config.o \
./boinc/boinc/client/app_control.o \
./boinc/boinc/client/app_graphics.o \
./boinc/boinc/client/app_start.o \
./boinc/boinc/client/async_file.o \
./boinc/boinc/client/auto_update.o \
./boinc/boinc/client/boinc_cmd.o \
./boinc/boinc/client/boinc_log.o \
./boinc/boinc/client/check_security.o \
./boinc/boinc/client/check_state.o \
./boinc/boinc/client/client_msgs.o \
./boinc/boinc/client/client_state.o \
./boinc/boinc/client/client_types.o \
./boinc/boinc/client/coproc_sched.o \
./boinc/boinc/client/cpu_sched.o \
./boinc/boinc/client/cs_account.o \
./boinc/boinc/client/cs_apps.o \
./boinc/boinc/client/cs_benchmark.o \
./boinc/boinc/client/cs_cmdline.o \
./boinc/boinc/client/cs_files.o \
./boinc/boinc/client/cs_notice.o \
./boinc/boinc/client/cs_platforms.o \
./boinc/boinc/client/cs_prefs.o \
./boinc/boinc/client/cs_proxy.o \
./boinc/boinc/client/cs_scheduler.o \
./boinc/boinc/client/cs_statefile.o \
./boinc/boinc/client/cs_trickle.o \
./boinc/boinc/client/current_version.o \
./boinc/boinc/client/dhrystone.o \
./boinc/boinc/client/dhrystone2.o \
./boinc/boinc/client/file_names.o \
./boinc/boinc/client/file_xfer.o \
./boinc/boinc/client/gpu_amd.o \
./boinc/boinc/client/gpu_detect.o \
./boinc/boinc/client/gpu_intel.o \
./boinc/boinc/client/gpu_nvidia.o \
./boinc/boinc/client/gpu_opencl.o \
./boinc/boinc/client/gui_http.o \
./boinc/boinc/client/gui_rpc_server.o \
./boinc/boinc/client/gui_rpc_server_ops.o \
./boinc/boinc/client/hostinfo_network.o \
./boinc/boinc/client/hostinfo_unix.o \
./boinc/boinc/client/hostinfo_unix_test.o \
./boinc/boinc/client/hostinfo_win.o \
./boinc/boinc/client/http_curl.o \
./boinc/boinc/client/log_flags.o \
./boinc/boinc/client/mac_address.o \
./boinc/boinc/client/main.o \
./boinc/boinc/client/net_stats.o \
./boinc/boinc/client/pers_file_xfer.o \
./boinc/boinc/client/project.o \
./boinc/boinc/client/result.o \
./boinc/boinc/client/rr_sim.o \
./boinc/boinc/client/rrsim_test.o \
./boinc/boinc/client/sandbox.o \
./boinc/boinc/client/scheduler_op.o \
./boinc/boinc/client/setprojectgrp.o \
./boinc/boinc/client/sim.o \
./boinc/boinc/client/sim_util.o \
./boinc/boinc/client/switcher.o \
./boinc/boinc/client/sysmon_win.o \
./boinc/boinc/client/thread.o \
./boinc/boinc/client/time_stats.o \
./boinc/boinc/client/whetstone.o \
./boinc/boinc/client/work_fetch.o 

CPP_DEPS += \
./boinc/boinc/client/acct_mgr.d \
./boinc/boinc/client/acct_setup.d \
./boinc/boinc/client/app.d \
./boinc/boinc/client/app_config.d \
./boinc/boinc/client/app_control.d \
./boinc/boinc/client/app_graphics.d \
./boinc/boinc/client/app_start.d \
./boinc/boinc/client/async_file.d \
./boinc/boinc/client/auto_update.d \
./boinc/boinc/client/boinc_cmd.d \
./boinc/boinc/client/boinc_log.d \
./boinc/boinc/client/check_security.d \
./boinc/boinc/client/check_state.d \
./boinc/boinc/client/client_msgs.d \
./boinc/boinc/client/client_state.d \
./boinc/boinc/client/client_types.d \
./boinc/boinc/client/coproc_sched.d \
./boinc/boinc/client/cpu_sched.d \
./boinc/boinc/client/cs_account.d \
./boinc/boinc/client/cs_apps.d \
./boinc/boinc/client/cs_benchmark.d \
./boinc/boinc/client/cs_cmdline.d \
./boinc/boinc/client/cs_files.d \
./boinc/boinc/client/cs_notice.d \
./boinc/boinc/client/cs_platforms.d \
./boinc/boinc/client/cs_prefs.d \
./boinc/boinc/client/cs_proxy.d \
./boinc/boinc/client/cs_scheduler.d \
./boinc/boinc/client/cs_statefile.d \
./boinc/boinc/client/cs_trickle.d \
./boinc/boinc/client/current_version.d \
./boinc/boinc/client/dhrystone.d \
./boinc/boinc/client/dhrystone2.d \
./boinc/boinc/client/file_names.d \
./boinc/boinc/client/file_xfer.d \
./boinc/boinc/client/gpu_amd.d \
./boinc/boinc/client/gpu_detect.d \
./boinc/boinc/client/gpu_intel.d \
./boinc/boinc/client/gpu_nvidia.d \
./boinc/boinc/client/gpu_opencl.d \
./boinc/boinc/client/gui_http.d \
./boinc/boinc/client/gui_rpc_server.d \
./boinc/boinc/client/gui_rpc_server_ops.d \
./boinc/boinc/client/hostinfo_network.d \
./boinc/boinc/client/hostinfo_unix.d \
./boinc/boinc/client/hostinfo_unix_test.d \
./boinc/boinc/client/hostinfo_win.d \
./boinc/boinc/client/http_curl.d \
./boinc/boinc/client/log_flags.d \
./boinc/boinc/client/mac_address.d \
./boinc/boinc/client/main.d \
./boinc/boinc/client/net_stats.d \
./boinc/boinc/client/pers_file_xfer.d \
./boinc/boinc/client/project.d \
./boinc/boinc/client/result.d \
./boinc/boinc/client/rr_sim.d \
./boinc/boinc/client/rrsim_test.d \
./boinc/boinc/client/sandbox.d \
./boinc/boinc/client/scheduler_op.d \
./boinc/boinc/client/setprojectgrp.d \
./boinc/boinc/client/sim.d \
./boinc/boinc/client/sim_util.d \
./boinc/boinc/client/switcher.d \
./boinc/boinc/client/sysmon_win.d \
./boinc/boinc/client/thread.d \
./boinc/boinc/client/time_stats.d \
./boinc/boinc/client/whetstone.d \
./boinc/boinc/client/work_fetch.d 


# Each subdirectory must supply rules for building sources it contributes
boinc/boinc/client/%.o: ../boinc/boinc/client/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


