################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../milkyway/src/milkyway_boinc_util.cc 

C_SRCS += \
../milkyway/src/milkyway_alloc.c \
../milkyway/src/milkyway_asprintf.c \
../milkyway/src/milkyway_cl_device.c \
../milkyway/src/milkyway_cl_program.c \
../milkyway/src/milkyway_cl_setup.c \
../milkyway/src/milkyway_cl_show_types.c \
../milkyway/src/milkyway_cl_util.c \
../milkyway/src/milkyway_cpuid.c \
../milkyway/src/milkyway_lua_dsfmt.c \
../milkyway/src/milkyway_lua_marshal.c \
../milkyway/src/milkyway_lua_math.c \
../milkyway/src/milkyway_lua_types.c \
../milkyway/src/milkyway_lua_util.c \
../milkyway/src/milkyway_lua_vector.c \
../milkyway/src/milkyway_rename.c \
../milkyway/src/milkyway_show.c \
../milkyway/src/milkyway_timing.c \
../milkyway/src/milkyway_util.c 

CC_DEPS += \
./milkyway/src/milkyway_boinc_util.d 

OBJS += \
./milkyway/src/milkyway_alloc.o \
./milkyway/src/milkyway_asprintf.o \
./milkyway/src/milkyway_boinc_util.o \
./milkyway/src/milkyway_cl_device.o \
./milkyway/src/milkyway_cl_program.o \
./milkyway/src/milkyway_cl_setup.o \
./milkyway/src/milkyway_cl_show_types.o \
./milkyway/src/milkyway_cl_util.o \
./milkyway/src/milkyway_cpuid.o \
./milkyway/src/milkyway_lua_dsfmt.o \
./milkyway/src/milkyway_lua_marshal.o \
./milkyway/src/milkyway_lua_math.o \
./milkyway/src/milkyway_lua_types.o \
./milkyway/src/milkyway_lua_util.o \
./milkyway/src/milkyway_lua_vector.o \
./milkyway/src/milkyway_rename.o \
./milkyway/src/milkyway_show.o \
./milkyway/src/milkyway_timing.o \
./milkyway/src/milkyway_util.o 

C_DEPS += \
./milkyway/src/milkyway_alloc.d \
./milkyway/src/milkyway_asprintf.d \
./milkyway/src/milkyway_cl_device.d \
./milkyway/src/milkyway_cl_program.d \
./milkyway/src/milkyway_cl_setup.d \
./milkyway/src/milkyway_cl_show_types.d \
./milkyway/src/milkyway_cl_util.d \
./milkyway/src/milkyway_cpuid.d \
./milkyway/src/milkyway_lua_dsfmt.d \
./milkyway/src/milkyway_lua_marshal.d \
./milkyway/src/milkyway_lua_math.d \
./milkyway/src/milkyway_lua_types.d \
./milkyway/src/milkyway_lua_util.d \
./milkyway/src/milkyway_lua_vector.d \
./milkyway/src/milkyway_rename.d \
./milkyway/src/milkyway_show.d \
./milkyway/src/milkyway_timing.d \
./milkyway/src/milkyway_util.d 


# Each subdirectory must supply rules for building sources it contributes
milkyway/src/%.o: ../milkyway/src/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

milkyway/src/%.o: ../milkyway/src/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


